# This software is Copyright ©2026. The Regents of the University of California
# ("Regents"). All Rights Reserved.
#
# See LICENSE.txt for license details.

import os
import sys
import subprocess
import tempfile
import shutil
import pysam
from Bio.Seq import Seq
from Bio.SeqRecord import SeqRecord
from Bio import SeqIO
from . import config
from .utils import run_quiet

def filter_vcf_gene(input_vcf, gene, filter_region, symbolic_vcf, pass_vcf, fail_vcf, sv_overlap_vcf, pass_unphased, filtered_vcf, platform, genotyper, force_include_unphased=False):
	# Extract region
	base = os.path.basename(filtered_vcf)
	prefix = base.replace("_PASS_phased.vcf.gz", "")
	region_vcf = os.path.join(os.path.dirname(filtered_vcf), f"{prefix}.vcf.gz")

	cmd = f"bcftools view -r {filter_region} {input_vcf} -Oz -o {region_vcf}"
	run_quiet(cmd)
	run_quiet(f"bcftools index -f {region_vcf}")

	# ========== FIRST PASS: Collect PASS pbsv SV regions with haplotype info ==========
	sv_regions = []  # list of (start, end, affected_haplotypes)
	vf_sv_pass = pysam.VariantFile(region_vcf)

	for rec in vf_sv_pass:
		# Check if this is a pbsv SV call (ID starts with "pbsv.") or a Sniffles/Sniffles2 call
		rec_id = rec.id or ""
		rec_id_l = rec.id.lower() if rec.id else ""
		if rec_id_l.startswith("pbsv.") or rec_id_l.startswith("sniffles"):
			# Only include PASS SVs
			filter_keys = list(rec.filter.keys())
			if filter_keys == [] or filter_keys == ["PASS"]:
				# Skip symbolic ALTs (BND breakends, <DUP>, <INV>, etc.)
				# These have unreliable coordinate spans that could incorrectly suppress real variants
				if rec.alts is None:
					continue
				is_symbolic = False
				for alt in rec.alts:
					alt_str = str(alt)
					# Check for symbolic notation <TYPE> or BND bracket notation
					if alt_str.startswith("<") or "]" in alt_str or "[" in alt_str:
						is_symbolic = True
						break
				if is_symbolic:
					continue

				sample = list(rec.samples.values())[0]
				gt = sample.get("GT")

				# Only use phased or homozygous SVs for overlap suppression
				# Unphased het SVs have ambiguous haplotype assignment
				if gt is None or None in gt:
					continue

				is_het = len(set(gt)) > 1
				if is_het and not sample.phased:
					# Unphased het SV - don't use for overlap suppression
					continue

				start = rec.pos
				end = rec.info.get("END", rec.pos + len(rec.ref) - 1)

				# Which haplotypes carry the SV? (indices where allele > 0)
				affected_haps = set()
				for i, allele in enumerate(gt):
					if allele is not None and allele > 0:
						affected_haps.add(i)

				if affected_haps:
					sv_regions.append((start, end, affected_haps))

	vf_sv_pass.close()

	if config.VERBOSE:
		print(f"[SV-OVERLAP] {gene}: Collected {len(sv_regions)} PASS pbsv SV regions for overlap checking")
		for sv_start, sv_end, sv_haps in sv_regions:
			print(f"[SV-OVERLAP]   SV region: {sv_start}-{sv_end}, haplotypes: {sv_haps}")

	# ========== FIRST PASS (part 2): Collect TRGT tandem repeat regions ==========
	# TRGT records have explicit REF/ALT allele sequences. Any bcftools-called
	# variants that fall inside a TRGT span are redundant (the TRGT allele already
	# captures the correct sequence for the entire repeat region).
	tr_regions = []  # list of (start, end) — TRGT regions to protect
	vf_tr_pass = pysam.VariantFile(region_vcf)
	for rec in vf_tr_pass:
		if "TRID" in rec.info and rec.info["TRID"] not in (None, "", "."):
			# TRGT record: use END if available, otherwise POS + len(REF)
			tr_start = rec.pos
			tr_end = rec.info.get("END", rec.pos + len(rec.ref))
			tr_regions.append((tr_start, tr_end))
	vf_tr_pass.close()

	if config.VERBOSE:
		print(f"[TR-OVERLAP] {gene}: Collected {len(tr_regions)} TRGT regions for overlap suppression")
		for tr_start, tr_end in tr_regions:
			print(f"[TR-OVERLAP]   TR region: {tr_start}-{tr_end}")

	# ========== FIRST PASS (part 3): Collect homozygous deletion spans ==========
	# A base inside a homozygous deletion is absent from BOTH haplotypes, so any other
	# variant called there (typically a low-depth het SNP from reads that misaligned across
	# the deletion) is spurious. pbsv/TRGT spans are handled above; this catches small
	# homozygous deletions from the small-variant caller (too small for pbsv, not in a TRGT
	# catalog), whose deleted bases would otherwise carry junk het calls into vcf2fasta.
	hom_del_spans = []  # (del_start, del_end) inclusive -- bases deleted on both haplotypes
	vf_hd = pysam.VariantFile(region_vcf)
	for rec in vf_hd:
		if rec.alts is None:
			continue
		hd_sample = list(rec.samples.values())[0]
		hd_gt = hd_sample.get("GT")
		if hd_gt is None or None in hd_gt or len(set(hd_gt)) != 1 or hd_gt[0] == 0:
			continue  # require a homozygous ALT genotype
		hd_alt = rec.alleles[hd_gt[0]]
		if len(hd_alt) < len(rec.ref):
			hom_del_spans.append((rec.pos + 1, rec.pos + len(rec.ref) - 1))
	vf_hd.close()

	if config.VERBOSE:
		print(f"[HOMDEL-OVERLAP] {gene}: Collected {len(hom_del_spans)} homozygous deletion spans")
		for hd_start, hd_end in hom_del_spans:
			print(f"[HOMDEL-OVERLAP]   deleted span: {hd_start}-{hd_end}")

	# Helper function: is a position inside any homozygous deletion span?
	def inside_homozygous_deletion(pos):
		for hd_start, hd_end in hom_del_spans:
			if hd_start <= pos <= hd_end:
				return True
		return False

	# Helper function: check if variant overlaps a PASS SV on the same haplotype
	def overlaps_sv_same_haplotype(pos, ref_len, var_haplotypes, indel_size=0):
		var_end = pos + ref_len - 1
		for sv_start, sv_end, sv_haplotypes in sv_regions:
			if not (var_haplotypes & sv_haplotypes):
				continue
			# Strict positional overlap
			if pos <= sv_end and var_end >= sv_start:
				return True
			# Large indels (>=50bp) within 10bp of an SV are almost certainly
			# the same event with different left-alignment by bcftools vs pbsv
			if indel_size >= 50 and abs(pos - sv_start) <= 10:
				return True
		return False

	# Helper function: check if a non-TRGT variant falls inside a TRGT region
	def overlaps_tr_region(pos, ref_len):
		var_end = pos + ref_len
		for tr_start, tr_end in tr_regions:
			if pos >= tr_start and var_end <= tr_end:
				return True
		return False

	vf = pysam.VariantFile(region_vcf)
	sym_out = pysam.VariantFile(symbolic_vcf, "wz", header=vf.header)
	pass_out = pysam.VariantFile(pass_vcf, "wz", header=vf.header)
	fail_out = pysam.VariantFile(fail_vcf, "wz", header=vf.header)
	sv_overlap_out = pysam.VariantFile(sv_overlap_vcf, "wz", header=vf.header)
	sv_overlap_count = 0  # Counter for suppressed variants

	# ========== PASS/FAIL CLASSIFICATION ==========
	for rec in vf:
		sample = list(rec.samples.values())[0]
		gt = sample.get("GT")

		# HARD EXCLUDE — must not exist downstream
		if gt is None or None in gt:
			continue

		# symbolic — exclude truly symbolic alleles (BND, <DEL>, etc.)
		# TRGT records with explicit REF/ALT sequences are NOT symbolic
		if (
			rec.alts is None or
			any(str(a).startswith("<") for a in rec.alts) or
			any("]" in str(a) or "[" in str(a) for a in rec.alts) or
			any(set(str(a)) - set("ACGTN") for a in rec.alts)
		):
			sym_out.write(rec)
			continue

		# pbsv SV (ID starts with "pbsv.") or Sniffles/Sniffles2 SVs
		rec_id = rec.id or ""
		rec_id_l = rec_id.lower()
		if rec_id_l.startswith("pbsv.") or rec_id_l.startswith("sniffles"):
			if rec.filter.keys() == ["PASS"] or rec.filter.keys() == []:
				pass_out.write(rec)
			else:
				fail_out.write(rec)
			continue

		# ========== HOMOZYGOUS DELETION OVERLAP CHECK ==========
		# A heterozygous call at a homozygously-deleted base cannot be real (the base is absent
		# from both haplotypes). Suppress it so vcf2fasta does not apply a spurious variant.
		# Homozygous variants are left alone -- they define the sequence.
		if hom_del_spans and len(set(gt)) > 1 and inside_homozygous_deletion(rec.pos):
			if config.VERBOSE:
				print(f"[HOMDEL-OVERLAP]   Suppressed: {rec.chrom}:{rec.pos} {rec.ref[:20]}->{rec.alts[0][:20] if rec.alts else '.'} (het inside homozygous deletion)")
			sv_overlap_out.write(rec)
			sv_overlap_count += 1
			continue

		# ========== TR OVERLAP CHECK ==========
		# Suppress non-TRGT variants that fall entirely inside a TRGT region.
		# The TRGT record already encodes the correct allele for the whole repeat;
		# bcftools-called SNVs/indels inside that span are redundant and would
		# cause double-counting in vcf2fasta.
		is_trgt = "TRID" in rec.info and rec.info["TRID"] not in (None, "", ".")
		if not is_trgt and tr_regions and overlaps_tr_region(rec.pos, len(rec.ref)):
			if config.VERBOSE:
				print(f"[TR-OVERLAP]   Suppressed: {rec.chrom}:{rec.pos} {rec.ref[:20]}->{rec.alts[0][:20] if rec.alts else '.'} (inside TRGT span)")
			sv_overlap_out.write(rec)
			sv_overlap_count += 1
			continue

		# TRGT records with explicit alleles bypass quality filters (no DP/GQ)
		if is_trgt:
			if config.VERBOSE:
				sample_gt = sample.get("GT")
				gt_str = "|".join(str(a) for a in sample_gt) if sample_gt else "."
				alts_str = ",".join(str(a)[:30] for a in rec.alts) if rec.alts else "."
				print(f"[TR-PASS]      {gene} {rec.chrom}:{rec.pos} REF={rec.ref[:30]} ALT={alts_str} GT={gt_str} TRID={rec.info.get('TRID','?')}")
			pass_out.write(rec)
			continue

		# ========== SV OVERLAP CHECK ==========
		# Skip non-SV variants that overlap a PASS SV on the same haplotype
		# Only apply haplotype-aware suppression for phased indels ≥30bp
		# Never suppress SNPs - they are real variants near SV breakpoints
		# Size and type come from the carried alleles, not rec.alts[0]. DeepVariant
		# emits multiallelic records holding both a SNP and an indel alt (T -> A,TG).
		ref_len = len(rec.ref)
		gt_alleles = [rec.alleles[a] for a in gt if a is not None and a < len(rec.alleles)]
		if not gt_alleles:
			gt_alleles = [a for a in rec.alts if a is not None]
		indel_size = max((abs(ref_len - len(a)) for a in gt_alleles), default=0)
		is_snp = ref_len == 1 and all(len(a) == 1 for a in gt_alleles)

		is_homozygous = gt is not None and len(set(a for a in gt if a is not None)) == 1
		if (sample.phased or is_homozygous) and not is_snp:
			var_haplotypes = set()
			for i, allele in enumerate(gt):
				if allele is not None and allele > 0:
					var_haplotypes.add(i)

			if var_haplotypes and overlaps_sv_same_haplotype(rec.pos, len(rec.ref), var_haplotypes, indel_size):
				sv_overlap_count += 1
				if config.VERBOSE:
					print(f"[SV-OVERLAP]   Suppressed indel ({indel_size}bp): {rec.chrom}:{rec.pos} {rec.ref[:20]}...>{rec.alts[0][:20]}... GT={gt} haps={var_haplotypes}")
				sv_overlap_out.write(rec)
				continue

		# DeepVariant/Clair3 specific filtering; for hybrid, apply PASS-only filter to indels from DeepVariant
		is_snp = len(rec.ref) == 1 and all(len(a) == 1 for a in rec.alts if a is not None)
		if genotyper in ("deepvariant", "clair3") or (genotyper == "hybrid" and not is_snp):
			if rec.filter.keys() == ["PASS"] or rec.filter.keys() == []:
				pass_out.write(rec)
			else:
				fail_out.write(rec)
			continue

		# small variants (QUAL/GQ/DP already filtered upstream by caller)
		pass_out.write(rec)

	sym_out.close()
	pass_out.close()
	fail_out.close()
	sv_overlap_out.close()

	if config.VERBOSE:
		print(f"[SV-OVERLAP] {gene}: Total small variants suppressed due to SV overlap: {sv_overlap_count}")

	run_quiet(f"bcftools index -f {symbolic_vcf}")
	run_quiet(f"bcftools index -f {pass_vcf}")
	run_quiet(f"bcftools index -f {fail_vcf}")
	run_quiet(f"bcftools index -f {sv_overlap_vcf}")

	# ========== WHITELIST LOGIC (RESTORED) ==========
	def is_unphased_het(rec):
		"""True for a heterozygous genotype left unphased by HiPhase/LongPhase.

		Tested on the parsed genotype so multiallelic hets (0/2, 1/3, ...) are
		caught, not just the common 0/1 and 1/2 forms.
		"""
		sample = list(rec.samples.values())[0]
		gt = sample.get("GT")
		if gt is None or None in gt:
			return False
		if len(set(gt)) < 2:
			return False
		return not sample.phased

	het_sites = []
	unphased_hets = []

	with pysam.VariantFile(pass_vcf) as pass_vf:
		for rec in pass_vf:
			sample = list(rec.samples.values())[0]
			gt = sample.get("GT")
			if gt is None or None in gt:
				continue
			if len(set(gt)) >= 2:  # heterozygous
				het_sites.append(rec)
				if not sample.phased:
					unphased_hets.append(rec)

	if config.VERBOSE:
		print(f"{gene}: het={len(het_sites)}, unphased={len(unphased_hets)}")
	allow_single_unphased = (len(het_sites) == 1 and len(unphased_hets) == 1)

	# Unphased hets are kept for CDS-rescued genes, and for genes whose only het is
	# that one site, where haplotype assignment is arbitrary but the pair is correct.
	drop_unphased_hets = not (force_include_unphased or allow_single_unphased)

	if config.VERBOSE:
		if force_include_unphased:
			print(f"[CDS-RESCUE] {gene}: force_include_unphased=True, keeping all variants")
		elif allow_single_unphased:
			print(f"{gene}: single unphased het whitelisted, keeping all variants")

	# Split the QC-pass records: unphased hets to pass_unphased (reported, not
	# reconstructed), everything else to the filtered VCF that feeds vcf2fasta.
	n_routed = 0
	with pysam.VariantFile(pass_vcf) as vf_in:
		with pysam.VariantFile(pass_unphased, "wz", header=vf_in.header) as unph_out, \
			 pysam.VariantFile(filtered_vcf, "wz", header=vf_in.header) as filt_out:
			for rec in vf_in:
				if drop_unphased_hets and is_unphased_het(rec):
					unph_out.write(rec)
					n_routed += 1
				else:
					filt_out.write(rec)

	run_quiet(f"bcftools index -f {pass_unphased}")
	run_quiet(f"bcftools index -f {filtered_vcf}")

	if config.VERBOSE and n_routed != len(unphased_hets) and drop_unphased_hets:
		print(f"WARNING: {gene}: routed {n_routed} unphased hets but counted {len(unphased_hets)}")

	# ========== UNPHASED PASS SUMMARY ==========
	if unphased_hets:
		print(f"{gene}: {len(unphased_hets)} unphased PASS het(s) → {pass_unphased}\n")

	if config.VERBOSE:
		if force_include_unphased and unphased_hets:
			print(f"\nUnphased PASS variants in {gene}:\n")
			for rec in unphased_hets:
				print(str(rec).strip())
			print("\n")
		else:
			unph = pysam.VariantFile(pass_unphased)
			records = [rec for rec in unph]

			if records:
				print(f"\nUnphased PASS variants in {gene}:\n")
				for rec in records:
					print(str(rec).strip())
				print("\n")

def read_gene_gff_cols(gff_dir, gene_lower):
	"""Return the 9 tab-separated columns of the 'gene' feature line from
	{gene_lower}_gene.gff3 (the same file used for full-gene vcf2fasta runs)."""
	path = os.path.join(gff_dir, f"{gene_lower}_gene.gff3")
	with open(path) as f:
		for line in f:
			if line.startswith("#"):
				continue
			cols = line.rstrip("\n").split("\t")
			if len(cols) >= 8 and cols[2] == "gene":
				return cols
	raise RuntimeError(f"No 'gene' feature line found in {path}")


def variant_ref_spans(vcf_path):
	"""Reference footprints an interval boundary must not land inside. vcf2fasta emits a
	feature record TWICE (a 2x "doubling") whenever the feature start or stop falls inside
	the REF span of a record that covers more than one reference base (a deletion or MNP).
	This is independent of genotype: hom-ref (0/0) records trigger it too, because the
	doubling is a parser artifact of the overlapping REF span, not of the applied allele.
	Returns the full inclusive footprint (pos, pos+len(ref)-1) of every such record."""
	spans = []
	vcf = pysam.VariantFile(vcf_path)
	for rec in vcf:
		if len(rec.ref) > 1:
			spans.append((rec.pos, rec.pos + len(rec.ref) - 1))
	vcf.close()
	return spans


def snap_left_out_of_variant_span(pos, spans):
	"""Push a left/start boundary forward to the first reference base outside every
	variant REF span. A feature start inside a multi-base REF span makes vcf2fasta emit
	the record twice (doubling); it also mis-indexed the old clamp."""
	moved = True
	while moved:
		moved = False
		for span_start, span_end in spans:
			if span_start <= pos <= span_end:
				pos = span_end + 1
				moved = True
	return pos


def snap_right_out_of_variant_span(pos, spans):
	"""Push a right/stop boundary back to the last reference base outside every variant REF span."""
	moved = True
	while moved:
		moved = False
		for span_start, span_end in spans:
			if span_start <= pos <= span_end:
				pos = span_start - 1
				moved = True
	return pos


def hap_net_indel(vcf_path, haplotype, start, stop):
	"""Net length change (sum of len(alt)-len(ref)) of variants on this haplotype
	whose POS falls within [start, stop]. Used to compute the expected length of a
	clean interval extraction for the doubling guard."""
	net = 0
	vcf = pysam.VariantFile(vcf_path)
	for rec in vcf:
		if rec.pos < start or rec.pos > stop:
			continue
		gt = list(rec.samples.values())[0].get("GT")
		if gt is None or None in gt:
			continue
		allele_idx = gt[haplotype]
		if allele_idx == 0:
			continue
		net += len(rec.alleles[allele_idx]) - len(rec.ref)
	vcf.close()
	return net


def extract_interval_vcf2fasta(gene, gene_lower, clamped_start, clamped_stop, gff_dir, reference_genome, vcf_path, vcf2fasta_output_dir):
	"""Rebuild a sub-interval of both haplotypes by re-running vcf2fasta on a
	one-feature GFF for [clamped_start, clamped_stop].

	Boundaries are snapped out of the REF span of any multi-base record, since a
	feature edge inside such a span makes vcf2fasta emit the record twice.

	Returns (allele_1, allele_2), strand-oriented by vcf2fasta, or ("", "") if the
	snapped interval collapses.
	"""
	left = min(clamped_start, clamped_stop)
	right = max(clamped_start, clamped_stop)

	spans = variant_ref_spans(vcf_path)
	snapped_left = snap_left_out_of_variant_span(left, spans)
	snapped_right = snap_right_out_of_variant_span(right, spans)
	if snapped_left != left or snapped_right != right:
		print(f"{gene}: snapped interval boundary out of indel "
			f"[{left}-{right}] -> [{snapped_left}-{snapped_right}]")
	# Write the contracted extraction into the persistent vcf2fasta_out location, replacing the
	# full-gene first-round output, so that dir holds the sequence actually used for matching
	# (the ARS-spanning interval) rather than the untrusted full gene. vcf2fasta appends "_gene"
	# to -o, so out_base -> real_out = out_base + "_gene".
	out_base = os.path.join(vcf2fasta_output_dir, gene)
	real_out = out_base + "_gene"
	shutil.rmtree(real_out, ignore_errors=True)

	if snapped_left > snapped_right:
		return "", ""

	cols = read_gene_gff_cols(gff_dir, gene_lower)
	cols[3] = str(snapped_left)
	cols[4] = str(snapped_right)

	# The interval GFF is scratch; only the vcf2fasta output (in real_out) is kept.
	gff_workdir = tempfile.mkdtemp(prefix=f"{gene_lower}_interval_gff_")
	try:
		gff_path = os.path.join(gff_workdir, f"{gene_lower}_interval.gff3")
		with open(gff_path, "w") as f:
			f.write("##gff-version 3\n")
			f.write("\t".join(cols) + "\n")

		run_vcf2fasta(
			input_vcf=vcf_path,
			input_gff=gff_path,
			reference_genome=reference_genome,
			output_dir=out_base,
			gene=gene,
			feature="gene",
		)

		find_cmd = f"find {real_out} -type f"
		result = subprocess.run(find_cmd, shell=True, capture_output=True, text=True)
		fasta_files = [x for x in result.stdout.split() if not x.endswith(".gff3")]
		if not fasta_files:
			raise FileNotFoundError(
				f"{gene}: no vcf2fasta output for interval {snapped_left}-{snapped_right}")

		with open(fasta_files[0]) as f:
			records = f.read().split(">")[1:]
		haplotypes = ["".join(rec.split("\n")[1:]).replace("-", "").strip() for rec in records]
		if len(haplotypes) < 2:
			raise ValueError(
				f"{gene}: expected 2 haplotype records from interval extraction, got {len(haplotypes)}")

		ref_len = snapped_right - snapped_left + 1
		for hap in (0, 1):
			expected = ref_len + hap_net_indel(vcf_path, hap, snapped_left, snapped_right)
			if len(haplotypes[hap]) > 1.5 * max(expected, 1):
				raise RuntimeError(
					f"{gene}: vcf2fasta interval extraction doubled "
					f"(hap{hap} len={len(haplotypes[hap])}, expected ~{expected}, "
					f"interval {snapped_left}-{snapped_right}). Boundary may sit inside an indel.")

		return haplotypes[0], haplotypes[1]
	finally:
		shutil.rmtree(gff_workdir, ignore_errors=True)


def extract_cds_subset_vcf2fasta(gene, gene_lower, keep_ranges, gff_dir, reference_genome, vcf_path, vcf2fasta_output_dir):
	"""Reconstruct a subset of a gene's CDS by clipping the CDS GFF to keep_ranges
	and re-running vcf2fasta.

	keep_ranges is [(haploblock_start, haploblock_end)] for the unphased trim, or
	the ARS CDS exon ranges for the ars_only tier. Output replaces the full-CDS
	first round in vcf2fasta_out/{gene}_CDS. Returns ("", "") if nothing is kept.
	"""
	cds_gff = os.path.join(gff_dir, f"{gene_lower}_cds_sorted.gff3")
	clipped = []  # (start, end, gff_line) for CDS records clipped to keep_ranges
	for raw in open(cds_gff):
		if raw.startswith("#"):
			continue
		cols = raw.rstrip("\n").split("\t")
		if len(cols) < 9 or cols[2] != "CDS":
			continue
		s, e = int(cols[3]), int(cols[4])
		for ks, ke in keep_ranges:
			cs, ce = max(s, ks), min(e, ke)
			if cs <= ce:
				rc = list(cols)
				rc[3], rc[4] = str(cs), str(ce)
				clipped.append((cs, ce, "\t".join(rc)))

	out_base = os.path.join(vcf2fasta_output_dir, gene)
	real_out = out_base + "_CDS"
	shutil.rmtree(real_out, ignore_errors=True)  # drop the full-CDS first-round output
	if not clipped:
		return "", ""
	clipped.sort()

	gff_workdir = tempfile.mkdtemp(prefix=f"{gene_lower}_cdssubset_gff_")
	try:
		gff_path = os.path.join(gff_workdir, f"{gene_lower}_cds_subset.gff3")
		with open(gff_path, "w") as f:
			f.write("##gff-version 3\n")
			for _, _, line in clipped:
				f.write(line + "\n")

		run_vcf2fasta(
			input_vcf=vcf_path,
			input_gff=gff_path,
			reference_genome=reference_genome,
			output_dir=out_base,
			gene=gene,
			feature="CDS",
		)

		find = subprocess.run(f"find {real_out} -type f", shell=True, capture_output=True, text=True)
		fasta_files = [x for x in find.stdout.split() if not x.endswith(".gff3")]
		if not fasta_files:
			raise FileNotFoundError(f"{gene}: no vcf2fasta CDS output for subset")

		with open(fasta_files[0]) as f:
			records = f.read().split(">")[1:]
		haplotypes = ["".join(rec.split("\n")[1:]).replace("-", "").strip() for rec in records]
		if len(haplotypes) < 2:
			raise ValueError(
				f"{gene}: expected 2 CDS haplotype records from subset extraction, got {len(haplotypes)}")

		expected = sum(ce - cs + 1 for cs, ce, _ in clipped)
		for hap in haplotypes:
			if len(hap) > 1.8 * max(expected, 1):
				raise RuntimeError(
					f"{gene}: vcf2fasta CDS subset extraction doubled (len={len(hap)}, "
					f"expected ~{expected}). A clipped CDS boundary may sit inside an indel.")

		return haplotypes[0], haplotypes[1]
	finally:
		shutil.rmtree(gff_workdir, ignore_errors=True)


def run_vcf2fasta(input_vcf, input_gff, reference_genome, output_dir, gene, feature):
	if feature == "CDS":
		vcf2fasta_cmd = f"vcf2fasta --fasta {reference_genome} --vcf {input_vcf} --gff {input_gff} -o {output_dir} --feat CDS --blend --force"
	elif feature == "gene":
		vcf2fasta_cmd = f"vcf2fasta --fasta {reference_genome} --vcf {input_vcf} --gff {input_gff} -o {output_dir} --feat gene --force"

	try:
		subprocess.run(vcf2fasta_cmd, shell=True, check=True, capture_output=True, text=True)
	except subprocess.CalledProcessError as e:
		if e.stdout:
			print(e.stdout, end="")
		if e.stderr:
			print(e.stderr, end="", file=sys.stderr)
		raise

def extract_cds_padded_and_trim(gene, gene_lower, gff_dir, reference_genome, vcf_path,
								vcf2fasta_output_dir, stop_codons, pad=5, accept_window=5):
	"""Repair a CDS that reconstructs without a stop codon.

	Some alleles (DPB1*13:01, *27:01, *107:01, *135:01) drop a base from a poly-A
	run at the CDS 3' terminus, so the fixed-coordinate window ends one base short
	of the shifted stop. Re-runs vcf2fasta on a CDS GFF whose 3'-terminal exon is
	extended by `pad` bases, then trims to the first in-frame stop.

	Returns (hap1, hap2). A haplotype is None if no stop was recovered within
	accept_window of the annotated CDS length, and the caller keeps the original.
	"""
	cds_gff = os.path.join(gff_dir, f"{gene_lower}_cds_sorted.gff3")
	recs = []
	strand = None
	for raw in open(cds_gff):
		if raw.startswith("#"):
			continue
		cols = raw.rstrip("\n").split("\t")
		if len(cols) < 9 or cols[2] != "CDS":
			continue
		recs.append((int(cols[3]), int(cols[4]), cols))
		strand = cols[6]
	if not recs or strand not in ("+", "-"):
		return None, None
	ann_len = sum(e - s + 1 for s, e, _ in recs)

	# Extend the transcript's 3'-terminal exon by `pad` reference bases: the largest-end
	# exon on the + strand, the smallest-start exon on the - strand.
	if strand == "+":
		idx = max(range(len(recs)), key=lambda i: recs[i][1])
		recs[idx][2][4] = str(recs[idx][1] + pad)
	else:
		idx = min(range(len(recs)), key=lambda i: recs[i][0])
		recs[idx][2][3] = str(max(1, recs[idx][0] - pad))

	gff_workdir = tempfile.mkdtemp(prefix=f"{gene_lower}_cdspad_gff_")
	pad_out = tempfile.mkdtemp(prefix=f"{gene_lower}_cdspad_out_")
	try:
		gff_path = os.path.join(gff_workdir, f"{gene_lower}_cds_padded.gff3")
		with open(gff_path, "w") as f:
			f.write("##gff-version 3\n")
			for _, _, cols in sorted(recs, key=lambda r: (r[0], r[1])):
				f.write("\t".join(cols) + "\n")

		run_vcf2fasta(
			input_vcf=vcf_path,
			input_gff=gff_path,
			reference_genome=reference_genome,
			output_dir=os.path.join(pad_out, gene),
			gene=gene,
			feature="CDS",
		)
		real_out = os.path.join(pad_out, gene) + "_CDS"
		find = subprocess.run(f"find {real_out} -type f", shell=True, capture_output=True, text=True)
		fasta_files = [x for x in find.stdout.split() if not x.endswith(".gff3")]
		if not fasta_files:
			return None, None
		with open(fasta_files[0]) as f:
			records = f.read().split(">")[1:]
		haps = ["".join(rec.split("\n")[1:]).replace("-", "").strip() for rec in records]

		def trim(seq):
			# First in-frame stop codon; accept only within +/-accept_window of ann_len,
			# else abort (a premature stop or nothing near the terminus -> keep original).
			for p in range(0, len(seq) - 2, 3):
				if seq[p:p + 3] in stop_codons:
					end = p + 3
					if ann_len - accept_window <= end <= ann_len + accept_window:
						return seq[:end]
					return None
			return None

		h1 = trim(haps[0]) if len(haps) > 0 else None
		h2 = trim(haps[1]) if len(haps) > 1 else None
		return h1, h2
	finally:
		shutil.rmtree(gff_workdir, ignore_errors=True)
		shutil.rmtree(pad_out, ignore_errors=True)

def parse_fastas(sample_ID, vcf2fasta_output_dir, outfile_gene, outfile_CDS, DNA_bases, stop_codons, unphased_genes=None, gene_dict=None, CDS_dict=None, gff_dir=None, cds_rescued_genes=None, ARS_dict=None, CLASS_I_GENES=None, gene_filtered_vcfs=None, reference_genome=None):
	# Use subprocess.run with capture_output to avoid race conditions with temporary files
	find_cmd = f"find {vcf2fasta_output_dir} -type f"
	result = subprocess.run(find_cmd, shell=True, check=True, capture_output=True, text=True)
	fasta_files = result.stdout.strip().split('\n') if result.stdout.strip() else []
	
	# Check if any FASTA files were found
	if not fasta_files:
		print(f"ERROR: No FASTA files found in {vcf2fasta_output_dir}")

		raise FileNotFoundError(f"No FASTA files found in {vcf2fasta_output_dir}")

	fasta_dict = dict()

	logging_strings = []
	for file in fasta_files:
		# Validate file exists and is readable
		if not os.path.exists(file):
			print(f"WARNING: File {file} does not exist, skipping")
			continue
		if not os.access(file, os.R_OK):
			print(f"WARNING: File {file} is not readable, skipping")
			continue
		
		if "_gene" in file:
			feat = "gene"
		elif "_CDS" in file:
			feat = "CDS"
		else:
			print(f"WARNING: File {file} does not contain '_gene' or '_CDS', skipping")
			continue
			
		gene = file.split("/")[-2].split(f"_{feat}")[0].upper().replace("_", "-")
		
		try:
			with open(file, "r") as f:
				lines = f.read().split(">")
		except (IOError, OSError) as e:
			print(f"ERROR: Failed to read file {file}: {e}")
			continue
		
		# Validate FASTA structure - should have at least 3 parts (empty, allele1, allele2)
		if len(lines) < 3:
			print(f"WARNING: File {file} has insufficient FASTA entries (expected 2 alleles, found {len(lines)-1}), skipping")
			continue
		
		# Old code
		# Remove deletion characters
		#allele_1 = lines[1].split("\n")[1].strip().replace("-","").strip()
		#allele_2 = lines[2].split("\n")[1].strip().replace("-","").strip()
		# Concatenate all lines after the header for each allele
		try:
			allele_1 = "".join(lines[1].split("\n")[1:]).replace("-", "").strip()
			allele_2 = "".join(lines[2].split("\n")[1:]).replace("-", "").strip()
		except IndexError as e:
			print(f"ERROR: Failed to parse FASTA alleles from {file}: {e}")
			continue

		if cds_rescued_genes and gene in cds_rescued_genes:
			rescue_info = cds_rescued_genes[gene]
			tier = rescue_info["tier"]
			all_hets = rescue_info["all_het_positions"]
			ars_start = rescue_info["ars_start"]
			ars_stop = rescue_info["ars_stop"]
			gene_start_coord, gene_stop_coord = gene_dict[gene]

			if tier == "cds_full":
				# Branch 3b-i: CDS hets <= 1
				if feat == "CDS":
					# No clamping — write full concatenated exon output
					logging_strings.append(f"{sample_ID} {gene} CDS rescue tier=cds_full: writing full CDS output")

				elif feat == "gene":
					cds_het_positions = rescue_info["cds_het_positions"]
					gene_lower = gene.lower().replace("-", "_")
					gene_coords_file = os.path.join(gff_dir, f"{gene_lower}_gene_coords.txt")
					gene_coords = [int(item) for item in open(gene_coords_file).read().splitlines()]

					interval = None
					if len(cds_het_positions) == 1:
						# Anchor on CDS het, extend to flanking hets
						cds_het = cds_het_positions[0]
						prev_hets = [h for h in all_hets if h < cds_het]
						next_hets = [h for h in all_hets if h > cds_het]
						left = max(prev_hets) + 1 if prev_hets else gene_start_coord
						right = min(next_hets) - 1 if next_hets else gene_stop_coord
						# Check if interval fully contains ARS
						if left <= ars_start and right >= ars_stop:
							interval = (left, right)
							logging_strings.append(f"{sample_ID} {gene} CDS rescue: CDS het at {cds_het}, interval chr6:{left}-{right} contains ARS")
						else:
							logging_strings.append(f"{sample_ID} {gene} CDS rescue: CDS het at {cds_het}, interval chr6:{left}-{right} does NOT contain ARS — skipping gene record")
					else:
						# CDS hets = 0: find largest ARS-overlapping 1-het interval
						best_interval = None
						best_size = 0
						for i in range(len(all_hets)):
							left = (all_hets[i-1] + 1) if i > 0 else gene_start_coord
							right = (all_hets[i+1] - 1) if i < len(all_hets) - 1 else gene_stop_coord
							# Check ARS overlap
							if left <= ars_stop and right >= ars_start:
								size = right - left
								if size > best_size:
									best_size = size
									best_interval = (left, right)
						if best_interval:
							interval = best_interval
							logging_strings.append(f"{sample_ID} {gene} CDS rescue: 0 CDS hets, best 1-het interval chr6:{interval[0]}-{interval[1]}")
						else:
							logging_strings.append(f"{sample_ID} {gene} CDS rescue: 0 CDS hets, no ARS-overlapping interval found — skipping gene record")

					if interval:
						clamped_start = max(interval[0], gene_start_coord)
						clamped_stop = min(interval[1], gene_stop_coord)
						vcf_path = gene_filtered_vcfs.get(gene) if gene_filtered_vcfs else None
						if vcf_path:
							allele_1, allele_2 = extract_interval_vcf2fasta(
								gene, gene_lower, clamped_start, clamped_stop, gff_dir, reference_genome, vcf_path, vcf2fasta_output_dir)
						else:
							idx1 = gene_coords.index(clamped_start)
							idx2 = gene_coords.index(clamped_stop)
							fasta_start, fasta_stop = sorted((idx1, idx2))
							allele_1 = allele_1[fasta_start:fasta_stop+1]
							allele_2 = allele_2[fasta_start:fasta_stop+1]
					else:
						allele_1, allele_2 = "", ""

			elif tier == "ars_only":
				# Branch 3b-ii-A: CDS hets > 1, ARS CDS hets <= 1
				if feat == "CDS":
					# Reconstruct only the ARS CDS exons via vcf2fasta on a clipped CDS GFF
					# (reference-anchored; replaces the cds_coords.index() slice that was off by
					# the in-CDS indel shift).
					gene_lower = gene.lower().replace("-", "_")
					ars_cds_ranges = rescue_info["ars_cds_ranges"]
					vcf_path = gene_filtered_vcfs.get(gene) if gene_filtered_vcfs else None
					if vcf_path and ars_cds_ranges:
						allele_1, allele_2 = extract_cds_subset_vcf2fasta(
							gene, gene_lower, ars_cds_ranges, gff_dir, reference_genome, vcf_path, vcf2fasta_output_dir)
						logging_strings.append(f"{sample_ID} {gene} CDS rescue tier=ars_only: extracted ARS CDS ({len(ars_cds_ranges)} exon range(s))")
					else:
						allele_1, allele_2 = "", ""
						logging_strings.append(f"{sample_ID} {gene} CDS rescue tier=ars_only: no ARS CDS ranges or VCF")

				elif feat == "gene":
					gene_lower = gene.lower().replace("-", "_")
					gene_coords_file = os.path.join(gff_dir, f"{gene_lower}_gene_coords.txt")
					gene_coords = [int(item) for item in open(gene_coords_file).read().splitlines()]

					skip_gene = False
					# Class I intron 2 check
					if CLASS_I_GENES and gene in CLASS_I_GENES:
						ars_cds_ranges = rescue_info["ars_cds_ranges"]
						if len(ars_cds_ranges) == 2:
							intron2_start = ars_cds_ranges[0][1] + 1
							intron2_stop = ars_cds_ranges[1][0] - 1
							intron2_hets = [h for h in all_hets if intron2_start <= h <= intron2_stop]
							if intron2_hets:
								skip_gene = True
								logging_strings.append(f"{sample_ID} {gene} CDS rescue tier=ars_only: het in intron 2 ({intron2_start}-{intron2_stop}) — skipping gene record")

					if skip_gene:
						allele_1, allele_2 = "", ""
					else:
						# Extend ARS outward to nearest het in each direction
						prev_hets = [h for h in all_hets if h < ars_start]
						next_hets = [h for h in all_hets if h > ars_stop]
						left = max(prev_hets) + 1 if prev_hets else gene_start_coord
						right = min(next_hets) - 1 if next_hets else gene_stop_coord

						clamped_start = max(left, gene_start_coord)
						clamped_stop = min(right, gene_stop_coord)
						vcf_path = gene_filtered_vcfs.get(gene) if gene_filtered_vcfs else None
						if vcf_path:
							allele_1, allele_2 = extract_interval_vcf2fasta(
								gene, gene_lower, clamped_start, clamped_stop, gff_dir, reference_genome, vcf_path, vcf2fasta_output_dir)
						else:
							idx1 = gene_coords.index(clamped_start)
							idx2 = gene_coords.index(clamped_stop)
							fasta_start, fasta_stop = sorted((idx1, idx2))
							allele_1 = allele_1[fasta_start:fasta_stop+1]
							allele_2 = allele_2[fasta_start:fasta_stop+1]
						logging_strings.append(f"{sample_ID} {gene} CDS rescue tier=ars_only: gene interval chr6:{clamped_start}-{clamped_stop}")

		elif unphased_genes and gene in unphased_genes:
			best_haploblock_start = unphased_genes[gene]["haploblock"][0]
			best_haploblock_end   = unphased_genes[gene]["haploblock"][1]

			if feat == "gene":
				# Load gene coords (1-based genomic positions)
				gene_lower = gene.lower().replace("-", "_")
				gene_coords_file = os.path.join(gff_dir, f"{gene_lower}_gene_coords.txt")
				gene_coords = [int(item) for item in open(gene_coords_file).read().splitlines()]

				clamped_start = max(best_haploblock_start, gene_dict[gene][0])
				clamped_stop  = min(best_haploblock_end,   gene_dict[gene][1])

				vcf_path = gene_filtered_vcfs.get(gene) if gene_filtered_vcfs else None
				if vcf_path:
					allele_1, allele_2 = extract_interval_vcf2fasta(
						gene, gene_lower, clamped_start, clamped_stop, gff_dir, reference_genome, vcf_path, vcf2fasta_output_dir)
				else:
					idx1 = gene_coords.index(clamped_start)
					idx2 = gene_coords.index(clamped_stop)
					fasta_start, fasta_stop = sorted((idx1, idx2))
					allele_1 = allele_1[fasta_start:fasta_stop+1]
					allele_2 = allele_2[fasta_start:fasta_stop+1]

			elif feat == "CDS":
				gene_lower = gene.lower().replace("-", "_")
				# If the CDS is effectively homozygous (<=1 CDS het), use the full untrimmed CDS
				# (mirrors the cds_full tier) so fields 1-3 are unaffected by the partial phasing.
				# Otherwise reconstruct the CDS clipped to the haploblock via vcf2fasta on a clipped
				# CDS GFF (reference-anchored; replaces the cds_coords.index() slice that was off by
				# the in-CDS indel shift).
				if unphased_genes[gene]["cds_hets"] > 1:
					vcf_path = gene_filtered_vcfs.get(gene) if gene_filtered_vcfs else None
					if vcf_path:
						allele_1, allele_2 = extract_cds_subset_vcf2fasta(
							gene, gene_lower, [(best_haploblock_start, best_haploblock_end)],
							gff_dir, reference_genome, vcf_path, vcf2fasta_output_dir)
					else:
						allele_1, allele_2 = "", ""

				pass_cds_counter = 0
				for cds_start, cds_stop in CDS_dict[gene]:
					if cds_stop < best_haploblock_start or cds_start > best_haploblock_end:
						status = "outside haploblock"
					elif cds_start >= best_haploblock_start and cds_stop <= best_haploblock_end:
						status = "fully contained"
						pass_cds_counter += 1
					else:
						status = "partially overlapping"
					logging_strings.append(f"{sample_ID} {gene} CDS {cds_start}-{cds_stop} is {status}")
				logging_strings.append(f"{sample_ID} {gene}: {pass_cds_counter} CDS fully contained in haploblock")


		if len(allele_1) == 0 or len(allele_2) == 0:
			print(f"File {file} has no sequence -- removing untrusted vcf2fasta output")
			try:
				os.remove(file)
			except OSError:
				pass
			continue
			
		if not set(allele_1).issubset(DNA_bases):
			print(f"{file} has invalid characters!")

		if not set(allele_2).issubset(DNA_bases):
			print(f"{file} has invalid characters!")

		# DPB1 3'-terminal homopolymer-contraction repair. A fully-phased DPB1 CDS can
		# reconstruct full-length but WITHOUT a stop because the allele (e.g. DPB1*13:01,
		# *27:01, *107:01, *135:01) drops 1bp of a poly-A at the CDS 3' terminus vs GRCh38,
		# sliding the stop's last base one position past the annotated CDS boundary. Re-extract
		# with a padded terminal exon and trim to the recovered in-frame stop. Scoped to DPB1;
		# skips rescued/partial records; no-op on records that already end in a stop.
		if (feat == "CDS" and gene == "HLA-DPB1"
				and (not cds_rescued_genes or gene not in cds_rescued_genes)
				and (not unphased_genes or gene not in unphased_genes)):
			vcf_path = gene_filtered_vcfs.get(gene) if gene_filtered_vcfs else None
			needs_repair = (allele_1[-3:] not in stop_codons) or (allele_2[-3:] not in stop_codons)
			if needs_repair and vcf_path and gff_dir and reference_genome:
				fixed_1, fixed_2 = extract_cds_padded_and_trim(
					gene, gene.lower().replace("-", "_"), gff_dir, reference_genome,
					vcf_path, vcf2fasta_output_dir, stop_codons)
				if allele_1[-3:] not in stop_codons and fixed_1:
					allele_1 = fixed_1
					print(f"{sample_ID} {gene} CDS terminal stop recovered via padded re-extraction (hap1)")
				if allele_2[-3:] not in stop_codons and fixed_2:
					allele_2 = fixed_2
					print(f"{sample_ID} {gene} CDS terminal stop recovered via padded re-extraction (hap2)")

		if feat == "CDS":
			if allele_1[0:3] != "ATG" or allele_2[0:3] != "ATG":
				print(f"{sample_ID} {gene} CDS sequence does not begin with start codon!\n")

			if not allele_1[-3:] in stop_codons or not allele_2[-3:] in stop_codons:
				print(f"{sample_ID} {gene} CDS sequence does not end with stop codon!\n")
		
		if feat not in fasta_dict:
			fasta_dict[feat] = {}
		if gene not in fasta_dict[feat]:
			fasta_dict[feat][gene] = []

		fasta_dict[feat][gene].append(allele_1)
		fasta_dict[feat][gene].append(allele_2)

	if config.VERBOSE:
		print("\n")
		print("Sanity check of partially phased genes")
		for string in logging_strings:
			print(string)
		print("\n")

	gene_records = []
	cds_records = []

	for feat, genes in fasta_dict.items():
		for gene, haplotypes in genes.items():
			hap1_name = f"{sample_ID}_{gene}_1"
			hap1_seq = haplotypes[0]
			hap2_name = f"{sample_ID}_{gene}_2"
			hap2_seq = haplotypes[1]
			#if gene in unphased_genes:
				#hap1_name = f"{sample_ID}_{gene}_1_incomplete"
				#hap2_name = f"{sample_ID}_{gene}_2_incomplete"

			if feat == "gene":
				gene_records.append(SeqRecord(Seq(hap1_seq), id=hap1_name, description = ""))
				gene_records.append(SeqRecord(Seq(hap2_seq), id=hap2_name, description = ""))
			elif feat == "CDS":
				cds_records.append(SeqRecord(Seq(hap1_seq), id=hap1_name, description = ""))
				cds_records.append(SeqRecord(Seq(hap2_seq), id=hap2_name, description = ""))

	print("FASTA output:")
	SeqIO.write(gene_records, outfile_gene, "fasta")
	print(f"  {len(gene_records)} records written to: {outfile_gene}")
	SeqIO.write(cds_records, outfile_CDS, "fasta")
	print(f"  {len(cds_records)} records written to: {outfile_CDS}")
	print("\n")
