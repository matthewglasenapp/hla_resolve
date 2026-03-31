# This software is Copyright ©2026. The Regents of the University of California
# ("Regents"). All Rights Reserved.
#
# See LICENSE.txt for license details.

import os
import subprocess
import pysam
from Bio.Seq import Seq
from Bio.SeqRecord import SeqRecord
from Bio import SeqIO

def filter_vcf_gene(input_vcf, gene, filter_region, symbolic_vcf, pass_vcf, fail_vcf, sv_overlap_vcf, pass_unphased, filtered_vcf, platform, genotyper, force_include_unphased=False):
	# Extract region
	base = os.path.basename(filtered_vcf)
	prefix = base.replace("_PASS_phased.vcf.gz", "")
	region_vcf = os.path.join(os.path.dirname(filtered_vcf), f"{prefix}.vcf.gz")

	cmd = f"bcftools view -r {filter_region} {input_vcf} -Oz -o {region_vcf}"
	subprocess.run(cmd, shell=True, check=True)
	subprocess.run(f"bcftools index -f {region_vcf}", shell=True, check=True)

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

	# Debug: Print collected SV regions
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

	print(f"[TR-OVERLAP] {gene}: Collected {len(tr_regions)} TRGT regions for overlap suppression")
	for tr_start, tr_end in tr_regions:
		print(f"[TR-OVERLAP]   TR region: {tr_start}-{tr_end}")

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

		# ========== TR OVERLAP CHECK ==========
		# Suppress non-TRGT variants that fall entirely inside a TRGT region.
		# The TRGT record already encodes the correct allele for the whole repeat;
		# bcftools-called SNVs/indels inside that span are redundant and would
		# cause double-counting in vcf2fasta.
		is_trgt = "TRID" in rec.info and rec.info["TRID"] not in (None, "", ".")
		if not is_trgt and tr_regions and overlaps_tr_region(rec.pos, len(rec.ref)):
			print(f"[TR-OVERLAP]   Suppressed: {rec.chrom}:{rec.pos} {rec.ref[:20]}->{rec.alts[0][:20] if rec.alts else '.'} (inside TRGT span)")
			sv_overlap_out.write(rec)
			sv_overlap_count += 1
			continue

		# TRGT records with explicit alleles bypass quality filters (no DP/GQ)
		if is_trgt:
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
		ref_len = len(rec.ref)
		alt_len = len(rec.alts[0]) if rec.alts else 0
		indel_size = abs(ref_len - alt_len)
		is_snp = (ref_len == 1 and alt_len == 1)

		is_homozygous = gt is not None and len(set(a for a in gt if a is not None)) == 1
		if (sample.phased or is_homozygous) and not is_snp:
			var_haplotypes = set()
			for i, allele in enumerate(gt):
				if allele is not None and allele > 0:
					var_haplotypes.add(i)

			if var_haplotypes and overlaps_sv_same_haplotype(rec.pos, len(rec.ref), var_haplotypes, indel_size):
				# This indel overlaps a PASS SV on the same haplotype - write to sv_overlap file
				sv_overlap_count += 1
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

	print(f"[SV-OVERLAP] {gene}: Total small variants suppressed due to SV overlap: {sv_overlap_count}")

	subprocess.run(f"bcftools index -f {symbolic_vcf}", shell=True, check=True)
	subprocess.run(f"bcftools index -f {pass_vcf}", shell=True, check=True)
	subprocess.run(f"bcftools index -f {fail_vcf}", shell=True, check=True)
	subprocess.run(f"bcftools index -f {sv_overlap_vcf}", shell=True, check=True)

	# ========== WHITELIST LOGIC (RESTORED) ==========
	het_sites = []
	unphased_hets = []

	pass_vf = pysam.VariantFile(pass_vcf)
	for rec in pass_vf:
		sample = list(rec.samples.values())[0]
		gt = sample.get("GT")
		if gt is None or None in gt:
			continue
		if len(set(gt)) == 2:  # heterozygous
			het_sites.append(rec)
			if not sample.phased:
				unphased_hets.append(rec)

	print(f"{gene}: het={len(het_sites)}, unphased={len(unphased_hets)}")
	allow_single_unphased = (len(het_sites) == 1 and len(unphased_hets) == 1)

	het_clauses = [
		'GT="0/1"', 'GT="1/0"', 'GT="1/2"',
		'GT="2/1"', 'GT="2/3"', 'GT="3/2"'
	]

	if force_include_unphased:
		# CDS rescue: include ALL variants (phased + unphased) in filtered VCF
		print(f"[CDS-RESCUE] {gene}: force_include_unphased=True, keeping all variants")
		unphased_expr = 'GT="9/9"'     # matches nothing → nothing sent to pass_unphased
		keep_expr = 'GT!="9/9"'        # matches everything → all variants kept

	elif allow_single_unphased:
		# one heterozygous site, unphased → treat as fully phased
		chosen = unphased_hets[0]
		chrom = chosen.chrom
		pos   = chosen.pos

		# NEGATED form for "keep all non-hets"
		negated = " && ".join([f'{c.replace("=", "!=")}' for c in het_clauses])

		# whitelist the one unphased site so it remains in phased VCF
		whitelist = f'(CHROM="{chrom}" && POS={pos})'

		keep_expr = f'({negated}) || {whitelist}'

		# IMPORTANT: prevent *anything* from being written to pass_unphased
		unphased_expr = 'GT="9/9"'     # matches nothing

	else:
		# Normal case: send all heterozygous unphased variants to pass_unphased
		unphased_expr = " || ".join(het_clauses)

		# phased variants = everything NOT matching the het genotypes
		keep_expr = " && ".join([f'{c.replace("=", "!=")}' for c in het_clauses])

	# extract unphased PASS
	cmd = f"bcftools view -i '{unphased_expr}' {pass_vcf} -Oz -o {pass_unphased}"
	subprocess.run(cmd, shell=True, check=True)
	subprocess.run(f"bcftools index -f {pass_unphased}", shell=True, check=True)

	# extract phased PASS (final filtered VCF)
	cmd = f"bcftools view -i '{keep_expr}' {pass_vcf} -Oz -o {filtered_vcf}"
	subprocess.run(cmd, shell=True, check=True)
	subprocess.run(f"bcftools index -f {filtered_vcf}", shell=True, check=True)

	# ========== PRINT UNPHASED RECORDS NEATLY ==========
	if force_include_unphased and unphased_hets:
		print(f"\nUnphased PASS variants in {gene}:\n")
		for rec in unphased_hets:
			print(str(rec).strip())
		print()
	else:
		unph = pysam.VariantFile(pass_unphased)
		records = [rec for rec in unph]

		if records:
			print(f"\nUnphased PASS variants in {gene}:\n")
			for rec in records:
				print(str(rec).strip())
			print()

def compute_indel_offset(vcf_path, haplotype, target_pos):
	"""
	Compute cumulative indel offset for a haplotype from all variants
	with pos < target_pos. Returns the net bp shift (positive = net insertion,
	negative = net deletion) in forward-strand coordinate space.
	"""
	offset = 0
	vcf = pysam.VariantFile(vcf_path)
	for rec in vcf:
		if rec.pos >= target_pos:
			break
		sample = list(rec.samples.values())[0]
		gt = sample.get('GT')
		if gt is None or None in gt:
			continue
		allele_idx = gt[haplotype]
		if allele_idx == 0:
			continue
		ref = rec.ref
		alt = rec.alleles[allele_idx]
		offset += len(alt) - len(ref)
	vcf.close()
	return offset


def compute_cds_indel_offset(vcf_path, haplotype, target_pos, cds_exon_ranges):
	"""
	Compute cumulative indel offset for a haplotype from variants within CDS
	exons that have pos < target_pos. Unlike compute_indel_offset, this skips
	intronic variants since vcf2fasta --blend only applies variants per-exon.

	Args:
		vcf_path: Path to the filtered VCF
		haplotype: 0 or 1
		target_pos: Genomic position cutoff
		cds_exon_ranges: List of (start, stop) tuples for CDS exons (1-based, inclusive)
	"""
	offset = 0
	vcf = pysam.VariantFile(vcf_path)
	for rec in vcf:
		if rec.pos >= target_pos:
			break
		# Check if variant falls within any CDS exon
		in_cds = any(start <= rec.pos <= stop for start, stop in cds_exon_ranges)
		if not in_cds:
			continue
		sample = list(rec.samples.values())[0]
		gt = sample.get('GT')
		if gt is None or None in gt:
			continue
		allele_idx = gt[haplotype]
		if allele_idx == 0:
			continue
		ref = rec.ref
		alt = rec.alleles[allele_idx]
		offset += len(alt) - len(ref)
	vcf.close()
	return offset


def clamp_cds_fasta_sequence(allele_seq, haplotype, clamped_start, clamped_stop, cds_coords, cds_exon_ranges, is_minus_strand, vcf_path):
	"""
	Clamp a vcf2fasta CDS (--blend) haplotype sequence to a genomic interval,
	accounting for indel-induced coordinate shifts within CDS exons.

	vcf2fasta --blend concatenates per-exon sequences (each with variants applied),
	then reverse complements for minus-strand genes. The cds_coords file provides
	the reference-based mapping, but indels within exons shift the actual FASTA
	indices.

	Args:
		allele_seq: Full CDS haplotype sequence from vcf2fasta --blend
		haplotype: 0 or 1
		clamped_start: Low genomic coordinate of desired CDS interval
		clamped_stop: High genomic coordinate of desired CDS interval
		cds_coords: List of all CDS base positions (flattened across exons)
		cds_exon_ranges: List of (start, stop) tuples for CDS exons
		is_minus_strand: Whether the gene is on the minus strand
		vcf_path: Path to the filtered VCF used by vcf2fasta
	"""
	total_len = len(allele_seq)

	# Reference-based CDS index (position in concatenated exon list)
	ref_idx_start = cds_coords.index(clamped_start)
	ref_idx_stop = cds_coords.index(clamped_stop)

	# Compute indel offsets only from variants within CDS exons
	offset_at_start = compute_cds_indel_offset(vcf_path, haplotype, clamped_start, cds_exon_ranges)
	offset_at_stop = compute_cds_indel_offset(vcf_path, haplotype, clamped_stop + 1, cds_exon_ranges)

	# Forward-strand CDS indices (before reverse complement)
	fwd_low = ref_idx_start + offset_at_start
	fwd_high = ref_idx_stop + offset_at_stop

	if is_minus_strand:
		fasta_start = total_len - 1 - fwd_high
		fasta_stop = total_len - 1 - fwd_low
	else:
		fasta_start = fwd_low
		fasta_stop = fwd_high

	return allele_seq[fasta_start:fasta_stop + 1]


def clamp_fasta_sequence(allele_seq, haplotype, clamped_start, clamped_stop, gene_start, is_minus_strand, vcf_path):
	"""
	Clamp a vcf2fasta haplotype sequence to a genomic interval, accounting for
	indel-induced coordinate shifts and strand orientation.

	vcf2fasta applies variants to the forward-strand reference, then reverse
	complements for minus-strand genes. Reference-based indices don't account
	for insertions/deletions that shift the sequence, so we compute per-haplotype
	offsets from the VCF to find the correct FASTA indices.

	Args:
		allele_seq: Full haplotype sequence from vcf2fasta
		haplotype: 0 or 1 (which haplotype to compute offsets for)
		clamped_start: Low genomic coordinate of desired interval
		clamped_stop: High genomic coordinate of desired interval
		gene_start: Lowest genomic coordinate of the gene (GFF start)
		is_minus_strand: Whether the gene is on the minus strand
		vcf_path: Path to the filtered VCF used by vcf2fasta
	"""
	total_len = len(allele_seq)

	offset_at_start = compute_indel_offset(vcf_path, haplotype, clamped_start)
	offset_at_stop = compute_indel_offset(vcf_path, haplotype, clamped_stop + 1)

	# Forward-strand intermediate indices (before reverse complement)
	fwd_low = (clamped_start - gene_start) + offset_at_start
	fwd_high = (clamped_stop - gene_start) + offset_at_stop

	if is_minus_strand:
		fasta_start = total_len - 1 - fwd_high
		fasta_stop = total_len - 1 - fwd_low
	else:
		fasta_start = fwd_low
		fasta_stop = fwd_high

	return allele_seq[fasta_start:fasta_stop + 1]


def run_vcf2fasta(input_vcf, input_gff, reference_genome, output_dir, gene, feature):
	if feature == "CDS":
		vcf2fasta_cmd = f"vcf2fasta --fasta {reference_genome} --vcf {input_vcf} --gff {input_gff} -o {output_dir} --feat CDS --blend --force"
	elif feature == "gene":
		vcf2fasta_cmd = f"vcf2fasta --fasta {reference_genome} --vcf {input_vcf} --gff {input_gff} -o {output_dir} --feat gene --force"

	subprocess.run(vcf2fasta_cmd, shell = True, check = True)

def parse_fastas(sample_ID, vcf2fasta_output_dir, outfile_gene, outfile_CDS, DNA_bases, stop_codons, unphased_genes=None, gene_dict=None, CDS_dict=None, gff_dir=None, cds_rescued_genes=None, ARS_dict=None, CLASS_I_GENES=None, gene_filtered_vcfs=None):
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
						is_minus = gene_coords[0] > gene_coords[-1]
						gff_gene_start = min(gene_start_coord, gene_stop_coord)
						vcf_path = gene_filtered_vcfs.get(gene) if gene_filtered_vcfs else None
						if vcf_path:
							allele_1 = clamp_fasta_sequence(allele_1, 0, clamped_start, clamped_stop, gff_gene_start, is_minus, vcf_path)
							allele_2 = clamp_fasta_sequence(allele_2, 1, clamped_start, clamped_stop, gff_gene_start, is_minus, vcf_path)
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
					# Extract only ARS CDS positions from CDS sequence
					gene_lower = gene.lower().replace("-", "_")
					cds_coords_file = os.path.join(gff_dir, f"{gene_lower}_cds_sorted_coords.txt")
					cds_coords = [int(item) for item in open(cds_coords_file).read().splitlines()]

					ars_cds_ranges = rescue_info["ars_cds_ranges"]
					ars_cds_positions = [pos for pos in cds_coords
										if any(s <= pos <= e for s, e in ars_cds_ranges)]

					if ars_cds_positions:
						idx1 = cds_coords.index(ars_cds_positions[0])
						idx2 = cds_coords.index(ars_cds_positions[-1])
						cds_fasta_start, cds_fasta_stop = sorted((idx1, idx2))
						allele_1 = allele_1[cds_fasta_start:cds_fasta_stop+1]
						allele_2 = allele_2[cds_fasta_start:cds_fasta_stop+1]
						logging_strings.append(f"{sample_ID} {gene} CDS rescue tier=ars_only: extracting ARS CDS ({len(ars_cds_positions)} positions)")
					else:
						allele_1, allele_2 = "", ""
						logging_strings.append(f"{sample_ID} {gene} CDS rescue tier=ars_only: no ARS CDS positions found")

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
						is_minus = gene_coords[0] > gene_coords[-1]
						gff_gene_start = min(gene_start_coord, gene_stop_coord)
						vcf_path = gene_filtered_vcfs.get(gene) if gene_filtered_vcfs else None
						if vcf_path:
							allele_1 = clamp_fasta_sequence(allele_1, 0, clamped_start, clamped_stop, gff_gene_start, is_minus, vcf_path)
							allele_2 = clamp_fasta_sequence(allele_2, 1, clamped_start, clamped_stop, gff_gene_start, is_minus, vcf_path)
						else:
							idx1 = gene_coords.index(clamped_start)
							idx2 = gene_coords.index(clamped_stop)
							fasta_start, fasta_stop = sorted((idx1, idx2))
							allele_1 = allele_1[fasta_start:fasta_stop+1]
							allele_2 = allele_2[fasta_start:fasta_stop+1]
						logging_strings.append(f"{sample_ID} {gene} CDS rescue tier=ars_only: gene interval chr6:{clamped_start}-{clamped_stop}")

		elif unphased_genes and gene in unphased_genes:
			best_haploblock_start = unphased_genes[gene][0]
			best_haploblock_end   = unphased_genes[gene][1]

			if feat == "gene":
				# Load gene coords (1-based genomic positions)
				gene_lower = gene.lower().replace("-", "_")
				gene_coords_file = os.path.join(gff_dir, f"{gene_lower}_gene_coords.txt")
				gene_coords = [int(item) for item in open(gene_coords_file).read().splitlines()]

				clamped_start = max(best_haploblock_start, gene_dict[gene][0])
				clamped_stop  = min(best_haploblock_end,   gene_dict[gene][1])

				is_minus = gene_coords[0] > gene_coords[-1]
				gff_gene_start = min(gene_dict[gene][0], gene_dict[gene][1])
				vcf_path = gene_filtered_vcfs.get(gene) if gene_filtered_vcfs else None
				if vcf_path:
					allele_1 = clamp_fasta_sequence(allele_1, 0, clamped_start, clamped_stop, gff_gene_start, is_minus, vcf_path)
					allele_2 = clamp_fasta_sequence(allele_2, 1, clamped_start, clamped_stop, gff_gene_start, is_minus, vcf_path)
				else:
					idx1 = gene_coords.index(clamped_start)
					idx2 = gene_coords.index(clamped_stop)
					fasta_start, fasta_stop = sorted((idx1, idx2))
					allele_1 = allele_1[fasta_start:fasta_stop+1]
					allele_2 = allele_2[fasta_start:fasta_stop+1]

			elif feat == "CDS":
				# Load CDS coords (flattened genomic positions from all CDS exons)
				gene_lower = gene.lower().replace("-", "_")
				cds_coords_file = os.path.join(gff_dir, f"{gene_lower}_cds_sorted_coords.txt")
				cds_coords = [int(item) for item in open(cds_coords_file).read().splitlines()]

				# Collect overlap of haploblock with CDS
				cds_overlap = [pos for pos in cds_coords if best_haploblock_start <= pos <= best_haploblock_end]

				if cds_overlap:
					idx1 = cds_coords.index(cds_overlap[0])
					idx2 = cds_coords.index(cds_overlap[-1])
					cds_fasta_start, cds_fasta_stop = sorted((idx1, idx2))
					allele_1 = allele_1[cds_fasta_start:cds_fasta_stop+1]
					allele_2 = allele_2[cds_fasta_start:cds_fasta_stop+1]
				else:
					# no overlap between haploblock and CDS, wipe to empty
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
			print(f"File {file} has no sequence!")
			continue
			
		if not set(allele_1).issubset(DNA_bases):
			print(f"{file} has invalid characters!")

		if not set(allele_2).issubset(DNA_bases):
			print(f"{file} has invalid characters!")
		
		if feat == "CDS":
			if allele_1[0:3] != "ATG" or allele_2[0:3] != "ATG":
				print(f"{sample_ID} {gene} CDS sequence does not begin with start codon!")

			if not allele_1[-3:] in stop_codons or not allele_2[-3:] in stop_codons:
				print(f"{sample_ID} {gene} CDS sequence does not end with stop codon!")
		
		if feat not in fasta_dict:
			fasta_dict[feat] = {}
		if gene not in fasta_dict[feat]:
			fasta_dict[feat][gene] = []

		fasta_dict[feat][gene].append(allele_1)
		fasta_dict[feat][gene].append(allele_2)

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

	SeqIO.write(gene_records, outfile_gene, "fasta")
	print(f"Wrote {len(gene_records)} records to {outfile_gene}")
	SeqIO.write(cds_records, outfile_CDS, "fasta")
	print(f"Wrote {len(cds_records)} records to {outfile_CDS}")
	print("\n")
