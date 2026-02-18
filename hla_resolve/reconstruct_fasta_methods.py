import os
import subprocess
import pysam
from Bio.Seq import Seq
from Bio.SeqRecord import SeqRecord
from Bio import SeqIO

def filter_vcf_gene(input_vcf, gene, filter_region, symbolic_vcf, pass_vcf, fail_vcf, pass_unphased, filtered_vcf, platform, genotyper, hla_genes_regions_file):
	# Extract region
	base = os.path.basename(filtered_vcf)
	prefix = base.replace("_PASS_phased.vcf.gz", "")
	region_vcf = os.path.join(os.path.dirname(filtered_vcf), f"{prefix}.vcf.gz")

	cmd = f"bcftools view -r {filter_region} {input_vcf} -Oz -o {region_vcf}"
	subprocess.run(cmd, shell=True, check=True)
	subprocess.run(f"bcftools index -f {region_vcf}", shell=True, check=True)

	vf = pysam.VariantFile(region_vcf)
	sym_out = pysam.VariantFile(symbolic_vcf, "wz", header=vf.header)
	pass_out = pysam.VariantFile(pass_vcf, "wz", header=vf.header)
	fail_out = pysam.VariantFile(fail_vcf, "wz", header=vf.header)

	# ========== PASS/FAIL CLASSIFICATION ==========
	for rec in vf:
		sample = list(rec.samples.values())[0]
		gt = sample.get("GT")

		# HARD EXCLUDE — must not exist downstream
		if gt is None or None in gt:
			continue

		# symbolic
		if (
			("TRID" in rec.info and rec.info["TRID"] not in (None, "", ".")) or
			rec.alts is None or
			any(str(a).startswith("<") for a in rec.alts) or
			any("]" in str(a) or "[" in str(a) for a in rec.alts) or
			any(set(str(a)) - set("ACGTN") for a in rec.alts)
		):
			sym_out.write(rec)
			continue

		# pbsv SV (ID starts with "pbsv.") or Sniffles(Sniffles2) SVs from ONT
		rec_id = rec.id or ""
		rec_id_l = rec_id.lower()
		if rec_id_l.startswith("pbsv.") or rec_id_l.startswith("sniffles"):
			if rec.filter.keys() == ["PASS"] or rec.filter.keys() == []:
				pass_out.write(rec)
			else:
				fail_out.write(rec)
			continue

		# DeepVariant specific filtering
		if genotyper == "deepvariant":
			if rec.filter.keys() == ["PASS"] or rec.filter.keys() == []:
				pass_out.write(rec)
			else:
				fail_out.write(rec)
			continue

		# small variants
		sample = list(rec.samples.values())[0]
		dp = sample.get("DP")
		gq = sample.get("GQ")
		qual = rec.qual or 0
		ref = rec.ref
		alt = rec.alts[0]

		# DP filter
		if dp is None or dp < 2:
			fail_out.write(rec)
			continue

		# SNP
		if len(ref) == 1 and len(alt) == 1:
			if gq not in (None, ".") and gq < 20:
				fail_out.write(rec); continue
			if qual < 10:
				fail_out.write(rec); continue
			pass_out.write(rec); continue

		# INDEL
		if gq not in (None, ".") and gq < 10:
			fail_out.write(rec); continue
		pass_out.write(rec)

	sym_out.close()
	pass_out.close()
	fail_out.close()

	subprocess.run(f"bcftools index -f {symbolic_vcf}", shell=True, check=True)
	subprocess.run(f"bcftools index -f {pass_vcf}", shell=True, check=True)
	subprocess.run(f"bcftools index -f {fail_vcf}", shell=True, check=True)

	# ========== WHITELIST LOGIC ==========
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

	print(f"[DEBUG] {gene}: het={len(het_sites)}, unphased={len(unphased_hets)}")
	allow_single_unphased = (len(het_sites) == 1 and len(unphased_hets) == 1)

	het_clauses = [
		'GT="0/1"', 'GT="1/0"', 'GT="1/2"',
		'GT="2/1"', 'GT="2/3"', 'GT="3/2"'
	]

	if allow_single_unphased:
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
	unph = pysam.VariantFile(pass_unphased)
	records = [rec for rec in unph]

	if records:
		print(f"\nUnphased PASS variants in {gene}:\n")
		for rec in records:
			print(str(rec).strip())
		print()


def filter_vcf_gene_test(input_vcf, gene, filter_region, symbolic_vcf, pass_vcf, fail_vcf, sv_overlap_vcf, pass_unphased, filtered_vcf, platform, genotyper, hla_genes_regions_file):
	"""Version with SV overlap suppression logic - kept for testing"""
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

	# Helper function: check if variant overlaps a PASS SV on the same haplotype
	def overlaps_sv_same_haplotype(pos, ref_len, var_haplotypes):
		var_end = pos + ref_len - 1
		for sv_start, sv_end, sv_haplotypes in sv_regions:
			if pos <= sv_end and var_end >= sv_start:  # position overlap
				if var_haplotypes & sv_haplotypes:      # haplotype overlap (set intersection)
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

		# symbolic
		if (
			("TRID" in rec.info and rec.info["TRID"] not in (None, "", ".")) or
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

		# ========== SV OVERLAP CHECK ==========
		# Skip non-SV variants that overlap a PASS SV on the same haplotype
		# Only apply haplotype-aware suppression for phased indels ≥30bp
		# Never suppress SNPs - they are real variants near SV breakpoints
		ref_len = len(rec.ref)
		alt_len = len(rec.alts[0]) if rec.alts else 0
		indel_size = abs(ref_len - alt_len)
		is_snp = (ref_len == 1 and alt_len == 1)

		if sample.phased and not is_snp and indel_size >= 10:
			var_haplotypes = set()
			for i, allele in enumerate(gt):
				if allele is not None and allele > 0:
					var_haplotypes.add(i)

			if var_haplotypes and overlaps_sv_same_haplotype(rec.pos, len(rec.ref), var_haplotypes):
				# This indel overlaps a PASS SV on the same haplotype - write to sv_overlap file
				sv_overlap_count += 1
				print(f"[SV-OVERLAP]   Suppressed indel ({indel_size}bp): {rec.chrom}:{rec.pos} {rec.ref[:20]}...>{rec.alts[0][:20]}... GT={gt} haps={var_haplotypes}")
				sv_overlap_out.write(rec)
				continue

		# DeepVariant specific filtering
		if genotyper == "deepvariant":
			if rec.filter.keys() == ["PASS"] or rec.filter.keys() == []:
				pass_out.write(rec)
			else:
				fail_out.write(rec)
			continue

		# small variants
		sample = list(rec.samples.values())[0]
		dp = sample.get("DP")
		gq = sample.get("GQ")
		qual = rec.qual or 0
		ref = rec.ref
		alt = rec.alts[0]

		# DP filter
		if dp is None or dp < 2:
			fail_out.write(rec)
			continue

		# SNP
		if len(ref) == 1 and len(alt) == 1:
			if gq not in (None, ".") and gq < 20:
				fail_out.write(rec); continue
			if qual < 10:
				fail_out.write(rec); continue
			pass_out.write(rec); continue

		# INDEL
		if gq not in (None, ".") and gq < 10:
			fail_out.write(rec); continue
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

	print(f"[DEBUG] {gene}: het={len(het_sites)}, unphased={len(unphased_hets)}")
	allow_single_unphased = (len(het_sites) == 1 and len(unphased_hets) == 1)

	het_clauses = [
		'GT="0/1"', 'GT="1/0"', 'GT="1/2"',
		'GT="2/1"', 'GT="2/3"', 'GT="3/2"'
	]

	if allow_single_unphased:
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
	unph = pysam.VariantFile(pass_unphased)
	records = [rec for rec in unph]

	if records:
		print(f"\nUnphased PASS variants in {gene}:\n")
		for rec in records:
			print(str(rec).strip())
		print()

def run_vcf2fasta(vcf2fasta, input_vcf, input_gff, reference_genome, output_dir, gene, feature):
	gene_id = gene.lower().replace("-", "_")
	
	if feature == "CDS":
		vcf2fasta_cmd = f"python3 {vcf2fasta} --fasta {reference_genome} --vcf {input_vcf} --gff {input_gff} -o {output_dir} --feat CDS --blend"
	elif feature == "gene":
		vcf2fasta_cmd = f"python3 {vcf2fasta} --fasta {reference_genome} --vcf {input_vcf} --gff {input_gff} -o {output_dir} --feat gene"
	
	subprocess.run(vcf2fasta_cmd, shell = True, check = True)

def parse_fastas(sample_ID, vcf2fasta_output_dir, outfile_gene, outfile_CDS, DNA_bases, stop_codons, unphased_genes=None, gene_dict=None, CDS_dict=None, gff_dir=None):
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

		if unphased_genes and gene in unphased_genes:
			best_haploblock_start = unphased_genes[gene][0]
			best_haploblock_end   = unphased_genes[gene][1]

			if feat == "gene":
				# Load gene coords (1-based genomic positions)
				gene_lower = gene.lower().replace("-", "_")
				gene_coords_file = os.path.join(gff_dir, f"{gene_lower}_gene_coords.txt")
				gene_coords = [int(item) for item in open(gene_coords_file).read().splitlines()]

				clamped_start = max(best_haploblock_start, gene_dict[gene][0])
				clamped_stop  = min(best_haploblock_end,   gene_dict[gene][1])

				# Map haploblock genomic coords into gene FASTA indices
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
				print(f"File {file} does not begin with start codon!")
		
			if not allele_1[-3:] in stop_codons or not allele_2[-3:] in stop_codons:
				print(f"File {file} does not end with stop codon!")
		
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
