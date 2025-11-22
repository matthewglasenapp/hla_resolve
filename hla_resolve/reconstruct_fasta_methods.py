import os
import subprocess
import pysam
from Bio.Seq import Seq
from Bio.SeqRecord import SeqRecord
from Bio import SeqIO

def filter_vcf_full(input_vcf, pass_vcf, fail_vcf, pass_unphased, filtered_vcf, unphased_overlap_tsv, platform, genotyper, hla_genes_regions_file):
	# bcftools call does not have a FILTER=PASS annotation, so drop that from the filtering expression
	# Include homozygous ALT (both biallelic and multiallelic) and phased heterozygous
	# Exclude symbolic variants 
	# Note: both pbsv and sniffles have the FILTER=PASS annotation and normal GT fields for INS/DEL. DUP,BND,INV are symbolic and will be excluded but reported
	if platform == "ONT":
		# For ONT data, skip TRID-related filtering entirely
		if genotyper == "bcftools":
			keep_expr = (
				'(GT="1/1" || GT="2/2" || GT="3/3" || GT="4/4" || GT="5/5" || '
				'GT="0|1" || GT="1|0" || GT="1|2" || GT="2|1" || GT="2|3" || GT="3|2") && '
				'(ALT!~"^<") && '
				'(FMT/DP>=10) && '
				'((REF~"^[ACGT]$" && ALT~"^[ACGT]$" && (GQ="." || (GQ>=60 && QUAL>=30))) || '
				'((REF!~"^[ACGT]$" || ALT!~"^[ACGT]$") && (GQ="." || GQ>=40)))'
			)
		else:
			keep_expr = (
				'(FILTER="PASS") && '
				'(ALT!~"^<") && '
				'(GT="1/1" || GT="2/2" || GT="3/3" || GT="4/4" || GT="5/5" || '
				'GT="0|1" || GT="1|0" || GT="1|2" || GT="2|1" || GT="2|3" || GT="3|2")'
			)

		# For unphased_expr, skip TRID check as well
		if genotyper == "bcftools":
			unphased_expr = (
				'(ALT~"^<" || '
				'((GT="0/1" || GT="1/0" || GT="1/2" || GT="2/1" || GT="2/3" || GT="3/2") && '
				'(FMT/DP>=10) && '
				'((REF~"^[ACGT]$" && ALT~"^[ACGT]$" && (GQ="." || (GQ>=60 && QUAL>=30))) || '
				'((REF!~"^[ACGT]$" || ALT!~"^[ACGT]$") && (GQ="." || GQ>=40)))))'
			)
		else:
			unphased_expr = (
				'(FILTER="PASS") && '
				'(GT="0/1" || GT="1/0" || GT="1/2" || GT="2/1" || GT="2/3" || GT="3/2")'
			)

	else:
		# Your existing TRID-aware expressions go here
		if genotyper == "bcftools":
			keep_expr = (
				'TRID=="" && '
				'(GT="1/1" || GT="2/2" || GT="3/3" || GT="4/4" || GT="5/5" || '
				'GT="0|1" || GT="1|0" || GT="1|2" || GT="2|1" || GT="2|3" || GT="3|2") && '
				'(ALT!~"^<") && '
				'((REF~"^[ACGT]$" && ALT~"^[ACGT]$" && (GQ="." || (GQ>=20 && QUAL>=10))) || '
				'((REF!~"^[ACGT]$" || ALT!~"^[ACGT]$") && (GQ="." || GQ>=10)))'
			)
		else:
			keep_expr = (
				'(FILTER="PASS") && '
				'(ALT!~"^<") && '
				'(GT="1/1" || GT="2/2" || GT="3/3" || GT="4/4" || GT="5/5" || '
				'GT="0|1" || GT="1|0" || GT="1|2" || GT="2|1" || GT="2|3" || GT="3|2")'
			)

		if genotyper == "bcftools":
			unphased_expr = (
				'(TRID!="" || ALT~"^<" || '
				'((GT="0/1" || GT="1/0" || GT="1/2" || GT="2/1" || GT="2/3" || GT="3/2") && '
				'((REF~"^[ACGT]$" && ALT~"^[ACGT]$" && (GQ="." || (GQ>=20 && QUAL>=10))) || '
				'((REF!~"^[ACGT]$" || ALT!~"^[ACGT]$") && (GQ="." || GQ>=10)))))'
			)
		else:
			unphased_expr = (
				'((FILTER="PASS") && (GT="0/1" || GT="1/0" || '
				'GT="1/2" || GT="2/1" || GT="2/3" || GT="3/2")) || TRID!=""'
			)

	# Create new pass-filter VCF
	subprocess.run(
		f'bcftools view -i \'{keep_expr}\' {input_vcf} -Oz -o {pass_vcf}',
		shell=True, check=True
	)
	subprocess.run(f"bcftools index -f {pass_vcf}", shell=True, check=True)

	# Create new fail-filter vcf
	subprocess.run(
		f'bcftools view -e \'{keep_expr}\' {input_vcf} -Oz -o {fail_vcf}',
		shell=True, check=True
	)
	subprocess.run(f"bcftools index -f {fail_vcf}", shell=True, check=True)

	subprocess.run(
		f'bcftools view -i \'{unphased_expr}\' {fail_vcf} -Oz -o {pass_unphased}',
		shell=True, check=True
	)
	subprocess.run(f"bcftools index -f {pass_unphased}", shell=True, check=True)

	# Intersect unphased PASS heterozygous genotypes from the fail-filter vcf with HLA BED to get overlapping variants that could not bee applied
	subprocess.run(
		f'bedtools intersect -a {pass_unphased} -b {hla_genes_regions_file} -wa -wb -header > {unphased_overlap_tsv}',
		shell=True, check=True
	)
	
	# Legacy code from before switching to single-gene VCFs as input to vcf2fasta. 
	# Report overlapping variants that could not be applied 
	# with open(unphased_overlap_tsv, "r") as f:
	# 	overlap_lines = f.read().splitlines()

	# counter = 0
	# header = ""
	# overlap_dict = dict()
	# for line in overlap_lines:
	# 	if line.startswith("#CHROM"):
	# 		header = line
	# 	if not line.startswith("chr6"):
	# 		continue
	# 	counter += 1
	# 	gene = line.split("\t")[-1].split("_")[0]
	# 	record = line.split("\t")[:-4]
	# 	new_record = "\t".join(record)
	# 	if not gene in overlap_dict:
	# 		overlap_dict[gene] = [new_record]
	# 	else:
	# 		overlap_dict[gene].append(new_record)

	# print(f"The following {counter} unphased variants overlapped with an HLA gene and could not be applied")
	# print("\n")
	# for gene, records in overlap_dict.items():
	# 	print(gene)
	# 	print(header)
	# 	for record in records:
	# 		print(record)
	# 	print("\n")

	# Find first phased variant
	# vcf = pysam.VariantFile(pass_vcf)
	# first_phased_pos = None
	# chrom = None

	# for record in vcf:
	# 	for sample_data in record.samples.values():
	# 		if sample_data.phased:
	# 			first_phased_pos = record.pos
	# 			chrom = record.chrom
	# 			break
	# 	if first_phased_pos:
	# 		break

	# if not first_phased_pos:
	# 	raise ValueError(f"No phased variants found in {pass_vcf}")

	# region = f"{chrom}:{first_phased_pos}-"

	# # Filter VCF from the first phased variant onward
	# subprocess.run(
	# 	f"bcftools view -r {region} {pass_vcf} -Oz -o {filtered_vcf}",
	# 	shell=True, check=True
	# )
	# subprocess.run(f"bcftools index {filtered_vcf}", shell=True, check=True)

def filter_vcf_gene(input_vcf, gene, filter_region, pass_vcf, fail_vcf, pass_unphased, filtered_vcf, unphased_overlap_tsv, platform, genotyper, hla_genes_regions_file):
	region_vcf = filtered_vcf.replace(".vcf.gz", f".{gene}.region.vcf.gz")
	subprocess.run(
		f"bcftools view -r {filter_region} {input_vcf} -Oz -o {region_vcf}",
		shell=True, check=True
	)
	subprocess.run(f"bcftools index -f {region_vcf}", shell=True, check=True)

	# Address edge-case of one unphased heteterozygous genotype in the gene region
	vcf_region = pysam.VariantFile(region_vcf)

	het_sites = []
	unphased_hets = []
	valid_dna = set("ACGTN")

	for rec in vcf_region:
		if "TRID" in rec.info:
			continue
		
		if any(set(str(a)) - valid_dna for a in rec.alts):
			continue
		
		sample = list(rec.samples.values())[0]
		gt = sample.get('GT')

		if gt is None or None in gt:
			continue
	
		# Check if genotype is heterozygous
		if len(set(gt)) == 2:
			het_sites.append(rec)

			# Check if heterozygous genotype is unphased
			if not sample.phased:
				unphased_hets.append(rec)

	# If exactly ONE heterozygous site AND it is unphased, allow it to pass
	allow_single_unphased = (len(het_sites) == 1 and len(unphased_hets) == 1)

	if platform == "ONT":
		# For ONT data, skip TRID-related filtering entirely
		if genotyper == "bcftools":
			keep_expr = (
				'(GT="1/1" || GT="2/2" || GT="3/3" || GT="4/4" || GT="5/5" || '
				'GT="0|1" || GT="1|0" || GT="1|2" || GT="2|1" || GT="2|3" || GT="3|2") && '
				'(ALT!~"^<") && '
				'(FMT/DP>=10) && '
				'((REF~"^[ACGT]$" && ALT~"^[ACGT]$" && (GQ="." || (GQ>=60 && QUAL>=30))) || '
				'((REF!~"^[ACGT]$" || ALT!~"^[ACGT]$") && (GQ="." || GQ>=40)))'
			)
		else:
			keep_expr = (
				'(FILTER="PASS") && '
				'(ALT!~"^<") && '
				'(GT="1/1" || GT="2/2" || GT="3/3" || GT="4/4" || GT="5/5" || '
				'GT="0|1" || GT="1|0" || GT="1|2" || GT="2|1" || GT="2|3" || GT="3|2")'
			)

		# For unphased_expr, skip TRID check as well
		if genotyper == "bcftools":
			unphased_expr = (
				'(ALT~"^<" || '
				'((GT="0/1" || GT="1/0" || GT="1/2" || GT="2/1" || GT="2/3" || GT="3/2") && '
				'(FMT/DP>=10) && '
				'((REF~"^[ACGT]$" && ALT~"^[ACGT]$" && (GQ="." || (GQ>=60 && QUAL>=30))) || '
				'((REF!~"^[ACGT]$" || ALT!~"^[ACGT]$") && (GQ="." || GQ>=40)))))'
			)
		else:
			unphased_expr = (
				'(FILTER="PASS") && '
				'(GT="0/1" || GT="1/0" || GT="1/2" || GT="2/1" || GT="2/3" || GT="3/2")'
			)

	else:
		# PacBio TRID-aware expressions
		if genotyper == "bcftools":
			keep_expr = (
				'TRID=="" && '
				'(GT="1/1" || GT="2/2" || GT="3/3" || GT="4/4" || GT="5/5" || '
				'GT="0|1" || GT="1|0" || GT="1|2" || GT="2|1" || GT="2|3" || GT="3|2") && '
				'(ALT!~"^<") && '
				'((REF~"^[ACGT]$" && ALT~"^[ACGT]$" && (GQ="." || (GQ>=20 && QUAL>=10))) || '
				'((REF!~"^[ACGT]$" || ALT!~"^[ACGT]$") && (GQ="." || GQ>=10)))'
			)
		else:
			keep_expr = (
				'(FILTER="PASS") && '
				'(ALT!~"^<") && '
				'(GT="1/1" || GT="2/2" || GT="3/3" || GT="4/4" || GT="5/5" || '
				'GT="0|1" || GT="1|0" || GT="1|2" || GT="2|1" || GT="2|3" || GT="3|2")'
			)

		if genotyper == "bcftools":
			unphased_expr = (
				'(TRID!="" || ALT~"^<" || '
				'((GT="0/1" || GT="1/0" || GT="1/2" || GT="2/1" || GT="2/3" || GT="3/2") && '
				'((REF~"^[ACGT]$" && ALT~"^[ACGT]$" && (GQ="." || (GQ>=20 && QUAL>=10))) || '
				'((REF!~"^[ACGT]$" || ALT!~"^[ACGT]$") && (GQ="." || GQ>=10)))))'
			)
		else:
			unphased_expr = (
				'((FILTER="PASS") && (GT="0/1" || GT="1/0" || '
				'GT="1/2" || GT="2/1" || GT="2/3" || GT="3/2")) || TRID!=""'
			)

	# Filter expression override if exactly one unphased heterozygous genotype exists
	if allow_single_unphased:
		rec = unphased_hets[0]
		chrom = rec.chrom
		pos = rec.pos

		print(f"Allowing single unphased heterozygous site at {gene} {chrom}:{pos}" + "\n")

		# Build whitelist condition for this one variant only
		whitelist_expr = f'(CHROM="{chrom}" && POS={pos})'

		# Modify keep_expr to include the whitelist_expr
		keep_expr = f'(({keep_expr}) || {whitelist_expr})'

		# Make unphased_expr match nothing
		unphased_expr = 'GT="0/0"'
	
	# Create new pass-filter VCF
	subprocess.run(
		f'bcftools view -i \'{keep_expr}\' {region_vcf} -Oz -o {pass_vcf}',
		shell=True, check=True
	)
	subprocess.run(f"bcftools index -f {pass_vcf}", shell=True, check=True)

	# Create new fail-filter vcf
	subprocess.run(
		f'bcftools view -e \'{keep_expr}\' {region_vcf} -Oz -o {fail_vcf}',
		shell=True, check=True
	)
	subprocess.run(f"bcftools index -f {fail_vcf}", shell=True, check=True)

	# Created new pass-filter VCF for unphased heterozygous genotypes
	subprocess.run(
		f'bcftools view -i \'{unphased_expr}\' {fail_vcf} -Oz -o {pass_unphased}',
		shell=True, check=True
	)
	subprocess.run(f"bcftools index -f {pass_unphased}", shell=True, check=True)

	# Intersect unphased PASS heterozygous genotypes with gene coordinates to get overlapping variants that could not be applied
	subprocess.run(
		f'bedtools intersect -a {pass_unphased} -b {hla_genes_regions_file} -wa -wb -header > {unphased_overlap_tsv}',
		shell=True, check=True
	)
	
	# Report overlapping variants that could not be applied 
	with open(unphased_overlap_tsv, "r") as f:
		overlap_lines = f.read().splitlines()

	header = None
	records_for_this_gene = []

	for line in overlap_lines:
		if line.startswith("#CHROM"):
			header = line
			continue

		if not line.startswith("chr6"):
			continue

		gene_name = line.split("\t")[-1].split("_")[0]

		# Only keep variants for the gene currently being filtered
		if gene_name != gene:
			continue

		record = line.split("\t")[:-4]
		new_record = "\t".join(record)
		records_for_this_gene.append(new_record)

	if len(records_for_this_gene) > 0:
		print(f"The following {len(records_for_this_gene)} unphased variants overlapped with {gene} and could not be applied\n")
		print(gene)
		print(header)
		for record in records_for_this_gene:
			print(record)
		print("\n")

def run_vcf2fasta(vcf2fasta, input_vcf, input_gff, reference_genome, output_dir, gene, feature):
	gene_id = gene.lower().replace("-", "_")
	
	if feature == "CDS":
		vcf2fasta_cmd = f"python3 {vcf2fasta} --fasta {reference_genome} --vcf {input_vcf} --gff {input_gff} -o {output_dir} --feat CDS --blend"
	elif feature == "gene":
		vcf2fasta_cmd = f"python3 {vcf2fasta} --fasta {reference_genome} --vcf {input_vcf} --gff {input_gff} -o {output_dir} --feat gene"
	
	subprocess.run(vcf2fasta_cmd, shell = True, check = True)

def parse_fastas(sample_ID, vcf2fasta_output_dir, outfile_gene, outfile_CDS, DNA_bases, stop_codons, unphased_genes=None, gene_dict=None, CDS_dict=None, gff_dir=None):
	print("\n")
	print(f"DEBUG: Current working directory: {os.getcwd()}")
	print(f"DEBUG: Looking for FASTA files in: {vcf2fasta_output_dir}")
	
	# Use subprocess.run with capture_output to avoid race conditions with temporary files
	find_cmd = f"find {vcf2fasta_output_dir} -type f"
	print(f"DEBUG: Running command: {find_cmd}")
	result = subprocess.run(find_cmd, shell=True, check=True, capture_output=True, text=True)
	fasta_files = result.stdout.strip().split('\n') if result.stdout.strip() else []
	print(f"DEBUG: Found {len(fasta_files)} FASTA files")
	
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
