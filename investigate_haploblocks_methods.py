import os
import csv
import json
import pysam

genes_bed = "/hb/groups/cornejo_lab/matt/hla_capture/input_data/reference/parse_haploblocks_bed.bed"
genes_of_interest = ("HLA-A", "HLA-B", "HLA-C", "HLA-DRB1", "HLA-DQA1", "HLA-DQA2", "HLA-DQB1", "HLA-DQB2", "HLA-DPA1", "HLA-DPB1")

# Extended MHC coordinates
mhc_start = 29555628
mhc_stop = 33409896

# Get list of haploblock intervals for MHC
def parse_haploblocks(self):
	sample_name = self.sample_ID
	heterozygous_sites = []
	haploblock_list = []

	if self.platform == "PACBIO":
		vcf_file = os.path.join(self.phased_vcf_dir, self.sample_ID + ".hiphase.joint.vcf.gz")
		haploblock_file = os.path.join(self.phased_vcf_dir, self.sample_ID + ".phased.blocks.txt")
	elif self.platform == "ONT":
		vcf_file = os.path.join(self.phased_vcf_dir, self.sample_ID + ".longphase.vcf.gz")
		haploblock_file = os.path.join(self.phased_vcf_dir, self.sample_ID + ".phased.haploblocks.txt")
	
	vcf = pysam.VariantFile(vcf_file)

	for record in vcf:
		if record.chrom != "chr6":
			continue
		if record.pos < mhc_start or record.pos > mhc_stop:
			continue

		# Safety check in case sample_name is missing in the VCF
		if sample_name not in record.samples:
			raise ValueError(f"Sample '{sample_name}' not found in {vcf_file}")

		genotype = record.samples[sample_name]["GT"]
		if genotype in [(0, 1), (1, 0)]:
			heterozygous_sites.append(record.pos)

	print(f"Sample {sample_name} has {len(heterozygous_sites)} heterozygous extended MHC genotypes")

	print(f"Parsing {sample_name} haploblock file: {haploblock_file}")

	with open(haploblock_file, "r") as f:
		haploblocks = f.read().splitlines()

	for line in haploblocks[1:]:
		fields = line.split("\t")
		
		# HiPhase and Longphase have slightly differently formatted haploblock tsv files. 
		if self.platform == "PACBIO":
			chromosome, start, stop = "chr6", int(fields[4]) - 1, int(fields[5])
		elif self.platform == "ONT":
			chromosome, start, stop = "chr6", int(fields[3]) - 1, int(fields[4])

		if chromosome == "chr6" and stop > mhc_start:
			haploblock_list.append([start,stop])

	return heterozygous_sites, haploblock_list

# Check whether each captured MHC gene is completely spanned by a haploblock
def evaluate_gene_haploblocks(self, het_sites, haploblocks):
	# List of fully phased genes
	haploblocks.sort()
	gene_list = []
	
	# List of genes with partially overlapping haploblock
	incomplete_data = []
	
	genes = open(genes_bed, "r").read().splitlines()
	for line in genes:
		fields = line.split("\t")
		name = fields[3].split("_")[0]
		gene = name
		gene_start = int(fields[1])
		gene_stop = int(fields[2])
		gene_length = gene_stop - gene_start

		gene_het_sites = [site for site in het_sites if site >= gene_start and site <= gene_stop]

		# If the gene has 0 or 1 heterozygous sites, it is effectively fully phased
		if len(gene_het_sites) <= 1:
			gene_list.append(gene)
			continue

		# If the gene is completely spanned by a single haploblock, it is fully phased 
		fully_phased = False
		for block_start, block_stop in haploblocks:
			if block_start <= gene_start and block_stop >= gene_stop:
				gene_list.append(gene)
				fully_phased = True
				break

			# Check to see if unphased genes become fully phased when extending haploblocks through homozygous regions
			# Find first heterozygous site upstream of block start. Do not extend if no heterozygous sites. 
			if block_start <= gene_stop and block_stop >= gene_start:
				extended_start = max([h for h in het_sites if h < block_start], default=block_start)
				extended_start = max(extended_start, gene_start)

				# Find first heterozgyous site downstream of block stop. Do not extend if no heterozygous sites
				extended_stop = min([h for h in het_sites if h > block_stop], default=block_stop)
				extended_stop = min(extended_stop, gene_stop)

				# Consider gene fully phased if haploblock extension through homozygous bases fully spans gene coordinates. 
				if extended_start <= gene_start and extended_stop >= gene_stop:
					gene_list.append(gene)
					fully_phased = True
					break

		# Get details on haploblock overlap for genes of interest (HLA Class I/II) that were not fully phased 
		if not fully_phased and gene in genes_of_interest:
			overlapping_haploblocks = []
			upstream_block = None
			downstream_block = None

			for block_start, block_stop in haploblocks:
				if block_stop < gene_start:
					upstream_block = (block_start, block_stop)
				elif block_start > gene_stop:
					downstream_block = (block_start, block_stop)
					break
				elif block_stop >= gene_start and block_start <= gene_stop:
					overlapping_haploblocks.append((block_start, block_stop))

			num_pre_merge_haploblocks = len(overlapping_haploblocks)  # Track count before extension & merging
			print(f"Processing {self.sample_ID} {gene}")
			print(f"Overlapping unextended haploblocks: {len(overlapping_haploblocks)}")

			if upstream_block:
				overlapping_haploblocks.insert(0, upstream_block)
			if downstream_block:
				overlapping_haploblocks.append(downstream_block)

			# Step 1: Extend each haploblock independently
			extended_haploblocks = []
			for block_start, block_stop in overlapping_haploblocks:
				extended_start = max([h for h in het_sites if h < block_start], default=block_start)
				extended_start = max(extended_start, gene_start)

				extended_stop = min([h for h in het_sites if h > block_stop], default=block_stop)
				extended_stop = min(extended_stop, gene_stop)

				extended_haploblocks.append((extended_start, extended_stop))

			print(f"Extended haploblocks: {len(extended_haploblocks)}")

			# Step 2: Merge overlapping extended haploblocks
			# merged_intervals = []
			# extended_haploblocks.sort()

			# for start, stop in extended_haploblocks:
			# 	if not merged_intervals or start > merged_intervals[-1][1]:
			# 		merged_intervals.append((start, stop))
			# 	else:
			# 		last_start, last_stop = merged_intervals[-1]
			# 		new_stop = max(last_stop, stop)
			# 		merged_intervals[-1] = (last_start, new_stop)

			# print(f"Merged haploblocks: {len(merged_intervals)}")

			# # Step 3: Compute overlap & percentage
			# total_overlap = 0
			# for start, stop in merged_intervals:
			# 	overlap_start = max(start, gene_start)
			# 	overlap_stop = min(stop, gene_stop)
			# 	overlap_length = max(0, overlap_stop - overlap_start)
			# 	total_overlap += overlap_length
			# prop_overlap = total_overlap / gene_length
			# prop_phased_string = f"{prop_overlap * 100:.2f}%"

			# Compute the proportion of the gene spanned by the largest extended haploblock
			max_overlap = 0

			for start, stop in extended_haploblocks:
				overlap_start = max(start, gene_start)
				overlap_stop = min(stop, gene_stop)
				overlap_length = max(0, overlap_stop - overlap_start)
				max_overlap = max(max_overlap, overlap_length)
			prop_overlap = max_overlap / gene_length
			largest_overlap_string = f"{prop_overlap*100:.2f}%"

			# Use the pre-merge haploblock count, not merged count!
			incomplete_data.append([self.sample_ID, gene, num_pre_merge_haploblocks, largest_overlap_string])
			print(f"{self.sample_ID} {gene}, Pre-Merge Haploblocks: {num_pre_merge_haploblocks}")
			print(f"Proportion of gene contained in largest overlapping haploblock: {largest_overlap_string}")

	with open(os.path.join(self.parsed_haploblock_dir, f"phased_genes.tsv"), "w", newline="") as csv_file:
		writer = csv.writer(csv_file, delimiter="\t")
		writer.writerow(["sample", "num_genes", "genes"])
		writer.writerow([self.sample_ID, len(gene_list), ",".join(gene_list)])


	if incomplete_data:
		with open(os.path.join(self.parsed_haploblock_dir, f"incomplete.csv"), "w", newline="") as csvfile:
			csv_writer = csv.writer(csvfile)
			csv_writer.writerow(["sample", "gene", "num_haploblocks", "largest_haploblock"])
			csv_writer.writerows(incomplete_data)

	return gene_list