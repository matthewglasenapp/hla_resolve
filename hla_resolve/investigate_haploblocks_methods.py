import os
import csv
import json
import pysam
import pandas as pd

# Get list of haploblock intervals for MHC
def parse_haploblocks(input_vcf, input_haploblock_file, platform,sample_ID, mhc_start, mhc_stop):
	sample_name = sample_ID
	heterozygous_sites = []
	haploblock_list = []
	
	vcf = pysam.VariantFile(input_vcf)

	for record in vcf:
		if record.chrom != "chr6":
			continue
		if record.pos < mhc_start or record.pos > mhc_stop:
			continue

		# Safety check in case sample_name is missing in the VCF
		if sample_name not in record.samples:
			raise ValueError(f"Sample '{sample_name}' not found in {input_vcf}")

		genotype = record.samples[sample_name]["GT"]
		if genotype in [(0, 1), (1, 0)]:
			heterozygous_sites.append(record.pos)

	print(f"Sample {sample_name} has {len(heterozygous_sites)} heterozygous extended MHC genotypes")

	print(f"Parsing {sample_name} haploblock file: {input_haploblock_file}")
	print("\n")

	with open(input_haploblock_file, "r") as f:
		haploblocks = f.read().splitlines()

	for line in haploblocks[1:]:
		fields = line.split("\t")
		
		# HiPhase and Longphase have slightly differently formatted haploblock tsv files. 
		if platform == "PACBIO":
			chromosome, start, stop = "chr6", int(fields[4]) - 1, int(fields[5])
		elif platform == "ONT":
			chromosome, start, stop = "chr6", int(fields[3]) - 1, int(fields[4])

		if chromosome == "chr6" and stop > mhc_start:
			haploblock_list.append([start,stop])

	return heterozygous_sites, haploblock_list

# Check whether each captured MHC gene is completely spanned by a haploblock
def evaluate_gene_haploblocks(output_file, incomplete_file, sample_ID, genes_bed, genes_of_interest, het_sites, haploblocks, ARS_dict=None):
	# List of fully phased genes
	haploblocks.sort()
	gene_list = []
	unphased_genes = dict()
	
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
			print(f"{sample_ID} {gene} has {len(gene_het_sites)} heterozygous sites")
			print(f"Treating as fully phased")
			continue

		# If the gene is completely spanned by a single haploblock, it is fully phased 
		fully_phased = False
		for block_start, block_stop in haploblocks:
			if block_start <= gene_start and block_stop >= gene_stop:
				gene_list.append(gene)
				fully_phased = True
				break

			# Check to see if unphased genes become fully phased when extending haploblocks through homozygous regions
			if block_start <= gene_stop and block_stop >= gene_start:
				# Find closest het site upstream
				upstream_het_sites = [h for h in het_sites if h < block_start]
				if upstream_het_sites:
					extended_start = max(upstream_het_sites) + 1
				else:
					extended_start = block_start
				# removed clamping to gene_start

				# Find closest het site downstream
				downstream_het_sites = [h for h in het_sites if h > block_stop]
				if downstream_het_sites:
					extended_stop = min(downstream_het_sites) - 1
				else:
					extended_stop = block_stop
				# removed clamping to gene_stop

				if extended_start <= gene_start and extended_stop >= gene_stop:
					gene_list.append(gene)
					fully_phased = True
					break

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

			num_pre_merge_haploblocks = len(overlapping_haploblocks)
			print(f"Processing {sample_ID} {gene}")
			print(f"Overlapping unextended haploblocks: {len(overlapping_haploblocks)}")
			print(overlapping_haploblocks)

			if upstream_block:
				print(f"Upstream block: {upstream_block}")
				overlapping_haploblocks.insert(0, upstream_block)
			if downstream_block:
				print(f"Downstream block: {downstream_block}")
				overlapping_haploblocks.append(downstream_block)

			# Step 1: Extend each haploblock independently
			extended_haploblocks = []
			for block_start, block_stop in overlapping_haploblocks:
				upstream_het_sites = [h for h in het_sites if h < block_start]
				if upstream_het_sites:
					extended_start = max(upstream_het_sites) + 1
				else:
					extended_start = block_start
				# removed clamping to gene_start

				downstream_het_sites = [h for h in het_sites if h > block_stop]
				if downstream_het_sites:
					extended_stop = min(downstream_het_sites) - 1
				else:
					extended_stop = block_stop
				# removed clamping to gene_stop

				extended_haploblocks.append((extended_start, extended_stop))

			print(f"Extended haploblocks: {len(extended_haploblocks)}")
			print(extended_haploblocks)
			
			# Filter out extended haploblocks that have 0 base overlap with the gene
			overlapping_extended_haploblocks = []
			for start, stop in extended_haploblocks:
				overlap_start = max(start, gene_start)
				overlap_stop = min(stop, gene_stop)
				overlap_length = max(0, overlap_stop - overlap_start)
				if overlap_length > 0:
					overlapping_extended_haploblocks.append((start, stop))
			
			print(f"Extended haploblocks with gene overlap: {len(overlapping_extended_haploblocks)}")
			print(overlapping_extended_haploblocks)

			# Step 2: Merge overlapping extended haploblocks
			# (left unchanged)

			# Step 3: Best haploblock selection (unchanged)
			best_haploblock = None
			if gene in ARS_dict:
				ARS_start, ARS_stop = map(int, ARS_dict[gene].split(":")[1].split("-"))
				ars_spanning_blocks = [(start, stop) for start, stop in overlapping_extended_haploblocks 
									if start <= ARS_start and stop >= ARS_stop]
				largest_ars_spanning_block = max(ars_spanning_blocks, key=lambda x: x[1] - x[0]) if ars_spanning_blocks else None
				
				if largest_ars_spanning_block:
					best_haploblock = largest_ars_spanning_block
					print(f"{sample_ID} {gene} largest ARS-spanning haploblock: chr6:{largest_ars_spanning_block[0]}-{largest_ars_spanning_block[1]}")
				else:
					print(f"{sample_ID} {gene} no haploblocks span the ARS")
			
			if not best_haploblock:
				max_overlap = 0
				for start, stop in overlapping_extended_haploblocks:
					overlap_start = max(start, gene_start)
					overlap_stop = min(stop, gene_stop)
					overlap_length = max(0, overlap_stop - overlap_start)
					if overlap_length > max_overlap:
						max_overlap = overlap_length
						best_haploblock = (start, stop)

			if best_haploblock:
				overlap_start = max(best_haploblock[0], gene_start)
				overlap_stop = min(best_haploblock[1], gene_stop)
				overlap_length = max(0, overlap_stop - overlap_start)
				prop_overlap = overlap_length / gene_length
				largest_overlap_string = f"{prop_overlap*100:.2f}%"
			else:
				largest_overlap_string = "0.00%"

			incomplete_data.append([sample_ID, gene, num_pre_merge_haploblocks, largest_overlap_string])
			print(f"Proportion of gene contained in largest overlapping haploblock: {largest_overlap_string}")

	with open(output_file, "w", newline="") as csv_file:
		writer = csv.writer(csv_file, delimiter="\t")
		writer.writerow(["sample", "num_genes", "genes"])
		writer.writerow([sample_ID, len(gene_list), ",".join(gene_list)])


	if incomplete_data:
		with open(incomplete_file, "w", newline="") as csvfile:
			csv_writer = csv.writer(csvfile)
			csv_writer.writerow(["sample", "gene", "num_haploblocks", "largest_haploblock"])
			csv_writer.writerows(incomplete_data)

	return gene_list, unphased_genes
