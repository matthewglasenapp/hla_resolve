# This software is Copyright ©2026. The Regents of the University of California
# ("Regents"). All Rights Reserved.
#
# See LICENSE.txt for license details.

import os
import csv
import pysam

# Get list of haploblock intervals for MHC
def parse_haploblocks(input_vcf, input_haploblock_file, platform,sample_ID, mhc_start, mhc_stop):
	heterozygous_sites = []
	haploblock_list = []

	vcf = pysam.VariantFile(input_vcf)

	# Get the sample name from the VCF header rather than relying on sample_ID,
	# since WGS/WES BAMs may have a different SM tag than the user-provided sample_ID
	vcf_samples = list(vcf.header.samples)
	if len(vcf_samples) != 1:
		raise ValueError(f"Expected exactly 1 sample in {input_vcf}, found {len(vcf_samples)}: {vcf_samples}")
	sample_name = vcf_samples[0]

	for record in vcf:
		if record.chrom != "chr6":
			continue
		if record.pos < mhc_start or record.pos > mhc_stop:
			continue

		# Skip non-PASS records from DV/Clair3 (e.g. RefCall).
		# bcftools/freebayes records have FILTER=. (empty keys) and pass through;
		# their quality filtering is already applied upstream by the caller command.
		if record.filter.keys() != ["PASS"] and record.filter.keys() != []:
			continue

		gt = record.samples[sample_name]["GT"]
		if gt is not None and len(gt) == 2 and gt[0] is not None and gt[1] is not None and gt[0] != gt[1]:
			heterozygous_sites.append(record.pos)

	#print(f"Sample {sample_name} has {len(heterozygous_sites)} heterozygous extended MHC genotypes")

	print(f"Parsing {sample_name} haploblock file: {input_haploblock_file}")
	print("\n")

	with open(input_haploblock_file, "r") as f:
		haploblocks = f.read().splitlines()

	for line in haploblocks[1:]:
		fields = line.split("\t")
		
		# HiPhase and Longphase have slightly differently formatted haploblock tsv files. 
		if platform == "PACBIO":
			chromosome, start, stop = "chr6", int(fields[4]), int(fields[5])
		elif platform == "ONT":
			chromosome, start, stop = "chr6", int(fields[3]), int(fields[4])

		if chromosome == "chr6" and stop > mhc_start:
			haploblock_list.append([start,stop])

	return heterozygous_sites, haploblock_list

# Check whether each captured MHC gene is completely spanned by a haploblock
def evaluate_gene_haploblocks(output_file, incomplete_file, sample_ID, genes_bed, genes_of_interest, het_sites, haploblocks, ARS_dict=None, CDS_dict=None):
	# List of fully phased genes
	haploblocks.sort()
	gene_list = []
	unphased_genes = dict()
	do_not_type_genes = []
	cds_rescued_genes = dict()

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
			if gene in genes_of_interest:
				if len(gene_het_sites) == 0:
					print(f"{sample_ID} {gene} has 0 QC-Pass heterozygous genotypes")
					print(f"Treating as fully phased" + "\n")
				elif len(gene_het_sites) == 1:
					print(f"{sample_ID} {gene} has 1 QC-Pass heterozygous genotype")
					print(f"Treating as fully phased" + "\n")
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

				# Find closest het site downstream
				downstream_het_sites = [h for h in het_sites if h > block_stop]
				if downstream_het_sites:
					extended_stop = min(downstream_het_sites) - 1
				else:
					extended_stop = block_stop

				if extended_start <= gene_start and extended_stop >= gene_stop:
					gene_list.append(gene)
					fully_phased = True
					break

		if not fully_phased and gene in genes_of_interest:
			print(f"{sample_ID} {gene} is not fully phased")
			print(f"Entering haplotype rescue mode for {gene}!")
			print(f"Searching for the largest CDS-overlapped haplotype block!")
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

			num_raw_haploblocks = len(overlapping_haploblocks)
			print(f"{sample_ID} {gene} has {num_raw_haploblocks} overlapping unextended haploblocks")
			print(f"{sample_ID} {gene} overlapping unextended haploblocks: {overlapping_haploblocks}")

			if upstream_block:
				#print(f"Upstream block: {upstream_block}")
				overlapping_haploblocks.insert(0, upstream_block)
			if downstream_block:
				#print(f"Downstream block: {downstream_block}")
				overlapping_haploblocks.append(downstream_block)

			# Step 1: Extend each haploblock independently
			print("Extending haploblocks through homozygous sites")
			extended_haploblocks = []
			for block_start, block_stop in overlapping_haploblocks:
				upstream_het_sites = [h for h in het_sites if h < block_start]
				if upstream_het_sites:
					extended_start = max(upstream_het_sites) + 1
				else:
					extended_start = block_start

				downstream_het_sites = [h for h in het_sites if h > block_stop]
				if downstream_het_sites:
					extended_stop = min(downstream_het_sites) - 1
				else:
					extended_stop = block_stop

				extended_haploblocks.append((extended_start, extended_stop))

			#print(f"Extended haploblocks: {len(extended_haploblocks)}")
			#print(extended_haploblocks)
			
			# Filter out extended haploblocks that have 0 base overlap with the gene
			overlapping_extended_haploblocks = []
			for start, stop in extended_haploblocks:
				overlap_start = max(start, gene_start)
				overlap_stop = min(stop, gene_stop)
				overlap_length = max(0, overlap_stop - overlap_start)
				if overlap_length > 0:
					overlapping_extended_haploblocks.append((start, stop))
			
			print(f"{sample_ID} {gene} has {len(overlapping_extended_haploblocks)} extended haploblocks with gene overlap")
			print(f"{sample_ID} {gene} extended haploblocks with gene overlap: {overlapping_extended_haploblocks}")

			# Step 3: Best haploblock selection
			best_haploblock = None
			if gene in ARS_dict:
				ARS_start, ARS_stop = map(int, ARS_dict[gene].split(":")[1].split("-"))
				ars_spanning_blocks = [(start, stop) for start, stop in overlapping_extended_haploblocks 
									if start <= ARS_start and stop >= ARS_stop]
				largest_ars_spanning_block = max(ars_spanning_blocks, key=lambda x: x[1] - x[0]) if ars_spanning_blocks else None
				
				if largest_ars_spanning_block:
					best_haploblock = largest_ars_spanning_block
					unphased_genes[gene] = best_haploblock
					print(f"{sample_ID} {gene} has a haploblock that fully spans the antigen recognition sequence")
					print(f"{sample_ID} {gene} largest ARS-spanning haploblock: chr6:{largest_ars_spanning_block[0]}-{largest_ars_spanning_block[1]}")
				else:
					print(f"{sample_ID} {gene} no haploblocks span the ARS")

					# ========== CDS-AWARE RESCUE (branch 3b) ==========
					# Count ALL QC-pass het sites falling in CDS exons
					cds_hets = [h for h in gene_het_sites
								if any(s <= h <= e for s, e in CDS_dict.get(gene, []))] if CDS_dict else []

					# Count ALL QC-pass het sites falling in ARS CDS exons
					ars_cds_ranges = [(s, e) for s, e in CDS_dict.get(gene, [])
									if s >= ARS_start and e <= ARS_stop] if CDS_dict else []
					ars_cds_hets = [h for h in gene_het_sites
								if any(s <= h <= e for s, e in ars_cds_ranges)]

					print(f"{sample_ID} {gene} CDS hets: {len(cds_hets)}, ARS CDS hets: {len(ars_cds_hets)}")

					if CDS_dict and len(ars_cds_hets) <= 1:
						if len(cds_hets) <= 1:
							# Branch 3b-i: CDS rescue "cds_full"
							tier = "cds_full"
							print(f"{sample_ID} {gene} CDS-aware rescue tier: {tier} (total CDS hets={len(cds_hets)})")
						else:
							# Branch 3b-ii-A: CDS rescue "ars_only"
							tier = "ars_only"
							print(f"{sample_ID} {gene} CDS-aware rescue tier: {tier} (total CDS hets={len(cds_hets)}, ARS CDS hets={len(ars_cds_hets)})")

						cds_rescued_genes[gene] = {
							"tier": tier,
							"n_cds_hets": len(cds_hets),
							"n_ars_cds_hets": len(ars_cds_hets),
							"all_het_positions": sorted(gene_het_sites),
							"cds_het_positions": sorted(cds_hets),
							"ars_start": ARS_start,
							"ars_stop": ARS_stop,
							"ars_cds_ranges": ars_cds_ranges,
						}
					else:
						# Branch 3b-ii-B: ARS CDS hets > 1 — no typing
						print(f"{sample_ID} {gene} typing will not be performed (ARS CDS hets > 1)" + "\n")
						do_not_type_genes.append(gene)

					# Fall back to largest overlap for incomplete_data reporting
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

			incomplete_data.append([sample_ID, gene, num_raw_haploblocks, largest_overlap_string])
			if largest_ars_spanning_block:
				print(f"Proportion of gene contained in largest ARS-spanning haploblock: {largest_overlap_string}" + "\n")

	with open(output_file, "w", newline="") as csv_file:
		writer = csv.writer(csv_file, delimiter="\t")
		writer.writerow(["sample", "num_genes", "genes"])
		writer.writerow([sample_ID, len(gene_list), ",".join(gene_list)])


	if incomplete_data:
		with open(incomplete_file, "w", newline="") as csvfile:
			csv_writer = csv.writer(csvfile)
			csv_writer.writerow(["sample", "gene", "num_haploblocks", "largest_haploblock"])
			csv_writer.writerows(incomplete_data)

	#print(f"Unphased genes: {unphased_genes}")
	return gene_list, unphased_genes, do_not_type_genes, cds_rescued_genes
