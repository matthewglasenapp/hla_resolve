# This software is Copyright ©2026. The Regents of the University of California
# ("Regents"). All Rights Reserved.
#
# Licensed under the UC Santa Cruz Noncommercial License (the "License").
# You may not use this file except in compliance with the License.
# A copy of the License is included in this repository as LICENSE.txt.

import os
import shutil
from .preprocess_methods import (
	run_mosdepth,
	parse_mosdepth
)
from .investigate_haploblocks_methods import (
	parse_haploblocks,
	evaluate_gene_haploblocks
)
from .reconstruct_fasta_methods import (
	filter_vcf_gene,
	run_vcf2fasta,
	parse_fastas
)
from .hla_typer import main as classify_hla_alleles

def convert_gene_name_for_gff(gene_name):
	"""
	Convert gene name from config format (HLA-A) to GFF filename format (hla_a)
	"""
	return gene_name.lower().replace('-', '_')

def print_results(config):
	"""
	Print the HLA typing results from the output file.
	"""
	results_file = os.path.join(config['hla_typing_dir'], "allele_output.csv")
	if os.path.exists(results_file):
		with open(results_file, "r") as f:
			lines = f.read().splitlines()
		
		# Check if file has at least 2 lines (header + data)
		if len(lines) < 2:
			print(f"Warning: Results file {results_file} exists but contains insufficient data (only {len(lines)} lines)")
			return
		
		# Check if the data line has enough columns
		data_line = lines[1]
		data_parts = data_line.split(",")
		if len(data_parts) < 2:
			print(f"Warning: Results file {results_file} data line has insufficient columns: {data_line}")
			return
		
		results = data_parts[1:]
		print(f"{config['sample_ID']} HLA Star Allele Calls:")
		for item in results:
			print(f"  {item}")
	else:
		print(f"Warning: Results file not found at {results_file}")

def resolve_alleles(config):
	"""
	Main orchestration function that runs the complete allele resolution workflow:
	1. Coverage analysis
	2. FASTA sequence reconstruction 
	3. HLA typing
	4. Print results
	"""	
	print("Running coverage analysis...")
	run_mosdepth(
		input_file=config['hg38_rmdup_chr6_bam'],
		output_dir=config['mosdepth_dir'],
		sample_ID=config['sample_ID'],
		regions_file=config['hla_genes_regions_file'],
		threads=config['threads']
	)
	
	sufficient_coverage_genes = parse_mosdepth(
		regions_file=config['mosdepth_regions'],
		thresholds_file=config['mosdepth_thresholds'], 
		depth_thresh=config['depth_thresh'],
		prop_20x_thresh=config['prop_20x_thresh'],
		prop_30x_thresh=config['prop_30x_thresh']
	)
	
	print("Reconstructing FASTA sequences")
	if config['platform'] == "PACBIO":
		phased_vcf = config['hiphase_joint_vcf']
		haploblock_file = config['phased_blocks']
	elif config['platform'] == "ONT":
		phased_vcf = config['longphase_vcf']
		haploblock_file = config['phased_haploblocks']
	
	heterozygous_sites, haploblock_list = parse_haploblocks(
		input_vcf=phased_vcf,
		input_haploblock_file=haploblock_file,
		platform=config['platform'],
		sample_ID=config['sample_ID'],
		mhc_start=config['mhc_start'],
		mhc_stop=config['mhc_stop']
	)

	phased_genes, unphased_genes, do_not_type_genes, cds_rescued_genes = evaluate_gene_haploblocks(
		output_file=config['phased_genes_tsv'],
		incomplete_file=config['incomplete_genes_csv'],
		sample_ID=config['sample_ID'],
		genes_bed=config['genes_bed'],
		genes_of_interest=config['genes_of_interest'],
		het_sites=heterozygous_sites,
		haploblocks=haploblock_list,
		ARS_dict=config.get('ARS_dict', None),
		CDS_dict=config.get('CDS_dict', None))
	
	# Print which genes were successfully phased
	print("Fully Phased Genes:")
	for gene in config['genes_of_interest']:
		if gene in phased_genes:
			print(f"  {gene}")
	print("\n")

	print("Partially Phased Genes:")
	for gene in config['genes_of_interest']:
		if gene in unphased_genes:
			print(f"  {gene}")
	print("\n")

	print("CDS-Rescued Genes:")
	for gene in config['genes_of_interest']:
		if gene in cds_rescued_genes:
			tier = cds_rescued_genes[gene]["tier"]
			print(f"  {gene} (tier: {tier})")
	print("\n")
	
	if config['platform'] == "PACBIO":
		input_vcf = config['hiphase_joint_vcf']
	elif config['platform'] == "ONT":
		input_vcf = config['longphase_merged_vcf']

	# Filter phased VCF by gene region
	gene_filtered_vcfs = {}
	for gene in config['genes_of_interest']:
		gff_gene_name = convert_gene_name_for_gff(gene)
		gff_file = os.path.join(config['gff_dir'], gff_gene_name + "_gene.gff3")
		gff_lines = [line.split("\t") for line in open(gff_file, "r").read().splitlines() if not line.startswith("#")]
		gff_record = gff_lines[0]
		chromosome, start, stop = gff_record[0], gff_record[3], gff_record[4]
		filter_region = f"{chromosome}:{start}-{stop}"

		# Create gene-specific output file paths in single_gene_vcf subdirectory
		gene_symbolic_vcf = os.path.join(config['single_gene_vcf_dir'], f"{config['sample_ID']}_{gene}.symbolic.vcf.gz")
		gene_pass_vcf = os.path.join(config['single_gene_vcf_dir'], f"{config['sample_ID']}_{gene}_PASS.vcf.gz")
		gene_fail_vcf = os.path.join(config['single_gene_vcf_dir'], f"{config['sample_ID']}_{gene}_FAIL.vcf.gz")
		gene_sv_overlap_vcf = os.path.join(config['single_gene_vcf_dir'], f"{config['sample_ID']}_{gene}_SV_OVERLAP.vcf.gz")
		gene_pass_unphased_vcf = os.path.join(config['single_gene_vcf_dir'], f"{config['sample_ID']}_{gene}_PASS_UNPHASED.vcf.gz")
		gene_filtered_vcf = os.path.join(config['single_gene_vcf_dir'], f"{config['sample_ID']}_{gene}_PASS_phased.vcf.gz")

		# Derive genotyper label from split caller config
		if config['snp_caller'] == config['indel_caller']:
			genotyper = config['snp_caller']
		else:
			genotyper = "hybrid"

		filter_vcf_gene(
			input_vcf=input_vcf,
			gene=gene,
			filter_region=filter_region,
			symbolic_vcf=gene_symbolic_vcf,
			pass_vcf=gene_pass_vcf,
			fail_vcf=gene_fail_vcf,
			sv_overlap_vcf=gene_sv_overlap_vcf,
			pass_unphased=gene_pass_unphased_vcf,
			filtered_vcf=gene_filtered_vcf,
			platform=config['platform'],
			genotyper=genotyper,
			hla_genes_regions_file=config['hla_genes_regions_file'],
			force_include_unphased=(gene in cds_rescued_genes)
		)
		
		gene_filtered_vcfs[gene] = gene_filtered_vcf
	
	# Reset vcf2fasta_out_dir for sequential runs 
	if any(os.scandir(config['vcf2fasta_out_dir'])):
		shutil.rmtree(config['vcf2fasta_out_dir'])
		os.makedirs(config['vcf2fasta_out_dir'], exist_ok=True)

	for gene in config['genes_of_interest']:
		if gene in sufficient_coverage_genes and gene not in do_not_type_genes:
			gene_filtered_vcf = gene_filtered_vcfs.get(gene)
				
			gff_gene_name = convert_gene_name_for_gff(gene)
			run_vcf2fasta(
				vcf2fasta=config['vcf2fasta_script'],
				input_vcf=gene_filtered_vcf,
				output_dir=os.path.join(config['vcf2fasta_out_dir'], gene),
				input_gff=os.path.join(config['gff_dir'], gff_gene_name + "_gene.gff3"),
				reference_genome=config['reference_genome'],
				gene=gene, 
				feature="gene")
			
			run_vcf2fasta(
				vcf2fasta=config['vcf2fasta_script'],
				input_vcf=gene_filtered_vcf,
				output_dir=os.path.join(config['vcf2fasta_out_dir'], gene),
				input_gff=os.path.join(config['gff_dir'], gff_gene_name + "_cds_sorted.gff3"),
				reference_genome=config['reference_genome'],
				gene=gene, 
				feature="CDS")
	
	parse_fastas(
		sample_ID=config['sample_ID'],
		vcf2fasta_output_dir=config['vcf2fasta_out_dir'],
		outfile_gene=config['hla_gene_fasta'],
		outfile_CDS=config['hla_cds_fasta'],
		DNA_bases=config['DNA_bases'],
		stop_codons=config['stop_codons'],
		unphased_genes=unphased_genes,
		gene_dict=config['gene_dict'],
		CDS_dict=config['CDS_dict'],
		gff_dir=config['gff_dir'],
		cds_rescued_genes=cds_rescued_genes,
		ARS_dict=config.get('ARS_dict', None),
		CLASS_I_GENES=config.get('CLASS_I_GENES', None)
	)
	
	# Step 3: HLA typing
	print("Typing HLA Alleles!")
	original_dir = os.getcwd()
	os.chdir(config['hla_typing_dir'])

	try:
		classify_hla_alleles(
			reference_xml_file=config['IMGT_XML'],
			hla_fasta_dir=config['hla_fasta_dir'],
			sample_ID=config['sample_ID'],
			generate_query_ref_comp=True
		)
	finally:
		os.chdir(original_dir)
	
	# Step 4: Print results
	print("Step 4: Printing HLA typing results...")
	print_results(config)
	
	print("HLA allele resolution workflow completed!")
