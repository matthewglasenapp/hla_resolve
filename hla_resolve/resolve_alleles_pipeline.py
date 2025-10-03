import os
import shutil
from preprocess_methods import (
	run_mosdepth,
	parse_mosdepth
)
from investigate_haploblocks_methods import (
	parse_haploblocks,
	evaluate_gene_haploblocks
)
from reconstruct_fasta_methods import (
	filter_vcf,
	run_vcf2fasta,
	parse_fastas
)
from hla_typer import main as classify_hla_alleles

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
			results = f.read().splitlines()[1].split(",")[1:]
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
	# Step 1: Coverage analysis
	print("Running coverage analysis...")
	run_mosdepth(
		input_file=config['hg38_rmdup_chr6_bam'],
		output_dir=config['mosdepth_dir'],
		sample_ID=config['sample_ID'],
		regions_file=config['mosdepth_regions_file'],
		threads=config['threads']
	)
	
	sufficient_coverage_genes = parse_mosdepth(
		regions_file=config['mosdepth_regions'],
		thresholds_file=config['mosdepth_thresholds'], 
		depth_thresh=config['depth_thresh'],
		prop_20x_thresh=config['prop_20x_thresh'],
		prop_30x_thresh=config['prop_30x_thresh']
	)
	
	# Step 2: FASTA reconstruction
	print("Step 2: Reconstructing FASTA sequences...")
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

	phased_genes, unphased_genes = evaluate_gene_haploblocks(
		output_file=config['phased_genes_tsv'],
		incomplete_file=config['incomplete_genes_csv'],
		sample_ID=config['sample_ID'],
		genes_bed=config['genes_bed'],  
		genes_of_interest=config['genes_of_interest'],
		het_sites=heterozygous_sites, 
		haploblocks=haploblock_list,
		ARS_dict=config.get('ARS_dict', None))
	
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
	
	if config['platform'] == "PACBIO":
		input_vcf = config['hiphase_joint_vcf']
	elif config['platform'] == "ONT":
		input_vcf = config['longphase_merged_vcf']
	
	filter_vcf(
		input_vcf=input_vcf,
		pass_vcf=config['pass_vcf'],
		fail_vcf=config['fail_vcf'],
		pass_unphased=config['pass_unphased_vcf'],
		filtered_vcf=config['filtered_vcf'],
		unphased_overlap_tsv=config['unphased_tsv'],
		platform=config['platform'],
		genotyper=config['genotyper'],
		hla_genes_regions_file=config['hla_genes_regions_file']
	)
	
	# Reset vcf2fasta_out_dir for sequential runs 
	if any(os.scandir(config['vcf2fasta_out_dir'])):
		shutil.rmtree(config['vcf2fasta_out_dir'])
		os.makedirs(config['vcf2fasta_out_dir'], exist_ok=True)

	for gene in config['genes_of_interest']:
		if gene in sufficient_coverage_genes or gene in unphased_genes:
			gff_gene_name = convert_gene_name_for_gff(gene)
			run_vcf2fasta(
				vcf2fasta=config['vcf2fasta_script'],
				input_vcf=config['filtered_vcf'],
				output_dir=os.path.join(config['vcf2fasta_out_dir'], gene),
				input_gff=os.path.join(config['gff_dir'], gff_gene_name + "_gene.gff3"),
				reference_genome=config['reference_genome'],
				gene=gene, 
				feature="gene")
			
			run_vcf2fasta(
				vcf2fasta=config['vcf2fasta_script'],
				input_vcf=config['filtered_vcf'],
				output_dir=os.path.join(config['vcf2fasta_out_dir'], gene),
				input_gff=os.path.join(config['gff_dir'], gff_gene_name + "_cds_sorted.gff3"),
				reference_genome=config['reference_genome'],
				gene=gene, 
				feature="CDS")
		
		elif gene in sufficient_coverage_genes and gene in unphased_genes:
			print(f"Gene {gene} was not fully phased. Recovering antigen recognition sequence and largest haplotype block.")
			print(f"Sample {config['sample_ID']} {gene} best haploblock: {unphased_genes[gene]}")
			#Modify vcf2fasta output
	
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
		gff_dir=config['gff_dir']
	)
	
	# Step 3: HLA typing
	print("Step 3: Typing HLA alleles...")
	original_dir = os.getcwd()
	os.chdir(config['hla_typing_dir'])
	
	classify_hla_alleles(
		reference_xml_file=config['IMGT_XML'], 
		hla_fasta_dir=config['hla_fasta_dir'], 
		sample_ID=config['sample_ID']
	)
	
	os.chdir(original_dir)
	
	# Step 4: Print results
	print("Step 4: Printing HLA typing results...")
	print_results(config)
	
	print("HLA allele resolution workflow completed!")

