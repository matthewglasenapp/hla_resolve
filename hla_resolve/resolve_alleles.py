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

def resolve_alleles(
	config,
	mosdepth_regions_file,
	depth_thresh,
	prop_20x_thresh,
	prop_30x_thresh,
	mhc_start,
	mhc_stop,
	genes_bed,
	genes_of_interest,
	hla_genes_regions_file,
	vcf2fasta_script,
	reference_genome,
	DNA_bases,
	stop_codons,
	IMGT_XML
):
	"""
	Main orchestration function that runs the complete allele resolution workflow:
	1. Coverage analysis
	2. FASTA sequence reconstruction 
	3. HLA typing
	"""
	print("Starting HLA allele resolution workflow...")
	
	# Step 1: Coverage analysis
	print("Step 1: Running coverage analysis...")
	sufficient_coverage_genes = run_coverage_analysis(
		config=config,
		mosdepth_regions_file=mosdepth_regions_file,
		depth_thresh=depth_thresh,
		prop_20x_thresh=prop_20x_thresh,
		prop_30x_thresh=prop_30x_thresh
	)
	
	# Step 2: FASTA reconstruction
	print("Step 2: Reconstructing FASTA sequences...")
	phased_genes = reconstruct_fasta_sequences(
		config=config,
		sufficient_coverage_genes=sufficient_coverage_genes,
		mhc_start=mhc_start,
		mhc_stop=mhc_stop,
		genes_bed=genes_bed,
		genes_of_interest=genes_of_interest,
		hla_genes_regions_file=hla_genes_regions_file,
		vcf2fasta_script=vcf2fasta_script,
		reference_genome=reference_genome,
		DNA_bases=DNA_bases,
		stop_codons=stop_codons
	)
	
	# Step 3: HLA typing
	print("Step 3: Typing HLA alleles...")
	type_hla_alleles(
		config=config,
		phased_genes=phased_genes,
		IMGT_XML=IMGT_XML
	)
	
	print("HLA allele resolution workflow completed!")
	return phased_genes

def run_coverage_analysis(
	config,
	mosdepth_regions_file,
	depth_thresh,
	prop_20x_thresh,
	prop_30x_thresh
):
	run_mosdepth(
		input_file=os.path.join(config['mapped_bam_dir'], config['sample_ID'] +".hg38.rmdup.chr6.bam"),
		output_dir=config['mosdepth_dir'],
		sample_ID=config['sample_ID'],
		regions_file=mosdepth_regions_file,
		threads=config['threads']
	)
	
	sufficient_coverage_genes = parse_mosdepth(
		regions_file=os.path.join(config['mosdepth_dir'], config['sample_ID'] + ".regions.bed.gz"),
		thresholds_file=os.path.join(config['mosdepth_dir'], config['sample_ID'] + ".thresholds.bed.gz"), 
		depth_thresh=depth_thresh,
		prop_20x_thresh=prop_20x_thresh,
		prop_30x_thresh=prop_30x_thresh
	)
	return sufficient_coverage_genes

def reconstruct_fasta_sequences(
	config,
	sufficient_coverage_genes,
	mhc_start,
	mhc_stop,
	genes_bed,
	genes_of_interest,
	hla_genes_regions_file,
	vcf2fasta_script,
	reference_genome,
	DNA_bases,
	stop_codons
):
	if config['platform'] == "PACBIO":
		phased_vcf = os.path.join(config['phased_vcf_dir'], config['sample_ID'] + ".hiphase.joint.vcf.gz")
		haploblock_file = os.path.join(config['phased_vcf_dir'], config['sample_ID'] + ".phased.blocks.txt")
	elif config['platform'] == "ONT":
		phased_vcf = os.path.join(config['phased_vcf_dir'], config['sample_ID'] + ".longphase.vcf.gz")
		haploblock_file = os.path.join(config['phased_vcf_dir'], config['sample_ID'] + ".phased.haploblocks.txt")
	
	heterozygous_sites, haploblock_list = parse_haploblocks(
		phased_vcf=phased_vcf,
		haploblock_file=haploblock_file,
		sample_ID=config['sample_ID'],
		platform=config['platform'],
		mhc_start=mhc_start,
		mhc_stop=mhc_stop
	)

	phased_genes = evaluate_gene_haploblocks(
		phased_genes_file=os.path.join(config['parsed_haploblock_dir'], f"phased_genes.tsv"),
		incomplete_genes_file=os.path.join(config['parsed_haploblock_dir'], f"incomplete.csv"),
		sample_ID=config['sample_ID'],
		genes_bed=genes_bed,  
		genes_of_interest=genes_of_interest,
		heterozygous_sites=heterozygous_sites, 
		haploblock_list=haploblock_list)
	
	if config['platform'] == "PACBIO":
		input_vcf = os.path.join(config['phased_vcf_dir'], config['sample_ID'] + ".hiphase.joint.vcf.gz")
	elif config['platform'] == "ONT":
		input_vcf = os.path.join(config['phased_vcf_dir'], config['sample_ID'] + ".longphase.merged.vcf.gz")
	
	filter_vcf(
		input_vcf=input_vcf,
		pass_vcf=os.path.join(config['phased_vcf_dir'], f"{config['sample_ID']}_PASS.vcf.gz"),
		fail_vcf=os.path.join(config['phased_vcf_dir'], f"{config['sample_ID']}_FAIL.vcf.gz"),
		pass_unphased_vcf=os.path.join(config['phased_vcf_dir'], f"{config['sample_ID']}_PASS_UNPHASED.vcf.gz"),
		filtered_vcf=os.path.join(config['filtered_vcf_dir'], f"{config['sample_ID']}_filtered.vcf.gz"),
		unphased_tsv=os.path.join(config['phased_vcf_dir'], config['sample_ID'] + ".unphased.tsv"),
		platform=config['platform'],
		genotyper=config['genotyper'],
		hla_genes_regions_file=hla_genes_regions_file
	)
	
	# Reset vcf2fasta_out_dir for sequential runs 
	if any(os.scandir(config['vcf2fasta_out_dir'])):
		shutil.rmtree(config['vcf2fasta_out_dir'])
		os.makedirs(config['vcf2fasta_out_dir'], exist_ok=True)

	for gene in phased_genes:
		if gene in genes_of_interest and gene in sufficient_coverage_genes:
			run_vcf2fasta(
				vcf2fasta=vcf2fasta_script,
				input_vcf=os.path.join(config['filtered_vcf_dir'], f"{config['sample_ID']}_filtered.vcf.gz"),
				output_dir=os.path.join(config['vcf2fasta_out_dir'], gene),
				input_gff=os.path.join(config['gff_dir'], gene + "_cds_sorted.gff3"),
				reference_genome=reference_genome,
				gene=gene, 
				feature="gene")
			
			run_vcf2fasta(
				vcf2fasta=vcf2fasta_script,
				input_vcf=os.path.join(config['filtered_vcf_dir'], f"{config['sample_ID']}_filtered.vcf.gz"),
				output_dir=os.path.join(config['vcf2fasta_out_dir'], gene),
				input_gff=os.path.join(config['gff_dir'], gene + "_gene.gff3"),
				reference_genome=reference_genome,
				gene=gene, 
				feature="CDS")
	
	parse_fastas(
		vcf2fasta_out_dir=config['vcf2fasta_out_dir'],
		output_gene_fasta=os.path.join(config['hla_fasta_dir'], config['sample_ID'] + "_HLA_haplotypes_gene.fasta"),
		output_cds_fasta=os.path.join(config['hla_fasta_dir'], config['sample_ID'] + "_HLA_haplotypes_CDS.fasta"),
		DNA_bases=DNA_bases,
		stop_codons=stop_codons
	)

	return phased_genes

def type_hla_alleles(
	config,
	phased_genes,
	IMGT_XML
):
	original_dir = os.getcwd()
	os.chdir(config['hla_typing_dir'])
	# Lazy import to avoid overhead when not using HLA typing
	from hla_resolve.hla_typer import main as classify_hla_alleles
	classify_hla_alleles(
		imgt_xml=IMGT_XML, 
		hla_fasta_dir=config['hla_fasta_dir'], 
		sample_ID=config['sample_ID']
	)
	# Note: print_results() and clean_up functionality would need to be handled differently
	# since we no longer have access to the sample object

	if config.get('clean_up', False):
		# Define directories to clean up - this would need to be passed in config
		combined_dirs = [
			config['fastq_raw_dir'],
			config['fastq_trimmed_dir'],
			config['mapped_bam_dir'],
			config['genotypes_dir'],
			config['sv_dir'],
			config['phased_vcf_dir'],
			config['mosdepth_dir'],
			config['parsed_haploblock_dir'],
			config['filtered_vcf_dir'],
			config['vcf2fasta_out_dir']
		]
		for directory in combined_dirs:
			if os.path.exists(directory) and directory != config['hla_typing_dir']:
				shutil.rmtree(directory)
	
	os.chdir(original_dir)

