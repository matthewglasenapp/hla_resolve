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

def run_coverage_analysis(
	sample,
	mosdepth_regions_file,
	depth_thresh,
	prop_20x_thresh,
	prop_30x_thresh
):
	run_mosdepth(
		input_file=os.path.join(sample.mapped_bam_dir, sample.sample_ID +".hg38.rmdup.chr6.bam"),
		output_dir=sample.mosdepth_dir,
		sample_ID=sample.sample_ID,
		regions_file=mosdepth_regions_file,
		threads=sample.threads
	)
	
	sufficient_coverage_genes = parse_mosdepth(
		regions_file=os.path.join(sample.mosdepth_dir, sample.sample_ID + ".regions.bed.gz"),
		thresholds_file=os.path.join(sample.mosdepth_dir, sample.sample_ID + ".thresholds.bed.gz"), 
		depth_thresh=depth_thresh,
		prop_20x_thresh=prop_20x_thresh,
		prop_30x_thresh=prop_30x_thresh
	)
	return sufficient_coverage_genes

def reconstruct_fasta_sequences(
	sample,
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
	if sample.platform == "PACBIO":
		phased_vcf = os.path.join(sample.phased_vcf_dir, sample.sample_ID + ".hiphase.joint.vcf.gz")
		haploblock_file = os.path.join(sample.phased_vcf_dir, sample.sample_ID + ".phased.blocks.txt")
	elif sample.platform == "ONT":
		phased_vcf = os.path.join(sample.phased_vcf_dir, sample.sample_ID + ".longphase.vcf.gz")
		haploblock_file = os.path.join(sample.phased_vcf_dir, sample.sample_ID + ".phased.haploblocks.txt")
	
	heterozygous_sites, haploblock_list = parse_haploblocks(
		phased_vcf=phased_vcf,
		haploblock_file=haploblock_file,
		sample_ID=sample.sample_ID,
		platform=sample.platform,
		mhc_start=mhc_start,
		mhc_stop=mhc_stop
	)

	phased_genes = evaluate_gene_haploblocks(
		phased_genes_file=os.path.join(sample.parsed_haploblock_dir, f"phased_genes.tsv"),
		incomplete_genes_file=os.path.join(sample.parsed_haploblock_dir, f"incomplete.csv"),
		sample_ID=sample.sample_ID,
		genes_bed=genes_bed,  
		genes_of_interest=genes_of_interest,
		heterozygous_sites=heterozygous_sites, 
		haploblock_list=haploblock_list)
	
	if sample.platform == "PACBIO":
		input_vcf = os.path.join(sample.phased_vcf_dir, sample.sample_ID + ".hiphase.joint.vcf.gz")
	elif sample.platform == "ONT":
		input_vcf = os.path.join(sample.phased_vcf_dir, sample.sample_ID + ".longphase.merged.vcf.gz")
	
	filter_vcf(
		input_vcf=input_vcf,
		pass_vcf=os.path.join(sample.phased_vcf_dir, f"{sample.sample_ID}_PASS.vcf.gz"),
		fail_vcf=os.path.join(sample.phased_vcf_dir, f"{sample.sample_ID}_FAIL.vcf.gz"),
		pass_unphased_vcf=os.path.join(sample.phased_vcf_dir, f"{sample.sample_ID}_PASS_UNPHASED.vcf.gz"),
		filtered_vcf=os.path.join(sample.filtered_vcf_dir, f"{sample.sample_ID}_filtered.vcf.gz"),
		unphased_tsv=os.path.join(sample.phased_vcf_dir, sample.sample_ID + ".unphased.tsv"),
		platform=sample.platform,
		genotyper=sample.genotyper,
		hla_genes_regions_file=hla_genes_regions_file
	)
	
	# Reset self.vcf2fasta_out_dir for sequential runs 
	if any(os.scandir(sample.vcf2fasta_out_dir)):
		shutil.rmtree(sample.vcf2fasta_out_dir)
		os.makedirs(sample.vcf2fasta_out_dir, exist_ok=True)

	for gene in phased_genes:
		if gene in genes_of_interest and gene in sufficient_coverage_genes:
			run_vcf2fasta(
				vcf2fasta=vcf2fasta_script,
				input_vcf=os.path.join(sample.filtered_vcf_dir, f"{sample.sample_ID}_filtered.vcf.gz"),
				output_dir=os.path.join(sample.vcf2fasta_out_dir, gene),
				input_gff=os.path.join(sample.gff_dir, gene + "_cds_sorted.gff3"),
				reference_genome=reference_genome,
				gene=gene, 
				feature="gene")
			
			run_vcf2fasta(
				vcf2fasta=vcf2fasta_script,
				input_vcf=os.path.join(sample.filtered_vcf_dir, f"{sample.sample_ID}_filtered.vcf.gz"),
				output_dir=os.path.join(sample.vcf2fasta_out_dir, gene),
				input_gff=os.path.join(sample.gff_dir, gene + "_gene.gff3"),
				reference_genome=reference_genome,
				gene=gene, 
				feature="CDS")
	
	parse_fastas(
		vcf2fasta_out_dir=sample.vcf2fasta_out_dir,
		output_gene_fasta=os.path.join(sample.hla_fasta_dir, sample.sample_ID + "_HLA_haplotypes_gene.fasta"),
		output_cds_fasta=os.path.join(sample.hla_fasta_dir, sample.sample_ID + "_HLA_haplotypes_CDS.fasta"),
		DNA_bases=DNA_bases,
		stop_codons=stop_codons
	)

	return phased_genes

def type_hla_alleles(
	sample,
	phased_genes,
	IMGT_XML
):
	original_dir = os.getcwd()
	os.chdir(sample.hla_typing_dir)
	# Lazy import to avoid overhead when not using HLA typing
	from hla_resolve.hla_typer import main as classify_hla_alleles
	classify_hla_alleles(
		imgt_xml=IMGT_XML, 
		hla_fasta_dir=sample.hla_fasta_dir, 
		sample_ID=sample.sample_ID
	)
	sample.print_results()

	if sample.clean_up:
		for directory in sample.combined_dirs:
			if os.path.exists(directory) and directory != sample.hla_typing_dir:
				shutil.rmtree(directory)
	
	os.chdir(original_dir)