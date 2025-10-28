import os
import subprocess
from .preprocess_methods import (
	trim_adapters,
	align_to_reference_minimap,
	bait_DRB_paralogs,
	align_to_reference_vg,
	reassign_mapq,
	mark_duplicates_picard,
	filter_reads,
	call_variants_bcftools,
	call_variants_deepvariant,
	call_variants_clair3,
	call_structural_variants_sniffles,
	phase_genotypes_longphase,
	merge_longphase_vcfs
)
from .config import min_reads_sample

def preprocess_ont_sample(config):
	# trim_adapters(
	# 	adapters=config['adapters'],
	# 	input_file=config['raw_fastq'],
	# 	output_file=config['trimmed_fastq'],
	# 	sample_ID=config['sample_ID'],
	# 	threads=config['threads'],
	# 	adapter_file=config['adapter_file'],
	# 	five_prime_adapter=config['five_prime_adapter'],
	# 	three_prime_adapter=config['three_prime_adapter']
	# )
	
	# align_to_reference_minimap(
	# 	input_file=config['trimmed_fastq'],
	# 	output_file=config['hg38_bam'],
	# 	read_group_string=config['read_group_string'],
	# 	reference_fasta=config['reference_genome'],
	# 	platform=config['platform'],
	# 	threads=config['threads']
	# )

	# bait_DRB_paralogs(
	# 	input_file=config['trimmed_fastq'],
	# 	output_file=config['hg38_bam_drb'],
	# 	DRB34_reads_file=config['DRB34_reads_file'],
	# 	read_group_string=config['read_group_string'],
	# 	reference_fasta=config['dummy_reference'],
	# 	platform=config['platform'],
	# 	threads=config['threads']
	# )
	
	# vg mapping discontinued
	# if config['aligner'] == "vg":
	# 	align_to_reference_vg(
	# 		vg=config['vg'],
	# 		input_file=config['trimmed_fastq'],
	# 		output_file=config['pangenome_bam'],  
	# 		sample_ID=config['sample_ID'],
	# 		read_group_string=config['read_group_string'],
	# 		reference_gbz=config['reference_gbz'],
	# 		ref_paths=config['ref_paths'],
	# 		platform=config['platform'],
	# 		threads=config['threads']
	# 		)
	
	# 	reassign_mapq(
	# 		bam_hg38=config['hg38_bam'],
	# 		bam_pg=config['pg_bam'],
	# 		reassigned_pg=config['pg_mapq_reassign_bam']
	# 	)
		
	# 	mark_duplicates_picard(
	# 		input_file=config['pg_mapq_reassign_bam'],
	# 		output_file=config['pg_mapq_reassign_mrkdup_bam'],
	# 		metrics_file=config['pg_mapq_reassign_mrkdup_metrics'],
	# 		temp_dir=os.path.join(config['mapped_bam_dir'], "mark_duplicates"),
	# 		picard=config['picard']
	# 	)

	# 	chr6_read_count = filter_reads(
	# 		input_file=config['pg_mapq_reassign_mrkdup_bam'],
	# 		output_file=config['hg38_rmdup_chr6_bam'],
	# 		threads=config['threads']
	# 	)

	# filter_reads(
	# 	input_file=config['hg38_bam'],
	# 	output_file=config['hg38_chr6_bam'],
	# 	DRB34_reads_file=config['DRB34_reads_file'],
	# 	threads=config['threads']
	# )
	
	# mark_duplicates_picard(
	# 	input_file=config['hg38_chr6_bam'],
	# 	output_file=config['hg38_rmdup_chr6_bam'],
	# 	metrics_file=config['hg38_mrkdup_metrics'],
	# 	temp_dir=os.path.join(config['mapped_bam_dir'], "mark_duplicates"),
	# 	picard=config['picard']
	# )

	chr6_read_count = int(subprocess.check_output(f"samtools view -c {config['hg38_rmdup_chr6_bam']}", shell=True).strip())
	if chr6_read_count >= min_reads_sample:
		# if config['genotyper'] == "bcftools":
		# 	call_variants_bcftools(
		# 		input_file=config['hg38_rmdup_chr6_bam'],
		# 		output_file=config['snv_vcf'],
		# 		reference_fasta=config['reference_genome'],
		# 		platform=config['platform'],
		# 		threads=config['threads']
		# 	)

		# elif config['genotyper'] == "deepvariant":
		# 	call_variants_deepvariant(
		# 		input_bam=config['hg38_rmdup_chr6_bam'],
		# 		output_vcf=config['snv_vcf'],
		# 		output_gvcf=config['snv_gvcf'],
		# 		platform=config['platform'],
		# 		deepvariant_sif=config['deepvariant_sif'],
		# 		reference_fasta=config['reference_genome'],
		# 		genotypes_dir=config['genotypes_dir'],
		# 		mapped_bam_dir=config['mapped_bam_dir'],
		# 		sample_ID=config['sample_ID']
		# 	)
		# elif config['genotyper'] == "clair3":
		call_variants_clair3(
			input_bam=config['hg38_rmdup_chr6_bam'],
			output_vcf=config['snv_vcf'],
			platform=config['platform'],
			reference_fasta=config['reference_genome'],
			threads=config['threads'],
			chr6_bed=config['chr6_bed'],
			clair3_ont_model_path=config['clair3_ont_model_path'],
			clair3_hifi_model_path=config['clair3_hifi_model_path'],
			genotypes_dir=config['genotypes_dir'],
			sample_ID=config['sample_ID']
		)
		# call_structural_variants_sniffles(
		# 	input_bam=config['hg38_rmdup_chr6_bam'],
		# 	output_vcf=config['sv_vcf'],
		# 	threads=config['threads'],
		# 	reference_fasta=config['reference_genome'],
		# 	chr6_bed=config['chr6_bed'],
		# 	tandem_repeat_bed=config['tandem_repeat_bed']
		# )
		# phase_genotypes_longphase(
		# 	input_bam=config['hg38_rmdup_chr6_bam'],
		# 	input_SNV_vcf=config['snv_vcf'],
		# 	input_SV_vcf=config['sv_vcf'],
		# 	output_blocks_file=config['phased_haploblocks'],
		# 	output_gtf_file=config['phased_haploblocks_gtf'],
		# 	phased_vcf=config['longphase_vcf'],
		# 	phased_SV_vcf=config['longphase_sv_vcf'],
		# 	haplotagged_bam=config['hg38_rmdup_chr6_haplotag_bam'],
		# 	longphase=config['longphase'],
		# 	reference_fasta=config['reference_genome'],
		# 	threads=config['threads'],
		# 	phased_vcf_dir=config['phased_vcf_dir'],
		# 	sample_ID=config['sample_ID']
		# )
		# merge_longphase_vcfs(
		# 	phased_vcf=config['longphase_vcf'],
		# 	phased_SV_vcf=config['longphase_sv_vcf'],
		# 	merged_vcf=config['longphase_merged_vcf'],
		# 	reference_fasta=config['reference_genome'],
		# 	phased_vcf_dir=config['phased_vcf_dir'],
		# 	sample_ID=config['sample_ID']
		# )
	
	else:
		print("Insufficient reads for variant calling")
		print("Sample {} had {} reads!".format(config['sample_ID'], chr6_read_count))