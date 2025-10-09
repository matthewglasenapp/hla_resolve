import os
from preprocess_methods import (
	trim_adapters,
	run_fastqc,
	mark_duplicates_pbmarkdup,
	align_to_reference_minimap,
	align_to_reference_vg,
	reassign_mapq,
	filter_reads,
	call_variants_bcftools,
	call_variants_deepvariant,
	call_variants_clair3,
	call_structural_variants_sawfish,
	genotype_tandem_repeats,
	phase_genotypes_hiphase,
	merge_hiphase_vcfs
)
from config import min_reads_sample

def preprocess_pacbio_sample(config):
	# run_fastqc(
	# 	input_file=config['raw_fastq']
	# )
	# )
	
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

	# run_fastqc(
	# 	input_file=config['trimmed_fastq']
	# )
	
	# mark_duplicates_pbmarkdup(
	# 	input_file=config['trimmed_fastq'],
	# 	output_file=config['trimmed_pbmarkdup_fastq'],
	# 	threads=config['threads']
	# )
	
	# align_to_reference_minimap(
	# 	input_file=config['trimmed_pbmarkdup_fastq_gz'],
	# 	output_file=config['hg38_bam'],
	# 	read_group_string=config['read_group_string'],
	# 	reference_fasta=config['reference_genome'],
	# 	platform=config['platform'],
	# 	threads=config['threads'],
	# )
	
	if config['aligner'] == "vg":
		align_to_reference_vg(
			vg=config['vg'],
			input_file=config['trimmed_pbmarkdup_fastq_gz'],
			output_file=config['pangenome_bam'],
			reheader_bam=config['pg_reheader_bam'],
			sample_ID=config['sample_ID'],
			read_group_string=config['read_group_string'],
			reference_gbz=config['reference_gbz'],
			ref_paths=config['ref_paths'],
			platform=config['platform'],
			threads=config['threads']
		)
		
		reassign_mapq(
			bam_hg38=config['hg38_bam'],
			bam_pg=config['pg_reheader_bam'],
			reassigned_pg=config['pg_mapq_reassign_bam']
		)
	
		chr6_read_count = filter_reads(
			input_file=config['pg_mapq_reassign_bam'],
			output_file=config['hg38_rmdup_chr6_bam'],
			threads=config['threads']
		)
	
	# else:
	# 	chr6_read_count = filter_reads(
	# 		input_file=config['hg38_bam'],
	# 		output_file=config['hg38_rmdup_chr6_bam'],
	# 		threads=config['threads']
	# 	)

	if chr6_read_count >= min_reads_sample:
		if config['genotyper'] == "bcftools":
			call_variants_bcftools(
				input_file=config['hg38_rmdup_chr6_bam'],
				output_file=config['snv_vcf'],
				reference_fasta=config['reference_genome'],
				threads=config['threads'],
				platform=config['platform']
			)
		
		elif config['genotyper'] == "deepvariant":
			call_variants_deepvariant(
				input_bam=config['hg38_rmdup_chr6_bam'],
				output_vcf=config['snv_vcf'],
				output_gvcf=config['snv_gvcf'],
				platform=config['platform'],
				deepvariant_sif=config['deepvariant_sif'],
				reference_fasta=config['reference_genome'],
				genotypes_dir=config['genotypes_dir'],
				mapped_bam_dir=config['mapped_bam_dir'],
				sample_ID=config['sample_ID']
			)
		
		elif config['genotyper'] == "clair3":
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
		
		# call_structural_variants_pbsv(sample)
		
		call_structural_variants_sawfish(
			input_bam=config['hg38_rmdup_chr6_bam'],
			small_variant_calls=config['snv_vcf'],
			output_vcf=config['sv_vcf'],
			sv_dir=config['sv_dir'],
			sawfish=config['sawfish'],
			reference_fasta=config['reference_genome']
		)
		
		genotype_tandem_repeats(
			input_bam=config['hg38_rmdup_chr6_bam'],
			output_vcf=config['tr_vcf'],
			pbtrgt_dir=config['pbtrgt_dir'],
			threads=config['threads'],
			reference_fasta=config['reference_genome'],
			pbtrgt_repeat_file=config['pbtrgt_repeat_file'],
			original_cwd=config['ORIGINAL_CWD']
		)
		
		phase_genotypes_hiphase(
			input_bam=config['hg38_rmdup_chr6_bam'],
			input_snv=config['snv_vcf'],
			input_SV=config['sv_vcf'],
			input_TR=config['tr_vcf'],
			output_bam=config['hg38_rmdup_chr6_haplotag_bam'],
			output_snv=config['hiphase_snv_vcf'],
			output_SV=config['hiphase_sv_vcf'],
			output_TR=config['hiphase_tr_vcf'],
			output_summary_file=config['phased_summary'],
			output_blocks_file=config['phased_blocks'],
			output_stats_file=config['phased_stats'],
			threads=config['threads'],
			reference_fasta=config['reference_genome'],
			phased_vcf_dir=config['phased_vcf_dir'],
			sample_ID=config['sample_ID']
		)
		
		merge_hiphase_vcfs(
			input_snv=config['hiphase_snv_vcf'],
			input_SV=config['hiphase_sv_vcf'],
			input_TR=config['hiphase_tr_vcf'],
			output_vcf=config['hiphase_joint_vcf'],
			reference_fasta=config['reference_genome']
		)
	
	else:
		print("Insufficient reads for variant calling")
		print("Sample {} had {} reads!".format(config['sample_ID'], chr6_read_count))