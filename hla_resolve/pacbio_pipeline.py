# This software is Copyright ©2026. The Regents of the University of California
# ("Regents"). All Rights Reserved.
#
# Licensed under the UC Santa Cruz Noncommercial License (the "License").
# You may not use this file except in compliance with the License.
# A copy of the License is included in this repository as LICENSE.txt.

import os
import subprocess
from .preprocess_methods import (
	trim_adapters,
	mark_duplicates_pbmarkdup,
	align_to_reference_minimap,
	mark_duplicates_picard,
	filter_reads,
	classify_DRB_reads,
	call_variants_bcftools,
	call_variants_deepvariant,
	call_variants_clair3,
	call_variants_freebayes,
	merge_hybrid_vcfs,
	rescue_refcalls,
	call_structural_variants_pbsv,
	genotype_tandem_repeats,
	phase_genotypes_hiphase,
	merge_hiphase_vcfs
)
from .config import min_reads_sample

def preprocess_pacbio_sample(config):
	trim_adapters(
		adapters=config['adapters'],
		input_file=config['raw_fastq'],
		output_file=config['trimmed_fastq'],
		sample_ID=config['sample_ID'],
		threads=config['threads'],
		adapter_file=config['adapter_file'],
		five_prime_adapter=config['five_prime_adapter'],
		three_prime_adapter=config['three_prime_adapter']
	)

	if config['scheme'] == "targeted":
		mark_duplicates_pbmarkdup(
			input_file=config['trimmed_fastq'],
			output_file=config['trimmed_pbmarkdup_fastq'],
			threads=config['threads']
		)
		
		align_to_reference_minimap(
			input_file=config['trimmed_pbmarkdup_fastq_gz'],
			output_file=config['hg38_bam'],
			read_group_string=config['read_group_string'],
			reference_fasta=config['reference_genome'],
			platform=config['platform'],
			threads=config['threads'],
		)

		classify_DRB_reads(
			input_file=config['trimmed_pbmarkdup_fastq_gz'],
			output_file=config['hg38_bam_drb'],
			DRB34_reads_file=config['DRB34_reads_file'],
			read_group_string=config['read_group_string'],
			reference_fasta=config['drb_multiallele_reference'],
			platform=config['platform'],
			threads=config['threads']
		)

		chr6_read_count = filter_reads(
			input_file=config['hg38_bam'],
			output_file=config['hg38_rmdup_chr6_bam'],
			DRB34_reads_file=config['DRB34_reads_file'],
			threads=config['threads']
		)
	
	elif config['scheme'] == "WGS" or config['scheme'] == "WES":
		align_to_reference_minimap(
			input_file=config['trimmed_fastq'],
			output_file=config['hg38_bam'],
			read_group_string=config['read_group_string'],
			reference_fasta=config['reference_genome'],
			platform=config['platform'],
			threads=config['threads'],
		)

		classify_DRB_reads(
			input_file=config['trimmed_fastq'],
			output_file=config['hg38_bam_drb'],
			DRB34_reads_file=config['DRB34_reads_file'],
			read_group_string=config['read_group_string'],
			reference_fasta=config['drb_multiallele_reference'],
			platform=config['platform'],
			threads=config['threads']
		)

		filter_reads(
			input_file=config['hg38_bam'],
			output_file=config['hg38_chr6_bam'],
			DRB34_reads_file=config['DRB34_reads_file'],
			threads=config['threads']
		)
		
		mark_duplicates_picard(
			input_file=config['hg38_chr6_bam'],
			output_file=config['hg38_rmdup_chr6_bam'],
			metrics_file=config['hg38_mrkdup_metrics'],
			temp_dir=os.path.join(config['mapped_bam_dir'], "mark_duplicates"),
			picard=config['picard']
		)

		chr6_read_count = int(subprocess.check_output(f"samtools view -c {config['hg38_rmdup_chr6_bam']}", shell=True).strip())

	if chr6_read_count >= min_reads_sample:
		snp_caller = config['snp_caller']
		indel_caller = config['indel_caller']

		# Map each caller to its intermediate output file (used in hybrid mode)
		caller_vcf_map = {
			"bcftools":    config['bcftools_snp_vcf'],
			"deepvariant": config['dv_full_vcf'],
			"clair3":      config['clair3_full_vcf'],
		}

		if snp_caller == indel_caller:
			# Single caller for both SNPs and indels
			if snp_caller == "bcftools":
				call_variants_bcftools(
					input_file=config['hg38_rmdup_chr6_bam'],
					output_file=config['snv_vcf'],
					reference_fasta=config['reference_genome'],
					threads=config['threads'],
					platform=config['platform']
				)
			elif snp_caller == "deepvariant":
				dv_output = config['dv_full_vcf'] if config['rescue_refcalls'] else config['snv_vcf']
				call_variants_deepvariant(
					input_bam=config['hg38_rmdup_chr6_bam'],
					output_vcf=dv_output,
					output_gvcf=config['snv_gvcf'],
					platform=config['platform'],
					deepvariant_sif=config['deepvariant_sif'],
					reference_fasta=config['reference_genome'],
					genotypes_dir=config['genotypes_dir'],
					mapped_bam_dir=config['mapped_bam_dir'],
					sample_ID=config['sample_ID'],
					threads=config['threads']
				)
				if config['rescue_refcalls']:
					rescue_refcalls(
						input_vcf=config['dv_full_vcf'],
						output_vcf=config['snv_vcf']
					)
			elif snp_caller == "clair3":
				call_variants_clair3(
					input_bam=config['hg38_rmdup_chr6_bam'],
					output_vcf=config['snv_vcf'],
					platform=config['platform'],
					clair3_sif=config['clair3_sif'],
					reference_fasta=config['reference_genome'],
					threads=config['threads'],
					genotypes_dir=config['genotypes_dir'],
					mapped_bam_dir=config['mapped_bam_dir'],
					sample_ID=config['sample_ID'],
					clair3_model=config['clair3_model']
				)
			elif snp_caller == "freebayes":
				call_variants_freebayes(
					input_bam=config['hg38_rmdup_chr6_bam'],
					output_vcf=config['snv_vcf'],
					reference_fasta=config['reference_genome']
				)
		else:
			# Hybrid: different callers for SNPs and indels
			snp_intermediate = caller_vcf_map[snp_caller]
			indel_intermediate = caller_vcf_map[indel_caller]

			if snp_caller == "bcftools":
				call_variants_bcftools(
					input_file=config['hg38_rmdup_chr6_bam'],
					output_file=snp_intermediate,
					reference_fasta=config['reference_genome'],
					threads=config['threads'],
					platform=config['platform']
				)
			elif snp_caller == "deepvariant":
				call_variants_deepvariant(
					input_bam=config['hg38_rmdup_chr6_bam'],
					output_vcf=snp_intermediate,
					output_gvcf=config['snv_gvcf'],
					platform=config['platform'],
					deepvariant_sif=config['deepvariant_sif'],
					reference_fasta=config['reference_genome'],
					genotypes_dir=config['genotypes_dir'],
					mapped_bam_dir=config['mapped_bam_dir'],
					sample_ID=config['sample_ID'],
					threads=config['threads']
				)
			elif snp_caller == "clair3":
				call_variants_clair3(
					input_bam=config['hg38_rmdup_chr6_bam'],
					output_vcf=snp_intermediate,
					platform=config['platform'],
					clair3_sif=config['clair3_sif'],
					reference_fasta=config['reference_genome'],
					threads=config['threads'],
					genotypes_dir=config['genotypes_dir'],
					mapped_bam_dir=config['mapped_bam_dir'],
					sample_ID=config['sample_ID'],
					clair3_model=config['clair3_model']
				)

			if indel_caller == "bcftools":
				call_variants_bcftools(
					input_file=config['hg38_rmdup_chr6_bam'],
					output_file=indel_intermediate,
					reference_fasta=config['reference_genome'],
					threads=config['threads'],
					platform=config['platform']
				)
			elif indel_caller == "deepvariant":
				call_variants_deepvariant(
					input_bam=config['hg38_rmdup_chr6_bam'],
					output_vcf=indel_intermediate,
					output_gvcf=config['snv_gvcf'],
					platform=config['platform'],
					deepvariant_sif=config['deepvariant_sif'],
					reference_fasta=config['reference_genome'],
					genotypes_dir=config['genotypes_dir'],
					mapped_bam_dir=config['mapped_bam_dir'],
					sample_ID=config['sample_ID'],
					threads=config['threads']
				)
			elif indel_caller == "clair3":
				call_variants_clair3(
					input_bam=config['hg38_rmdup_chr6_bam'],
					output_vcf=indel_intermediate,
					platform=config['platform'],
					clair3_sif=config['clair3_sif'],
					reference_fasta=config['reference_genome'],
					threads=config['threads'],
					genotypes_dir=config['genotypes_dir'],
					mapped_bam_dir=config['mapped_bam_dir'],
					sample_ID=config['sample_ID'],
					clair3_model=config['clair3_model']
				)

			if config['rescue_refcalls'] and indel_caller == "deepvariant":
				rescue_refcalls(
					input_vcf=indel_intermediate,
					output_vcf=config['dv_rescued_vcf'],
					indels_only=True
				)
				indel_intermediate = config['dv_rescued_vcf']

			merge_hybrid_vcfs(
				snp_vcf=snp_intermediate,
				indel_vcf=indel_intermediate,
				indel_only_vcf=config['hybrid_indel_vcf'],
				merged_vcf=config['snv_vcf'],
				filter_indel_pass=indel_caller in ("deepvariant", "clair3")
			)

		call_structural_variants_pbsv(
			input_bam=config['hg38_rmdup_chr6_bam'],
			output_svsig=config['sv_svsig'],
			output_vcf=config['sv_vcf'].replace('.vcf.gz', '.vcf'),
			threads=config['threads'],
			tandem_repeat_bed=config['tandem_repeat_bed'],
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