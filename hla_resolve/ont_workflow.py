import os
from preprocess_methods import (
	trim_adapters,
	align_to_reference_minimap,
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

def preprocess_ont_sample(
	config,
	reference_fasta,
	vg,
	reference_gbz,
	ref_paths,
	deepvariant_sif,
	chr6_bed,
	clair3_ont_model_path,
	clair3_hifi_model_path,
	longphase,
	tandem_repeat_bed
):
	trim_adapters(
	adapters=config['adapters'],
	input_file=os.path.join(config['fastq_raw_dir'], config['sample_ID'] + ".fastq.gz"),
	output_file=os.path.join(config['fastq_trimmed_dir'], config['sample_ID'] + ".trimmed.fastq.gz"),
	sample_ID=config['sample_ID'],
	threads=config['threads'],
	adapter_file=config['adapter_file']
	)
	
	align_to_reference_minimap(
		input_file=os.path.join(config['fastq_trimmed_dir'], config['sample_ID'] + ".trimmed.fastq.gz"),
		output_file=os.path.join(config['mapped_bam_dir'], config['sample_ID'] + ".hg38.bam"),
		read_group_string=config['read_group_string'],
		reference_fasta=reference_fasta,
		platform=config['platform'],
		threads=config['threads']
	)
	
	if config['aligner'] == "vg":
		align_to_reference_vg(
			vg=vg,
			input_file=os.path.join(config['fastq_trimmed_dir'], config['sample_ID'] + ".trimmed.fastq.gz"),
			output_file=os.path.join(config['mapped_bam_dir'], config['sample_ID'] + ".pangenome.bam"),  
			sample_ID=config['sample_ID'],
			read_group_string=config['read_group_string'],
			reference_gbz=reference_gbz,
			ref_paths=ref_paths,
			platform=config['platform'],
			threads=config['threads']
			)
	
		reassign_mapq(
			bam_hg38=os.path.join(config['mapped_bam_dir'], config['sample_ID'] + ".hg38.bam"),
			bam_pg=os.path.join(config['mapped_bam_dir'], config['sample_ID'] + ".pangenome.bam"),
			reassigned_pg=os.path.join(config['mapped_bam_dir'], config['sample_ID'] + ".pg.mapq_reassign.bam")
		)
		
		mark_duplicates_picard(
			input_file=os.path.join(config['mapped_bam_dir'], config['sample_ID'] + ".pg.mapq_reassign.bam"),
			output_file=os.path.join(config['mapped_bam_dir'], config['sample_ID'] + ".pg.mapq_reassign.mrkdup.bam"),
			metrics_file=os.path.join(config['mapped_bam_dir'], config['sample_ID'] + ".pg.mapq_reassign.mrkdup.metrics.txt"),
			temp_dir=os.path.join(config['mapped_bam_dir'], "mark_duplicates")
		)

		filter_reads(
			input_file=os.path.join(config['mapped_bam_dir'], config['sample_ID'] + ".pg.mapq_reassign.mrkdup.bam"),
			output_file=os.path.join(config['mapped_bam_dir'], config['sample_ID'] + ".hg38.rmdup.chr6.bam"),
			threads=config['threads']
		)

	else:
		mark_duplicates_picard(
			input_file=os.path.join(config['mapped_bam_dir'], config['sample_ID'] + ".hg38.bam"),
			output_file=os.path.join(config['mapped_bam_dir'], config['sample_ID'] + ".hg38.mrkdup.bam"),
			metrics_file=os.path.join(config['mapped_bam_dir'], config['sample_ID'] + ".hg38.mrkdup.metrics.txt"),
			temp_dir=os.path.join(config['mapped_bam_dir'], "mark_duplicates")
		)

		filter_reads(
			input_file=os.path.join(config['mapped_bam_dir'], config['sample_ID'] + ".hg38.mrkdup.bam"),
			output_file=os.path.join(config['mapped_bam_dir'], config['sample_ID'] + ".hg38.rmdup.chr6.bam"),
			threads=config['threads']
		)

	if config['genotyper'] == "bcftools":
		call_variants_bcftools(
			input_file=os.path.join(config['mapped_bam_dir'], config['sample_ID'] +".hg38.rmdup.chr6.bam"),
			output_file=os.path.join(config['genotypes_dir'], config['sample_ID'] + ".vcf.gz"),
			reference_fasta=reference_fasta,
			platform=config['platform'],
			threads=config['threads']
		)

	elif config['genotyper'] == "deepvariant":
		call_variants_deepvariant(
			input_bam=os.path.join(config['mapped_bam_dir'], config['sample_ID'] + ".hg38.rmdup.chr6.bam"),
			output_vcf=os.path.join(config['genotypes_dir'], config['sample_ID'] + ".vcf.gz"),
			output_gvcf=os.path.join(config['genotypes_dir'], config['sample_ID'] + ".g.vcf.gz"),
			platform=config['platform'],
			deepvariant_sif=deepvariant_sif,
			reference_fasta=reference_fasta,
			genotypes_dir=config['genotypes_dir'],
			mapped_bam_dir=config['mapped_bam_dir'],
			sample_ID=config['sample_ID']
		)
	elif config['genotyper'] == "clair3":
		call_variants_clair3(
			input_bam=os.path.join(config['mapped_bam_dir'], config['sample_ID'] + ".hg38.rmdup.chr6.bam"),
			output_vcf=os.path.join(config['genotypes_dir'], config['sample_ID'] + ".vcf.gz"),
			platform=config['platform'],
			reference_fasta=reference_fasta,
			threads=config['threads'],
			chr6_bed=chr6_bed,
			clair3_ont_model_path=clair3_ont_model_path,
			clair3_hifi_model_path=clair3_hifi_model_path,
			genotypes_dir=config['genotypes_dir'],
			sample_ID=config['sample_ID']
		)
	call_structural_variants_sniffles(
		input_bam=os.path.join(config['mapped_bam_dir'], config['sample_ID'] + ".hg38.rmdup.chr6.bam"),
		output_vcf=os.path.join(config['sv_dir'], config['sample_ID'] + ".SV.vcf.gz"),
		threads=config['threads'],
		reference_fasta=reference_fasta,
		chr6_bed=chr6_bed,
		tandem_repeat_bed=tandem_repeat_bed
	)
	phase_genotypes_longphase(
		input_bam=os.path.join(config['mapped_bam_dir'], config['sample_ID'] + ".hg38.rmdup.chr6.bam"),
		input_SNV_vcf=os.path.join(config['genotypes_dir'], config['sample_ID'] + ".vcf.gz"),
		input_SV_vcf=os.path.join(config['sv_dir'], config['sample_ID'] + ".SV.vcf.gz"),
		output_blocks_file=os.path.join(config['phased_vcf_dir'], config['sample_ID'] + ".phased.haploblocks.txt"),
		output_gtf_file=os.path.join(config['phased_vcf_dir'], config['sample_ID'] + ".phased.haploblocks.gtf"),
		phased_vcf=os.path.join(config['phased_vcf_dir'], config['sample_ID'] + ".longphase.vcf.gz"),
		phased_SV_vcf=os.path.join(config['phased_vcf_dir'], config['sample_ID'] + ".longphase_SV.vcf.gz"),
		haplotagged_bam=os.path.join(config['mapped_bam_dir'], config['sample_ID'] + ".hg38.rmdup.chr6.haplotag.bam"),
		longphase=longphase,
		reference_fasta=reference_fasta,
		threads=config['threads'],
		phased_vcf_dir=config['phased_vcf_dir'],
		sample_ID=config['sample_ID']
	)
	merge_longphase_vcfs(
		phased_vcf=os.path.join(config['phased_vcf_dir'], config['sample_ID'] + ".longphase.vcf.gz"),
		phased_SV_vcf=os.path.join(config['phased_vcf_dir'], config['sample_ID'] + ".longphase_SV.vcf.gz"),
		merged_vcf=os.path.join(config['phased_vcf_dir'], config['sample_ID'] + ".longphase.merged.vcf.gz"),
		reference_fasta=reference_fasta,
		phased_vcf_dir=config['phased_vcf_dir'],
		sample_ID=config['sample_ID']
	)