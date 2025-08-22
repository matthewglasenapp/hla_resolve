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

def preprocess_pacbio_sample(
	config,
	reference_fasta,
	vg,
	reference_gbz,
	ref_paths,
	deepvariant_sif,
	chr6_bed,
	clair3_ont_model_path,
	clair3_hifi_model_path,
	sawfish,
	pbtrgt_repeat_file
):
	trim_adapters(
		adapters=config['adapters'],
		input_file=os.path.join(config['fastq_raw_dir'], config['sample_ID'] + ".fastq.gz"),
		output_file=os.path.join(config['fastq_trimmed_dir'], config['sample_ID'] + ".trimmed.fastq.gz"),
		sample_ID=config['sample_ID'],
		threads=config['threads'],
		adapter_file=config['adapter_file']
	)

	run_fastqc(
		input_file=os.path.join(config['fastq_trimmed_dir'], config['sample_ID'] + ".trimmed.fastq.gz")
	)
	
	mark_duplicates_pbmarkdup(
		input_file=os.path.join(config['fastq_trimmed_dir'], config['sample_ID'] + ".trimmed.fastq.gz"),
		output_file=os.path.join(config['fastq_trimmed_dir'], config['sample_ID'] + ".trimmed.pbmarkdup.fastq"),
		threads=config['threads']
	)
	
	run_fastqc(
		input_file=os.path.join(config['fastq_trimmed_dir'], config['sample_ID'] + ".trimmed.pbmarkdup.fastq.gz")
	)
	
	align_to_reference_minimap(
		input_file=os.path.join(config['fastq_trimmed_dir'], config['sample_ID'] + ".trimmed.pbmarkdup.fastq.gz"),
		output_file=os.path.join(config['mapped_bam_dir'], config['sample_ID'] + ".hg38.bam"),
		read_group_string=config['read_group_string'],
		reference_fasta=reference_fasta,
		platform=config['platform'],
		threads=config['threads'],
	)
	
	if config['aligner'] == "vg":
		align_to_reference_vg(
			vg=vg,
			input_file=os.path.join(config['fastq_trimmed_dir'], config['sample_ID'] + ".trimmed.pbmarkdup.fastq.gz"),
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
	
		filter_reads(
			input_file=os.path.join(config['mapped_bam_dir'], config['sample_ID'] + ".pg.mapq_reassign.bam"),
			output_file=os.path.join(config['mapped_bam_dir'], config['sample_ID'] + ".hg38.rmdup.chr6.bam"),
			threads=config['threads']
		)
	
	else:
		filter_reads(
			input_file=os.path.join(config['mapped_bam_dir'], config['sample_ID'] + ".hg38.bam"),
			output_file=os.path.join(config['mapped_bam_dir'], config['sample_ID'] + ".hg38.rmdup.chr6.bam"),
			threads=config['threads']
		)

	if config['genotyper'] == "bcftools":
		call_variants_bcftools(
			input_file=os.path.join(config['mapped_bam_dir'], config['sample_ID'] +".hg38.rmdup.chr6.bam"),
			output_file=os.path.join(config['genotypes_dir'], config['sample_ID'] + ".vcf.gz"),
			reference_fasta=reference_fasta,
			threads=config['threads'],
			platform=config['platform']
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
	
	# call_structural_variants_pbsv(sample)
	call_structural_variants_sawfish(
		input_bam=os.path.join(config['mapped_bam_dir'], config['sample_ID'] + ".hg38.rmdup.chr6.bam"),
		small_variant_calls=os.path.join(config['genotypes_dir'], config['sample_ID'] + ".vcf.gz"),
		output_vcf=os.path.join(config['sv_dir'], config['sample_ID'] + ".SV.vcf.gz"),
		sv_dir=config['sv_dir'],
		sawfish=sawfish,
		reference_fasta=reference_fasta
	)
	
	genotype_tandem_repeats(
		input_bam=os.path.join(config['mapped_bam_dir'], config['sample_ID'] + ".hg38.rmdup.chr6.bam"),
		output_vcf=os.path.join(config['pbtrgt_dir'], config['sample_ID'] + ".TR.vcf.gz"),
		pbtrgt_dir=config['pbtrgt_dir'],
		threads=config['threads'],
		reference_fasta=reference_fasta,
		pbtrgt_repeat_file=pbtrgt_repeat_file,
		original_cwd=config['ORIGINAL_CWD']
	)
	
	phase_genotypes_hiphase(
		input_bam=os.path.join(config['mapped_bam_dir'], config['sample_ID'] + ".hg38.rmdup.chr6.bam"),
		input_snv=os.path.join(config['genotypes_dir'], config['sample_ID'] + ".vcf.gz"),
		input_SV=os.path.join(config['sv_dir'], config['sample_ID'] + ".SV.vcf.gz"),
		input_TR=os.path.join(config['pbtrgt_dir'], config['sample_ID'] + ".TR.vcf.gz"),
		output_bam=os.path.join(config['mapped_bam_dir'], config['sample_ID'] + ".hg38.rmdup.chr6.haplotag.bam"),
		output_snv=os.path.join(config['phased_vcf_dir'], config['sample_ID'] + ".hiphase.vcf.gz"),
		output_SV=os.path.join(config['phased_vcf_dir'], config['sample_ID'] + ".hiphase.SV.vcf.gz"),
		output_TR=os.path.join(config['phased_vcf_dir'], config['sample_ID'] + ".hiphase.TR.vcf.gz"),
		output_summary_file=os.path.join(config['phased_vcf_dir'], config['sample_ID'] + ".phased.summary.txt"),
		output_blocks_file=os.path.join(config['phased_vcf_dir'], config['sample_ID'] + ".phased.blocks.txt"),
		threads=config['threads'],
		reference_fasta=reference_fasta,
		phased_vcf_dir=config['phased_vcf_dir'],
		sample_ID=config['sample_ID']
	)
	
	merge_hiphase_vcfs(
		input_snv=os.path.join(config['phased_vcf_dir'], config['sample_ID'] + ".hiphase.vcf.gz"),
		input_SV=os.path.join(config['phased_vcf_dir'], config['sample_ID'] + ".hiphase.SV.vcf.gz"),
		input_TR=os.path.join(config['phased_vcf_dir'], config['sample_ID'] + ".hiphase.TR.vcf.gz"),
		output_vcf=os.path.join(config['phased_vcf_dir'], config['sample_ID'] + ".hiphase.joint.vcf.gz"),
		reference_fasta=reference_fasta
	)