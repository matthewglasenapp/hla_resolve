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
	sample,
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
		adapters=sample.adapters,
		input_file=os.path.join(sample.fastq_raw_dir, sample.sample_ID + ".fastq.gz"),
		output_file=os.path.join(sample.fastq_trimmed_dir, sample.sample_ID + ".trimmed.fastq.gz"),
		sample_ID=sample.sample_ID,
		threads=sample.threads,
		adapter_file=sample.adapter_file
	)

	run_fastqc(
		input_file=os.path.join(sample.fastq_trimmed_dir, sample.sample_ID + ".trimmed.fastq.gz")
	)
	
	mark_duplicates_pbmarkdup(
		input_file=os.path.join(sample.fastq_trimmed_dir, sample.sample_ID + ".trimmed.fastq.gz"),
		output_file=os.path.join(sample.fastq_trimmed_dir, sample.sample_ID + ".trimmed.pbmarkdup.fastq"),
		threads=sample.threads
	)
	
	run_fastqc(
		input_file=os.path.join(sample.fastq_trimmed_dir, sample.sample_ID + ".trimmed.pbmarkdup.fastq.gz")
	)
	
	align_to_reference_minimap(
		input_file=os.path.join(sample.fastq_trimmed_dir, sample.sample_ID + ".trimmed.pbmarkdup.fastq.gz"),
		output_file=os.path.join(sample.mapped_bam_dir, sample.sample_ID + ".hg38.bam"),
		read_group_string=sample.read_group_string,
		reference_fasta=reference_fasta,
		platform=sample.platform,
		threads=sample.threads,
	)
	
	if sample.aligner == "vg":
		align_to_reference_vg(
			vg=vg,
			input_file=os.path.join(sample.fastq_trimmed_dir, sample.sample_ID + ".trimmed.pbmarkdup.fastq.gz"),
			output_file=os.path.join(sample.mapped_bam_dir, sample.sample_ID + ".pangenome.bam"),
			sample_ID=sample.sample_ID,
			read_group_string=sample.read_group_string,
			reference_gbz=reference_gbz,
			ref_paths=ref_paths,
			platform=sample.platform,
			threads=sample.threads
		)
		
		reassign_mapq(
			bam_hg38=os.path.join(sample.mapped_bam_dir, sample.sample_ID + ".hg38.bam"),
			bam_pg=os.path.join(sample.mapped_bam_dir, sample.sample_ID + ".pangenome.bam"),
			reassigned_pg=os.path.join(sample.mapped_bam_dir, sample.sample_ID + ".pg.mapq_reassign.bam")
		)
	
		filter_reads(
			input_file=os.path.join(sample.mapped_bam_dir, sample.sample_ID + ".pg.mapq_reassign.bam"),
			output_file=os.path.join(sample.mapped_bam_dir, sample.sample_ID + ".hg38.rmdup.chr6.bam"),
			threads=sample.threads
		)
	
	else:
		filter_reads(
			input_file=os.path.join(sample.mapped_bam_dir, sample.sample_ID + ".hg38.bam"),
			output_file=os.path.join(sample.mapped_bam_dir, sample.sample_ID + ".hg38.rmdup.chr6.bam"),
			threads=sample.threads
		)

	if sample.genotyper == "bcftools":
		call_variants_bcftools(
			input_file=os.path.join(sample.mapped_bam_dir, sample.sample_ID +".hg38.rmdup.chr6.bam"),
			output_file=os.path.join(sample.genotypes_dir, sample.sample_ID + ".vcf.gz"),
			reference_fasta=reference_fasta,
			threads=sample.threads,
			platform=sample.platform
		)
	
	elif sample.genotyper == "deepvariant":
		call_variants_deepvariant(
			input_bam=os.path.join(sample.mapped_bam_dir, sample.sample_ID + ".hg38.rmdup.chr6.bam"),
			output_vcf=os.path.join(sample.genotypes_dir, sample.sample_ID + ".vcf.gz"),
			output_gvcf=os.path.join(sample.genotypes_dir, sample.sample_ID + ".g.vcf.gz"),
			platform=sample.platform,
			deepvariant_sif=deepvariant_sif,
			reference_fasta=reference_fasta,
			genotypes_dir=sample.genotypes_dir,
			mapped_bam_dir=sample.mapped_bam_dir,
			sample_ID=sample.sample_ID
		)
	
	elif sample.genotyper == "clair3":
		call_variants_clair3(
			input_bam=os.path.join(sample.mapped_bam_dir, sample.sample_ID + ".hg38.rmdup.chr6.bam"),
			output_vcf=os.path.join(sample.genotypes_dir, sample.sample_ID + ".vcf.gz"),
			platform=sample.platform,
			reference_fasta=reference_fasta,
			threads=sample.threads,
			chr6_bed=chr6_bed,
			clair3_ont_model_path=clair3_ont_model_path,
			clair3_hifi_model_path=clair3_hifi_model_path,
			genotypes_dir=sample.genotypes_dir,
			sample_ID=sample.sample_ID
		)
	
	# call_structural_variants_pbsv(sample)
	call_structural_variants_sawfish(
		input_bam=os.path.join(sample.mapped_bam_dir, sample.sample_ID + ".hg38.rmdup.chr6.bam"),
		small_variant_calls=os.path.join(sample.genotypes_dir, sample.sample_ID + ".vcf.gz"),
		output_vcf=os.path.join(sample.sv_dir, sample.sample_ID + ".SV.vcf.gz"),
		sv_dir=sample.sv_dir,
		sawfish=sawfish,
		reference_fasta=reference_fasta
	)
	
	genotype_tandem_repeats(
		input_bam=os.path.join(sample.mapped_bam_dir, sample.sample_ID + ".hg38.rmdup.chr6.bam"),
		output_vcf=os.path.join(sample.pbtrgt_dir, sample.sample_ID + ".TR.vcf.gz"),
		pbtrgt_dir=sample.pbtrgt_dir,
		threads=sample.threads,
		reference_fasta=reference_fasta,
		pbtrgt_repeat_file=pbtrgt_repeat_file,
		original_cwd=sample.ORIGINAL_CWD
	)
	
	phase_genotypes_hiphase(
		input_bam=os.path.join(sample.mapped_bam_dir, sample.sample_ID + ".hg38.rmdup.chr6.bam"),
		input_snv=os.path.join(sample.genotypes_dir, sample.sample_ID + ".vcf.gz"),
		input_SV=os.path.join(sample.sv_dir, sample.sample_ID + ".SV.vcf.gz"),
		input_TR=os.path.join(sample.pbtrgt_dir, sample.sample_ID + ".TR.vcf.gz"),
		output_bam=os.path.join(sample.mapped_bam_dir, sample.sample_ID + ".hg38.rmdup.chr6.haplotag.bam"),
		output_snv=os.path.join(sample.phased_vcf_dir, sample.sample_ID + ".hiphase.vcf.gz"),
		output_SV=os.path.join(sample.phased_vcf_dir, sample.sample_ID + ".hiphase.SV.vcf.gz"),
		output_TR=os.path.join(sample.phased_vcf_dir, sample.sample_ID + ".hiphase.TR.vcf.gz"),
		output_summary_file=os.path.join(sample.phased_vcf_dir, sample.sample_ID + ".phased.summary.txt"),
		output_blocks_file=os.path.join(sample.phased_vcf_dir, sample.sample_ID + ".phased.blocks.txt"),
		threads=sample.threads,
		reference_fasta=reference_fasta,
		phased_vcf_dir=sample.phased_vcf_dir,
		sample_ID=sample.sample_ID
	)
	
	merge_hiphase_vcfs(
		input_snv=os.path.join(sample.phased_vcf_dir, sample.sample_ID + ".hiphase.vcf.gz"),
		input_SV=os.path.join(sample.phased_vcf_dir, sample.sample_ID + ".hiphase.SV.vcf.gz"),
		input_TR=os.path.join(sample.phased_vcf_dir, sample.sample_ID + ".hiphase.TR.vcf.gz"),
		output_vcf=os.path.join(sample.phased_vcf_dir, sample.sample_ID + ".hiphase.joint.vcf.gz"),
		reference_fasta=reference_fasta
	)