def preprocess_ont_sample(sample):
    trim_adapters(
    adapters=sample.adapters,
    input_file=os.path.join(sample.fastq_raw_dir, sample.sample_ID + ".fastq.gz"),
    output_file=os.path.join(sample.fastq_trimmed_dir, sample.sample_ID + ".trimmed.fastq.gz"),
    sample_ID=sample.sample_ID,
    threads=sample.threads,
    adapter_file=sample.adapter_file
    )
    
    align_to_reference_minimap(
        input_file=os.path.join(sample.fastq_trimmed_dir, sample.sample_ID + ".trimmed.fastq.gz"),
        output_file=os.path.join(sample.mapped_bam_dir, sample.sample_ID + ".hg38.bam"),
        read_group_string=sample.read_group_string,
        reference_fasta=Samples.reference_fasta,
        platform=sample.platform,
        threads=sample.threads
    )
    
    if sample.aligner == "vg":
        align_to_reference_vg(
            vg=Samples.vg,
            input_file=os.path.join(sample.fastq_trimmed_dir, sample.sample_ID + ".trimmed.fastq.gz"),
            output_file=os.path.join(sample.mapped_bam_dir, sample.sample_ID + ".pangenome.bam"),  
            sample_ID=sample.sample_ID,
            read_group_string=sample.read_group_string,
            reference_gbz=Samples.reference_gbz,
            ref_paths=Samples.ref_paths,
            platform=sample.platform,
            threads=sample.threads
            )
    
        reassign_mapq(
            bam_hg38=os.path.join(sample.mapped_bam_dir, sample.sample_ID + ".hg38.bam"),
            bam_pg=os.path.join(sample.mapped_bam_dir, sample.sample_ID + ".pangenome.bam"),
            reassigned_pg=os.path.join(sample.mapped_bam_dir, sample.sample_ID + ".pg.mapq_reassign.bam")
        )
        
        mark_duplicates_picard(
            input_file=os.path.join(sample.mapped_bam_dir, sample.sample_ID + ".pg.mapq_reassign.bam"),
            output_file=os.path.join(sample.mapped_bam_dir, sample.sample_ID + ".pg.mapq_reassign.mrkdup.bam"),
            metrics_file=os.path.join(sample.mapped_bam_dir, sample.sample_ID + ".pg.mapq_reassign.mrkdup.metrics.txt"),
            temp_dir=os.path.join(sample.mapped_bam_dir, "mark_duplicates")
        )

        filter_reads(
            input_file=os.path.join(sample.mapped_bam_dir, sample.sample_ID + ".pg.mapq_reassign.mrkdup.bam"),
            output_file=os.path.join(sample.mapped_bam_dir, sample.sample_ID + ".hg38.rmdup.chr6.bam"),
            threads=sample.threads
        )

    else:
        mark_duplicates_picard(
            input_file=os.path.join(sample.mapped_bam_dir, sample.sample_ID + ".hg38.bam"),
            output_file=os.path.join(sample.mapped_bam_dir, sample.sample_ID + ".hg38.mrkdup.bam"),
            metrics_file=os.path.join(sample.mapped_bam_dir, sample.sample_ID + ".hg38.mrkdup.metrics.txt"),
            temp_dir=os.path.join(sample.mapped_bam_dir, "mark_duplicates")
        )

        filter_reads(
            input_file=os.path.join(sample.mapped_bam_dir, sample.sample_ID + ".hg38.mrkdup.bam"),
            output_file=os.path.join(sample.mapped_bam_dir, sample.sample_ID + ".hg38.rmdup.chr6.bam"),
            threads=sample.threads
        )

    if sample.genotyper == "bcftools":
        call_variants_bcftools(
            input_file=os.path.join(sample.mapped_bam_dir, sample.sample_ID +".hg38.rmdup.chr6.bam"),
            output_file=os.path.join(sample.genotypes_dir, sample.sample_ID + ".vcf.gz"),
            reference_fasta=Samples.reference_fasta,
            platform=sample.platform,
            threads=sample.threads
        )

    elif sample.genotyper == "deepvariant":
        call_variants_deepvariant(
            input_bam=os.path.join(sample.mapped_bam_dir, sample.sample_ID + ".hg38.rmdup.chr6.bam"),
            output_vcf=os.path.join(sample.genotypes_dir, sample.sample_ID + ".vcf.gz"),
            output_gvcf=os.path.join(sample.genotypes_dir, sample.sample_ID + ".g.vcf.gz"),
            platform=sample.platform,
            deepvariant_sif=Samples.deepvariant_sif,
            reference_fasta=Samples.reference_fasta,
            genotypes_dir=sample.genotypes_dir,
            mapped_bam_dir=sample.mapped_bam_dir,
            sample_ID=sample.sample_ID
        )
    elif sample.genotyper == "clair3":
        call_variants_clair3(
            input_bam=os.path.join(sample.mapped_bam_dir, sample.sample_ID + ".hg38.rmdup.chr6.bam"),
            output_vcf=os.path.join(sample.genotypes_dir, sample.sample_ID + ".vcf.gz"),
            platform=sample.platform,
            reference_fasta=Samples.reference_fasta,
            threads=sample.threads,
            chr6_bed=Samples.chr6_bed,
            clair3_ont_model_path=Samples.clair3_ont_model_path,
            clair3_hifi_model_path=Samples.clair3_hifi_model_path,
            genotypes_dir=sample.genotypes_dir,
            sample_ID=sample.sample_ID
        )
    call_structural_variants_sniffles(
        input_bam=os.path.join(sample.mapped_bam_dir, sample.sample_ID + ".hg38.rmdup.chr6.bam"),
        output_vcf=os.path.join(sample.sv_dir, sample.sample_ID + ".SV.vcf.gz"),
        threads=sample.threads,
        reference_fasta=Samples.reference_fasta,
        chr6_bed=Samples.chr6_bed,
        tandem_repeat_bed=Samples.tandem_repeat_bed
    )
    phase_genotypes_longphase(
        input_bam=os.path.join(sample.mapped_bam_dir, sample.sample_ID + ".hg38.rmdup.chr6.bam"),
        input_SNV_vcf=os.path.join(sample.genotypes_dir, sample.sample_ID + ".vcf.gz"),
        input_SV_vcf=os.path.join(sample.sv_dir, sample.sample_ID + ".SV.vcf.gz"),
        output_blocks_file=os.path.join(sample.phased_vcf_dir, sample.sample_ID + ".phased.haploblocks.txt"),
        output_gtf_file=os.path.join(sample.phased_vcf_dir, sample.sample_ID + ".phased.haploblocks.gtf"),
        phased_vcf=os.path.join(sample.phased_vcf_dir, sample.sample_ID + ".longphase.vcf.gz"),
        phased_SV_vcf=os.path.join(sample.phased_vcf_dir, sample.sample_ID + ".longphase_SV.vcf.gz"),
        haplotagged_bam=os.path.join(sample.mapped_bam_dir, sample.sample_ID + ".hg38.rmdup.chr6.haplotag.bam"),
        longphase=Samples.longphase,
        reference_fasta=Samples.reference_fasta,
        threads=sample.threads,
        phased_vcf_dir=sample.phased_vcf_dir,
        sample_ID=sample.sample_ID
    )
    merge_longphase_vcfs(
        phased_vcf=os.path.join(sample.phased_vcf_dir, sample.sample_ID + ".longphase.vcf.gz"),
        phased_SV_vcf=os.path.join(sample.phased_vcf_dir, sample.sample_ID + ".longphase_SV.vcf.gz"),
        merged_vcf=os.path.join(sample.phased_vcf_dir, sample.sample_ID + ".longphase.merged.vcf.gz"),
        reference_fasta=Samples.reference_fasta,
        phased_vcf_dir=sample.phased_vcf_dir,
        sample_ID=sample.sample_ID
    )