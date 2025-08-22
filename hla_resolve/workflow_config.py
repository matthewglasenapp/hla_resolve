def build_workflow_config(sample):
    """
    Build a configuration dictionary for workflows from a sample object.
    This avoids passing the entire sample object to workflow functions.
    """
    config = {
        # Core sample information
        'sample_ID': sample.sample_ID,
        'threads': sample.threads,
        'platform': sample.platform,
        'aligner': sample.aligner,
        'genotyper': sample.genotyper,
        
        # Adapter settings
        'adapters': sample.adapters,
        'adapter_file': sample.adapter_file,
        'read_group_string': sample.read_group_string,
        
        # Directory paths
        'fastq_raw_dir': sample.fastq_raw_dir,
        'fastq_trimmed_dir': sample.fastq_trimmed_dir,
        'mapped_bam_dir': sample.mapped_bam_dir,
        'genotypes_dir': sample.genotypes_dir,
        'sv_dir': sample.sv_dir,
        'phased_vcf_dir': sample.phased_vcf_dir,
        'mosdepth_dir': sample.mosdepth_dir,
        'parsed_haploblock_dir': sample.parsed_haploblock_dir,
        'filtered_vcf_dir': sample.filtered_vcf_dir,
        'vcf2fasta_out_dir': sample.vcf2fasta_out_dir,
        'hla_fasta_dir': sample.hla_fasta_dir,
        'hla_typing_dir': sample.hla_typing_dir,
        'gff_dir': sample.gff_dir,
        'clean_up': sample.clean_up,
    }
    
    # Platform-specific directories
    if sample.platform == "PACBIO":
        config['pbtrgt_dir'] = getattr(sample, 'pbtrgt_dir', None)
        config['ORIGINAL_CWD'] = getattr(sample, 'ORIGINAL_CWD', None)
    
    return config
