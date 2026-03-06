import os
import shutil
import subprocess
import json
import pysam
from Bio import SeqIO
from .preprocess_methods import convert_bam_to_fastq
from .config import (
	min_reads_sample, min_read_length, reference_genome_vg_gbz, reference_genome_vg_paths, vg, 
	longphase, sawfish, clair3_sif, clair3_ont_model, clair3_hifi_model,
	depth_thresh, prop_20x_thresh, prop_30x_thresh,
	mhc_start, mhc_stop, genes_bed, genes_of_interest, genes_of_interest_extended,
	hla_genes_regions_file, vcf2fasta_script, reference_genome_minimap2, reference_genome_vg,
	DNA_bases, stop_codons, IMGT_XML, gff_dir, ARS_dict, gene_dict, CDS_dict, CLASS_I_GENES, dummy_reference, drb_multiallele_reference,
	deepvariant_sif, tandem_repeat_bed, chr6_bed, pbtrgt_repeat_file, picard
)

class Samples:
    # Class variables for reference file paths (imported from config.py)
    deepvariant_sif = deepvariant_sif
    clair3_sif = clair3_sif
    tandem_repeat_bed = tandem_repeat_bed
    chr6_bed = chr6_bed
    pbtrgt_repeat_file = pbtrgt_repeat_file
    
    def __init__(self, input_file, sample_name, platform, output_dir,
                 aligner, snp_caller, indel_caller, trim_adapters=False, adapter_file=None,
                 threads=1, read_group_string=None, clean_up=False, scheme=None,
                 clair3_model=None, rescue_refcalls=False):
        
        # Original initialization code
        self.ORIGINAL_CWD = os.getcwd()
        self.input_file = os.path.realpath(os.path.abspath(input_file))
        self.scheme = scheme

        if not os.path.exists(self.input_file):
            raise FileNotFoundError(f"Input file not found: {self.input_file}")
        if not os.access(self.input_file, os.R_OK):
            raise OSError(f"Input file is not readable: {self.input_file}")

        self.sample_ID = str(sample_name).strip()
        if not self.sample_ID:
            raise ValueError("Sample name cannot be empty or only whitespace.")

        self.adapters = trim_adapters
        self.adapter_file = adapter_file

        # Validate adapter file if specified by user
        if self.adapters and self.adapter_file:
            if not os.path.exists(self.adapter_file):
                raise FileNotFoundError(f"Adapter file not found: {self.adapter_file}")
            if not os.access(self.adapter_file, os.R_OK):
                raise OSError(f"Adapter file is not readable: {self.adapter_file}")

            with open(self.adapter_file) as f:
                sequences = [str(record.seq).strip().upper() for record in SeqIO.parse(f, "fasta")]

            if len(sequences) != 2:
                raise ValueError(
                    f"Adapter file must contain exactly two sequences (forward, reverse), found {len(sequences)}."
                )

            self.five_prime_adapter, self.three_prime_adapter = sequences[:2]
        
        else:
            self.five_prime_adapter = None
            self.three_prime_adapter = None

        self.platform = platform.upper()
        self.threads = threads
        self.aligner = aligner
        self.snp_caller = snp_caller
        self.indel_caller = indel_caller
        self.rescue_refcalls = rescue_refcalls
        self.clean_up = clean_up
        self.clair3_model = clair3_model if clair3_model else (clair3_ont_model if self.platform == "ONT" else clair3_hifi_model)

        output_dir_abs = os.path.realpath(os.path.abspath(os.path.join(output_dir, self.sample_ID)))

        # Ensure input dir is not inside the output directory 
        try:
            inside = os.path.commonpath([self.input_file, output_dir_abs]) == output_dir_abs
        except ValueError:
            inside = False
        if inside:
            raise ValueError(
                f"Input file {self.input_file} is inside the output directory {output_dir_abs}. "
                "Please place the input file outside the output directory."
            )

        self.output_dir = output_dir_abs
        os.makedirs(self.output_dir, exist_ok=True)

        # Platform-agnostic output directories
        self.fastq_raw_dir = os.path.join(self.output_dir, "fastq_raw")
        self.fastq_trimmed_dir = os.path.join(self.output_dir, "fastq_trimmed")
        self.mapped_bam_dir = os.path.join(self.output_dir, "mapped_bam")
        self.parsed_haploblock_dir = os.path.join(self.output_dir, "haploblocks")
        self.genotypes_dir = os.path.join(self.output_dir, "genotype_calls")
        self.sv_dir = os.path.join(self.output_dir, "structural_variant_vcf")
        self.filtered_vcf_dir = os.path.join(self.output_dir, "filtered_vcf")
        self.single_gene_vcf_dir = os.path.join(self.filtered_vcf_dir, "single_gene_vcf")
        self.vcf2fasta_out_dir = os.path.join(self.output_dir, "vcf2fasta_out")
        self.hla_fasta_dir = os.path.join(self.output_dir, "hla_fasta_haplotypes")
        self.hla_typing_dir = os.path.join(self.output_dir, "hla_typing_results")
        self.mosdepth_dir = os.path.join(self.output_dir, "mosdepth")
        self.phased_vcf_dir = os.path.join(self.output_dir, "phased_vcf")

        platform_dirs = []

        # PacBio-specific directories
        if self.platform == "PACBIO":
            self.pbtrgt_dir = os.path.join(self.output_dir, "pbtrgt_vcf")
            platform_dirs.extend([self.pbtrgt_dir])

        self.combined_dirs = [
            self.fastq_raw_dir, self.fastq_trimmed_dir, 
            self.mapped_bam_dir, self.parsed_haploblock_dir, 
            self.genotypes_dir, self.sv_dir, 
            self.filtered_vcf_dir, self.single_gene_vcf_dir, self.vcf2fasta_out_dir, 
            self.hla_fasta_dir, self.hla_typing_dir,
            self.mosdepth_dir, self.phased_vcf_dir
        ] + platform_dirs

        for directory in self.combined_dirs:
            os.makedirs(directory, exist_ok=True)

        # Define all file paths as properties to avoid os.path.join() in workflow functions
        self._define_file_paths()

        parsed_rg = self.parse_input_file(self.input_file)

        if read_group_string is not None and str(read_group_string).strip():
            self.read_group_string = str(read_group_string).strip()
        else:
            self.read_group_string = parsed_rg

        print(f"Processing Sample {self.sample_ID}!\n\n")
        print(f"Sample ID: {self.sample_ID}")
        print(f"Read Group: {self.read_group_string}")
        print("\n\n")

        self.prepare_raw_fastq()

    def parse_input_file(self, input_path):
        if input_path.endswith(".bam"):
            self.format = "BAM"
            self.verify_bam_integrity(input_path)
            read_count = self.count_bam_reads(input_path)

            if read_count < min_reads_sample:
                raise ValueError(f"Input BAM file {input_path} contains too few reads: {read_count:,}")

            with pysam.AlignmentFile(input_path, "rb", check_sq=False) as bamfile:
                header = bamfile.header.to_dict()
                rg_list = header.get("RG", [])
                if not rg_list:
                    raise ValueError(f"No @RG entry found in BAM header for {input_path}")
                if len(rg_list) > 1:
                    raise ValueError(f"BAM file {input_path} contains multiple @RG entries. Only one read group is supported per file/sample.")

                rg = rg_list[0]
                rg_id = rg.get("ID") or f"{self.sample_ID}_RG"
                rg_pl = self.platform
                rg_sm = self.sample_ID
                rg_lb = rg.get("LB", self.sample_ID)
                rg_pu = rg.get("PU", f"{self.sample_ID}_PU")

        elif input_path.endswith((".fastq", ".fq", ".fastq.gz", ".fq.gz")):
            if input_path.endswith(".gz"):
                self.format = "FASTQ.GZ"
            else:
                self.format = "FASTQ"

            read_count, mean_read_length = self.run_fastplong(input_path)
            if read_count < min_reads_sample:
                raise ValueError(f"Input fastq file {input_path} contains too few reads: {read_count:,}")
            if mean_read_length < min_read_length:
                raise ValueError(f"Input fastq file {input_path} contains short-read data. Mean read length: {mean_read_length}")

            rg_id = f"{self.sample_ID}_RG"
            rg_pl = self.platform
            rg_sm = self.sample_ID
            rg_lb = self.sample_ID
            rg_pu = f"{self.sample_ID}_PU"

        else:
            raise ValueError(f"Unsupported input file format: {input_path}")

        rg_string = f"@RG\\tID:{rg_id}\\tSM:{rg_sm}\\tPL:{rg_pl}\\tLB:{rg_lb}\\tPU:{rg_pu}"

        return rg_string

    def verify_bam_integrity(self, bam_path):
        quickcheck_cmd = f"samtools quickcheck -u -v {bam_path}"
        result = subprocess.run(quickcheck_cmd, shell=True, capture_output=True, text=True)
        if result.returncode != 0:
            raise ValueError(f"Input BAM is corrupt. Samtools quickcheck failed for {bam_path}:\n{result.stderr}")

    def count_bam_reads(self, bam_path):
        count_cmd = f"samtools view -@ {self.threads} -c {bam_path}"
        result = subprocess.run(count_cmd, shell=True, capture_output=True, text=True, check=True)
        count = int(result.stdout.strip())
        print(f"Total BAM records in {bam_path}: {count:,}")
        return count

    def run_fastplong(self, fq_path):
        html_path = os.path.join(self.fastq_raw_dir, self.sample_ID + ".fastplong.html")
        json_path = os.path.join(self.fastq_raw_dir, self.sample_ID + ".fastplong.json")

        # -A: disable adapter trimming, -Q: disable quality filtering, -L: disable length filtering
        # -m 0: no mean quality threshold, -n 100000: discard reads with >100000 N bases (i.e. keep all)
        # QC-only mode: generate reports without any read filtering or trimming
        fastplong_cmd = f"fastplong -i {fq_path} -h {html_path} -j {json_path} -w {self.threads} -A -Q -L -m 0 -n 100000"

        result = subprocess.run(fastplong_cmd, shell=True, check=False, stdout=subprocess.DEVNULL, stderr=subprocess.PIPE, text=True)
        if result.returncode != 0:
            raise ValueError(f"Input fastq is corrupt. Fastplong failed for {fq_path}:\n{result.stderr}")

        with open(json_path) as f:
            data = json.load(f)
        total_reads = data["summary"]["after_filtering"]["total_reads"]
        mean_read_length = int(data["summary"]["after_filtering"]["read_mean_length"])
        print(f"Total FASTQ records in {fq_path}: {total_reads:,}")
        return total_reads, mean_read_length

    def prepare_raw_fastq(self):
        if self.format == "BAM":
            output_file = os.path.join(self.fastq_raw_dir, self.sample_ID + ".fastq.gz")
            convert_bam_to_fastq(self.input_file, output_file, self.platform, self.threads)
        elif self.format == "FASTQ":
            new_fq = os.path.join(self.fastq_raw_dir, self.sample_ID + ".fastq")
            shutil.copy(self.input_file, new_fq)
            pigz_cmd = f"pigz -f -p {self.threads} {new_fq}"
            subprocess.run(pigz_cmd, shell=True, check=True)
            expected_output = os.path.join(self.fastq_raw_dir, self.sample_ID + ".fastq.gz")
            if not os.path.exists(expected_output):
                raise RuntimeError(f"Compression failed: {expected_output} not found")

        elif self.format == "FASTQ.GZ":
            new_fq = os.path.join(self.fastq_raw_dir, self.sample_ID + ".fastq.gz")
            shutil.copy(self.input_file, new_fq)
            expected_output = os.path.join(self.fastq_raw_dir, self.sample_ID + ".fastq.gz")
            if not os.path.exists(expected_output):
                raise RuntimeError(f"Compression failed: {expected_output} not found")

    def _define_file_paths(self):
        """Define all file paths as properties to avoid path construction in workflow functions"""
        
        # Raw and trimmed FASTQ files
        self.raw_fastq = os.path.join(self.fastq_raw_dir, f"{self.sample_ID}.fastq.gz")
        self.trimmed_fastq = os.path.join(self.fastq_trimmed_dir, f"{self.sample_ID}.trimmed.fastq.gz")
        self.trimmed_pbmarkdup_fastq = os.path.join(self.fastq_trimmed_dir, f"{self.sample_ID}.trimmed.pbmarkdup.fastq")
        self.trimmed_pbmarkdup_fastq_gz = os.path.join(self.fastq_trimmed_dir, f"{self.sample_ID}.trimmed.pbmarkdup.fastq.gz")
        
        # BAM files
        self.hg38_bam = os.path.join(self.mapped_bam_dir, f"{self.sample_ID}.hg38.bam")
        self.hg38_bam_drb = os.path.join(self.mapped_bam_dir, f"{self.sample_ID}.hg38.bam.drb")
        self.hg38_chr6_bam = os.path.join(self.mapped_bam_dir, f"{self.sample_ID}.hg38.chr6.bam")
        self.pangenome_bam = os.path.join(self.mapped_bam_dir, f"{self.sample_ID}.pangenome.bam")
        self.pg_bam = os.path.join(self.mapped_bam_dir, f"{self.sample_ID}.pg.bam")
        self.pg_reheader_bam = os.path.join(self.mapped_bam_dir, f"{self.sample_ID}.pg.reheader.bam")
        self.pg_mapq_reassign_bam = os.path.join(self.mapped_bam_dir, f"{self.sample_ID}.pg.mapq_reassign.bam")
        self.pg_mapq_reassign_mrkdup_bam = os.path.join(self.mapped_bam_dir, f"{self.sample_ID}.pg.mapq_reassign.mrkdup.bam")
        self.hg38_mrkdup_bam = os.path.join(self.mapped_bam_dir, f"{self.sample_ID}.hg38.mrkdup.bam")
        self.hg38_rmdup_chr6_bam = os.path.join(self.mapped_bam_dir, f"{self.sample_ID}.hg38.rmdup.chr6.bam")
        self.hg38_rmdup_chr6_haplotag_bam = os.path.join(self.mapped_bam_dir, f"{self.sample_ID}.hg38.rmdup.chr6.haplotag.bam")
        
        # Metrics files
        self.pg_mapq_reassign_mrkdup_metrics = os.path.join(self.mapped_bam_dir, f"{self.sample_ID}.pg.mapq_reassign.mrkdup.metrics.txt")
        self.hg38_mrkdup_metrics = os.path.join(self.mapped_bam_dir, f"{self.sample_ID}.hg38.mrkdup.metrics.txt")
        
        # VCF files
        self.snv_vcf = os.path.join(self.genotypes_dir, f"{self.sample_ID}.vcf.gz")
        self.snv_gvcf = os.path.join(self.genotypes_dir, f"{self.sample_ID}.g.vcf.gz")
        self.bcftools_snp_vcf = os.path.join(self.genotypes_dir, f"{self.sample_ID}.bcftools.snps.vcf.gz")
        self.dv_full_vcf = os.path.join(self.genotypes_dir, f"{self.sample_ID}.dv.vcf.gz")
        self.dv_rescued_vcf = os.path.join(self.genotypes_dir, f"{self.sample_ID}.dv.rescued.vcf.gz")
        self.dv_indel_vcf = os.path.join(self.genotypes_dir, f"{self.sample_ID}.dv.indels.vcf.gz")
        self.clair3_full_vcf = os.path.join(self.genotypes_dir, f"{self.sample_ID}.clair3.vcf.gz")
        self.hybrid_indel_vcf = os.path.join(self.genotypes_dir, f"{self.sample_ID}.hybrid.indels.vcf.gz")
        self.sv_vcf = os.path.join(self.sv_dir, f"{self.sample_ID}.SV.vcf.gz")
        self.sv_svsig = os.path.join(self.sv_dir, f"{self.sample_ID}.svsig.gz")
        self.tr_vcf = os.path.join(self.pbtrgt_dir, f"{self.sample_ID}.TR.vcf.gz") if hasattr(self, 'pbtrgt_dir') else None
        
        # Phased VCF files
        self.hiphase_snv_vcf = os.path.join(self.phased_vcf_dir, f"{self.sample_ID}.hiphase.vcf.gz")
        self.hiphase_sv_vcf = os.path.join(self.phased_vcf_dir, f"{self.sample_ID}.hiphase.SV.vcf.gz")
        self.hiphase_tr_vcf = os.path.join(self.phased_vcf_dir, f"{self.sample_ID}.hiphase.TR.vcf.gz")
        self.hiphase_joint_vcf = os.path.join(self.phased_vcf_dir, f"{self.sample_ID}.hiphase.joint.vcf.gz")
        
        # LongPhase VCF files (for ONT)
        self.longphase_vcf = os.path.join(self.phased_vcf_dir, f"{self.sample_ID}.longphase.vcf.gz")
        self.longphase_sv_vcf = os.path.join(self.phased_vcf_dir, f"{self.sample_ID}.longphase_SV.vcf.gz")
        self.longphase_merged_vcf = os.path.join(self.phased_vcf_dir, f"{self.sample_ID}.longphase.merged.vcf.gz")
        
        # Phasing output files
        self.phased_summary = os.path.join(self.phased_vcf_dir, f"{self.sample_ID}.phased.summary.txt")
        self.phased_stats = os.path.join(self.phased_vcf_dir, f"{self.sample_ID}.phased.stats.txt")
        self.phased_blocks = os.path.join(self.phased_vcf_dir, f"{self.sample_ID}.phased.blocks.txt")
        self.phased_haploblocks = os.path.join(self.phased_vcf_dir, f"{self.sample_ID}.phased.haploblocks.txt")
        self.phased_haploblocks_gtf = os.path.join(self.phased_vcf_dir, f"{self.sample_ID}.phased.haploblocks.gtf")
        
        # Filtered VCF files
        self.pass_vcf = os.path.join(self.phased_vcf_dir, f"{self.sample_ID}_PASS.vcf.gz")
        self.fail_vcf = os.path.join(self.phased_vcf_dir, f"{self.sample_ID}_FAIL.vcf.gz")
        self.pass_unphased_vcf = os.path.join(self.phased_vcf_dir, f"{self.sample_ID}_PASS_UNPHASED.vcf.gz")
        self.filtered_vcf = os.path.join(self.filtered_vcf_dir, f"{self.sample_ID}_filtered.vcf.gz")
        self.unphased_tsv = os.path.join(self.phased_vcf_dir, f"{self.sample_ID}.unphased.tsv")
        
        # HLA typing output files
        self.hla_gene_fasta = os.path.join(self.hla_fasta_dir, f"{self.sample_ID}_HLA_haplotypes_gene.fasta")
        self.hla_cds_fasta = os.path.join(self.hla_fasta_dir, f"{self.sample_ID}_HLA_haplotypes_CDS.fasta")
        
        # Mosdepth files
        self.mosdepth_regions = os.path.join(self.mosdepth_dir, f"{self.sample_ID}.regions.bed.gz")
        self.mosdepth_thresholds = os.path.join(self.mosdepth_dir, f"{self.sample_ID}.thresholds.bed.gz")
        
        # Haploblock analysis files
        self.phased_genes_tsv = os.path.join(self.parsed_haploblock_dir, "phased_genes.tsv")
        self.incomplete_genes_csv = os.path.join(self.parsed_haploblock_dir, "incomplete.csv")

        # DRB34 reads file
        self.DRB34_reads_file = os.path.join(self.mapped_bam_dir, f"{self.sample_ID}.drb34_reads.txt")
        


def build_workflow_config(sample):
	"""
	Build a configuration dictionary for workflows from a sample object.
	This avoids passing the entire sample object to workflow functions.
	"""
	# Set reference genome based on aligner
	if sample.aligner == "minimap2":
		reference_genome = reference_genome_minimap2
	elif sample.aligner == "vg":
		reference_genome = reference_genome_vg
	else:
		raise ValueError(f"Unknown aligner: {sample.aligner}")
	
	config = {
		# Core sample information
		'sample_ID': sample.sample_ID,
		'threads': sample.threads,
		'platform': sample.platform,
		'aligner': sample.aligner,
		'snp_caller': sample.snp_caller,
		'indel_caller': sample.indel_caller,
		'rescue_refcalls': sample.rescue_refcalls,
		'scheme': sample.scheme,
		
		# Adapter settings
		'adapters': sample.adapters,
		'adapter_file': sample.adapter_file,
		'five_prime_adapter': sample.five_prime_adapter,
		'three_prime_adapter': sample.three_prime_adapter,
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
		'single_gene_vcf_dir': sample.single_gene_vcf_dir,
		'vcf2fasta_out_dir': sample.vcf2fasta_out_dir,
		'hla_fasta_dir': sample.hla_fasta_dir,
		'hla_typing_dir': sample.hla_typing_dir,
		'gff_dir': gff_dir,  # Use global GFF directory from config.py
		'clean_up': sample.clean_up,
		
		# File paths (pre-constructed to avoid os.path.join() in workflows)
		'raw_fastq': sample.raw_fastq,
		'trimmed_fastq': sample.trimmed_fastq,
		'trimmed_pbmarkdup_fastq': sample.trimmed_pbmarkdup_fastq,
		'trimmed_pbmarkdup_fastq_gz': sample.trimmed_pbmarkdup_fastq_gz,
		'hg38_bam': sample.hg38_bam,
		'hg38_bam_drb': sample.hg38_bam_drb,
		'hg38_chr6_bam': sample.hg38_chr6_bam,
		'DRB34_reads_file': sample.DRB34_reads_file,
		'pangenome_bam': sample.pangenome_bam,
		'pg_bam': sample.pg_bam,
		'pg_reheader_bam': sample.pg_reheader_bam,
		'pg_mapq_reassign_bam': sample.pg_mapq_reassign_bam,
		'pg_mapq_reassign_mrkdup_bam': sample.pg_mapq_reassign_mrkdup_bam,
		'hg38_mrkdup_bam': sample.hg38_mrkdup_bam,
		'hg38_mrkdup_metrics': sample.hg38_mrkdup_metrics,
		'hg38_rmdup_chr6_bam': sample.hg38_rmdup_chr6_bam,
		'hg38_rmdup_chr6_haplotag_bam': sample.hg38_rmdup_chr6_haplotag_bam,
		'snv_vcf': sample.snv_vcf,
		'snv_gvcf': sample.snv_gvcf,
		'bcftools_snp_vcf': sample.bcftools_snp_vcf,
		'dv_full_vcf': sample.dv_full_vcf,
		'dv_rescued_vcf': sample.dv_rescued_vcf,
		'dv_indel_vcf': sample.dv_indel_vcf,
		'clair3_full_vcf': sample.clair3_full_vcf,
		'hybrid_indel_vcf': sample.hybrid_indel_vcf,
		'sv_vcf': sample.sv_vcf,
		'sv_svsig': sample.sv_svsig,
		'tr_vcf': sample.tr_vcf,
		'hiphase_snv_vcf': sample.hiphase_snv_vcf,
		'hiphase_sv_vcf': sample.hiphase_sv_vcf,
		'hiphase_tr_vcf': sample.hiphase_tr_vcf,
		'hiphase_joint_vcf': sample.hiphase_joint_vcf,
		'longphase_vcf': sample.longphase_vcf,
		'longphase_sv_vcf': sample.longphase_sv_vcf,
		'longphase_merged_vcf': sample.longphase_merged_vcf,
		'phased_summary': sample.phased_summary,
		'phased_stats': sample.phased_stats,
		'phased_blocks': sample.phased_blocks,
		'phased_haploblocks': sample.phased_haploblocks,
		'phased_haploblocks_gtf': sample.phased_haploblocks_gtf,
		'pass_vcf': sample.pass_vcf,
		'fail_vcf': sample.fail_vcf,
		'pass_unphased_vcf': sample.pass_unphased_vcf,
		'filtered_vcf': sample.filtered_vcf,
		'unphased_tsv': sample.unphased_tsv,
		'hla_gene_fasta': sample.hla_gene_fasta,
		'hla_cds_fasta': sample.hla_cds_fasta,
		'mosdepth_regions': sample.mosdepth_regions,
		'mosdepth_thresholds': sample.mosdepth_thresholds,
		'phased_genes_tsv': sample.phased_genes_tsv,
		'incomplete_genes_csv': sample.incomplete_genes_csv,
		
		# Reference files and tool paths (from Samples class and config.py)
		'reference_genome': reference_genome,
		'dummy_reference': dummy_reference,
		'drb_multiallele_reference': drb_multiallele_reference,
		'deepvariant_sif': Samples.deepvariant_sif,
		'clair3_sif': Samples.clair3_sif,
		'clair3_model': sample.clair3_model,
		'chr6_bed': Samples.chr6_bed,
		'tandem_repeat_bed': Samples.tandem_repeat_bed,
		'pbtrgt_repeat_file': Samples.pbtrgt_repeat_file,
		'reference_gbz': reference_genome_vg_gbz,
		'ref_paths': reference_genome_vg_paths,
		'vg': vg,
		'longphase': longphase,
		'sawfish': sawfish,
		'picard': picard,
		'depth_thresh': depth_thresh,
		'prop_20x_thresh': prop_20x_thresh,
		'prop_30x_thresh': prop_30x_thresh,
		'mhc_start': mhc_start,
		'mhc_stop': mhc_stop,
		'genes_bed': genes_bed,
		'genes_of_interest': genes_of_interest,
		'genes_of_interest_extended': genes_of_interest_extended,
		'hla_genes_regions_file': hla_genes_regions_file,
		'vcf2fasta_script': vcf2fasta_script,
		'DNA_bases': DNA_bases,
		'stop_codons': stop_codons,
		'IMGT_XML': IMGT_XML,
		'ARS_dict': ARS_dict,
		'gene_dict': gene_dict,
		'CDS_dict': CDS_dict,
		'CLASS_I_GENES': CLASS_I_GENES,
	}
	
	# Platform-specific directories
	if sample.platform == "PACBIO":
		config['pbtrgt_dir'] = getattr(sample, 'pbtrgt_dir', None)
		config['ORIGINAL_CWD'] = getattr(sample, 'ORIGINAL_CWD', None)
	
	return config
