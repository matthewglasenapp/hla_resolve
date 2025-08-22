#!/usr/bin/env python3

import os
import shutil
import subprocess
import sys
import time
import pysam
import argparse
import json
import textwrap
from Bio import SeqIO

# Core preprocessing functions (always needed)
from hla_resolve.preprocess_methods import (
    convert_bam_to_fastq,
    trim_adapters,
    mark_duplicates_pbmarkdup,
    trim_reads,
    run_porechop_abi,
    run_fastqc,
    align_to_reference_minimap,
    align_to_reference_vg,
    reassign_mapq,
    filter_reads,
    run_mosdepth,
    parse_mosdepth,
    filter_reads,
    mark_duplicates_picard
)

# Workflow functions
from hla_resolve.ont_workflow import preprocess_ont_sample
from hla_resolve.pacbio_workflow import preprocess_pacbio_sample

# Variant calling functions (platform/genotyper dependent)
from hla_resolve.preprocess_methods import (
    call_variants_bcftools,
    call_variants_deepvariant,
    call_variants_clair3,
    genotype_tandem_repeats,
    call_structural_variants_pbsv,
    call_structural_variants_sawfish,
    call_structural_variants_sniffles
)

# Phasing functions (platform dependent)
from hla_resolve.preprocess_methods import (
    phase_genotypes_hiphase,
    phase_genotypes_longphase,
    merge_hiphase_vcfs,
    merge_longphase_vcfs
)

# Haploblock analysis functions
from hla_resolve.investigate_haploblocks_methods import (
    parse_haploblocks,
    evaluate_gene_haploblocks
)

# FASTA reconstruction functions
from hla_resolve.reconstruct_fasta_methods import (
    filter_vcf,
    run_vcf2fasta,
    parse_fastas
)

# Note: hla_typer.main is imported lazily in the type_hla_alleles function to avoid overhead

genes_of_interest = ("HLA-A", "HLA-B", "HLA-C", "HLA-DRB1", "HLA-DQA1", "HLA-DQB1", "HLA-DPA1", "HLA-DPB1")

# Minimum reads per sample
# DeepVariant is stalling and not exiting for samples with very few BAM records (e.g., HG01891: 35 mapped reads to chr6)
# Set mapped chr6 reads threshold at which variant calling should not proceed
min_reads_sample = 100

# This program is for long-read data only. 
# Require that mean read length is at least 300 bp or higher
min_read_length = 300

# IPD/IMGT HLA XML file is now a class variable in the Samples class

# Ensure all required tools are installed and executable
def check_required_commands():    
    print("Checking the installation status of the required bioinformatics tools!")

    required_commands = [
        "bam2fastq",
        "bcftools",
        "bgzip",
        "cutadapt",
        "fastplong",
        "fastqc",
        "gatk",
        "hiphase",
        "pbmarkdup",
        "pbmm2",
        "pbsv",
        "pigz",
        "samtools",
        "singularity",
        "sniffles",
        "tabix",
        "trgt",
    ]

    missing_commands = []
    for command in required_commands:
        if shutil.which(command) is None:
            missing_commands.append(command)
    if len(missing_commands) != 0:
        print(f"Error: Missing the following commands: {', '.join(missing_commands)}")
        sys.exit(1)
    else:
        print("All tools required are installed!")
        print("\n\n")

class Samples:
    # Class variables for reference file paths
    reference_fasta = None
    deepvariant_sif = None
    tandem_repeat_bed = None
    chr6_bed = None
    pbtrgt_repeat_file = None
    
    # Class variables for other hardcoded paths and constants
    # Pangenome graph reference info 
    reference_gbz = "/hb/scratch/ogarci12/hybridcapture_pangenome/ref/hprc-v1.0-mc-grch38-minaf.0.1.gbz"
    ref_paths = "/hb/scratch/ogarci12/hybridcapture_pangenome/ref/hprc-v1.0-mc-grch38-minaf.0.1.dict"
    vg = "/hb/scratch/ogarci12/hybridcapture_pangenome/vg"
    mosdepth_regions_file = "/hb/scratch/mglasena/test_hla_resolve/hla_resolve/hla_resolve/hla_genes.bed"
    # Transposase mosaic end binding sequence
    # The TE sequence (and its reverse complement) introduced during tagmentation still needs to be removed
    # Adapters and barcodes were removed by PacBio with lima
    me = "AGATGTGTATAAGAGACAG"
    me_rc = "CTGTCTCTTATACACATCT"
    longphase = "/hb/home/mglasena/software/longphase/longphase_linux-x64"
    prowler_trimmer = "/hb/home/mglasena/software/ProwlerTrimmer/TrimmerLarge.py"
    sawfish = "/hb/home/mglasena/software/sawfish-v2.0.3-x86_64-unknown-linux-gnu/bin/sawfish"
    clair3_ont_model_path = "/hb/home/mglasena/.conda/envs/clair3/bin/models/r941_prom_sup_g5014"
    clair3_hifi_model_path = "/hb/home/mglasena/.conda/envs/clair3/bin/models/hifi_revio"
    # Coverage Thresholds
    # Might have to relax if you didn't get the flanking regions of the gene (i.e., UTR)
    depth_thresh = 30
    prop_20x_thresh = 0.8
    prop_30x_thresh = 0.8
    
    # Additional paths from reconstruct_fasta_methods.py
    vcf2fasta_script = "/hb/scratch/mglasena/vcf2fasta/vcf2fasta.py"
    reference_genome = "/hb/groups/cornejo_lab/matt/hla_capture/input_data/reference/GCA_000001405.15_GRCh38_no_alt_analysis_set.fa"
    gff_dir = "/hb/groups/cornejo_lab/matt/hla_capture/input_data/hla_gff/"
    hla_genes_regions_file = "/hb/scratch/mglasena/test_hla_resolve/hla_resolve/hla_resolve/hla_genes.bed"
    
    # Additional paths and constants from investigate_haploblocks_methods.py
    genes_bed = "/hb/groups/cornejo_lab/matt/hla_capture/input_data/reference/parse_haploblocks_bed.bed"
    genes_of_interest = ("HLA-A", "HLA-B", "HLA-C", "HLA-DRB1", "HLA-DQA1", "HLA-DQA2", "HLA-DQB1", "HLA-DQB2", "HLA-DPA1", "HLA-DPB1")
    # Extended MHC coordinates
    mhc_start = 29555628
    mhc_stop = 33409896
    # DNA bases and stop codons
    DNA_bases = {"A", "T", "G", "C"}
    stop_codons = ["TAA", "TAG", "TGA"]
    
    # IPD/IMGT HLA XML file
    IMGT_XML = "/hb/scratch/mglasena/test_hla_resolve/hla_resolve/hla_resolve/hla.xml"
    
    def __init__(self, input_file, sample_name, platform, output_dir, 
                 aligner, genotyper, trim_adapters=False, adapter_file=None, 
                 threads=1, read_group_string=None, clean_up=False):
        
        # Set the class variables based on data_dir
        if not Samples.reference_fasta:
            data_dir = os.path.join(os.path.dirname(os.path.abspath(__file__)), "data")
            Samples.reference_fasta = os.path.join(data_dir, "reference/augmented_hg38.fa")
            Samples.deepvariant_sif = os.path.join(data_dir, "deepvariant_sif/deepvariant.sif")
            Samples.tandem_repeat_bed = os.path.join(data_dir, "repeats_bed/human_GRCh38_no_alt_analysis_set.trf.bed")
            Samples.chr6_bed = os.path.join(data_dir, "reference/chr6.bed")
            Samples.pbtrgt_repeat_file = os.path.join(data_dir, "repeats_bed/test_chr6_polymorphic_repeats.hg38.bed")
        
        # Class variables are now accessed directly via Samples.class_variable_name
        # No need to copy them to instance variables anymore
        
        # data_dir points directly to the data/ subdirectory
        self.data_dir = os.path.join(os.path.dirname(os.path.abspath(__file__)), "data")
        
        # Original initialization code
        self.ORIGINAL_CWD = os.getcwd()
        self.input_file = os.path.realpath(os.path.abspath(input_file))

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
        self.genotyper = genotyper
        self.clean_up = clean_up

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
            self.filtered_vcf_dir, self.vcf2fasta_out_dir, 
            self.hla_fasta_dir, self.hla_typing_dir,
            self.mosdepth_dir, self.phased_vcf_dir
        ] + platform_dirs

        for directory in self.combined_dirs:
            os.makedirs(directory, exist_ok=True)

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

    def print_results(self):
        results_file = os.path.join(self.hla_typing_dir, "refined_allele_output.csv")
        with open(results_file, "r") as f:
            results = f.read().splitlines()[1].split(",")[1:]
        print(f"{self.sample_ID} HLA Star Allele Calls")
        for item in results:
            print(item)

def main():
    parser = argparse.ArgumentParser(
    description="Run HLA-Resolve",
    formatter_class=argparse.ArgumentDefaultsHelpFormatter,
    epilog=textwrap.dedent("""\
        Examples:
          python3 -m hla_resolve.cli --input_file reads.bam --sample_name HG002 --platform pacbio --output_dir out --aligner minimap2 --genotyper deepvariant --threads 10
    """),
)
    parser.add_argument("--input_file", required=True, help="Path to the raw sequencing reads file")
    parser.add_argument("--sample_name", required=True, help="Override the parsed sample name", default=None)
    parser.add_argument("--platform", choices=["pacbio", "ont"], required=True, help="Specify sequencing platform (pacbio, ont)")
    parser.add_argument("--output_dir", required=True, help="Output Directory", default=None)
    parser.add_argument("--aligner", choices=["minimap2", "vg"], required=True, help="Tool for reference genome alignment", default=None)
    parser.add_argument("--genotyper", choices=["bcftools", "clair3", "deepvariant"], required=False, help="Tool for genotyping", default="deepvariant")
    parser.add_argument("--trim_adapters", action="store_true", help="Enable adapter trimming before processing")
    parser.add_argument("--adapter_file", type=str, required=False, default=None, help="Path to a file with custom adapter sequences (FASTA/FASTQ). If not provided, default adapters will be used.")
    parser.add_argument("--threads", type=int, required=False, help="Number of threads to use", default=6)
    parser.add_argument("--read_group_string", required=False, help="Override the parsed read group string", default=None)
    parser.add_argument("--clean-up", action="store_true", help="Remove intermediate files")
    
    # Show help and exit if no arguments were provided
    if len(sys.argv) == 1:
        parser.print_help()
        parser.exit()

    args = parser.parse_args()

    # Check that all required tools are installed
    check_required_commands()
    
    start_time = time.time()
    
    sample = Samples(input_file=args.input_file, sample_name=args.sample_name, platform=args.platform, output_dir=args.output_dir, aligner=args.aligner, genotyper=args.genotyper, trim_adapters=args.trim_adapters, adapter_file=args.adapter_file, threads=args.threads, read_group_string=args.read_group_string, clean_up=args.clean_up)

    if sample.platform == "PACBIO":
        preprocess_pacbio_sample(
            sample=sample,
            reference_fasta=Samples.reference_fasta,
            vg=Samples.vg,
            reference_gbz=Samples.reference_gbz,
            ref_paths=Samples.ref_paths,
            deepvariant_sif=Samples.deepvariant_sif,
            chr6_bed=Samples.chr6_bed,
            clair3_ont_model_path=Samples.clair3_ont_model_path,
            clair3_hifi_model_path=Samples.clair3_hifi_model_path,
            sawfish=Samples.sawfish,
            pbtrgt_repeat_file=Samples.pbtrgt_repeat_file
        )
    elif sample.platform == "ONT":
        preprocess_ont_sample(
            sample=sample,
            reference_fasta=Samples.reference_fasta,
            vg=Samples.vg,
            reference_gbz=Samples.reference_gbz,
            ref_paths=Samples.ref_paths,
            deepvariant_sif=Samples.deepvariant_sif,
            chr6_bed=Samples.chr6_bed,
            clair3_ont_model_path=Samples.clair3_ont_model_path,
            clair3_hifi_model_path=Samples.clair3_hifi_model_path,
            longphase=Samples.longphase,
            tandem_repeat_bed=Samples.tandem_repeat_bed
        )
    
    sufficient_coverage_genes = run_coverage_analysis(sample)
    
    phased_genes = reconstruct_fasta_sequences(sample, sufficient_coverage_genes)
    
    type_hla_alleles(sample, phased_genes)

    end_time = time.time()
    elapsed_time = end_time - start_time
    minutes, seconds = divmod(elapsed_time,60)
    print(f"Processed sampled in {int(minutes)}:{seconds:.2f}!")

if __name__ == "__main__":
    main()
