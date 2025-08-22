import os
import shutil
import subprocess
import json
import pysam
from Bio import SeqIO
from preprocess_methods import convert_bam_to_fastq
from config import min_reads_sample, min_read_length

class Samples:
    # Class variables for reference file paths
    reference_fasta = None
    deepvariant_sif = None
    tandem_repeat_bed = None
    chr6_bed = None
    pbtrgt_repeat_file = None
    
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
