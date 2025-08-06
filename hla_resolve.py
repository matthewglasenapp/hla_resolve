import os
import shutil
import subprocess
import sys
import time
import pysam
import argparse
from preprocess_methods import convert_bam_to_fastq, mark_duplicates_pbmarkdup, mark_duplicates_picard, trim_adapters, run_fastqc, trim_reads, align_to_reference_minimap, align_to_reference_vg, reassign_mapq, filter_reads, call_variants_deepvariant, call_variants_clair3, call_structural_variants_pbsv, call_structural_variants_sniffles, genotype_tandem_repeats, phase_genotypes_hiphase, merge_hiphase_vcfs, phase_genotypes_whatshap, phase_genotypes_longphase, merge_longphase_vcfs, run_porechop_abi
from investigate_haploblocks_methods import parse_haploblocks, evaluate_gene_haploblocks
from reconstruct_fasta_methods import filter_vcf, run_vcf2fasta, parse_fastas

# genes_of_interest = ("HLA-A", "HLA-B", "HLA-C", "HLA-DRB1", "HLA-DRB5", "HLA-DQA1", "HLA-DQA2", "HLA-DQB1", "HLA-DQB2", "HLA-DPA1", "HLA-DPB1")
genes_of_interest = ("HLA-A", "HLA-B", "HLA-C")

# Minimum reads per sample
# DeepVariant is stalling and not exiting for samples with very few BAM records (e.g., HG01891: 35 mapped reads to chr6)
# Set mapped chr6 reads threshold at which variant calling should not proceed
min_reads_sample = 100

# Ensure all required tools are installed and executable
def check_required_commands():    
	print("Checking the installation status of the required bioinformatics tools!")

	required_commands = [
		"bam2fastq",
		"bcftools",
		"bgzip",
		"cutadapt",
		"fastqc",
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
		"whatshap"
	]

	missing_commands = []
	for command in required_commands:
		if shutil.which(command) is None:
			missing_commands.append(command)
	if len(missing_commands) != 0:
		print("Error: Missing the following commands: {}".format(", ".join(missing_commands)))
		sys.exit(1)
	else:
		print("All tools required are installed!")
		print("\n\n")

class Samples:
	def __init__(self, unmapped_bam, sample_name, platform, output_dir, threads, read_group_string=None):
		self.unmapped_bam = unmapped_bam
		self.sample_ID = sample_name
		self.platform = platform.upper()
		self.threads = threads

		if read_group_string:
			self.read_group_string = read_group_string
		else:
			self.read_group_string = self.parse_unmapped_bam(unmapped_bam)

		self.output_dir = os.path.join(output_dir, self.sample_ID)
		os.makedirs(self.output_dir, exist_ok=True)

		# Platform-agnostic output directories
		self.fastq_raw_dir = os.path.join(self.output_dir, "fastq_raw")
		self.mapped_bam_dir = os.path.join(self.output_dir, "mapped_bam")
		self.parsed_haploblock_dir = os.path.join(self.output_dir, "haploblocks")
		self.whatshap_phased_vcf_dir = os.path.join(self.output_dir, "phased_vcf_whatshap")
		self.filtered_vcf_dir = os.path.join(self.output_dir, "filtered_vcf")
		self.vcf2fasta_out_dir = os.path.join(self.output_dir, "vcf2fasta_out")
		self.hla_fasta_dir = os.path.join(self.output_dir, "hla_fasta_haplotypes")

		platform_dirs = []

		# PacBio-specific directories
		if self.platform == "PACBIO":
			self.fastq_rmdup_dir = os.path.join(self.output_dir, "fastq_rmdup")
			self.fastq_rmdup_cutadapt_dir = os.path.join(self.output_dir, "fastq_rmdup_cutadapt")
			self.deepvariant_dir = os.path.join(self.output_dir, "deepvariant_vcf")
			self.pbsv_dir = os.path.join(self.output_dir, "pbsv_vcf")
			self.pbtrgt_dir = os.path.join(self.output_dir, "pbtrgt_vcf")
			self.hiphase_phased_vcf_dir = os.path.join(self.output_dir, "phased_vcf_hiphase")

			platform_dirs.extend([
				self.fastq_rmdup_dir, self.fastq_rmdup_cutadapt_dir,
				self.deepvariant_dir, self.pbsv_dir, self.pbtrgt_dir,
				self.hiphase_phased_vcf_dir
			])

		# ONT-specific directories
		elif self.platform == "ONT":
			self.fastq_porechop_dir = os.path.join(self.output_dir, "fastq_porechop")
			self.fastq_prowler_dir = os.path.join(self.output_dir, "fastq_prowler")
			self.clair3_dir = os.path.join(self.output_dir, "clair3_vcf")
			self.sniffles_dir = os.path.join(self.output_dir, "sniffles")
			self.merged_vcf_dir = os.path.join(self.output_dir, "merged_vcf")
			self.longphase_phased_vcf_dir = os.path.join(self.output_dir, "phased_vcf_longphase")

			platform_dirs.extend([
				self.fastq_porechop_dir, self.fastq_prowler_dir, self.clair3_dir, self.sniffles_dir,
				self.merged_vcf_dir, self.longphase_phased_vcf_dir
			])

		combined_dirs = [
			self.fastq_raw_dir, self.mapped_bam_dir,
			self.parsed_haploblock_dir, self.whatshap_phased_vcf_dir,
			self.filtered_vcf_dir, self.vcf2fasta_out_dir, self.hla_fasta_dir
		] + platform_dirs

		for directory in combined_dirs:
			os.makedirs(directory, exist_ok=True)

		print(f"Processing Sample {self.sample_ID}!\n\n")
		print(f"Sample ID: {self.sample_ID}")
		print(f"Read Group: {self.read_group_string}")

	def parse_unmapped_bam(self, bam_path):
		import pysam

		with pysam.AlignmentFile(bam_path, "rb", check_sq=False) as bamfile:
			header = bamfile.header.to_dict()
			rg_list = header.get("RG", [])
			if not rg_list:
				raise ValueError(f"No @RG entry found in BAM header for {bam_path}")

			rg = rg_list[0]
			rg_id = rg.get("ID") or f"{self.sample_ID}_RG"
			rg_pl = self.platform
			rg_sm = self.sample_ID
			rg_lb = rg.get("LB", self.sample_ID)
			rg_pu = rg.get("PU", f"{self.sample_ID}_PU")

			rg_string = f"@RG\\tID:{rg_id}\\tSM:{rg_sm}\\tPL:{rg_pl}\\tLB:{rg_lb}\\tPU:{rg_pu}"
			return rg_string

Samples.convert_bam_to_fastq = convert_bam_to_fastq
Samples.mark_duplicates_pbmarkdup = mark_duplicates_pbmarkdup
Samples.mark_duplicates_picard = mark_duplicates_picard
Samples.trim_adapters = trim_adapters
Samples.run_porechop_abi = run_porechop_abi
Samples.trim_reads = trim_reads
Samples.run_fastqc = run_fastqc
Samples.align_to_reference_minimap = align_to_reference_minimap
Samples.align_to_reference_vg = align_to_reference_vg
Samples.reassign_mapq = reassign_mapq
Samples.filter_reads = filter_reads
Samples.call_variants_deepvariant = call_variants_deepvariant
Samples.call_variants_clair3 = call_variants_clair3
Samples.call_structural_variants_pbsv = call_structural_variants_pbsv
Samples.call_structural_variants_sniffles = call_structural_variants_sniffles
Samples.genotype_tandem_repeats = genotype_tandem_repeats
Samples.phase_genotypes_hiphase = phase_genotypes_hiphase
Samples.merge_hiphase_vcfs = merge_hiphase_vcfs
Samples.phase_genotypes_whatshap = phase_genotypes_whatshap
Samples.phase_genotypes_longphase = phase_genotypes_longphase
Samples.merge_longphase_vcfs = merge_longphase_vcfs
Samples.parse_haploblocks = parse_haploblocks
Samples.evaluate_gene_haploblocks = evaluate_gene_haploblocks
Samples.filter_vcf = filter_vcf
Samples.run_vcf2fasta = run_vcf2fasta
Samples.parse_fastas = parse_fastas

def main():
	parser = argparse.ArgumentParser(description="Run HLA-Resolve")
	parser.add_argument("--input_bam", required=True, help="Path to the unmapped HiFi BAM file")
	parser.add_argument("--sample_name", required=True, help="Override the parsed sample name", default=None)
	parser.add_argument("--platform", choices=["pacbio", "ont"], required=True, help="Specify sequencing platform (pacbio, ont)")
	parser.add_argument("--output_dir", required=True, help="Output Directory", default=None)
	parser.add_argument("--threads", type=int, required=False, help="Number of threads to use", default=6)
	parser.add_argument("--read_group_string", required=False, help="Override the parsed read group string", default=None)
	args = parser.parse_args()

	print("=============================")
	print("         HLA-RESOLVE         ")
	print("=============================")

	# Check that all required tools are installed
	check_required_commands()
	start_time = time.time()
	sample = Samples(unmapped_bam=args.input_bam, sample_name=args.sample_name, platform =args.platform, output_dir=args.output_dir, threads=args.threads, read_group_string=args.read_group_string)

	if args.platform.upper() == "PACBIO":
		sample.convert_bam_to_fastq()
		sample.mark_duplicates_pbmarkdup()
		sample.run_fastqc(os.path.join(sample.fastq_rmdup_dir, sample.sample_ID + ".dedup.fastq.gz"))
		sample.trim_adapters()
		sample.run_fastqc(os.path.join(sample.fastq_rmdup_cutadapt_dir, sample.sample_ID + ".dedup.trimmed.fastq.gz"))
		sample.align_to_reference_minimap()
		sample.align_to_reference_vg()
		sample.reassign_mapq()

		chr6_reads = sample.filter_reads()
		if chr6_reads > min_reads_sample:
			sample.call_variants_deepvariant()
			sample.call_structural_variants_pbsv()
			sample.genotype_tandem_repeats()
			sample.phase_genotypes_hiphase()
			sample.merge_hiphase_vcfs()
		else:
			print("Insufficient reads for variant calling")
			print("Sample {sample_id} had {num_reads} reads!".format(sample_id = sample.sample_ID, num_reads = chr6_reads))

	elif args.platform.upper() == "ONT":
		# sample.convert_bam_to_fastq()
		# sample.run_porechop_abi()
		sample.trim_reads()
		sample.align_to_reference_minimap()
		sample.align_to_reference_vg()
		sample.reassign_mapq()
		
		chr6_reads = sample.filter_reads()
		if chr6_reads > min_reads_sample:
			sample.call_variants_clair3()
			sample.call_structural_variants_sniffles()
			sample.phase_genotypes_whatshap()
			sample.phase_genotypes_longphase()
			sample.merge_longphase_vcfs()
		else:
			print("Insufficient reads for variant calling")
			print("Sample {sample_id} had {num_reads} reads!".format(sample_id = sample.sample_ID, num_reads = chr6_reads))
			
	heterozygous_sites, haploblock_list = sample.parse_haploblocks()
	phased_genes = sample.evaluate_gene_haploblocks(heterozygous_sites, haploblock_list)
	sample.filter_vcf()
	for gene in phased_genes:
		if gene in genes_of_interest:
			sample.run_vcf2fasta(gene, "gene")
			sample.run_vcf2fasta(gene, "CDS")
	sample.parse_fastas()
	
	end_time = time.time()
	elapsed_time = end_time - start_time
	minutes, seconds = divmod(elapsed_time,60)
	print("Processed sampled in {}:{:.2f}!".format(int(minutes), seconds))

if __name__ == "__main__":
	main()
