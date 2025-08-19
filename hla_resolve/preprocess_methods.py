import os
import subprocess
import pysam
import gzip
import shutil
from Bio import SeqIO

# Absolute path to the directory containing this script
SCRIPT_DIR = os.path.dirname(os.path.abspath(__file__))

# Path to the data directory bundled with the script
DATA_DIR = os.path.join(SCRIPT_DIR, "data")

# Input file paths
# Pangenome graph reference info 
reference_gbz="/hb/scratch/ogarci12/hybridcapture_pangenome/ref/hprc-v1.0-mc-grch38-minaf.0.1.gbz"
ref_paths="/hb/scratch/ogarci12/hybridcapture_pangenome/ref/hprc-v1.0-mc-grch38-minaf.0.1.dict"

# Use reference fasta with no alternate contigs.
# reference_fasta = os.path.join(DATA_DIR, "reference/GCA_000001405.15_GRCh38_no_alt_analysis_set.fa")
# reference_fasta = os.path.join(DATA_DIR, "reference/GRCh38_primary_only.fa")
# Rename fasta headers
# sed 's/^>GRCh38\.chr/>chr/' /hb/scratch/ogarci12/hybridcapture_pangenome/ref/hprc-v1.0-mc-grch38-minaf.0.1.fa > /hb/groups/cornejo_lab/matt/hla_capture/input_data/reference/hprc-v1.0-chr-renamed.fa

# Reference to use after mapping to graph and surjecting to GRCh38
# reference_fasta = os.path.join(DATA_DIR, "reference/hprc-v1.0-chr-renamed.fa")

# Referencew with added Y scaffold!
reference_fasta = os.path.join(DATA_DIR, "reference/augmented_hg38.fa")

# DeepVariant sif file
deepvariant_sif = os.path.join(DATA_DIR, "deepvariant_sif/deepvariant.sif")

# Path to vg install
vg = "/hb/scratch/ogarci12/hybridcapture_pangenome/vg"

# GRCh38 tandem repeat mask file for pbsv
# Downloaded from https://github.com/PacificBiosciences/pbsv/blob/master/annotations/human_GRCh38_no_alt_analysis_set.trf.bed
tandem_repeat_bed = os.path.join(DATA_DIR, "repeats_bed/human_GRCh38_no_alt_analysis_set.trf.bed")

chr6_bed = os.path.join(DATA_DIR, "reference/chr6.bed")

mosdepth_regions_file = "/hb/scratch/mglasena/test_hla_resolve/hla_resolve/hla_resolve/hla_genes.bed"

# GRCh38 tandem repeat definition file for pbtrgt
# Downloaded from https://zenodo.org/records/8329210
# pbtrgt_repeat_file = os.path.join(DATA_DIR, "repeats_bed/polymorphic_repeats.hg38.bed")
pbtrgt_repeat_file = os.path.join(DATA_DIR, "repeats_bed/test_chr6_polymorphic_repeats.hg38.bed")

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

# Convert BAM file of unmapped HiFi (ccs) reads to FASTQ format for marking duplicates and trimming adapters
def convert_bam_to_fastq(sample):
	if sample.platform == "PACBIO":
		print("Converting HiFi ccs reads to fastq format using pbtk bam2fastq!")
		print(f"bam2fastq input file: {sample.input_file}")
		
		os.chdir(sample.fastq_raw_dir)
		
		bam2fastq_cmd = f"bam2fastq -j {sample.threads} {sample.input_file} -o {sample.sample_ID}"
		
		subprocess.run(bam2fastq_cmd, shell=True, check=True)
				
		print(f"Raw fastq reads written to: {os.path.join(sample.fastq_raw_dir, sample.sample_ID + '.fastq.gz')}")
		print("\n\n")

		os.chdir(sample.ORIGINAL_CWD)
	
	elif sample.platform == "ONT":
		print("Converting ONT raw reads to fastq format using Samtools fastq!")
		print(f"Samtools fastq input file: {sample.input_file}")

		output_fastq = os.path.join(sample.fastq_raw_dir, sample.sample_ID + ".fastq")

		samtools_fq_cmd = f"samtools fastq -@ {sample.threads} {sample.input_file} > {output_fastq}"
		
		subprocess.run(samtools_fq_cmd, shell=True, check=True)

# Adapter/barcode trimming for ONT data
def run_porechop_abi(sample):
	print("Removing adapters and barcodes with porechop_abi!")
	
	input_fastq = os.path.join(sample.fastq_raw_dir, sample.sample_ID + ".fastq")
	
	print(f"porechop_abi input file: {input_fastq}")

	output_fastq = os.path.join(sample.fastq_porechop_dir, sample.sample_ID + ".porechop.fastq")
	
	porechop_cmd = f"porechop_abi --ab_initio -v 1 -i {input_fastq} -t {sample.threads} -o {output_fastq} --format fastq"

	subprocess.run(porechop_cmd, shell=True, check=True)

# Quality trimming step for ONT data 
def trim_reads(sample):
	print("Trimiming reads with ProwlerTrimmer!")
	
	input_fastq_dir = sample.fastq_porechop_dir + "/"
	input_fastq_file = sample.sample_ID + ".porechop.fastq"
	output_dir = sample.fastq_prowler_dir + "/"

	prowler_trimmer_cmd = f'python3 {prowler_trimmer} -i {input_fastq_dir} -f {input_fastq_file} -o {output_dir} -m "D" -q 20'

	subprocess.run(prowler_trimmer_cmd, shell=True, check=True, stdout=subprocess.DEVNULL, stderr=subprocess.DEVNULL)

	trimmed_fastq = os.path.join(sample.fastq_prowler_dir, sample.sample_ID + ".porechopTrimLT-U0-D20W100L100R0.fastq")
	renamed_trimmed_fastq = os.path.join(sample.fastq_prowler_dir, sample.sample_ID + ".porechop.prowler.fastq")
	os.rename(trimmed_fastq, renamed_trimmed_fastq)
	pigz_cmd = f"pigz -f -p {sample.threads} {renamed_trimmed_fastq}"
	subprocess.run(pigz_cmd, shell=True, check=True)

# Run fastqc on fastq file
def run_fastqc(self, fastqc_file):
	print(f"Running fastqc on {fastqc_file}")
	
	fastqc_cmd = f"fastqc --memory 5000 {fastqc_file}"
	
	subprocess.run(fastqc_cmd, shell=True, check=True)
	
	print("\n\n")

# Trim Adapter Sequences and polyA tails
# Customize this command based on the fastqc output
def legacy_trim_adapters(sample):
	print("Trimming adapter sequences with cutadapt!")
	
	input_file = os.path.join(sample.fastq_rmdup_dir, sample.sample_ID + ".pbmarkdup.fastq.gz")
	
	print(f"cutadapt input file: {input_file}")

	output_file = os.path.join(sample.fastq_rmdup_cutadapt_dir, sample.sample_ID + ".pbmarkdup.cutadapt.fastq.gz")
	
	# The following argument is used to trim polyA tails: -a 'A{{10}}N{{90}}. The double brackets are for python string interpolation
	# me and me_rc are the transposase mosaic end binding sequences defined at the top of the script
	#cutadapt_cmd = "cutadapt -j {threads} --quiet -n {num_cuts_allowed} -g {five_prime_adapter} -a {three_prime_adapter} -a 'A{{10}}N{{90}}' -o {output_fastq} {input_fastq}".format(threads = sample.threads, num_cuts_allowed = 3, five_prime_adapter = me, three_prime_adapter = me_rc, output_fastq = output_file, input_fastq = input_file)
	cutadapt_cmd = f"cutadapt -j {sample.threads} --quiet -n 2 -g {me} -a {me_rc} -o {output_file} {input_file}"

	subprocess.run(cutadapt_cmd, shell=True, check=True)
	
	print(f"Deduplicated trimmed reads written to: {output_file}")
	print("\n\n")

def trim_adapters(sample):
	input_file = os.path.join(sample.fastq_raw_dir, sample.sample_ID + ".fastq.gz")
	output_file = os.path.join(sample.fastq_trimmed_dir, sample.sample_ID + ".trimmed.fastq.gz")

	if sample.adapters:

		if sample.adapter_file:
			print("Trimming adapter sequences with cutadapt!")
			print(f"cutadapt input file: {input_file}")

			cutadapt_cmd = f"cutadapt -j {sample.threads} --minimum-length 21 --quiet -n 2 -g {sample.five_prime_adapter} -a {sample.three_prime_adapter} -o {output_file} {input_file}"

			subprocess.run(cutadapt_cmd, shell=True, check=True)
	
			print(f"Trimmed reads written to: {output_file}")
			print("\n\n")

		else:
			print("Trimming adapter sequences with fastplong in AUTO mode!")
			print(f"fastplong input file: {input_file}")

			html_path = os.path.join(sample.fastq_trimmed_dir, sample.sample_ID + ".fastplong.html")
			json_path = os.path.join(sample.fastq_trimmed_dir, sample.sample_ID + ".fastplong.json")

			fastplong_cmd = f"fastplong -i {input_file} -h {html_path} -j {json_path} -w {sample.threads} -n 100000 -o {output_file}"
			
			subprocess.run(fastplong_cmd, shell=True, check=False)

			print(f"Trimmed reads written to: {output_file}")

	else:
		print("Users specified no adapters present")
		print(f"Skipping adapter trimming and transfering raw fastq file {input_file} to {output_file}")
		shutil.copy(input_file, output_file)

# Mark PCR duplicates with pbmarkdup
def mark_duplicates_pbmarkdup(sample):
	print("Removing PCR duplicates using pbmarkdup!")
	
	input_fastq = os.path.join(sample.fastq_trimmed_dir, sample.sample_ID + ".trimmed.fastq.gz")
	
	print(f"pbmarkdup input file: {input_fastq}")

	output_fastq = os.path.join(sample.fastq_trimmed_dir, sample.sample_ID + ".trimmed.pbmarkdup.fastq")

	pbmarkdup_cmd = f"pbmarkdup -j {sample.threads} --rmdup {input_fastq} {output_fastq}"
	
	subprocess.run(pbmarkdup_cmd, shell=True, check=True)

	gzip_cmd = f"pigz -f -p {sample.threads} {output_fastq}"
	subprocess.run(gzip_cmd, shell=True, check=True)
	
	print(f"De-duplicated reads written to: {output_fastq}.gz!")
	print("\n\n")

# Align to GRCh38 reference genome with pbmm2
def align_to_reference_minimap(sample):
	print("Aligning reads to GRCh38 reference genome with minimap2!")
	
	if sample.platform == "PACBIO":
		input_fastq = os.path.join(sample.fastq_trimmed_dir, sample.sample_ID + ".trimmed.pbmarkdup.fastq.gz")
		platform_string = "map-hifi"
	elif sample.platform == "ONT":
		input_fastq = os.path.join(sample.fastq_trimmed_dir, sample.sample_ID + ".trimmed.fastq.gz")
		platform_string = "map-ont"

	output_bam = os.path.join(sample.mapped_bam_dir, sample.sample_ID + ".hg38.bam")

	print(f"minimap2 input file: {input_fastq}")

	# Deprecate pbmm2. Moving forward with minimap2 for ease of comparison with ONT data.
	# print("Aligning reads to GRCh38 reference genome with pbmm2!")
	# print("pbmm2 input file: {}".format(input_fastq))
	# pbmm2_cmd = "pbmm2 align -j {threads} {reference_genome} {input_file} {output_file} --sort --log-level INFO --unmapped --bam-index BAI --rg '{rg_string}'".format(threads = sample.threads, reference_genome = reference_fasta, input_file = input_fastq, output_file = output_bam, rg_string = sample.read_group_string)
	# subprocess.run(pbmm2_cmd, shell=True, check=True)

	minimap_threads = int(sample.threads * 2 / 3)
	samtools_threads = sample.threads - minimap_threads
	minimap_rg_string = "'{}'".format(sample.read_group_string.replace("\t", "\\t"))

	minimap2_cmd = f"minimap2 -Y -t {minimap_threads} -ax {platform_string} {reference_fasta} {input_fastq} -R {minimap_rg_string} | samtools sort -@ {samtools_threads} -o {output_bam}"
	index_bam = f"samtools index {output_bam}"
	
	subprocess.run(minimap2_cmd, shell=True, check=True)
	subprocess.run(index_bam, shell=True, check=True)

	print(f"Mapped bam written to: {output_bam}")
	print("\n\n")

def align_to_reference_vg(sample):
	print("Aligning reads to pangenome reference genome with vg giraffe!")
	
	if sample.platform == "PACBIO":
		input_fastq = os.path.join(sample.fastq_trimmed_dir, sample.sample_ID + ".trimmed.pbmarkdup.fastq.gz")
		parameter_preset = "hifi"

	elif sample.platform == "ONT":
		input_fastq = os.path.join(sample.fastq_prowler_dir, sample.sample_ID + ".trimmed.fastq.gz")
		parameter_preset = "r10"
	
	output_bam = os.path.join(sample.mapped_bam_dir, sample.sample_ID + ".pangenome.bam")
		
	print(f"vg giraffe input file: {input_fastq}")
	
	vg_threads = int(sample.threads * 2 / 3)
	samtools_threads = sample.threads - vg_threads

	rg_fields = dict(field.split(":", 1) for field in sample.read_group_string.split("\\t")[1:])
	read_group_id = rg_fields["ID"]

	vg_command = f"{vg} giraffe -b {parameter_preset} -Z {reference_gbz} -f {input_fastq} -p -P -o BAM --threads {vg_threads} --ref-paths {ref_paths} -R {read_group_id} -N {sample.sample_ID} | samtools sort -@ {samtools_threads} -o {output_bam}"
	index_bam = f"samtools index {output_bam}"

	subprocess.run(vg_command, shell=True, check=True)
	subprocess.run(index_bam, shell=True, check=True)

	# Reheader
	temp_sam_1 = os.path.join(sample.mapped_bam_dir, sample.sample_ID + "_temp1.sam")
	temp_sam_2 = os.path.join(sample.mapped_bam_dir, sample.sample_ID + "_temp2.sam")
	# converted_bam = output_bam.replace(".bam", ".reaheader.bam")
	converted_bam = output_bam.replace(".pangenome.bam", ".pg.bam")
	convert_to_sam = f"samtools view -h {output_bam} > {temp_sam_1}"
	rename_records = f"sed \"s/\\<GRCh38\\.chr/chr/g\" {temp_sam_1} > {temp_sam_2}"
	convert_to_bam = f"samtools view -b -o {converted_bam} {temp_sam_2}"
	index_new_bam = f"samtools index {converted_bam}"
	subprocess.run(convert_to_sam, shell=True, check=True)
	subprocess.run(rename_records, shell=True, check=True)
	subprocess.run(convert_to_bam, shell=True, check=True)
	subprocess.run(index_new_bam, shell=True, check=True)
	
	# Add read group 
	# final_bam = output_bam.replace(".pangenome.bam", ".pg.bam")
	
	# rg_fields = dict(field.split(":", 1) for field in sample.read_group_string.split("\t")[1:])
	# rg_id = rg_fields["ID"]

	# if rg_fields.get("SM") and rg_fields["SM"] != sample.sample_ID:
	# 	print(f"[WARNING] Overriding SM: {rg_fields['SM']} → {sample.sample_ID}")
	
	# add_rg_cmd = "samtools addreplacerg -r ID:{rg_id} -r SM:{SM} -o {final_bam} {converted_bam}".format(rg_id = rg_id, SM = sample.sample_ID, final_bam = final_bam, converted_bam = converted_bam)
	
	# subprocess.run(add_rg_cmd, shell=True, check=True)
	# index_final_bam = "samtools index {input_bam}".format(input_bam = final_bam)
	# subprocess.run(index_final_bam, shell=True, check=True)

	clean_up = f"rm {output_bam} {temp_sam_1} {temp_sam_2}"
	subprocess.run(clean_up, shell=True, check=True)

def reassign_mapq(sample):
	mapq_dict = dict()

	bam_hg38 = os.path.join(sample.mapped_bam_dir, sample.sample_ID + ".hg38.bam")
	bam_pg = os.path.join(sample.mapped_bam_dir, sample.sample_ID + ".pg.bam")
	reassigned_pg = os.path.join(sample.mapped_bam_dir, sample.sample_ID + ".pg.mapq_reassign.bam")

	with pysam.AlignmentFile(bam_hg38, "rb") as f:
		for read in f:
			mapq_dict[read.query_name] = read.mapping_quality

	missing_reads = []

	with pysam.AlignmentFile(bam_pg) as inbam, pysam.AlignmentFile(reassigned_pg, "wb", template=inbam) as outbam:
		for read in inbam:
			new_mapq = mapq_dict.get(read.query_name)
			if new_mapq is not None:
				read.mapping_quality = new_mapq
			else:
				missing_reads.append(read.query_name)
			outbam.write(read)

	index_bam = f"samtools index {reassigned_pg}"
	subprocess.run(index_bam, shell=True, check=True)

	if missing_reads:
		print(f"{len(missing_reads)} reads in {bam_pg} were missing from {bam_hg38}")

# Only marking duplicates for ONT data becuase pbmarkdup was used earlier for PacBio data
def mark_duplicates_picard(sample):
	if sample.aligner == "minimap2":
		input_bam = os.path.join(sample.mapped_bam_dir, sample.sample_ID + ".hg38.bam")
		output_bam = os.path.join(sample.mapped_bam_dir, sample.sample_ID + ".hg38.mrkdup.bam")
		metrics_file = os.path.join(sample.mapped_bam_dir, sample.sample_ID + ".hg38.mrkdup.metrics.txt")

	elif sample.aligner == "vg":
		input_bam = os.path.join(sample.mapped_bam_dir, sample.sample_ID + ".pg.mapq_reassign.bam")
		output_bam = os.path.join(sample.mapped_bam_dir, sample.sample_ID + ".pg.mapq_reassign.mrkdup.bam")
		metrics_file = os.path.join(sample.mapped_bam_dir, sample.sample_ID + ".pg.mapq_reassign.mrkdup.metrics.txt")
	
	temp_dir = os.path.join(sample.mapped_bam_dir, "mark_duplicates")
	os.makedirs(temp_dir, exist_ok=True)
	
	mark_duplicates_cmd = f"gatk MarkDuplicates -I {input_bam} -O {output_bam} --TMP_DIR {temp_dir} -M {metrics_file} --CREATE_INDEX true --VALIDATION_STRINGENCY LENIENT"
	
	subprocess.run(mark_duplicates_cmd, shell=True, check=True)

# Filer reads that did not map to chromosome 6
def filter_reads(sample):
	print("Excluding BAM records that don't map to chromosome 6!")

	output_bam = os.path.join(sample.mapped_bam_dir, sample.sample_ID +".hg38.rmdup.chr6.bam")

	if sample.aligner == "minimap2":
		if sample.platform == "PACBIO":
			input_bam = os.path.join(sample.mapped_bam_dir, sample.sample_ID + ".hg38.bam")
		elif sample.platform == "ONT":
			input_bam = os.path.join(sample.mapped_bam_dir, sample.sample_ID + ".hg38.mrkdup.bam")

	elif sample.aligner == "vg":
		if sample.platform == "PACBIO":
			input_bam = os.path.join(sample.mapped_bam_dir, sample.sample_ID + ".pg.mapq_reassign.bam")
		elif sample.platform == "ONT":
			input_bam = os.path.join(sample.mapped_bam_dir, sample.sample_ID + ".pg.mapq_reassign.mrkdup.bam")

	print(f"Samtools input file: {input_bam}")

	# Extract chromosome 6 and exclude secondary and supplementary alignments
	# samtools_cmd = "samtools view -@ {threads} -F 2304 -b {input_file} chr6 > '{output_file}'".format(threads = sample.threads, input_file = input_bam, output_file = output_bam)
	samtools_cmd = f"samtools view -F 1024 -@ {sample.threads} -b {input_bam} chr6 > '{output_bam}'"
	index_cmd = f"samtools index {output_bam}"

	subprocess.run(samtools_cmd, shell=True, check=True)
	subprocess.run(index_cmd, shell=True, check=True)

	count_reads_cmd = f"samtools view -c {output_bam}"

	read_count = int(subprocess.check_output(count_reads_cmd, shell=True).strip())
	
	print(f"Filtered BAM records written to: {output_bam}")
	print("\n\n")

	return read_count

def run_mosdepth(sample):
	input_file = os.path.join(sample.mapped_bam_dir, sample.sample_ID +".hg38.rmdup.chr6.bam")

	print(f"Running mosdepth on {input_file}")
	
	prefix = os.path.join(sample.mosdepth_dir, sample.sample_ID)
	
	# --flag 3328 excludes duplicates and secondary/supplementary alignments
	mosdepth = f"mosdepth --flag 3328 --by {mosdepth_regions_file} --thresholds 20,30 -t {sample.threads} {prefix} {input_file}"
	
	subprocess.run(mosdepth, shell=True, check=True)
	
	print("\n\n")

def parse_mosdepth(sample):
	print("Parsing mosdepth results!")
	print("\n\n")
	coverage_dict = dict()
	regions_file = os.path.join(sample.mosdepth_dir, sample.sample_ID + ".regions.bed.gz")
	thresholds_file = os.path.join(sample.mosdepth_dir, sample.sample_ID + ".thresholds.bed.gz")
	
	print("Gene,mean_depth,prop_20x,prop_30x")

	with gzip.open(regions_file, "rt") as f1, gzip.open(thresholds_file,"rt") as f2:
		regions = f1.read().splitlines()
		thresholds = f2.read().splitlines()[1:]

		for regions_line, thresholds_line in zip(regions,thresholds):
			regions_fields = regions_line.split("\t")
			gene = regions_fields[3].split("_")[0]
			start = regions_fields[1]
			stop = regions_fields[2]
			length = int(stop) - int(start)
			coverage_depth = float(regions_fields[4])
			threshold_fields = thresholds_line.split("\t")
			num_20x = int(threshold_fields[4])
			num_30x = int(threshold_fields[5])
			prop_20x = num_20x / length
			prop_30x = num_30x / length

			print(gene, coverage_depth, prop_20x, prop_30x)

			if coverage_depth >= depth_thresh and prop_20x >= prop_20x_thresh and prop_30x >= prop_30x_thresh:
				sample.sufficient_coverage_genes.append(gene)
			else:
				print(f"Gene {gene} has insufficient coverage for haplotyping and star allele calling")
	print("\n\n")

# Call SNV with DeepVariant
def call_variants_deepvariant(sample):
	input_file = os.path.join(sample.mapped_bam_dir, sample.sample_ID +".hg38.rmdup.chr6.bam")

	if sample.platform == "PACBIO":
		model_type = "PACBIO"
	elif sample.platform == "ONT":
		model_type = "ONT_R104"
	
	print("Calling SNVs and small INDELS with DeepVariant!")
	print(f"DeepVariant input file: {input_file}")
	
	output_vcf = os.path.join(sample.genotypes_dir, sample.sample_ID + ".vcf.gz")
	output_gvcf = os.path.join(sample.genotypes_dir, sample.sample_ID + ".g.vcf.gz")

	bind_paths = [
		f"{sample.genotypes_dir}:/data",
		f"{sample.mapped_bam_dir}:/input",
		f"{os.path.dirname(reference_fasta)}:/reference"
	]

	bind_flags = " ".join(f"--bind {path}" for path in bind_paths)

	deepvariant_cmd = f"""
		singularity exec {bind_flags} {deepvariant_sif} /opt/deepvariant/bin/run_deepvariant \
			--model_type={model_type} \
			--ref=/reference/{os.path.basename(reference_fasta)} \
			--reads=/input/{sample.sample_ID}.hg38.rmdup.chr6.bam \
			--output_vcf=/data/{sample.sample_ID}.vcf.gz \
			--output_gvcf=/data/{sample.sample_ID}.g.vcf.gz \
			--regions chr6 \
			--num_shards=8
		"""

	# Log DeepVariant in own output file so it doesn't clog up STDOUT
	deepvariant_log = os.path.join(sample.genotypes_dir, sample.sample_ID + ".deepvariant.log")
	print(f"Writing stdout to {deepvariant_log}")

	with open(deepvariant_log, "w") as log_file:
		subprocess.run(deepvariant_cmd, shell=True, check=True, stdout=log_file, stderr=log_file)

	print(f"VCF written to {output_vcf}")
	print(f"GVCF written to {output_gvcf}")
	print("\n\n")

def call_variants_clair3(sample):
	input_file = os.path.join(sample.mapped_bam_dir, sample.sample_ID +".hg38.rmdup.chr6.bam")

	if sample.platform == "ONT":
		platform = "ont"
		clair3_model_path = clair3_ont_model_path
	elif sample.platform == "PACBIO":
		platform = "hifi"
		clair3_model_path = clair3_hifi_model_path

	print("Calling SNVs and small INDELS with Clair3!")

	output_dir = os.path.join(sample.genotypes_dir, sample.sample_ID)
	os.makedirs(output_dir, exist_ok=True)

	clair3_cmd = f"run_clair3.sh --bam_fn={input_file} --ref_fn={reference_fasta} --platform={platform} --model_path={clair3_model_path} --output={output_dir} --threads={sample.threads} --sample_name={sample.sample_ID} --bed_fn={chr6_bed}"
	
	subprocess.run(clair3_cmd, shell=True, check=True)

	raw_genotypes_file = os.path.join(sample.genotypes_dir, sample.sample_ID, "merge_output.vcf.gz")
	relocated_genotypes_file = os.path.join(sample.genotypes_dir, sample.sample_ID + ".vcf.gz")
	shutil.copy(raw_genotypes_file, relocated_genotypes_file)
	subprocess.run(f"tabix -p vcf {relocated_genotypes_file}", shell=True, check=True)
	
	print(f"VCF written to {relocated_genotypes_file}")
	print("\n\n")

# Call SNV with bcftools
def call_variants_bcftools(sample):
	input_file = os.path.join(sample.mapped_bam_dir, sample.sample_ID +".hg38.rmdup.chr6.bam")

	if sample.platform == "PACBIO":
		config = "pacbio-ccs-1.20"
	elif sample.platform == "ONT":
		config = "ont-sup-1.20"

	print("Calling SNVs and small INDELS with bcftools!")

	print(f"Bcftools input file: {input_file}")
	
	output_vcf = os.path.join(sample.genotypes_dir, sample.sample_ID + ".vcf.gz")

	pileup_threads = str(sample.threads // 2)
	call_threads = str(sample.threads // 2)

	bcftools_command = (
		f"bcftools mpileup --config {config} --threads {pileup_threads} "
		f"-f {reference_fasta} -d 1000000 -r chr6:28000000-34000000 "
		f"-a FORMAT/DP,AD,ADF,ADR,SP {input_file} | "
		f"bcftools call -mv -f GQ --threads {call_threads} -Ou | "
		f"bcftools view -i '(TYPE=\"snp\" && GQ>=20 && QUAL>=10) || (TYPE=\"indel\" && GQ>=10)' "
		f"-Oz -o {output_vcf}")		

	subprocess.run(bcftools_command, shell=True, check=True)
	subprocess.run(f"tabix -p vcf {output_vcf}", shell=True, check=True)

	print(f"VCF written to {output_vcf}")
	print("\n\n")

# Run pbsv to call structural variants (SV)
def call_structural_variants_pbsv(sample):
	print("Calling structural variants with pbsv!")

	input_bam = os.path.join(sample.mapped_bam_dir, sample.sample_ID +".hg38.rmdup.chr6.bam")

	print(f"pbsv input file: {input_bam}")
	
	output_svsig = os.path.join(sample.sv_dir, sample.sample_ID + ".svsig.gz")
	output_vcf = os.path.join(sample.sv_dir, sample.sample_ID + ".SV.vcf")

	# -a Don't downsample
	# -q Don't filter by MapQ
	pbsv_discover_cmd = f"pbsv discover -a 0 -q 0 --region chr6 --tandem-repeats {tandem_repeat_bed} {input_bam} {output_svsig}"
	
	subprocess.run(pbsv_discover_cmd, shell=True, check=True)

	index_svsig_cmd = f"tabix -c '#' -s 3 -b 4 -e 4 {output_svsig}"
	
	subprocess.run(index_svsig_cmd, shell=True, check=True)

	pbsv_call_cmd = f"pbsv call -j {sample.threads} --min-sv-length 51 --region chr6 --hifi {reference_fasta} {output_svsig} {output_vcf}"
	
	subprocess.run(pbsv_call_cmd, shell=True, check=True)

	compress_cmd = f"bgzip -c {output_vcf} > {output_vcf}.gz"
	index_vcf_cmd = f"tabix -p vcf {output_vcf}.gz"
	
	subprocess.run(compress_cmd, shell=True, check=True)
	subprocess.run(index_vcf_cmd, shell=True, check=True)

	print(f"pbsv SV VCF written to: {output_vcf}")
	print("\n\n")

# Run sawfish to call structural variants (SV)
def call_structural_variants_sawfish(sample):
	print("Calling structural variants with sawfish!")

	# Hardcode threads because discover needs 8GB RAM per thread
	# 4 threads is plenty for target capture, sawfish is wicked fast
	threads = 4

	input_bam = os.path.join(sample.mapped_bam_dir, sample.sample_ID +".hg38.rmdup.chr6.bam")
	small_variant_calls = os.path.join(sample.genotypes_dir, sample.sample_ID + ".vcf.gz")
	discover_dir = os.path.join(sample.sv_dir, "discover")

	print(f"Sawfish input file: {input_bam}")

	sawfish_discover_cmd = f"{sawfish} discover --threads {threads} --ref {reference_fasta} --bam {input_bam} --disable-cnv --maf {small_variant_calls} --output-dir {discover_dir} --min-indel-size 51 --min-sv-mapq 0"
	
	# subprocess.run(sawfish_discover_cmd, shell=True, check=True, stdout=subprocess.DEVNULL, stderr=subprocess.DEVNULL)
	subprocess.run(sawfish_discover_cmd, shell=True, check=True)
	
	sawfish_call_cmd = f"{sawfish} joint-call --threads {threads} --sample {discover_dir} --output-dir {sample.sv_dir} --min-sv-mapq 0"
	
	# subprocess.run(sawfish_call_cmd, shell=True, check=True, stdout=subprocess.DEVNULL, stderr=subprocess.DEVNULL)
	subprocess.run(sawfish_call_cmd, shell=True, check=True)

	sawfish_output_file = os.path.join(sample.sv_dir, "genotyped.sv.vcf.gz")
	renamed_file = os.path.join(sample.sv_dir, sample.sample_ID + ".SV.vcf.gz")
	if os.path.exists(sawfish_output_file):
		shutil.move(sawfish_output_file, renamed_file)
		index_cmd = f"tabix -p vcf {renamed_file}"
		subprocess.run(index_cmd, shell=True, check=True)
		print(f"Sawfish SV VCF written to: {renamed_file}")
	else:
		print("Sawfish failed!")
	
	print("\n\n")

def call_structural_variants_sniffles(sample):
	print("Calling structural variants with Sniffles!")
	
	input_bam = os.path.join(sample.mapped_bam_dir, sample.sample_ID +".hg38.rmdup.chr6.bam")
	output_vcf = os.path.join(sample.sv_dir, sample.sample_ID + ".SV.vcf.gz")

	sniffles_cmd = f"sniffles --output-rnames --allow-overwrite -t {sample.threads} --reference {reference_fasta} --regions {chr6_bed} -i {input_bam} -v {output_vcf} --tandem-repeats {tandem_repeat_bed}"
	subprocess.run(sniffles_cmd, shell=True, check=True)

	print(f"Sniffles SV VCF written to: {output_vcf}")
	print("\n\n")

# Genotype tandem repeats with pbtrgt
def genotype_tandem_repeats(sample):
	print("Genotyping tandem repeats with pbtrgt!")

	input_bam = os.path.join(sample.mapped_bam_dir, sample.sample_ID +".hg38.rmdup.chr6.bam")

	print(f"trgt input file: {input_bam}")
	
	output_prefix = sample.sample_ID + ".TR"
	
	os.chdir(sample.pbtrgt_dir)
	
	trgt_cmd = f"trgt genotype --threads {sample.threads} --genome {reference_fasta} --reads {input_bam} --repeats {pbtrgt_repeat_file} --output-prefix {output_prefix} --preset targeted"
	
	subprocess.run(trgt_cmd, shell=True, check=True)
	
	sort_cmd = f"bcftools sort -O z -o {output_prefix + '.sorted.vcf.gz'} {output_prefix + '.vcf.gz'}"

	subprocess.run(sort_cmd, shell=True, check=True)

	os.rename(output_prefix + ".sorted.vcf.gz", output_prefix + ".vcf.gz")

	index_cmd = f"tabix -p vcf {output_prefix + '.vcf.gz'}"
	
	subprocess.run(index_cmd, shell=True, check=True)

	print(f"TR VCF written to {os.path.join(sample.pbtrgt_dir, output_prefix + '.vcf.gz')}")
	print("\n\n")

	os.chdir(sample.ORIGINAL_CWD)

# Phase genotypes with HiPhase
def phase_genotypes_hiphase(sample):
	print("Phasing Genotypes with HiPhase!")

	input_snv = os.path.join(sample.genotypes_dir, sample.sample_ID + ".vcf.gz")
	input_SV = os.path.join(sample.sv_dir, sample.sample_ID + ".SV.vcf.gz")
	input_TR = os.path.join(sample.pbtrgt_dir, sample.sample_ID + ".TR.vcf.gz")

	input_bam = os.path.join(sample.mapped_bam_dir, sample.sample_ID +".hg38.rmdup.chr6.bam")

	print(f"Input BAM: {input_bam}")
	print(f"Input SNV: {input_snv}")
	print(f"Input SV: {input_SV}")
	print(f"Input TR: {input_TR}")
	
	output_bam = os.path.join(sample.mapped_bam_dir, sample.sample_ID + ".hg38.rmdup.chr6.haplotag.bam")
	output_snv = os.path.join(sample.phased_vcf_dir, sample.sample_ID + ".hiphase.vcf.gz")
	output_SV = os.path.join(sample.phased_vcf_dir, sample.sample_ID + ".hiphase.SV.vcf.gz")
	output_TR = os.path.join(sample.phased_vcf_dir, sample.sample_ID + ".hiphase.TR.vcf.gz")

	output_summary_file = os.path.join(sample.phased_vcf_dir, sample.sample_ID + ".phased.summary.txt")
	output_blocks_file = os.path.join(sample.phased_vcf_dir, sample.sample_ID + ".phased.blocks.txt")
	output_stats_file = os.path.join(sample.phased_vcf_dir, sample.sample_ID + ".phased.stats.txt")

	hiphase_cmd = f"hiphase --threads {sample.threads} --ignore-read-groups --reference {reference_fasta} --bam {input_bam} --output-bam {output_bam} --vcf {input_snv} --output-vcf {output_snv} --vcf {input_SV} --output-vcf {output_SV} --vcf {input_TR} --output-vcf {output_TR} --stats-file {output_stats_file} --blocks-file {output_blocks_file} --summary-file {output_summary_file}"
	
	# Log HiPhase in own output file so it doesn't clog up STDOUT
	hiphase_log = os.path.join(sample.phased_vcf_dir, sample.sample_ID + ".hiphase.log")
	print(f"Writing stdout to {hiphase_log}")

	with open(hiphase_log, "w") as log_file:
		subprocess.run(hiphase_cmd, shell=True, check=True, stdout=log_file, stderr=log_file)

	print(f"HiPhase phased SNV VCF: {output_snv}")
	print(f"HiPhase phased SV VCF: {output_SV}")
	print(f"HiPhase phased TR VCF: {output_TR}")
	print(f"HiPhase haplotagged BAM written to: {output_bam}")
	print(f"HiPhase phasing summary written to: {output_summary_file}")
	print(f"HiPhase phasing stats written to: {output_stats_file}")
	print(f"HiPhase phase blocks written to: {output_blocks_file}")
	print("\n\n")

# Merge phased SNV (DeepVariant), tandem repeat (TRGT), and structural variant (pbsv) VCFs with bcftools concat 
def merge_hiphase_vcfs(sample):
	print("Merging phased DeepVariant, pbsv, and TRGT VCF files!")

	input_snv = os.path.join(sample.phased_vcf_dir, sample.sample_ID + ".hiphase.vcf.gz")
	input_SV = os.path.join(sample.phased_vcf_dir, sample.sample_ID + ".hiphase.SV.vcf.gz")
	input_TR = os.path.join(sample.phased_vcf_dir, sample.sample_ID + ".hiphase.TR.vcf.gz")

	print(f"DeepVariant input file: {input_snv}")
	print(f"pbsv input file: {input_SV}")
	print(f"pbtrgt input file: {input_TR}")

	output_vcf = os.path.join(sample.phased_vcf_dir, sample.sample_ID + ".hiphase.joint.vcf.gz")
	
	# Copied directly from Holt et al. 2024 Supplementary Material Program 5. 
	concat_cmd = f"bcftools concat --allow-overlaps {input_snv} {input_SV} {input_TR} | grep -v 'chrX|chrY' | grep -v 'SVTYPE=BND|SVTYPE=INV|SVTYPE=DUP' | bcftools norm -d none --fasta-ref {reference_fasta} | bcftools sort | bgzip > {output_vcf}"
	
	subprocess.run(concat_cmd, shell=True, check=True)

	index_cmd = f"tabix {output_vcf}"
	
	subprocess.run(index_cmd, shell=True, check=True)

	print(f"Merged VCF written to: {output_vcf}")
	print("\n\n")

# Phase genotypes with WhatsHap
def phase_genotypes_whatshap(sample):
	print("Phasing Genotypes with WhatsHap!")

	input_bam = os.path.join(sample.mapped_bam_dir, sample.sample_ID +".hg38.rmdup.chr6.bam")
	input_snv = os.path.join(sample.genotypes_dir, sample.sample_ID + ".vcf.gz")
	haplotagged_bam = os.path.join(sample.mapped_bam_dir, sample.sample_ID + ".hg38.rumdup.chr6.whatshap.bam")
	phased_vcf = os.path.join(sample.phased_vcf_dir, sample.sample_ID + ".whatshap.vcf.gz")

	print(f"Input BAM: {input_bam}")
	print(f"Input VCF: {input_snv}")

	output_blocks_file = os.path.join(sample.whatshap_phased_vcf_dir, sample.sample_ID + ".phased.haploblocks.txt")
	output_gtf_file = os.path.join(sample.whatshap_phased_vcf_dir, sample.sample_ID + ".phased.haploblocks.gtf")

	whatshap_phase_cmd = f"whatshap phase --ignore-read-groups --output {phased_vcf} --reference {reference_fasta} {input_vcf} {input_bam}"

	index_cmd = f"bcftools index {phased_vcf}"
	tabix_cmd = f"tabix {phased_vcf}"

	whatshap_haplotag_cmd = f"whatshap haplotag --ignore-read-groups --output {haplotagged_bam} --reference {reference_fasta} {input_vcf} {input_bam}"

	whatshap_stats_cmd = f"whatshap stats --block-list={output_blocks_file} --gtf={output_gtf_file} {phased_vcf}"

	# Log WhatsHap in own output file so it doesn't clog up STDOUT
	whatshap_log = os.path.join(sample.whatshap_phased_vcf_dir, sample.sample_ID + ".whatshap.log")

	with open(whatshap_log, "w") as log_file:
		log_file.write("\n==== Running WhatsHap Phase ====\n")
		subprocess.run(whatshap_phase_cmd, shell=True, check=True, stdout=log_file, stderr=log_file)
		subprocess.run(index_cmd, shell=True, check=True, stdout=log_file, stderr=log_file)
		subprocess.run(tabix_cmd, shell=True, check=True, stdout=log_file, stderr=log_file)
		log_file.write("\n==== Running WhatsHap Haplotag ====\n")
		subprocess.run(whatshap_haplotag_cmd, shell=True, check=True, stdout=log_file, stderr=log_file)
		log_file.write("\n==== Running WhatsHap Stats ====\n")
		subprocess.run(whatshap_stats_cmd, shell=True, check=True, stdout=log_file, stderr=log_file)

	print(f"WhatsHap phased VCF written to: {phased_vcf}")
	print(f"WhatsHap haplotagged BAM written to: {haplotagged_bam}")
	print(f"WhatsHap phase block gtf written to: {output_gtf_file}")
	print(f"WhatsHap phase blocks written to: {output_blocks_file}")
	print("\n\n")

def phase_genotypes_longphase(sample):
	print("Phasing Genotypes with LongPhase!")

	input_bam = os.path.join(sample.mapped_bam_dir, sample.sample_ID +".hg38.rmdup.chr6.bam")
	input_SNV_vcf = os.path.join(sample.genotypes_dir, sample.sample_ID + ".vcf.gz")
	input_SV_vcf = os.path.join(sample.sv_dir, sample.sample_ID + ".SV.vcf.gz")

	print(f"Input BAM: {input_bam}")
	print(f"Input SNV VCF: {input_SNV_vcf}")
	print(f"Input SV VCF: {input_SV_vcf}")


	output_blocks_file = os.path.join(sample.phased_vcf_dir, sample.sample_ID + ".phased.haploblocks.txt")
	output_gtf_file = os.path.join(sample.phased_vcf_dir, sample.sample_ID + ".phased.haploblocks.gtf")

	phased_vcf = os.path.join(sample.phased_vcf_dir, sample.sample_ID + ".longphase.vcf.gz")
	phased_vcf_prefix = phased_vcf.split(".vcf.gz")[0]
	longphase_phase_cmd = f"{longphase} phase -s {input_SNV_vcf} --sv {input_SV_vcf} -b {input_bam} -r {reference_fasta} -t {sample.threads} -o {phased_vcf_prefix} --ont"

	# Compress and index SNV VCF
	compress_cmd = f"bgzip -f {phased_vcf_prefix}.vcf"
	index_cmd = f"bcftools index {phased_vcf_prefix + '.vcf.gz'}"
	tabix_cmd = f"tabix {phased_vcf_prefix + '.vcf.gz'}"

	# Compress and index SV VCF
	SV_prefix = phased_vcf_prefix + "_SV"
	phased_SV_vcf = os.path.join(sample.phased_vcf_dir, sample.sample_ID + ".longphase_SV.vcf.gz")
	compress_SV_cmd = f"bgzip -f {SV_prefix}.vcf"
	index_SV_cmd = f"bcftools index {phased_SV_vcf}"
	tabix_SV_cmd = f"tabix {phased_SV_vcf}"

	haplotagged_bam = os.path.join(sample.mapped_bam_dir, sample.sample_ID + ".hg38.rmdup.chr6.haplotag.bam")
	haplotagged_bam_prefix = haplotagged_bam.split(".bam")[0]
	longphase_haplotag_cmd = f"{longphase} haplotag -r {reference_fasta} -s {phased_vcf} --sv-file {input_SV_vcf} -b {input_bam} -t {sample.threads} -o {haplotagged_bam_prefix}"

	whatshap_stats_cmd = f"whatshap stats --block-list={output_blocks_file} --gtf={output_gtf_file} {phased_vcf}"

	# Log WhatsHap in own output file so it doesn't clog up STDOUT
	longphase_log = os.path.join(sample.phased_vcf_dir, sample.sample_ID + ".longphase.log")

	with open(longphase_log, "w") as log_file:
		log_file.write("\n==== Running LongPhase Phase ====\n")
		subprocess.run(longphase_phase_cmd, shell=True, check=True, stdout=log_file, stderr=log_file)
		subprocess.run(compress_cmd, shell=True, check=True, stdout=log_file, stderr=log_file)
		subprocess.run(index_cmd, shell=True, check=True, stdout=log_file, stderr=log_file)
		subprocess.run(tabix_cmd, shell=True, check=True, stdout=log_file, stderr=log_file)
		subprocess.run(compress_SV_cmd, shell=True, check=True, stdout=log_file, stderr=log_file)
		subprocess.run(index_SV_cmd, shell=True, check=True, stdout=log_file, stderr=log_file)
		subprocess.run(tabix_SV_cmd, shell=True, check=True, stdout=log_file, stderr=log_file)
		log_file.write("\n==== Running LongPhase Haplotag ====\n")
		subprocess.run(longphase_haplotag_cmd, shell=True, check=True, stdout=log_file, stderr=log_file)
		log_file.write("\n==== Running WhatsHap Stats ====\n")
		subprocess.run(whatshap_stats_cmd, shell=True, check=True, stdout=log_file, stderr=log_file)

	print(f"LongPhase phased VCF written to: {phased_vcf}")
	print(f"LongPhase haplotagged BAM written to: {haplotagged_bam}")
	print(f"LongPhase phase block gtf written to: {output_gtf_file}")
	print(f"LongPhase phase blocks written to: {output_blocks_file}")
	print("\n\n")

def merge_longphase_vcfs(sample):
	print("Merging longphase SNV and SV VCFs using bcftools...")

	phased_vcf = os.path.join(sample.phased_vcf_dir, sample.sample_ID + ".longphase.vcf.gz")
	phased_SV_vcf = os.path.join(sample.phased_vcf_dir, sample.sample_ID + ".longphase_SV.vcf.gz")
	merged_vcf = os.path.join(sample.phased_vcf_dir, sample.sample_ID + ".longphase.merged.vcf.gz")
	reheadered_SV_vcf = phased_SV_vcf.replace(".vcf.gz", ".reheader.vcf.gz")

	header_file = os.path.join(sample.phased_vcf_dir, sample.sample_ID + ".header.txt")
	merge_log = os.path.join(sample.phased_vcf_dir, sample.sample_ID + ".merge.log")

	with open(header_file, "w") as hf:
		hf.write(f"SAMPLE\t{sample.sample_ID}\n")

	reheader_cmd = f"bcftools reheader -s {header_file} {phased_SV_vcf} -o {reheadered_SV_vcf}"
	index_sv_cmd = f"bcftools index {reheadered_SV_vcf}"

	merge_cmd = (
		f"bcftools concat --allow-overlaps -a {phased_vcf} {reheadered_SV_vcf} | "
		f"bcftools norm -m -any -f {reference_fasta} | "
		f"bcftools sort -Oz -o {merged_vcf} -"
	)
	index_merged_cmd = f"bcftools index {merged_vcf}"

	with open(merge_log, "w") as log_file:
		log_file.write("==== Reheadering SV VCF ====\n")
		subprocess.run(reheader_cmd, shell=True, check=True, stdout=log_file, stderr=log_file)
		subprocess.run(index_sv_cmd, shell=True, check=True, stdout=log_file, stderr=log_file)

		log_file.write("\n==== Merging SNV and SV VCFs ====\n")
		subprocess.run(merge_cmd, shell=True, check=True, stdout=log_file, stderr=log_file)

		log_file.write("\n==== Indexing Merged VCF ====\n")
		subprocess.run(index_merged_cmd, shell=True, check=True, stdout=log_file, stderr=log_file)

	print(f" Merged VCF written to: {merged_vcf}")
	print("\n\n")
