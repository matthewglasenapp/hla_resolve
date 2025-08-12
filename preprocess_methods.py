import os
import subprocess
import pysam

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
reference_fasta = os.path.join(DATA_DIR, "reference/hprc-v1.0-chr-renamed.fa")

# Referencew with added Y scaffold!
# reference_fasta = os.path.join(DATA_DIR, "reference/augmented_hg38.fa")

# DeepVariant sif file
deepvariant_sif = os.path.join(DATA_DIR, "deepvariant_sif/deepvariant.sif")

# Path to vg install
vg = "/hb/scratch/ogarci12/hybridcapture_pangenome/vg"

# GRCh38 tandem repeat mask file for pbsv
# Downloaded from https://github.com/PacificBiosciences/pbsv/blob/master/annotations/human_GRCh38_no_alt_analysis_set.trf.bed
tandem_repeat_bed = os.path.join(DATA_DIR, "repeats_bed/human_GRCh38_no_alt_analysis_set.trf.bed")

chr6_bed = os.path.join(DATA_DIR, "reference/chr6.bed")

mosdepth_regions_file = "/hb/scratch/mglasena/hla_resolve/hla_genes.bed"

# GRCh38 tandem repeat definition file for pbtrgt
# Downloaded from https://zenodo.org/records/8329210
pbtrgt_repeat_file = os.path.join(DATA_DIR, "repeats_bed/polymorphic_repeats.hg38.bed")

# Transposase mosaic end binding sequence
# The TE sequence (and its reverse complement) introduced during tagmentation still needs to be removed
# Adapters and barcodes were removed by PacBio with lima
me = "AGATGTGTATAAGAGACAG"
me_rc = "CTGTCTCTTATACACATCT"

longphase = "/hb/home/mglasena/software/longphase/longphase_linux-x64"
prowler_trimmer = "/hb/home/mglasena/software/ProwlerTrimmer/TrimmerLarge.py"

clair3_model_path = "/hb/home/mglasena/.conda/envs/clair3/bin/models/r941_prom_sup_g5014"

# Convert BAM file of unmapped HiFi (ccs) reads to FASTQ format for marking duplicates and trimming adapters
def convert_bam_to_fastq(self):
	if self.platform == "PACBIO":
		print("Converting HiFi ccs reads to fastq format using pbtk bam2fastq!")
		print("bam2fastq input file: {}".format(self.input_file))
		
		os.chdir(self.fastq_raw_dir)
		
		bam2fastq_cmd = "bam2fastq -j {threads} {input_file} -o {output_prefix}".format(threads = self.threads, input_file = self.input_file, output_prefix = self.sample_ID)
		
		subprocess.run(bam2fastq_cmd, shell=True, check=True)
				
		print("Raw fastq reads written to: {}".format(os.path.join(self.fastq_raw_dir, self.sample_ID + ".fastq.gz")))
		print("\n\n")

		os.chdir(self.ORIGINAL_CWD)
	
	elif self.platform == "ONT":
		print("Converting ONT raw reads to fastq format using Samtools fastq!")
		print("Samtools fastq input file: {}".format(self.input_file))

		output_fastq = os.path.join(self.fastq_raw_dir, self.sample_ID + ".fastq")

		samtools_fq_cmd = "samtools fastq -@ {threads} {input_file} > {output_file}".format(threads = self.threads, input_file = self.input_file, output_file = output_fastq)
		
		subprocess.run(samtools_fq_cmd, shell=True, check=True)

# Adapter/barcode trimming for ONT data
def run_porechop_abi(self):
	print("Removing adapters and barcodes with porechop_abi!")
	
	input_fastq = os.path.join(self.fastq_raw_dir, self.sample_ID + ".fastq")
	
	print("porechop_abi input file: {}".format(input_fastq))

	output_fastq = os.path.join(self.fastq_porechop_dir, self.sample_ID + ".porechop.fastq")
	
	porechop_cmd = "porechop_abi --ab_initio -v 1 -i {input_file} -t {threads} -o {output_file} --format fastq".format(input_file = input_fastq, threads = self.threads, output_file = output_fastq)

	subprocess.run(porechop_cmd, shell=True, check=True)

# Quality trimming step for ONT data 
def trim_reads(self):
	print("Trimiming reads with ProwlerTrimmer!")
	
	input_fastq_dir = self.fastq_porechop_dir + "/"
	input_fastq_file = self.sample_ID + ".porechop.fastq"
	output_dir = self.fastq_prowler_dir + "/"

	prowler_trimmer_cmd = 'python3 {prowler_trimmer} -i {input_dir} -f {input_file} -o {output_dir} -m "D" -q 20'.format(prowler_trimmer = prowler_trimmer, input_dir = input_fastq_dir, input_file = input_fastq_file, output_dir = output_dir)

	subprocess.run(prowler_trimmer_cmd, shell=True, check=True, stdout=subprocess.DEVNULL, stderr=subprocess.DEVNULL)

	trimmed_fastq = os.path.join(self.fastq_prowler_dir, self.sample_ID + ".porechopTrimLT-U0-D20W100L100R0.fastq")
	pigz_cmd = "pigz -f -p {threads} {input_file}".format(threads = self.threads, input_file = trimmed_fastq)
	subprocess.run(pigz_cmd, shell=True, check=True)

# Mark PCR duplicates with pbmarkdup
def mark_duplicates_pbmarkdup(self):
	print("Removing PCR duplicates using pbmarkdup!")
	
	input_fastq = os.path.join(self.fastq_raw_dir, self.sample_ID + ".fastq.gz")
	
	print("pbmarkdup input file: {}".format(input_fastq))

	output_fastq = os.path.join(self.fastq_rmdup_dir, self.sample_ID + ".dedup.fastq")

	pbmarkdup_cmd = "pbmarkdup -j {threads} --rmdup {input_file} {output_file}".format(threads = self.threads, input_file = input_fastq, output_file = output_fastq)
	
	subprocess.run(pbmarkdup_cmd, shell=True, check=True)

	gzip_cmd = "pigz -f -p {threads} {input_file}".format(threads = self.threads, input_file = output_fastq)
	subprocess.run(gzip_cmd, shell=True, check=True)
	
	print("De-duplicated reads written to: {}!".format(output_fastq))
	print("\n\n")

# Run fastqc on fastq file
def run_fastqc(self, fastqc_file):
	print("Running fastqc on {}".format(fastqc_file))
	
	fastqc_cmd = "fastqc --memory {memory} {input_fastq}".format(memory = 5000, input_fastq = fastqc_file)
	
	subprocess.run(fastqc_cmd, shell=True, check=True)
	
	print("\n\n")

# Trim Adapter Sequences and polyA tails
# Customize this command based on the fastqc output
def trim_adapters(self):
	print("Trimming adapter sequences with cutadapt!")
	
	input_file = os.path.join(self.fastq_rmdup_dir, self.sample_ID + ".dedup.fastq.gz")
	
	print("cutadapt input file: {}".format(input_file))

	output_file = os.path.join(self.fastq_rmdup_cutadapt_dir, self.sample_ID + ".dedup.trimmed.fastq.gz")
	
	# The following argument is used to trim polyA tails: -a 'A{{10}}N{{90}}. The double brackets are for python string interpolation
	# me and me_rc are the transposase mosaic end binding sequences defined at the top of the script
	#cutadapt_cmd = "cutadapt -j {threads} --quiet -n {num_cuts_allowed} -g {five_prime_adapter} -a {three_prime_adapter} -a 'A{{10}}N{{90}}' -o {output_fastq} {input_fastq}".format(threads = self.threads, num_cuts_allowed = 3, five_prime_adapter = me, three_prime_adapter = me_rc, output_fastq = output_file, input_fastq = input_file)
	cutadapt_cmd = "cutadapt -j {threads} --quiet -n {num_cuts_allowed} -g {five_prime_adapter} -a {three_prime_adapter} -o {output_fastq} {input_fastq}".format(threads = self.threads, num_cuts_allowed = 2, five_prime_adapter = me, three_prime_adapter = me_rc, output_fastq = output_file, input_fastq = input_file)

	subprocess.run(cutadapt_cmd, shell=True, check=True)
	
	print("Deduplicated trimmed reads written to: {}".format(output_file))
	print("\n\n")

# Align to GRCh38 reference genome with pbmm2
def align_to_reference_minimap(self):
	print("Aligning reads to GRCh38 reference genome with minimap2!")
	
	if self.platform == "PACBIO":
		input_fastq = os.path.join(self.fastq_rmdup_cutadapt_dir, self.sample_ID + ".dedup.trimmed.fastq.gz")
		output_bam = os.path.join(self.mapped_bam_dir, self.sample_ID + ".dedup.trimmed.hg38.bam")
		platform_string = "map-hifi"
	elif self.platform == "ONT":
		input_fastq = os.path.join(self.fastq_prowler_dir, self.sample_ID + ".porechopTrimLT-U0-D20W100L100R0.fastq.gz")
		output_bam = os.path.join(self.mapped_bam_dir, self.sample_ID + ".porechop.trimmed.hg38.bam")
		platform_string = "map-ont"

	print("minimap2 input file: {}".format(input_fastq))

	# Deprecate pbmm2. Moving forward with minimap2 for ease of comparison with ONT data.
	# print("Aligning reads to GRCh38 reference genome with pbmm2!")
	# print("pbmm2 input file: {}".format(input_fastq))
	# pbmm2_cmd = "pbmm2 align -j {threads} {reference_genome} {input_file} {output_file} --sort --log-level INFO --unmapped --bam-index BAI --rg '{rg_string}'".format(threads = self.threads, reference_genome = reference_fasta, input_file = input_fastq, output_file = output_bam, rg_string = self.read_group_string)
	# subprocess.run(pbmm2_cmd, shell=True, check=True)

	minimap_threads = int(self.threads * 2 / 3)
	samtools_threads = self.threads - minimap_threads
	minimap_rg_string = "'{}'".format(self.read_group_string.replace("\t", "\\t"))

	minimap2_cmd = "minimap2 -t {minimap_threads} -ax {platform_string} {reference_genome} {input_file} -R {rg_string} | samtools sort -@ {samtools_threads} -o {output_file}".format(minimap_threads = minimap_threads, platform_string = platform_string, reference_genome = reference_fasta, input_file = input_fastq, rg_string = minimap_rg_string, samtools_threads = samtools_threads, output_file = output_bam)
	index_bam = "samtools index {input_file}".format(input_file = output_bam)
	
	subprocess.run(minimap2_cmd, shell=True, check=True)
	subprocess.run(index_bam, shell=True, check=True)

	print("Mapped bam written to: {}".format(output_bam))
	print("\n\n")

def mark_duplicates_picard(self):
	# input_bam = os.path.join(self.mapped_bam_dir, self.sample_ID + ".porechop.trimmed.hg38.bam")
	input_bam = os.path.join(self.mapped_bam_dir, self.sample_ID + ".porechop.trimmed.pg.mapq_reassign.bam")
	# output_bam = os.path.join(self.mapped_bam_dir, self.sample_ID + ".porechop.trimmed.hg38.mrkdup.bam")
	output_bam = os.path.join(self.mapped_bam_dir, self.sample_ID + ".porechop.trimmed.pg.mapq_reassign.mrkdup.bam")
	metrics_file = os.path.join(self.mapped_bam_dir, self.sample_ID + ".porechop.trimmed.pg.mapq_reassign.mrkdup.metrics.txt")
	temp_dir = os.path.join(self.mapped_bam_dir, "mark_duplicates")
	os.makedirs(temp_dir, exist_ok=True)
	
	mark_duplicates_cmd = "gatk MarkDuplicates -I {input_file} -O {output_file} --TMP_DIR {temp_dir} -M {metrics_file} --CREATE_INDEX true --VALIDATION_STRINGENCY LENIENT".format(input_file = input_bam, output_file = output_bam, temp_dir = temp_dir, metrics_file = metrics_file)
	
	subprocess.run(mark_duplicates_cmd, shell=True, check=True)

def align_to_reference_vg(self):
	print("Aligning reads to pangenome reference genome with vg giraffe!")
	
	if self.platform == "PACBIO":
		input_fastq = os.path.join(self.fastq_rmdup_cutadapt_dir, self.sample_ID + ".dedup.trimmed.fastq.gz")
		output_bam = os.path.join(self.mapped_bam_dir, self.sample_ID + ".dedup.trimmed.pangenome.bam")
		parameter_preset = "hifi"

	elif self.platform == "ONT":
		input_fastq = os.path.join(self.fastq_prowler_dir, self.sample_ID + ".porechopTrimLT-U0-D20W100L100R0.fastq.gz")
		output_bam = os.path.join(self.mapped_bam_dir, self.sample_ID + ".porechop.trimmed.pangenome.bam")
		parameter_preset = "r10"

	print("vg giraffe input file: {}".format(input_fastq))
	
	vg_threads = int(self.threads * 2 / 3)
	samtools_threads = self.threads - vg_threads

	rg_fields = dict(field.split(":", 1) for field in self.read_group_string.split("\\t")[1:])
	read_group_id = rg_fields["ID"]

	vg_command = "{vg} giraffe -b {platform} -Z {reference_gbz} -f {input_file} -p -P -o BAM --threads {vg_threads} --ref-paths {ref_paths} -R {read_group} -N {sample_name} | samtools sort -@ {samtools_threads} -o {output_file}".format(vg = vg, platform = parameter_preset, reference_gbz = reference_gbz, input_file = input_fastq, vg_threads = vg_threads, ref_paths = ref_paths, read_group = read_group_id, sample_name = self.sample_ID, samtools_threads = samtools_threads, output_file = output_bam)
	index_bam = "samtools index {input_file}".format(input_file = output_bam)

	subprocess.run(vg_command, shell=True, check=True)
	subprocess.run(index_bam, shell=True, check=True)

	# Reheader
	temp_sam_1 = os.path.join(self.mapped_bam_dir, self.sample_ID + "_temp1.sam")
	temp_sam_2 = os.path.join(self.mapped_bam_dir, self.sample_ID + "_temp2.sam")
	# converted_bam = output_bam.replace(".bam", ".reaheader.bam")
	converted_bam = output_bam.replace(".pangenome.bam", ".pg.bam")
	convert_to_sam = "samtools view -h {input_file} > {temp}".format(input_file = output_bam, temp = temp_sam_1)
	rename_records = "sed \"s/\\<GRCh38\\.chr/chr/g\" {temp1} > {temp2}".format(temp1=temp_sam_1, temp2=temp_sam_2)
	convert_to_bam = "samtools view -b -o {output_bam} {temp2}".format(output_bam = converted_bam, temp2 = temp_sam_2)
	index_new_bam = "samtools index {input_file}".format(input_file = converted_bam)
	subprocess.run(convert_to_sam, shell=True, check=True)
	subprocess.run(rename_records, shell=True, check=True)
	subprocess.run(convert_to_bam, shell=True, check=True)
	subprocess.run(index_new_bam, shell=True, check=True)
	
	# Add read group 
	# final_bam = output_bam.replace(".pangenome.bam", ".pg.bam")
	
	# rg_fields = dict(field.split(":", 1) for field in self.read_group_string.split("\t")[1:])
	# rg_id = rg_fields["ID"]

	# if rg_fields.get("SM") and rg_fields["SM"] != self.sample_ID:
	# 	print(f"[WARNING] Overriding SM: {rg_fields['SM']} → {self.sample_ID}")
	
	# add_rg_cmd = "samtools addreplacerg -r ID:{rg_id} -r SM:{SM} -o {final_bam} {converted_bam}".format(rg_id = rg_id, SM = self.sample_ID, final_bam = final_bam, converted_bam = converted_bam)
	
	# subprocess.run(add_rg_cmd, shell=True, check=True)
	# index_final_bam = "samtools index {input_bam}".format(input_bam = final_bam)
	# subprocess.run(index_final_bam, shell=True, check=True)

	clean_up = "rm {bam1} {temp1} {temp2}".format(bam1 = output_bam, temp1 = temp_sam_1, temp2 = temp_sam_2)
	subprocess.run(clean_up, shell=True, check=True)

def reassign_mapq(self):
	mapq_dict = dict()

	if self.platform == "PACBIO":
		bam_hg38 = os.path.join(self.mapped_bam_dir, self.sample_ID + ".dedup.trimmed.hg38.bam")
		bam_pg = os.path.join(self.mapped_bam_dir, self.sample_ID + ".dedup.trimmed.pg.bam")
		reassigned_pg = os.path.join(self.mapped_bam_dir, self.sample_ID + ".dedup.trimmed.pg.mapq_reassign.bam")
	elif self.platform == "ONT":
		bam_hg38 = os.path.join(self.mapped_bam_dir, self.sample_ID + ".porechop.trimmed.hg38.bam")
		bam_pg = os.path.join(self.mapped_bam_dir, self.sample_ID + ".porechop.trimmed.pg.bam")
		reassigned_pg = os.path.join(self.mapped_bam_dir, self.sample_ID + ".porechop.trimmed.pg.mapq_reassign.bam")

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

	index_bam = "samtools index {input_bam}".format(input_bam = reassigned_pg)
	subprocess.run(index_bam, shell=True, check=True)

	if missing_reads:
		print("{Count} reads in {pg_bam} were missing from {hg38_bam}".format(Count = len(missing_reads), pg_bam = bam_pg, hg38_bam = bam_hg38))

# Filer reads that did not map to chromosome 6
def filter_reads(self):
	print("Excluding BAM records that don't map to chromosome 6!")

	if self.platform == "PACBIO":
		input_bam = os.path.join(self.mapped_bam_dir, self.sample_ID + ".dedup.trimmed.pg.mapq_reassign.bam")
		output_bam = os.path.join(self.mapped_bam_dir, self.sample_ID + ".dedup.trimmed.hg38.chr6.bam")
	elif self.platform == "ONT":
		input_bam = os.path.join(self.mapped_bam_dir, self.sample_ID + ".porechop.trimmed.pg.mapq_reassign.mrkdup.bam")
		output_bam = os.path.join(self.mapped_bam_dir, self.sample_ID + ".porechop.trimmed.hg38.rmdup.chr6.bam")

	print("Samtools input file: {}".format(input_bam))

	# Extract chromosome 6 and exclude secondary and supplementary alignments
	#samtools_cmd = "samtools view -@ {threads} -F 2304 -b {input_file} chr6 > '{output_file}'".format(threads = self.threads, input_file = input_bam, output_file = output_bam)
	samtools_cmd = "samtools view -@ {threads} -F 2048 -b {input_file} chr6 > '{output_file}'".format(threads = self.threads, input_file = input_bam, output_file = output_bam)
	index_cmd = "samtools index {input_file}".format(input_file = output_bam)

	subprocess.run(samtools_cmd, shell=True, check=True)
	subprocess.run(index_cmd, shell=True, check=True)

	count_reads_cmd = "samtools view -c {input_file}".format(input_file = output_bam)

	read_count = int(subprocess.check_output(count_reads_cmd, shell=True).strip())
	
	print("Filtered BAM records written to: {}".format(output_bam))
	print("\n\n")

	return read_count

def run_mosdepth(self):
	input_file = os.path.join(self.mapped_bam_dir, self.sample_ID + ".dedup.trimmed.hg38.chr6.bam")
	prefix = os.path.join(self.mosdepth_dir, self.sample_ID)
	# --flag 3328 excludes duplicates and secondary/supplementary alignments
	mosdepth = "mosdepth --flag 3328 --by {regions_file} --thresholds 20,30 -t {threads} {prefix} {bam_file}".format(regions_file = mosdepth_regions_file, threads = self.threads, prefix = prefix, bam_file = input_file)
	subprocess.run(mosdepth, shell=True, check=True)

def parse_mosdepth(self):
	coverage_dict = dict()
	regions_file = os.path.join(self.mosdepth_dir, self.sample_ID + ".regions.bed.gz")
	thresholds_file = os.path.join(self.mosdepth_dir, self.sample_ID + "*.thresholds.bed.gz")
	
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

			if coverage_depth >= 30 and prop_20x >= 0.9 and prop_30x >= 0.9:
				self.sufficient_coverage_genes.append(gene)
			else:
				print("Gene {gene} has insufficient coverage for haplotyping and star allele calling".format(gene = gene))

# Call SNV with DeepVariant
def call_variants_deepvariant(self):
	print("Calling SNVs and small INDELS with DeepVariant!")

	input_bam = os.path.join(self.mapped_bam_dir, self.sample_ID + ".dedup.trimmed.hg38.chr6.bam")

	print("DeepVariant input file: {}".format(input_bam))
	
	output_vcf = os.path.join(self.deepvariant_dir, self.sample_ID + ".dedup.trimmed.hg38.chr6.SNV.vcf.gz")
	output_gvcf = os.path.join(self.deepvariant_dir, self.sample_ID + ".dedup.trimmed.hg38.chr6.SNV.g.vcf.gz")

	bind_paths = [
		f"{self.deepvariant_dir}:/data",
		f"{self.mapped_bam_dir}:/input",
		f"{os.path.dirname(reference_fasta)}:/reference"
	]

	bind_flags = " ".join("--bind {}".format(path) for path in bind_paths)

	deepvariant_cmd = """
		singularity exec {binds} {sif} /opt/deepvariant/bin/run_deepvariant \
			--model_type=PACBIO \
			--ref=/reference/{ref_filename} \
			--reads=/input/{sample}.dedup.trimmed.hg38.chr6.bam \
			--output_vcf=/data/{sample}.dedup.trimmed.hg38.chr6.SNV.vcf.gz \
			--output_gvcf=/data/{sample}.dedup.trimmed.hg38.chr6.SNV.g.vcf.gz \
			--regions chr6 \
			--num_shards=8
		""".format(
			binds=bind_flags,
			sif=deepvariant_sif,
			ref_filename=os.path.basename(reference_fasta),
			sample=self.sample_ID
			)

	# Log DeepVariant in own output file so it doesn't clog up STDOUT
	deepvariant_log = os.path.join(self.deepvariant_dir, self.sample_ID + ".deepvariant.log")

	with open(deepvariant_log, "w") as log_file:
		subprocess.run(deepvariant_cmd, shell=True, check=True, stdout=log_file, stderr=log_file)

	print("VCF written to {}".format(output_vcf))
	print("GVCF written to {}".format(output_gvcf))
	print("\n\n")

def call_variants_clair3(self):
	print("Calling SNVs and small INDELS with Clair3!")

	input_bam = os.path.join(self.mapped_bam_dir, self.sample_ID + ".porechop.trimmed.hg38.rmdup.chr6.bam")
	output_dir = os.path.join(self.clair3_dir, self.sample_ID)
	os.makedirs(output_dir, exist_ok=True)

	clair3_cmd = "run_clair3.sh --bam_fn={input_file} --ref_fn={reference_genome} --platform=ont --model_path={clair3_model} --output={output_dir} --threads={threads} --sample_name={sample} --bed_fn={bed_file}".format(input_file = input_bam, reference_genome = reference_fasta, clair3_model = clair3_model_path, output_dir = output_dir, threads = self.threads, sample = self.sample_ID, bed_file = chr6_bed)
	
	subprocess.run(clair3_cmd, shell=True, check=True)

	print("VCF written to {}".format(os.path.join(self.clair3_dir, self.sample_ID, "merge_output.vcf.gz")))
	print("\n\n")

# Call SNV with bcftools
def call_variants_bcftools(self):
	print("Calling SNVs and small INDELS with bcftools!")

	input_bam = os.path.join(self.mapped_bam_dir, self.sample_ID + ".dedup.trimmed.hg38.chr6.bam")

	print("Bcftools input file: {}".format(input_bam))
	
	output_vcf = os.path.join(self.deepvariant_dir, self.sample_ID + ".dedup.trimmed.hg38.chr6.SNV.vcf.gz")

	pileup_threads = str(self.threads // 2)
	call_threads = str(self.threads // 2)

	bcftools_command = (
		"bcftools mpileup --config pacbio-ccs-1.20 --threads {pileup_threads} "
		"-f {reference_genome} -d 1000000 -r chr6:28000000-34000000 "
		"-a FORMAT/DP,AD,ADF,ADR,SP {input_bam} | "
		"bcftools call -mv -f GQ --threads {call_threads} -Ou | "
		"bcftools view -i '(TYPE=\"snp\" && GQ>=20 && QUAL>=10) || (TYPE=\"indel\" && GQ>=10)' "
		"-Oz -o {output_vcf}").format(
		pileup_threads=pileup_threads,
		reference_genome=reference_fasta,
		input_bam=input_bam,
		call_threads=call_threads,
		output_vcf=output_vcf)		

	subprocess.run(bcftools_command, shell=True, check=True)
	subprocess.run(f"tabix -p vcf {output_vcf}", shell=True, check=True)

	print("VCF written to {}".format(output_vcf))
	print("\n\n")

# Run pbsv to call structural variants (SV)
def call_structural_variants_pbsv(self):
	print("Calling structural variants with pbsv!")

	input_bam = os.path.join(self.mapped_bam_dir, self.sample_ID + ".dedup.trimmed.hg38.chr6.bam")

	print("pbsv input file: {}".format(input_bam))
	
	output_svsig = os.path.join(self.pbsv_dir, self.sample_ID + ".dedup.trimmed.hg38.chr6.svsig.gz")
	output_vcf = os.path.join(self.pbsv_dir, self.sample_ID + ".dedup.trimmed.hg38.chr6.SV.vcf")

	# -a Don't downsample
	# -q Don't filter by MapQ
	pbsv_discover_cmd = "pbsv discover -a 0 -q 0 --region chr6 --tandem-repeats {tandem_repeat_file} {input_file} {output_file}".format(tandem_repeat_file = tandem_repeat_bed, input_file = input_bam, output_file = output_svsig)
	
	subprocess.run(pbsv_discover_cmd, shell=True, check=True)

	index_svsig_cmd = "tabix -c '#' -s 3 -b 4 -e 4 {svsig_file}".format(svsig_file = output_svsig)
	
	subprocess.run(index_svsig_cmd, shell=True, check=True)

	pbsv_call_cmd = "pbsv call -j {threads} --region chr6 --hifi {reference_genome} {input_file} {output_file}".format(threads = self.threads, reference_genome = reference_fasta, input_file = output_svsig, output_file = output_vcf)
	
	subprocess.run(pbsv_call_cmd, shell=True, check=True)

	compress_cmd = "bgzip -c {input_file} > {input_file}.gz".format(input_file = output_vcf)
	index_vcf_cmd = "tabix -p vcf {input_file}.gz".format(input_file = output_vcf)
	
	subprocess.run(compress_cmd, shell=True, check=True)
	subprocess.run(index_vcf_cmd, shell=True, check=True)

	print("pbsv SV VCF written to: {}".format(output_vcf))
	print("\n\n")

def call_structural_variants_sniffles(self):
	print("Calling structural variants with Sniffles!")
	
	input_bam = os.path.join(self.mapped_bam_dir, self.sample_ID + ".porechop.trimmed.hg38.rmdup.chr6.bam")
	output_vcf = os.path.join(self.sniffles_dir, self.sample_ID + ".porechop.trimmed.hg38.rmdup.chr6.sniffles.vcf.gz")

	sniffles_cmd = "sniffles --output-rnames --allow-overwrite -t {threads} --reference {reference_genome} --regions {bed_file} -i {input_bam} -v {output_vcf} --tandem-repeats {tandem_repeat_bed}".format(threads = self.threads, reference_genome = reference_fasta, bed_file = chr6_bed, input_bam = input_bam, output_vcf = output_vcf, tandem_repeat_bed = tandem_repeat_bed)
	subprocess.run(sniffles_cmd, shell=True, check=True)

	print("Sniffles SV VCF written to: {}".format(output_vcf))
	print("\n\n")

# Genotype tandem repeats with pbtrgt
def genotype_tandem_repeats(self):
	print("Genotyping tandem repeats with pbtrgt!")

	input_bam = os.path.join(self.mapped_bam_dir, self.sample_ID + ".dedup.trimmed.hg38.chr6.bam")

	print("trgt input file: {}".format(input_bam))
	
	output_prefix = self.sample_ID + ".dedup.trimmed.hg38.chr6.TR"
	
	os.chdir(self.pbtrgt_dir)
	
	trgt_cmd = "trgt genotype --threads {threads} --genome {reference_genome} --reads {input_file} --repeats {repeat_file} --output-prefix {output_prefix} --preset targeted".format(threads = self.threads, reference_genome = reference_fasta, input_file = input_bam, repeat_file = pbtrgt_repeat_file, output_prefix = output_prefix)
	
	subprocess.run(trgt_cmd, shell=True, check=True)
	
	sort_cmd = "bcftools sort -O z -o {output_file} {input_file}".format(output_file = output_prefix + ".sorted.vcf.gz", input_file = output_prefix + ".vcf.gz")

	subprocess.run(sort_cmd, shell=True, check=True)

	os.rename(output_prefix + ".sorted.vcf.gz", output_prefix + ".vcf.gz")

	index_cmd = "tabix -p vcf {input_file}".format(input_file = output_prefix + ".vcf.gz")
	
	subprocess.run(index_cmd, shell=True, check=True)

	print("TR VCF written to {}".format(os.path.join(self.pbtrgt_dir, output_prefix + ".vcf.gz")))
	print("\n\n")

	os.chdir(self.ORIGINAL_CWD)

# Phase genotypes with HiPhase
def phase_genotypes_hiphase(self):
	print("Phasing Genotypes with HiPhase!")

	input_snv = os.path.join(self.deepvariant_dir, self.sample_ID + ".dedup.trimmed.hg38.chr6.SNV.vcf.gz")
	input_SV = os.path.join(self.pbsv_dir, self.sample_ID + ".dedup.trimmed.hg38.chr6.SV.vcf.gz")
	input_TR = os.path.join(self.pbtrgt_dir, self.sample_ID + ".dedup.trimmed.hg38.chr6.TR.vcf.gz")

	input_bam = os.path.join(self.mapped_bam_dir, self.sample_ID + ".dedup.trimmed.hg38.chr6.bam")

	print("Input BAM: {}".format(input_bam))
	print("Input SNV: {}".format(input_snv))
	print("Input SV: {}".format(input_SV))
	print("Input TR: {}".format(input_TR))
	
	output_bam = os.path.join(self.mapped_bam_dir, self.sample_ID + ".dedup.trimmed.hg38.chr6.hiphase.haplotag.bam")
	output_snv = os.path.join(self.hiphase_phased_vcf_dir, self.sample_ID + ".dedup.trimmed.hg38.chr6.phased.SNV.vcf.gz")
	output_SV = os.path.join(self.hiphase_phased_vcf_dir, self.sample_ID + ".dedup.trimmed.hg38.chr6.phased.SV.vcf.gz")
	output_TR = os.path.join(self.hiphase_phased_vcf_dir, self.sample_ID + ".dedup.trimmed.hg38.chr6.phased.TR.vcf.gz")

	output_summary_file = os.path.join(self.hiphase_phased_vcf_dir, self.sample_ID + ".phased.summary.txt")
	output_blocks_file = os.path.join(self.hiphase_phased_vcf_dir, self.sample_ID + ".phased.blocks.txt")
	output_stats_file = os.path.join(self.hiphase_phased_vcf_dir, self.sample_ID + ".phased.stats.txt")

	hiphase_cmd = f"hiphase --threads {self.threads} --ignore-read-groups --reference {reference_fasta} --bam {input_bam} --output-bam {output_bam} --vcf {input_snv} --output-vcf {output_snv} --vcf {input_SV} --output-vcf {output_SV} --vcf {input_TR} --output-vcf {output_TR} --stats-file {output_stats_file} --blocks-file {output_blocks_file} --summary-file {output_summary_file}"
	
	# Log HiPhase in own output file so it doesn't clog up STDOUT
	hiphase_log = os.path.join(self.hiphase_phased_vcf_dir, self.sample_ID + ".hiphase.log")

	with open(hiphase_log, "w") as log_file:
		subprocess.run(hiphase_cmd, shell=True, check=True, stdout=log_file, stderr=log_file)

	print("HiPhase phased SNV VCF: {}".format(output_snv))
	print("HiPhase phased SV VCF: {}".format(output_SV))
	print("HiPhase phased TR VCF: {}".format(output_TR))
	print("HiPhase haplotagged BAM written to: {}".format(output_bam))
	print("HiPhase phasing summary written to: {}".format(output_summary_file))
	print("HiPhase phasing stats written to: {}".format(output_stats_file))
	print("HiPhase phase blocks written to: {}".format(output_blocks_file))
	print("\n\n")

# Merge phased SNV (DeepVariant), tandem repeat (TRGT), and structural variant (pbsv) VCFs with bcftools concat 
def merge_hiphase_vcfs(self):
	print("Merging phased DeepVariant, pbsv, and TRGT VCF files!")

	input_snv = os.path.join(self.hiphase_phased_vcf_dir, self.sample_ID + ".dedup.trimmed.hg38.chr6.phased.SNV.vcf.gz")
	input_SV = os.path.join(self.hiphase_phased_vcf_dir, self.sample_ID + ".dedup.trimmed.hg38.chr6.phased.SV.vcf.gz")
	input_TR = os.path.join(self.hiphase_phased_vcf_dir, self.sample_ID + ".dedup.trimmed.hg38.chr6.phased.TR.vcf.gz")

	print("DeepVariant input file: {}".format(input_snv))
	print("pbsv input file: {}".format(input_SV))
	print("pbtrgt input file: {}".format(input_TR))

	output_vcf = os.path.join(self.hiphase_phased_vcf_dir, self.sample_ID + ".dedup.trimmed.hg38.chr6.phased.joint.vcf.gz")
	
	# Copied directly from Holt et al. 2024 Supplementary Material Program 5. 
	concat_cmd = "bcftools concat --allow-overlaps {SNV_vcf} {SV_vcf} {TR_vcf} | grep -v 'chrX|chrY' | grep -v 'SVTYPE=BND|SVTYPE=INV|SVTYPE=DUP' | bcftools norm -d none --fasta-ref {reference_genome} | bcftools sort | bgzip > {output_file}".format(SNV_vcf = input_snv, SV_vcf = input_SV, TR_vcf = input_TR, reference_genome = reference_fasta, output_file = output_vcf)
	
	subprocess.run(concat_cmd, shell=True, check=True)

	index_cmd = "tabix {output_file}".format(output_file = output_vcf)
	
	subprocess.run(index_cmd, shell=True, check=True)

	print("Merged VCF written to: {}".format(output_vcf))
	print("\n\n")

# Phase genotypes with WhatsHap
def phase_genotypes_whatshap(self):
	print("Phasing Genotypes with WhatsHap!")

	if self.platform == "PACBIO":
		input_bam = os.path.join(self.mapped_bam_dir, self.sample_ID + ".dedup.trimmed.hg38.chr6.bam")
		input_vcf = os.path.join(self.deepvariant_dir, self.sample_ID + ".dedup.trimmed.hg38.chr6.SNV.vcf.gz")
		haplotagged_bam = os.path.join(self.mapped_bam_dir, self.sample_ID + ".dedup.trimmed.hg38.chr6.whatshap.haplotag.bam")
		phased_vcf = os.path.join(self.whatshap_phased_vcf_dir, self.sample_ID + ".dedup.trimmed.hg38.chr6.phased.vcf.gz")
	elif self.platform == "ONT":
		input_bam = os.path.join(self.mapped_bam_dir, self.sample_ID + ".porechop.trimmed.hg38.rmdup.chr6.bam")
		input_vcf = os.path.join(self.clair3_dir, self.sample_ID, "merge_output.vcf.gz")
		haplotagged_bam = os.path.join(self.mapped_bam_dir, self.sample_ID + ".porechop.trimmed.hg38.rmdup.chr6.whatshap.haplotag.bam")
		phased_vcf = os.path.join(self.whatshap_phased_vcf_dir, self.sample_ID + ".porechop.trimmed.hg38.rmdup.chr6.phased.vcf.gz")

	print("Input BAM: {}".format(input_bam))
	print("Input VCF: {}".format(input_vcf))

	output_blocks_file = os.path.join(self.whatshap_phased_vcf_dir, self.sample_ID + ".phased.haploblocks.txt")
	output_gtf_file = os.path.join(self.whatshap_phased_vcf_dir, self.sample_ID + ".phased.haploblocks.gtf")

	whatshap_phase_cmd = "whatshap phase --ignore-read-groups --output {output_file} --reference {reference_genome} {input_vcf} {input_bam}".format(output_file = phased_vcf, reference_genome = reference_fasta, input_vcf = input_vcf, input_bam = input_bam)

	index_cmd = "bcftools index {input_file}".format(input_file = phased_vcf)
	tabix_cmd = "tabix {input_file}".format(input_file = phased_vcf)

	whatshap_haplotag_cmd = "whatshap haplotag --ignore-read-groups --output {output_file} --reference {reference_genome} {input_vcf} {input_bam}".format(output_file = haplotagged_bam, reference_genome = reference_fasta, input_vcf = phased_vcf, input_bam = input_bam)

	whatshap_stats_cmd = "whatshap stats --block-list={block_list_file} --gtf={gtf_file} {input_vcf}".format(block_list_file = output_blocks_file, gtf_file = output_gtf_file, input_vcf = phased_vcf)

	# Log WhatsHap in own output file so it doesn't clog up STDOUT
	whatshap_log = os.path.join(self.whatshap_phased_vcf_dir, self.sample_ID + ".whatshap.log")

	with open(whatshap_log, "w") as log_file:
		log_file.write("\n==== Running WhatsHap Phase ====\n")
		subprocess.run(whatshap_phase_cmd, shell=True, check=True, stdout=log_file, stderr=log_file)
		subprocess.run(index_cmd, shell=True, check=True, stdout=log_file, stderr=log_file)
		subprocess.run(tabix_cmd, shell=True, check=True, stdout=log_file, stderr=log_file)
		log_file.write("\n==== Running WhatsHap Haplotag ====\n")
		subprocess.run(whatshap_haplotag_cmd, shell=True, check=True, stdout=log_file, stderr=log_file)
		log_file.write("\n==== Running WhatsHap Stats ====\n")
		subprocess.run(whatshap_stats_cmd, shell=True, check=True, stdout=log_file, stderr=log_file)

	print("WhatsHap phased VCF written to: {}".format(phased_vcf))
	print("WhatsHap haplotagged BAM written to: {}".format(haplotagged_bam))
	print("WhatsHap phase block gtf written to: {}".format(output_gtf_file))
	print("WhatsHap phase blocks written to: {}".format(output_blocks_file))
	print("\n\n")

def phase_genotypes_longphase(self):
	print("Phasing Genotypes with LongPhase!")

	input_bam = os.path.join(self.mapped_bam_dir, self.sample_ID + ".porechop.trimmed.hg38.rmdup.chr6.bam")
	input_SNV_vcf = os.path.join(self.clair3_dir, self.sample_ID, "merge_output.vcf.gz")
	input_SV_vcf = os.path.join(self.sniffles_dir, self.sample_ID + ".porechop.trimmed.hg38.rmdup.chr6.sniffles.vcf.gz")

	print("Input BAM: {}".format(input_bam))
	print("Input SNV VCF: {}".format(input_SNV_vcf))
	print("Input SV VCF: {}".format(input_SV_vcf))

	haplotagged_bam = os.path.join(self.mapped_bam_dir, self.sample_ID + ".porechop.trimmed.hg38.rmdup.chr6.longphase.haplotag.bam")
	phased_vcf = os.path.join(self.longphase_phased_vcf_dir, self.sample_ID + ".porechop.trimmed.hg38.rmdup.chr6.longphase.vcf.gz")
	output_blocks_file = os.path.join(self.longphase_phased_vcf_dir, self.sample_ID + ".phased.haploblocks.txt")
	output_gtf_file = os.path.join(self.longphase_phased_vcf_dir, self.sample_ID + ".phased.haploblocks.gtf")

	phased_vcf_prefix = phased_vcf.split(".vcf.gz")[0]
	longphase_phase_cmd = "{longphase} phase -s {input_snv_vcf} --sv {input_sv_vcf} -b {input_bam} -r {reference_genome} -t {threads} -o {phased_prefix} --ont".format(longphase = longphase, input_snv_vcf = input_SNV_vcf, input_sv_vcf = input_SV_vcf, input_bam = input_bam, reference_genome = reference_fasta, threads = self.threads, phased_prefix = phased_vcf_prefix)

	# Compress and index SNV VCF
	compress_cmd = "bgzip -f {prefix}.vcf".format(prefix=phased_vcf_prefix)
	index_cmd = "bcftools index {input_file}".format(input_file = phased_vcf)
	tabix_cmd = "tabix {input_file}".format(input_file = phased_vcf)

	# Compress and index SV VCF
	SV_prefix = phased_vcf_prefix + "_SV"
	phased_SV_vcf = os.path.join(self.longphase_phased_vcf_dir, self.sample_ID + ".porechop.trimmed.hg38.rmdup.chr6.longphase_SV.vcf.gz")
	compress_SV_cmd = "bgzip -f {prefix}.vcf".format(prefix=SV_prefix)
	index_SV_cmd = "bcftools index {input_file}".format(input_file = phased_SV_vcf)
	tabix_SV_cmd = "tabix {input_file}".format(input_file = phased_SV_vcf)

	haplotagged_bam_prefix = haplotagged_bam.split(".bam")[0]
	longphase_haplotag_cmd = "{longphase} haplotag -r {reference_genome} -s {input_snv_vcf} --sv-file {input_sv_vcf} -b {input_bam} -t {threads} -o {prefix}".format(longphase = longphase, reference_genome = reference_fasta, input_snv_vcf = phased_vcf, input_sv_vcf = input_SV_vcf, input_bam = input_bam, threads = self.threads, prefix = haplotagged_bam_prefix)

	whatshap_stats_cmd = "whatshap stats --block-list={block_list_file} --gtf={gtf_file} {input_vcf}".format(block_list_file = output_blocks_file, gtf_file = output_gtf_file, input_vcf = phased_vcf)

	# Log WhatsHap in own output file so it doesn't clog up STDOUT
	longphase_log = os.path.join(self.longphase_phased_vcf_dir, self.sample_ID + ".longphase.log")

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

	print("LongPhase phased VCF written to: {}".format(phased_vcf))
	print("LongPhase haplotagged BAM written to: {}".format(haplotagged_bam))
	print("LongPhase phase block gtf written to: {}".format(output_gtf_file))
	print("LongPhase phase blocks written to: {}".format(output_blocks_file))
	print("\n\n")

def merge_longphase_vcfs(self):
	print("Merging longphase SNV and SV VCFs using bcftools...")

	phased_vcf = os.path.join(self.longphase_phased_vcf_dir, self.sample_ID + ".porechop.trimmed.hg38.rmdup.chr6.longphase.vcf.gz")
	phased_SV_vcf = os.path.join(self.longphase_phased_vcf_dir, self.sample_ID + ".porechop.trimmed.hg38.rmdup.chr6.longphase_SV.vcf.gz")
	merged_vcf = os.path.join(self.longphase_phased_vcf_dir, self.sample_ID + ".porechop.trimmed.hg38.rmdup.chr6.longphase.merged.vcf.gz")
	reheadered_SV_vcf = phased_SV_vcf.replace(".vcf.gz", ".reheader.vcf.gz")

	header_file = os.path.join(self.longphase_phased_vcf_dir, self.sample_ID + ".header.txt")
	merge_log = os.path.join(self.longphase_phased_vcf_dir, self.sample_ID + ".merge.log")

	with open(header_file, "w") as hf:
		hf.write(f"SAMPLE\t{self.sample_ID}\n")

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
