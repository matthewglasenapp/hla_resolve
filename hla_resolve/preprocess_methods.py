# This software is Copyright ©2026. The Regents of the University of California
# ("Regents"). All Rights Reserved.
#
# See LICENSE.txt for license details.

import os
import sys
import subprocess
import pysam
import gzip
import shutil
import tempfile
from . import config
from .utils import run_quiet

# Convert BAM file of unmapped HiFi (ccs) reads to FASTQ format for marking duplicates and trimming adapters
def convert_bam_to_fastq(input_file, output_file, platform, threads):
	if platform == "PACBIO":
		print("Converting HiFi ccs reads to fastq format using pbtk bam2fastq!")
		print(f"bam2fastq input file: {input_file}")

		pbi_file = input_file + ".pbi"
		if not os.path.exists(pbi_file):
			print(f"PBI index not found for {input_file}. Indexing with pbindex...")
			run_quiet(f"pbindex {input_file}")

		output_prefix = output_file.split(".fastq.gz")[0]

		bam2fastq_cmd = f"bam2fastq -j {threads} {input_file} -o {output_prefix}"

		run_quiet(bam2fastq_cmd)

	elif platform == "ONT":
		print("Converting ONT raw reads to fastq format using samtools fastq!")
		print(f"Samtools fastq input file: {input_file}")
		samtools_threads = int(threads * 2 / 3)
		pigz_threads = threads - samtools_threads
		samtools_fq_cmd = f"samtools fastq -@ {samtools_threads} {input_file} | pigz -p {pigz_threads} > {output_file}"
		
		run_quiet(samtools_fq_cmd)
	
	print(f"Raw fastq reads written to: {output_file}")
	print("\n")

def trim_adapters(adapters, input_file, output_file, sample_ID, threads, adapter_file = None, five_prime_adapter = None, three_prime_adapter = None):
	if adapters:

		if adapter_file and five_prime_adapter and three_prime_adapter:
			print("Trimming adapter sequences with cutadapt!")
			print(f"cutadapt input file: {input_file}")
			print(f"5' adapter: {five_prime_adapter}")
			print(f"3' adapter: {three_prime_adapter}")

			# Use cutadapt with specific adapter sequences from the adapter file
			cutadapt_cmd = f"cutadapt -j {threads} --quiet -n 2 --minimum-length 100 -g {five_prime_adapter} -a {three_prime_adapter} -o {output_file} {input_file}"

			run_quiet(cutadapt_cmd)

			print(f"Trimmed reads written to: {output_file}")
			print("\n")

		else:
			print("Trimming adapter sequences with fastplong in AUTO mode!")
			print(f"fastplong input file: {input_file}")

			output_dir = os.path.dirname(output_file)
			html_path = os.path.join(output_dir, sample_ID + ".fastplong.html")
			json_path = os.path.join(output_dir, sample_ID + ".fastplong.json")

			fastplong_cmd = f"fastplong -i {input_file} -h {html_path} -j {json_path} -w {threads} -n 100000 -o {output_file}"

			subprocess.run(fastplong_cmd, shell=True, check=False)

			print(f"Trimmed reads written to: {output_file}")
			print("\n")

	else:
		print("Users specified no adapters present")
		print(f"Skipping adapter trimming and transfering raw fastq file {input_file} to {output_file}")
		shutil.copy(input_file, output_file)

# Mark PCR duplicates with pbmarkdup
def mark_duplicates_pbmarkdup(input_file, output_file, threads):
	print("Removing PCR duplicates using pbmarkdup!")

	print(f"pbmarkdup input file: {input_file}")

	pbmarkdup_cmd = f"pbmarkdup -j {threads} --rmdup --clobber {input_file} {output_file}"

	result = subprocess.run(pbmarkdup_cmd, shell=True, check=True, capture_output=True, text=True)
	if result.stdout:
		print(result.stdout, end="")
	if result.stderr:
		print(result.stderr, end="", file=sys.stderr)

	gzip_cmd = f"pigz -f -p {threads} {output_file}"
	run_quiet(gzip_cmd)
	
	print(f"De-duplicated reads written to: {output_file}.gz")
	print("\n")

# Align to GRCh38 reference genome with rammap
def align_to_reference_rammap(input_file, output_file, read_group_string, reference_fasta, platform, threads):
	print("Aligning reads to GRCh38 reference genome with rammap!")

	if platform == "PACBIO":
		platform_string = "map-hifi"
	elif platform == "ONT":
		platform_string = "map-ont"

	print(f"rammap input file: {input_file}")

	rammap_threads = int(threads * 2 / 3)
	samtools_threads = threads - rammap_threads
	rammap_rg_string = "'{}'".format(read_group_string.replace("\t", "\\t"))

	# rammap (minimap2) only accepts FASTA/FASTQ. PacBio WGS/WES input is an
	# unmapped BAM, so stream it through samtools fastq and read from stdin ("-")
	# rather than materializing a whole-genome FASTQ on disk.
	if input_file.endswith(".bam"):
		rammap_cmd = (
			f"samtools fastq -@ {samtools_threads} {input_file} | "
			f"{config.rammap} -Y -t {rammap_threads} -ax {platform_string} -R {rammap_rg_string} {reference_fasta} - | "
			f"samtools sort -@ {samtools_threads} -o {output_file}"
		)
	else:
		rammap_cmd = f"{config.rammap} -Y -t {rammap_threads} -ax {platform_string} -R {rammap_rg_string} {reference_fasta} {input_file} | samtools sort -@ {samtools_threads} -o {output_file}"
	index_bam = f"samtools index {output_file}"

	run_quiet(rammap_cmd)
	run_quiet(index_bam)

	print(f"Mapped bam written to: {output_file}")
	print("\n")

def _parse_drb_paralog_reads(output_file, drb_paralog_reads_file):
	"""
	Parse aligned BAM to identify DRB paralog reads (DRB3, DRB4, DRB5, DRB6, DRB7, DRB8, DRB9)
	based on primary alignment. Any read whose best-matching allele is NOT a DRB1
	allele is flagged for removal. Shared helper for classify_DRB_reads.
	"""
	drb_paralog_read_ids = set()

	with pysam.AlignmentFile(output_file, "rb") as bam:
		for read in bam:
			# Skip secondary, supplementary, and unmapped reads
			if read.is_secondary or read.is_supplementary or read.is_unmapped:
				continue

			# Flag any read whose best-matching allele is not DRB1
			# (IPD-style headers start with "DRB1*"; pseudogene contigs use names
			# like "DRB5_GRCh38" / "DRB6_GRCh38" / "DRB9_GRCh38".)
			reference_name = read.reference_name
			if not reference_name.startswith("DRB1*"):
				drb_paralog_read_ids.add(read.query_name)

	with open(drb_paralog_reads_file, "w") as f:
		for read_id in sorted(drb_paralog_read_ids):
			f.write(read_id + "\n")

	print(f"Classified {len(drb_paralog_read_ids)} reads as DRB paralog (DRB3/4/5/6/7/8/9)")
	print(f"DRB paralog read IDs written to: {drb_paralog_reads_file}")
	print("\n")

def classify_DRB_reads(input_file, output_file, drb_paralog_reads_file, read_group_string, reference_fasta, platform, threads, region=None):
	"""
	Identify reads originating from DRB paralogs (DRB3/4/5/6/7/8/9) using competitive
	mapping against a multi-allele reference containing representative genomic
	sequences from DRB1, DRB3, DRB4, DRB5, DRB6, DRB7, DRB8, and DRB9.

	For each read, rammap picks the best-matching allele as the primary
	alignment. If the best match is anything other than a DRB1 allele, that
	read ID is written to drb_paralog_reads_file for downstream removal by
	filter_reads().

	When `region` is given, `input_file` is treated as a coordinate-sorted,
	indexed GRCh38 BAM and only the primary reads overlapping that region (the
	DR cluster) are competitively mapped. This restricts paralog classification
	to reads already placed in the DR region instead of re-mapping the full read
	set. When `region` is None the whole `input_file` is mapped.
	"""
	print("Classifying DRB reads using multi-allele competitive mapping (rammap)!")

	if platform == "PACBIO":
		platform_string = "map-hifi"
	elif platform == "ONT":
		platform_string = "map-ont"

	rammap_threads = int(threads * 2 / 3)
	samtools_threads = threads - rammap_threads
	rammap_rg_string = "'{}'".format(read_group_string.replace("\t", "\\t"))

	# Restrict competitive mapping to primary reads already placed in the DR region.
	map_input = input_file
	if region is not None:
		map_input = output_file.replace(".bam", ".dr_region.fastq")
		run_quiet(f"samtools view -b -F 0x900 -@ {samtools_threads} {input_file} {region} | samtools fastq -@ {samtools_threads} - > {map_input}")

	# Map reads against the multi-allele DRB reference (DRB1, DRB3, DRB4 alleles)
	rammap_cmd = f"{config.rammap} -Y -t {rammap_threads} -ax {platform_string} -R {rammap_rg_string} {reference_fasta} {map_input} | samtools sort -@ {samtools_threads} -o {output_file}"
	index_bam = f"samtools index {output_file}"

	run_quiet(rammap_cmd)
	run_quiet(index_bam)

	_parse_drb_paralog_reads(output_file, drb_paralog_reads_file)

# Mark duplicates for ONT data or WGS PacBio data
# pbmarkdup used for hybrid-capture PacBio data but does not scale well for WGS data
def mark_duplicates_picard(input_file, output_file, metrics_file, temp_dir, picard):
	os.makedirs(temp_dir, exist_ok=True)
	
	mark_duplicates_cmd = f"java -jar {picard} MarkDuplicates -I {input_file} -O {output_file} --TMP_DIR {temp_dir} -M {metrics_file} --CREATE_INDEX true --REMOVE_DUPLICATES true --VALIDATION_STRINGENCY LENIENT"
	
	run_quiet(mark_duplicates_cmd)

	index_bam = f"samtools index {output_file}"
	run_quiet(index_bam)

# Filter reads that did not map to chromosome 6
def filter_reads(input_file, output_file, drb_paralog_reads_file, threads):
	print("Excluding BAM records that don't map to chromosome 6!")

	print(f"Samtools input file: {input_file}")

	# Exclude DRB paralog reads (DRB3/4/5/6/7/8/9, identified by competitive mapping),
	# restrict to the chr6 MHC window, and drop secondary/supplementary records.
	samtools_cmd = f"samtools view -h -F 2304 -@ {threads} {input_file} chr6:28000000-34000000 | grep -v -F -f {drb_paralog_reads_file} -- | samtools view -b -o {output_file}"

	index_cmd = f"samtools index {output_file}"

	run_quiet(samtools_cmd)
	run_quiet(index_cmd)

	count_reads_cmd = f"samtools view -c {output_file}"

	read_count = int(subprocess.check_output(count_reads_cmd, shell=True).strip())

	print(f"Filtered BAM records written to: {output_file}")
	print("\n")

	return read_count
	
# Call SNV with DeepVariant
def call_variants_deepvariant(input_bam, output_vcf, output_gvcf, platform, deepvariant_sif, reference_fasta, genotypes_dir, mapped_bam_dir, sample_ID, threads):
	if platform == "PACBIO":
		model_type = "PACBIO"
	elif platform == "ONT":
		model_type = "ONT_R104"
	
	print("Calling SNVs and small indels with DeepVariant!")
	print(f"DeepVariant input file: {input_bam}")
	
	bind_paths = [
		f"{genotypes_dir}:/data",
		f"{mapped_bam_dir}:/input",
		f"{os.path.dirname(reference_fasta)}:/reference"
	]

	bind_flags = " ".join(f"--bind {path}" for path in bind_paths)

	deepvariant_shards = min(threads, 8)

	deepvariant_cmd = f"""
		singularity exec {bind_flags} {deepvariant_sif} /opt/deepvariant/bin/run_deepvariant \
			--model_type={model_type} \
			--ref=/reference/{os.path.basename(reference_fasta)} \
			--reads=/input/{os.path.basename(input_bam)} \
			--output_vcf=/data/{os.path.basename(output_vcf)} \
			--output_gvcf=/data/{os.path.basename(output_gvcf)} \
			--regions chr6 \
			--num_shards={deepvariant_shards}
		"""

	# Log DeepVariant in own output file so it doesn't clog up STDOUT
	deepvariant_log = os.path.join(genotypes_dir, sample_ID + ".deepvariant.log")
	print(f"Writing stdout to {deepvariant_log}")

	with open(deepvariant_log, "w") as log_file:
		subprocess.run(deepvariant_cmd, shell=True, check=True, stdout=log_file, stderr=log_file)

	print(f"VCF written to: {output_vcf}")
	print(f"GVCF written to: {output_gvcf}")
	print("\n")

def call_variants_clair3(input_bam, output_vcf, platform, clair3_sif, reference_fasta, threads, genotypes_dir, mapped_bam_dir, sample_ID, clair3_model):
	if platform == "ONT":
		platform_type = "ont"
	elif platform == "PACBIO":
		platform_type = "hifi"

	print("Calling SNVs and small indels with Clair3!")
	print(f"Clair3 input file: {input_bam}")
	print(f"Clair3 model: {clair3_model}")

	output_dir = os.path.join(genotypes_dir, sample_ID)
	os.makedirs(output_dir, exist_ok=True)

	bind_paths = [
		f"{output_dir}:/output",
		f"{mapped_bam_dir}:/input",
		f"{os.path.dirname(reference_fasta)}:/reference"
	]

	bind_flags = " ".join(f"--bind {path}" for path in bind_paths)

	clair3_cmd = f"""
		singularity exec {bind_flags} {clair3_sif} /opt/bin/run_clair3.sh \
			--bam_fn=/input/{os.path.basename(input_bam)} \
			--ref_fn=/reference/{os.path.basename(reference_fasta)} \
			--threads={threads} \
			--platform={platform_type} \
			--model_path=/opt/models/{clair3_model} \
			--output=/output \
			--sample_name={sample_ID} \
			--ctg_name=chr6
		"""

	clair3_log = os.path.join(genotypes_dir, sample_ID + ".clair3.log")
	print(f"Writing stdout to {clair3_log}")

	with open(clair3_log, "w") as log_file:
		subprocess.run(clair3_cmd, shell=True, check=True, stdout=log_file, stderr=log_file)

	raw_genotypes_file = os.path.join(output_dir, "merge_output.vcf.gz")
	shutil.copy(raw_genotypes_file, output_vcf)
	run_quiet(f"tabix -p vcf {output_vcf}")

	print(f"VCF written to: {output_vcf}")
	print("\n")

# Call SNV with bcftools
def call_variants_bcftools(input_file, output_file, reference_fasta, platform, threads):
	if platform == "PACBIO":
		config = "pacbio-ccs-1.20"
	elif platform == "ONT":
		config = "ont-sup-1.20"

	print("Calling SNVs and small indels with bcftools!")

	print(f"Bcftools input file: {input_file}")
	
	pileup_threads = str(threads // 2)
	call_threads = str(threads // 2)

	bcftools_command = (
		f"bcftools mpileup --config {config} --threads {pileup_threads} "
		f"-f {reference_fasta} -d 1000000 -r chr6:28000000-34000000 "
		f"-a FORMAT/DP,AD,ADF,ADR,SP {input_file} | "
		f"bcftools call -mv -f GQ --threads {call_threads} -Ou | "
		f"bcftools view -i 'FORMAT/DP>=2 && ((TYPE=\"snp\" && GQ>=20 && QUAL>=10) || (TYPE=\"indel\" && GQ>=20 && QUAL>=10))' "
		f"-Oz -o {output_file}")

	run_quiet(bcftools_command)
	run_quiet(f"tabix -p vcf {output_file}")

	print(f"VCF written to: {output_file}")
	print("\n")

# Call SNVs and small indels with FreeBayes
def call_variants_freebayes(input_bam, output_vcf, reference_fasta):
	print("Calling SNVs and small indels with FreeBayes!")
	print(f"FreeBayes input file: {input_bam}")

	freebayes_cmd = (
		f"freebayes "
		f"-f {reference_fasta} "
		f"-r chr6:28000000-34000000 "
		f"--genotype-qualities "
		f"{input_bam} | "
		f"bcftools view -i 'FORMAT/DP>=2 && ((INFO/TYPE=\"snp\" && GQ>=20 && QUAL>=10) || (INFO/TYPE!=\"snp\" && GQ>=20 && QUAL>=10))' | "
		f"bgzip > {output_vcf}"
	)

	run_quiet(freebayes_cmd)
	run_quiet(f"tabix -p vcf {output_vcf}")

	print(f"VCF written to: {output_vcf}")
	print("\n")

# Rescue high-confidence variants from DeepVariant RefCall filter
# DeepVariant can conservatively filter real variants as RefCall (especially in HLA)
# This function overturns RefCall -> PASS for variants meeting quality thresholds
# When indels_only=True (hybrid mode), only indels are rescued; SNP alts are skipped
# When indels_only=False (pure deepvariant mode), both SNPs and indels are rescued
def rescue_refcalls(input_vcf, output_vcf, indels_only=False):
	mode_label = "indels" if indels_only else "SNPs and indels"
	print(f"Rescuing high-confidence RefCall {mode_label} from DeepVariant output!")
	print(f"Input VCF: {input_vcf}")

	vcf_in = pysam.VariantFile(input_vcf)

	tmp_fd, tmp_path = tempfile.mkstemp(suffix=".vcf.gz", dir=os.path.dirname(output_vcf))
	os.close(tmp_fd)
	vcf_out = pysam.VariantFile(tmp_path, 'wz', header=vcf_in.header)

	rescued = 0
	total_refcalls = 0

	for record in vcf_in:
		if record.alts and "RefCall" in record.filter:
			# In indels_only mode, skip records that have no indel alt alleles
			if indels_only:
				has_indel_alt = any(len(alt) != len(record.ref) for alt in record.alts)
				if not has_indel_alt:
					vcf_out.write(record)
					continue

			total_refcalls += 1
			sample = record.samples[0]

			gq = sample.get('GQ')
			dp = sample.get('DP')
			ad = sample.get('AD')
			vaf = sample.get('VAF')

			# To re-enable GQ filtering, add: gq >= 20 and
			if all(v is not None for v in (gq, dp, ad, vaf)):
				if dp >= 30:
					# Normalize VAF to tuple (pysam may return scalar for biallelic)
					if isinstance(vaf, (int, float)):
						vaf = (vaf,)

					# Evaluate each alt allele independently
					# Classify as het (0.3-0.7) or hom-alt (>0.7)
					passing_het = []
					passing_hom = []
					for i, alt in enumerate(record.alts):
						# In indels_only mode, skip SNP alts in multiallelic records
						if indels_only and len(alt) == len(record.ref):
							continue
						alt_ad = ad[i + 1]  # AD[0] is ref depth
						alt_vaf = vaf[i]
						if alt_ad is not None and alt_vaf is not None and alt_ad >= 10:
							if 0.3 <= alt_vaf <= 0.7:
								passing_het.append(i + 1)  # 1-indexed for VCF GT
							elif alt_vaf > 0.7:
								passing_hom.append(i + 1)

					if passing_hom or passing_het:
						record.filter.clear()
						record.filter.add('PASS')
						# Set genotype: hom-alt takes precedence over het
						if passing_hom:
							new_gt = (passing_hom[0], passing_hom[0])
							zygosity = "hom-alt"
						elif len(passing_het) >= 2:
							new_gt = (passing_het[0], passing_het[1])
							zygosity = "compound-het"
						else:
							new_gt = (0, passing_het[0])
							zygosity = "het"
						sample['GT'] = new_gt
						rescued += 1

						if config.VERBOSE:
							gt_str = "/".join(str(a) for a in new_gt)
							alts_str = ",".join(record.alts)
							ad_str = ",".join(str(a) for a in ad)
							vaf_str = ",".join(f"{v:.3f}" for v in vaf)
							is_snp = len(record.ref) == 1 and all(len(a) == 1 for a in record.alts)
							var_type = "SNP" if is_snp else "INDEL"
							print(f"  RESCUED {var_type}: {record.chrom}:{record.pos} {record.ref}>{alts_str} GT={gt_str} ({zygosity}) GQ={gq} DP={dp} AD={ad_str} VAF={vaf_str}")

		vcf_out.write(record)

	vcf_in.close()
	vcf_out.close()

	os.replace(tmp_path, output_vcf)
	run_quiet(f"tabix -p vcf {output_vcf}")

	print(f"Rescued {rescued}/{total_refcalls} RefCall {mode_label}")
	print(f"Output VCF: {output_vcf}")
	print("\n")

# Extract SNPs from snp_caller output and indels from indel_caller output, merge into a single VCF
_CALLER_DISPLAY = {
	"bcftools": "bcftools",
	"deepvariant": "DeepVariant",
	"clair3": "Clair3",
	"freebayes": "FreeBayes",
}

def merge_hybrid_vcfs(snp_vcf, indel_vcf, indel_only_vcf, merged_vcf, snp_caller, indel_caller, filter_indel_pass=True):
	snp_name = _CALLER_DISPLAY.get(snp_caller, snp_caller)
	indel_name = _CALLER_DISPLAY.get(indel_caller, indel_caller)
	print(f"Merging hybrid SNP ({snp_name}) and indel ({indel_name}) VCFs!")

	snp_only_vcf = snp_vcf.replace('.vcf.gz', '.snps_only.vcf.gz')
	snp_cmd = f"bcftools view -v snps {snp_vcf} -Oz -o {snp_only_vcf}"
	run_quiet(snp_cmd)
	run_quiet(f"tabix -p vcf {snp_only_vcf}")

	if filter_indel_pass:
		indel_cmd = f"bcftools view -v indels -f PASS {indel_vcf} -Oz -o {indel_only_vcf}"
	else:
		indel_cmd = f"bcftools view -v indels {indel_vcf} -Oz -o {indel_only_vcf}"
	run_quiet(indel_cmd)
	run_quiet(f"tabix -p vcf {indel_only_vcf}")

	merge_cmd = f"bcftools concat -a {snp_only_vcf} {indel_only_vcf} | bcftools sort | bgzip > {merged_vcf}"
	run_quiet(merge_cmd)
	run_quiet(f"tabix -p vcf {merged_vcf}")

	print(f"Merged VCF written to: {merged_vcf}")
	print("\n")

# Run pbsv to call structural variants (SV)
def call_structural_variants_pbsv(input_bam, output_svsig, output_vcf, threads, tandem_repeat_bed, reference_fasta):
	print("Calling structural variants with pbsv!")

	# Strip SA tags so pbsv cannot reconstruct split reads into large cross-gene
	# SVs. Supplementary records are already removed upstream (-F 2304 in
	# filter_reads), but the SA:Z: tags on retained primaries still let pbsv call
	# implausible multi-kb cross-gene events from 2-8 kb reads. SVs are taken from
	# CIGAR signatures of primary alignments only. The variant-calling/phasing BAM
	# is left untouched; only this pbsv input copy is stripped.
	sa_stripped_bam = input_bam.replace(".bam", ".pbsv_nosa.bam")
	run_quiet(f"samtools view -h -x SA -b -@ {threads} -o {sa_stripped_bam} {input_bam}")
	run_quiet(f"samtools index {sa_stripped_bam}")
	input_bam = sa_stripped_bam

	print(f"pbsv input file: {input_bam}")
	
	# -a Don't downsample
	# -q Don't filter by MapQ
	pbsv_discover_cmd = f"pbsv discover -a 0 -q 0 --region chr6 --tandem-repeats {tandem_repeat_bed} {input_bam} {output_svsig}"
	
	run_quiet(pbsv_discover_cmd)

	index_svsig_cmd = f"tabix -c '#' -s 3 -b 4 -e 4 {output_svsig}"
	
	run_quiet(index_svsig_cmd)

	pbsv_call_cmd = f"pbsv call -j {threads} --min-sv-length 20 --region chr6 --hifi {reference_fasta} {output_svsig} {output_vcf}"
	
	run_quiet(pbsv_call_cmd)

	compress_cmd = f"bgzip -c {output_vcf} > {output_vcf}.gz"
	index_vcf_cmd = f"tabix -p vcf {output_vcf}.gz"
	
	run_quiet(compress_cmd)
	run_quiet(index_vcf_cmd)

	print(f"pbsv SV VCF written to: {output_vcf}")
	print("\n")

def call_structural_variants_sniffles(input_bam, output_vcf, threads, reference_fasta, chr6_bed, tandem_repeat_bed):
	print("Calling structural variants with Sniffles2!")

	sniffles_cmd = f"sniffles --output-rnames --allow-overwrite -t 1 --reference {reference_fasta} --regions {chr6_bed} -i {input_bam} -v {output_vcf} --tandem-repeats {tandem_repeat_bed}"

	run_quiet(sniffles_cmd)

	print(f"Sniffles2 SV VCF written to: {output_vcf}")
	print("\n")

# Genotype tandem repeats with pbtrgt
def genotype_tandem_repeats(input_bam, output_vcf, pbtrgt_dir, threads, reference_fasta, pbtrgt_repeat_file, original_cwd, scheme):
	print("Genotyping tandem repeats with TRGT!")

	print(f"TRGT input file: {input_bam}")
	
	# Extract the base name without any extensions
	output_prefix = os.path.basename(output_vcf)
	if output_prefix.endswith('.vcf.gz'):
		output_prefix = output_prefix[:-7]  # Remove .vcf.gz
	elif output_prefix.endswith('.vcf'):
		output_prefix = output_prefix[:-4]  # Remove .vcf
	
	os.chdir(pbtrgt_dir)

	try:
		# TRGT's 'targeted' preset (disables rq>=0.98 read filtering, uses the
		# cluster genotyper, and targeted flank scoring) suits enriched capture/
		# amplicon data; WGS/WES use TRGT's default 'wgs' preset.
		trgt_preset = "targeted" if scheme in ("hybrid_capture", "amplicon") else "wgs"
		trgt_cmd = f"trgt genotype --threads {threads} --genome {reference_fasta} --reads {input_bam} --repeats {pbtrgt_repeat_file} --output-prefix {output_prefix} --preset {trgt_preset}"

		run_quiet(trgt_cmd)

		sort_cmd = f"bcftools sort -O z -o {output_prefix + '.sorted.vcf.gz'} {output_prefix + '.vcf.gz'}"

		run_quiet(sort_cmd)

		os.rename(output_prefix + ".sorted.vcf.gz", output_prefix + ".vcf.gz")

		index_cmd = f"tabix -p vcf {output_prefix + '.vcf.gz'}"

		run_quiet(index_cmd)

		print(f"TR VCF written to: {output_vcf}")
		print("\n")
	finally:
		os.chdir(original_cwd)

# Phase genotypes with HiPhase
def phase_genotypes_hiphase(input_bam, input_snv, input_SV, input_TR, output_bam, output_snv, output_SV, output_TR, output_summary_file, output_blocks_file, output_stats_file, threads, reference_fasta, phased_vcf_dir, sample_ID):
	print("Phasing genotypes with HiPhase!")

	print(f"Input BAM: {input_bam}")
	print(f"Input SNV: {input_snv}")
	print(f"Input SV: {input_SV}")
	print(f"Input TR: {input_TR}")
	
	# Check if input files exist
	missing_files = []
	if not os.path.exists(input_bam):
		missing_files.append(f"BAM: {input_bam}")
	if not os.path.exists(input_snv):
		missing_files.append(f"SNV VCF: {input_snv}")
	if not os.path.exists(input_SV):
		missing_files.append(f"SV VCF: {input_SV}")
	if not os.path.exists(input_TR):
		missing_files.append(f"TR VCF: {input_TR}")
	
	if missing_files:
		print(f"ERROR: Missing input files for HiPhase:")
		for file in missing_files:
			print(f"  - {file}")
		print("Skipping HiPhase phasing step.")
		return
	
	hiphase_cmd = f"hiphase --threads {threads} --ignore-read-groups --reference {reference_fasta} --bam {input_bam} --output-bam {output_bam} --vcf {input_snv} --output-vcf {output_snv} --vcf {input_SV} --output-vcf {output_SV} --vcf {input_TR} --output-vcf {output_TR} --stats-file {output_stats_file} --blocks-file {output_blocks_file} --summary-file {output_summary_file}"
	
	# Log HiPhase in own output file so it doesn't clog up STDOUT
	hiphase_log = os.path.join(phased_vcf_dir, sample_ID + ".hiphase.log")
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
	print("\n")

# Merge phased SNV (DeepVariant), tandem repeat (TRGT), and structural variant (pbsv) VCFs with bcftools concat
def merge_hiphase_vcfs(input_snv, input_SV, input_TR, output_vcf, reference_fasta):
	print("Merging phased small variant, pbsv, and TRGT VCF files!")

	print(f"Small variant input file: {input_snv}")
	print(f"pbsv input file: {input_SV}")
	print(f"TRGT input file: {input_TR}")

	# Merge SNV + SV with bcftools norm (normalizes indels, removes dups).
	# TR records are concatenated WITHOUT norm to preserve explicit allele sequences
	# that TRGT encodes in REF/ALT. bcftools norm destroys these by left-aligning
	# and trimming, making them unusable by vcf2fasta.
	output_dir = os.path.dirname(output_vcf)
	snv_sv_normed = os.path.join(output_dir, "snv_sv_normed.tmp.vcf.gz")

	# Step 1: Merge and normalize SNV + SV only
	norm_cmd = f"bcftools concat --allow-overlaps {input_snv} {input_SV} | grep -vE 'chrX|chrY' | grep -vE 'SVTYPE=BND|SVTYPE=INV|SVTYPE=DUP' | bcftools norm -d none --fasta-ref {reference_fasta} | bcftools sort | bgzip > {snv_sv_normed}"
	run_quiet(norm_cmd)
	run_quiet(f"tabix {snv_sv_normed}")

	# Step 2: Filter TR VCF (remove chrX/chrY) without norm
	tr_filtered = os.path.join(output_dir, "tr_filtered.tmp.vcf.gz")
	tr_cmd = f"bcftools view {input_TR} | grep -vE 'chrX|chrY' | bcftools sort | bgzip > {tr_filtered}"
	run_quiet(tr_cmd)
	run_quiet(f"tabix {tr_filtered}")

	# Step 3: Concat normed SNV/SV with un-normed TR, sort
	concat_cmd = f"bcftools concat --allow-overlaps {snv_sv_normed} {tr_filtered} | bcftools sort | bgzip > {output_vcf}"
	run_quiet(concat_cmd)
	run_quiet(f"tabix {output_vcf}")

	# Clean up temp files
	for tmp in [snv_sv_normed, snv_sv_normed + ".tbi", tr_filtered, tr_filtered + ".tbi"]:
		if os.path.exists(tmp):
			os.remove(tmp)

	print(f"Merged VCF written to: {output_vcf}")
	print("\n")

def phase_genotypes_longphase(input_bam, input_SNV_vcf, input_SV_vcf, output_blocks_file, output_gtf_file, phased_vcf, phased_SV_vcf, haplotagged_bam, longphase, reference_fasta, threads, phased_vcf_dir, sample_ID):
	print("Phasing genotypes with LongPhase!")

	print(f"Input BAM: {input_bam}")
	print(f"Input SNV VCF: {input_SNV_vcf}")
	print(f"Input SV VCF: {input_SV_vcf}")

	phased_vcf_prefix = phased_vcf.split(".vcf.gz")[0]
	longphase_phase_cmd = f"{longphase} phase -s {input_SNV_vcf} --sv-file {input_SV_vcf} -b {input_bam} -r {reference_fasta} -t {threads} -o {phased_vcf_prefix} --ont --indels"

	# Compress and index SNV VCF
	compress_cmd = f"bgzip -f {phased_vcf_prefix}.vcf"
	index_cmd = f"bcftools index {phased_vcf_prefix + '.vcf.gz'}"
	tabix_cmd = f"tabix {phased_vcf_prefix + '.vcf.gz'}"

	# Compress and index SV VCF
	SV_prefix = phased_vcf_prefix + "_SV"
	compress_SV_cmd = f"bgzip -f {SV_prefix}.vcf"
	index_SV_cmd = f"bcftools index {phased_SV_vcf}"
	tabix_SV_cmd = f"tabix {phased_SV_vcf}"

	haplotagged_bam_prefix = haplotagged_bam.rsplit(".bam", 1)[0]
	longphase_haplotag_cmd = f"{longphase} haplotag -r {reference_fasta} -s {phased_vcf} --sv-file {phased_SV_vcf} -b {input_bam} -t {threads} -o {haplotagged_bam_prefix}"

	whatshap_stats_cmd = f"whatshap stats --block-list={output_blocks_file} --gtf={output_gtf_file} {phased_vcf}"

	# Log WhatsHap in own output file so it doesn't clog up STDOUT
	longphase_log = os.path.join(phased_vcf_dir, sample_ID + ".longphase.log")

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
	print("\n")

def merge_longphase_vcfs(phased_vcf, phased_SV_vcf, merged_vcf, reference_fasta, phased_vcf_dir, sample_ID):
	print("Merging LongPhase SNV and SV VCFs using bcftools...")

	reheadered_SV_vcf = phased_SV_vcf.replace(".vcf.gz", ".reheader.vcf.gz")

	header_file = os.path.join(phased_vcf_dir, sample_ID + ".header.txt")
	merge_log = os.path.join(phased_vcf_dir, sample_ID + ".merge.log")

	with open(header_file, "w") as hf:
		hf.write(f"SAMPLE\t{sample_ID}\n")

	reheader_cmd = f"bcftools reheader -s {header_file} {phased_SV_vcf} -o {reheadered_SV_vcf}"
	index_sv_cmd = f"bcftools index {reheadered_SV_vcf}"

	merge_cmd = (
		f"bcftools concat --allow-overlaps -a {phased_vcf} {reheadered_SV_vcf} | "
		f"bcftools norm -d none -f {reference_fasta} | "
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
	print("\n")

def run_mosdepth(input_file, output_dir, sample_ID, regions_file, threads):
	print(f"Running mosdepth on {input_file}")
	
	prefix = os.path.join(output_dir, sample_ID)
	
	# --flag 3328 excludes duplicates and secondary/supplementary alignments
	mosdepth = f"mosdepth --flag 3328 --by {regions_file} --thresholds 10,20,30 -t {threads} {prefix} {input_file}"
	
	run_quiet(mosdepth)
	
	print("\n")

def parse_mosdepth(regions_file, thresholds_file, cds_depth_thresh, cds_prop_20x_thresh, cds_prop_30x_thresh,
					ars_depth_thresh, ars_prop_20x_thresh, ars_prop_30x_thresh):
	"""Aggregate mosdepth per-interval stats into per-gene full-gene + CDS + ARS gates.

	BED rows are expected to be named `<gene>_CDS_<n>` (one row per coding
	exon), `<gene>_ARS` (one row per gene), or `<gene>_<other>` (one row per
	gene for the full-gene interval — e.g. `HLA-A_ENSG00000206503`). CDS rows
	are pooled by gene using length-weighted sums, so `prop_30x` is the
	proportion of bases meeting the threshold across concatenated exons — not
	an average of per-exon proportions. Full-gene rows are reported but not
	used as a gate.
	"""
	# Per-gene accumulators for CDS rows
	cds_totals = {}       # gene -> {total_len, depth_sum, n_10x, n_20x, n_30x}
	ars_stats = {}        # gene -> {depth, prop_10x, prop_20x, prop_30x}
	gene_stats = {}       # gene -> {depth, prop_10x, prop_20x, prop_30x}

	with gzip.open(regions_file, "rt") as f1, gzip.open(thresholds_file, "rt") as f2:
		regions = f1.read().splitlines()
		thresholds = f2.read().splitlines()[1:]

		# regions.bed.gz cols: chrom, start, end, name, mean_depth
		# thresholds.bed.gz cols (with --thresholds 10,20,30): chrom, start, end, name, 10X, 20X, 30X
		for regions_line, thresholds_line in zip(regions, thresholds):
			rf = regions_line.split("\t")
			tf = thresholds_line.split("\t")
			name = rf[3]
			length = int(rf[2]) - int(rf[1])
			mean_depth = float(rf[4])
			n_10x, n_20x, n_30x = int(tf[4]), int(tf[5]), int(tf[6])

			if name.endswith("_ARS"):
				gene = name[:-len("_ARS")]
				ars_stats[gene] = {
					"depth": mean_depth,
					"prop_10x": n_10x / length,
					"prop_20x": n_20x / length,
					"prop_30x": n_30x / length,
				}
			elif "_CDS_" in name:
				gene = name.split("_CDS_")[0]
				b = cds_totals.setdefault(gene, {
					"total_len": 0, "depth_sum": 0.0,
					"n_10x": 0, "n_20x": 0, "n_30x": 0,
				})
				b["total_len"] += length
				b["depth_sum"] += mean_depth * length
				b["n_10x"] += n_10x
				b["n_20x"] += n_20x
				b["n_30x"] += n_30x
			else:
				# Full-gene row: `<gene>_<id>` (e.g. HLA-A_ENSG00000206503)
				gene = name.split("_")[0]
				gene_stats[gene] = {
					"depth": mean_depth,
					"prop_10x": n_10x / length,
					"prop_20x": n_20x / length,
					"prop_30x": n_30x / length,
				}

	# Per-gene aggregated CDS stats
	cds_stats = {}
	for gene, b in cds_totals.items():
		L = b["total_len"]
		cds_stats[gene] = {
			"depth": b["depth_sum"] / L if L else 0.0,
			"prop_10x": b["n_10x"] / L if L else 0.0,
			"prop_20x": b["n_20x"] / L if L else 0.0,
			"prop_30x": b["n_30x"] / L if L else 0.0,
		}

	# Human-readable per-gene report (gene, CDS aggregated, ARS as-is)
	print("Region, Mean Depth, % Bases 10X, % Bases 20X, % Bases 30X")
	for gene in sorted(set(list(cds_stats.keys()) + list(ars_stats.keys()) + list(gene_stats.keys()))):
		g = gene_stats.get(gene)
		if g is not None:
			print(f"{gene} gene", f"{g['depth']:.1f}",
				  f"{g['prop_10x']*100:.1f}%", f"{g['prop_20x']*100:.1f}%", f"{g['prop_30x']*100:.1f}%")
		c = cds_stats.get(gene)
		if c is not None:
			print(f"{gene} CDS", f"{c['depth']:.1f}",
				  f"{c['prop_10x']*100:.1f}%", f"{c['prop_20x']*100:.1f}%", f"{c['prop_30x']*100:.1f}%")
		a = ars_stats.get(gene)
		if a is not None:
			print(f"{gene} ARS", f"{a['depth']:.1f}",
				  f"{a['prop_10x']*100:.1f}%", f"{a['prop_20x']*100:.1f}%", f"{a['prop_30x']*100:.1f}%")

	# Apply gates
	cds_pass = {
		gene: (c["depth"] >= cds_depth_thresh
			   and c["prop_20x"] >= cds_prop_20x_thresh
			   and c["prop_30x"] >= cds_prop_30x_thresh)
		for gene, c in cds_stats.items()
	}
	ars_pass = {
		gene: (a["depth"] >= ars_depth_thresh
			   and a["prop_20x"] >= ars_prop_20x_thresh
			   and a["prop_30x"] >= ars_prop_30x_thresh)
		for gene, a in ars_stats.items()
	}

	sufficient_coverage_genes = []
	for gene in cds_pass:
		cds_ok = cds_pass.get(gene, False)
		ars_ok = ars_pass.get(gene, False)
		if cds_ok and ars_ok:
			sufficient_coverage_genes.append(gene)
		else:
			reasons = []
			if not cds_ok:
				c = cds_stats[gene]
				reasons.append(
					f"CDS [depth={c['depth']:.1f} prop_20x={c['prop_20x']:.2f} prop_30x={c['prop_30x']:.2f}]"
				)
			if not ars_ok:
				a = ars_stats.get(gene, {"depth": 0.0, "prop_20x": 0.0, "prop_30x": 0.0})
				reasons.append(
					f"ARS [depth={a['depth']:.1f} prop_20x={a['prop_20x']:.2f} prop_30x={a['prop_30x']:.2f}]"
				)
			print(f"Gene {gene} has insufficient {' and '.join(reasons)} coverage for haplotyping and star allele calling")

	print("\n")
	return sufficient_coverage_genes
