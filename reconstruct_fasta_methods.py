import os
import subprocess
import pysam
from Bio.Seq import Seq
from Bio.SeqRecord import SeqRecord
from Bio import SeqIO


vcf2fasta_script = "/hb/scratch/mglasena/vcf2fasta/vcf2fasta.py"
reference_genome = "/hb/groups/cornejo_lab/matt/hla_capture/input_data/reference/GCA_000001405.15_GRCh38_no_alt_analysis_set.fa"
gff_dir = "/hb/groups/cornejo_lab/matt/hla_capture/input_data/hla_gff/"

DNA_bases = {"A", "T", "G", "C"}
stop_codons = ["TAA", "TAG", "TGA"]

def filter_vcf(self):
	input_vcf = os.path.join(self.hiphase_phased_vcf_dir, self.sample_ID + ".dedup.trimmed.hg38.chr6.phased.joint.vcf.gz")
	pass_vcf = os.path.join(self.filtered_vcf_dir, f"{self.sample_ID}_PASS.vcf.gz")
	filtered_vcf = os.path.join(self.filtered_vcf_dir, f"{self.sample_ID}_filtered.vcf.gz")

	# Step 1: Filter for PASS variants and remove unphased hets and unsupported variant types
	#filter_expr = '(GT="hom" || GT~"\\|") && (TYPE="snp" || TYPE="indel" || SVTYPE="INS" || SVTYPE="DEL") && ALT!~"^<"'
	# New filter expression, deal with ./. genotypes
	filter_expr = '(GT="hom" || GT~"\\|") && GT!="./." && (TYPE="snp" || TYPE="indel" || SVTYPE="INS" || SVTYPE="DEL") && ALT!~"^<"'
	
	# Use for DeepVariant
	subprocess.run(f'bcftools view -f PASS -i \'{filter_expr}\' {input_vcf} -Oz -o {pass_vcf}', shell=True, check=True)
	
	# Use for bcftools 
	# subprocess.run(f'bcftools view -i \'{filter_expr}\' {input_vcf} -Oz -o {pass_vcf}', shell=True, check=True)

	subprocess.run(f"bcftools index {pass_vcf}", shell=True, check=True)

	# Step 2: Find first phased variant
	vcf = pysam.VariantFile(pass_vcf)
	first_phased_pos = None
	chrom = None

	for record in vcf:
		for sample_data in record.samples.values():
			if sample_data.phased:
				first_phased_pos = record.pos
				chrom = record.chrom
				break
		if first_phased_pos:
			break

	if not first_phased_pos:
		raise ValueError(f"No phased variants found in {pass_vcf}")

	region = f"{chrom}:{first_phased_pos}-"

	# Step 3: Filter VCF from the first phased variant onward
	subprocess.run(
		f"bcftools view -r {region} {pass_vcf} -Oz -o {filtered_vcf}",
		shell=True, check=True
	)
	subprocess.run(f"bcftools index {filtered_vcf}", shell=True, check=True)

def run_vcf2fasta(self, gene, feature):
	gene_id = gene.lower().replace("-", "_")
	input_vcf = os.path.join(self.filtered_vcf_dir, f"{self.sample_ID}_filtered.vcf.gz")
	out_dir = os.path.join(self.vcf2fasta_out_dir, gene_id)
	
	if feature == "CDS":
		input_gff = os.path.join(gff_dir, gene_id + "_cds_sorted.gff3")
		vcf2fasta_cmd = "python3 {vcf2fasta} --fasta {reference_fasta} --vcf {input_vcf} --gff {input_gff} -o {output_dir} --feat CDS --blend".format(vcf2fasta = vcf2fasta_script, reference_fasta = reference_genome, input_vcf = input_vcf, input_gff = input_gff, output_dir = out_dir)
	elif feature == "gene":
		input_gff = os.path.join(gff_dir, gene_id + "_gene.gff3")
		vcf2fasta_cmd = "python3 {vcf2fasta} --fasta {reference_fasta} --vcf {input_vcf} --gff {input_gff} -o {output_dir} --feat gene".format(
	vcf2fasta = vcf2fasta_script, reference_fasta = reference_genome, input_vcf = input_vcf, input_gff = input_gff, output_dir = out_dir)
	
	subprocess.run(vcf2fasta_cmd, shell = True, check = True)

def parse_fastas(self):
	find_cmd = "find {input_dir} -type f > fasta_files.txt".format(input_dir = self.vcf2fasta_out_dir)
	subprocess.run(find_cmd, shell = True, check = True)
	fasta_files = open("fasta_files.txt", "r").read().splitlines()
	os.remove("fasta_files.txt")

	fasta_dict = dict()

	for file in fasta_files:
		if "_gene" in file:
			feat = "gene"
		elif "_CDS" in file:
			feat = "CDS"
		gene = file.split("/")[-2].split(f"_{feat}")[0].upper().replace("_", "-")
		with open(file, "r") as f:
			lines = f.read().split(">")
		# Remove deletion characters
		allele_1 = lines[1].split("\n")[1].strip().replace("-","").strip()
		allele_2 = lines[2].split("\n")[1].strip().replace("-","").strip()

		if allele_1[0:3] != "ATG" or allele_2[0:3] != "ATG":
			print("File {} does not begin with start codon!".format(file))
		
		elif not allele_1[-3:] in stop_codons or not allele_2[-3:] in stop_codons:
			print("File {} does not end with stop codon!".format(file))

		if not set(allele_1).issubset(DNA_bases):
			print("{} has invalid characters!".format(file))

		if not set(allele_2).issubset(DNA_bases):
			print("{} has invalid characters!".format(file))

		if feat not in fasta_dict:
			fasta_dict[feat] = {}
		if gene not in fasta_dict[feat]:
			fasta_dict[feat][gene] = []

		fasta_dict[feat][gene].append(allele_1)
		fasta_dict[feat][gene].append(allele_2)

	outfile_gene = os.path.join(self.hla_fasta_dir, self.sample_ID + "_HLA_haplotypes_gene.fasta")
	outfile_CDS = os.path.join(self.hla_fasta_dir, self.sample_ID + "_HLA_haplotypes_CDS.fasta")

	gene_records = []
	cds_records = []

	for feat, genes in fasta_dict.items():
		for gene, haplotypes in genes.items():
			hap1_name = f"{self.sample_ID}_{gene}_1"
			hap1_seq = haplotypes[0]
			hap2_name = f"{self.sample_ID}_{gene}_2"
			hap2_seq = haplotypes[1]

			if feat == "gene":
				gene_records.append(SeqRecord(Seq(hap1_seq), id=hap1_name, description = ""))
				gene_records.append(SeqRecord(Seq(hap2_seq), id=hap2_name, description = ""))
			elif feat == "CDS":
				cds_records.append(SeqRecord(Seq(hap1_seq), id=hap1_name, description = ""))
				cds_records.append(SeqRecord(Seq(hap2_seq), id=hap2_name, description = ""))

	SeqIO.write(gene_records, outfile_gene, "fasta")
	print(f"Wrote {len(gene_records)} records to {outfile_gene}")
	SeqIO.write(cds_records, outfile_CDS, "fasta")
	print(f"Wrote {len(cds_records)} records to {outfile_CDS}")
