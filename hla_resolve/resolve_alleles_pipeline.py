# This software is Copyright ©2026. The Regents of the University of California
# ("Regents"). All Rights Reserved.
#
# See LICENSE.txt for license details.

import csv
import os
import shutil
from .preprocess_methods import (
	run_mosdepth,
	parse_mosdepth
)
from .investigate_haploblocks_methods import (
	parse_haploblocks,
	evaluate_gene_haploblocks
)
from .reconstruct_fasta_methods import (
	filter_vcf_gene,
	run_vcf2fasta,
	parse_fastas
)
from .hla_typer import main as classify_hla_alleles
from .utils import announce, detail, finish_stage, stage
from . import config as hla_config

def convert_gene_name_for_gff(gene_name):
	"""
	Convert gene name from config format (HLA-A) to GFF filename format (hla_a)
	"""
	return gene_name.lower().replace('-', '_')

NOT_TYPED = "not_typed"

def split_result_column(column):
	# HLA-A_1, or HLA-A_1_incomplete where the reconstruction fell short of the
	# whole gene. Returns (None, None) for the sample column.
	name = column[:-len("_incomplete")] if column.endswith("_incomplete") else column
	gene, sep, index = name.rpartition("_")
	if not sep or index not in ("1", "2"):
		return None, None
	return gene, index

def strip_hla_prefix(call):
	# The gene is already the row label, so HLA-A*01:01:01:01 prints as A*01:01:01:01.
	# G group names arrive from the IPD XML already in that short form.
	return call[len("HLA-"):] if call.startswith("HLA-") else call

# Every classification pass writes the same CSV shape, so one reader serves all
# three resolutions.
# (table heading, short label, best-guess file, genotype list string file)
RESULT_TABLES = (
	("G-group resolution", "G-group", "g_group_output.csv", "g_group_output_full.csv"),
	("Three-field resolution", "Three-field", "3_field_allele_output.csv", "3_field_allele_output_full.csv"),
	("Four-field resolution", "Four-field", "allele_output.csv", "allele_output_full.csv"),
)

def read_calls(results_file):
	"""Read one result CSV into (gene order, {gene: {"1": call, "2": call}}).
	Returns (None, None) when the file is missing, empty, or has no gene columns."""
	if not os.path.exists(results_file):
		return None, None

	with open(results_file, newline="") as f:
		rows = list(csv.DictReader(f))

	if not rows:
		return None, None

	row = rows[0]

	# Genes keep the column order of the CSV, which lines up with the deposited
	# haplotype FASTA records.
	calls = {}
	genes = []
	for column, value in row.items():
		gene, index = split_result_column(column)
		if gene is None:
			continue
		if gene not in calls:
			calls[gene] = {}
			genes.append(gene)
		call = (value or "").strip()
		calls[gene][index] = strip_hla_prefix(call) if call else NOT_TYPED

	if not genes:
		return None, None

	return genes, calls

def print_results(config):
	tables = []
	genes = []
	for label, _, filename, _full in RESULT_TABLES:
		table_genes, calls = read_calls(os.path.join(config['hla_typing_dir'], filename))
		if calls is None:
			announce(f"WARNING: no {label} calls found in {filename}")
			continue
		tables.append((label, calls))
		# A gene can reach one resolution and not another, so take the union.
		for gene in table_genes:
			if gene not in genes:
				genes.append(gene)

	if not tables:
		return

	# Pass 1 assigns a G group only on a perfect ARS match, so a blank means the
	# G group was never directly observed, either because the allele belongs to
	# none or because the reconstruction differs by a base or two. Reporting the
	# nearest group would assert a match that was not seen, so the three-field
	# call stands in. A G group name ends in G, so the two stay distinguishable.
	substituted = False
	by_label = dict(tables)
	g_calls = by_label.get("G-group resolution")
	three_calls = by_label.get("Three-field resolution")
	if g_calls is not None and three_calls is not None:
		for gene, haplotypes in three_calls.items():
			for index, call in haplotypes.items():
				if call == NOT_TYPED:
					continue
				if g_calls.setdefault(gene, {}).get(index, NOT_TYPED) == NOT_TYPED:
					g_calls[gene][index] = call
					substituted = True

	# One set of widths across all three blocks, so they line up under each other.
	gene_width = max([len("gene")] + [len(gene) for gene in genes])
	first_width = max([len("_1")] + [len(calls.get(gene, {}).get("1", NOT_TYPED))
	                                 for _, calls in tables for gene in genes])

	announce(f"Sample: {config['sample_ID']}")
	for label, calls in tables:
		announce()
		announce(label)
		announce(f"{'gene':<{gene_width}}   {'_1':<{first_width}}   _2")
		for gene in genes:
			first = calls.get(gene, {}).get("1", NOT_TYPED)
			second = calls.get(gene, {}).get("2", NOT_TYPED)
			announce(f"{gene:<{gene_width}}   {first:<{first_width}}   {second}")

	announce()
	announce("Note: Allele order within each gene is arbitrary and is not consistent between genes.")
	announce("Note: Ambiguous calls are written as genotype list strings to the ambiguities files below.")
	if substituted:
		announce("Note: Where no perfect G group match was found, the three-field call is shown.")

def print_output_files(config, phased_vcf):
	entries = [
		("Haplotagged BAM", config['hg38_rmdup_chr6_haplotag_bam']),
		("Phased VCF", phased_vcf),
		("Gene FASTA", config['hla_gene_fasta']),
		("CDS FASTA", config['hla_cds_fasta']),
	]
	for _, label, filename, ambiguities in RESULT_TABLES:
		entries.append((f"{label} calls", os.path.join(config['hla_typing_dir'], filename)))
		entries.append((f"{label} ambiguities", os.path.join(config['hla_typing_dir'], ambiguities)))
	label_width = max(len(label) for label, _ in entries) + 1

	announce("Output files:")
	for label, path in entries:
		announce(f"  {label + ':':<{label_width}}  {path}")

def resolve_alleles(config):
	"""
	Main orchestration function that runs the complete allele resolution workflow:
	1. Coverage analysis
	2. FASTA sequence reconstruction 
	3. HLA typing
	4. Print results
	"""	
	stage("Coverage assessment")
	run_mosdepth(
		input_file=config['hg38_rmdup_chr6_bam'],
		output_dir=config['mosdepth_dir'],
		sample_ID=config['sample_ID'],
		regions_file=config['hla_genes_regions_file'],
		threads=config['threads']
	)
	
	sufficient_coverage_genes, coverage_stats = parse_mosdepth(
		regions_file=config['mosdepth_regions'],
		thresholds_file=config['mosdepth_thresholds'],
		cds_depth_thresh=config['cds_depth_thresh'],
		cds_prop_20x_thresh=config['cds_prop_20x_thresh'],
		cds_prop_30x_thresh=config['cds_prop_30x_thresh'],
		ars_depth_thresh=config['ars_depth_thresh'],
		ars_prop_20x_thresh=config['ars_prop_20x_thresh'],
		ars_prop_30x_thresh=config['ars_prop_30x_thresh']
	)
	
	if config['platform'] == "PACBIO":
		phased_vcf = config['hiphase_joint_vcf']
		haploblock_file = config['phased_blocks']
	elif config['platform'] == "ONT":
		phased_vcf = config['longphase_vcf']
		haploblock_file = config['phased_haploblocks']
	
	heterozygous_sites, haploblock_list = parse_haploblocks(
		input_vcf=phased_vcf,
		input_haploblock_file=haploblock_file,
		platform=config['platform'],
		sample_ID=config['sample_ID'],
		mhc_start=config['mhc_start'],
		mhc_stop=config['mhc_stop']
	)

	stage("Haploblock evaluation")
	phased_genes, unphased_genes, do_not_type_genes, cds_rescued_genes = evaluate_gene_haploblocks(
		output_file=config['phased_genes_tsv'],
		incomplete_file=config['incomplete_genes_csv'],
		sample_ID=config['sample_ID'],
		genes_bed=config['genes_bed'],
		genes_of_interest=config['genes_of_interest'],
		het_sites=heterozygous_sites,
		haploblocks=haploblock_list,
		ARS_dict=config.get('ARS_dict', None),
		CDS_dict=config.get('CDS_dict', None))
	
	# Print which genes were successfully phased (skip empty sections)
	fully_phased = [g for g in config['genes_of_interest'] if g in phased_genes]
	if fully_phased:
		print(f"Fully phased genes: {', '.join(fully_phased)}")

	partially_phased = [g for g in config['genes_of_interest'] if g in unphased_genes]
	if partially_phased:
		print(f"Partially phased genes: {', '.join(partially_phased)}")

	rescued = [g for g in config['genes_of_interest'] if g in cds_rescued_genes]
	if rescued:
		tiers = [f"{gene} (tier: {cds_rescued_genes[gene]['tier']})" for gene in rescued]
		print(f"CDS-rescued genes: {', '.join(tiers)}")
	print()
	
	if config['platform'] == "PACBIO":
		input_vcf = config['hiphase_joint_vcf']
	elif config['platform'] == "ONT":
		input_vcf = config['longphase_merged_vcf']

	stage("Variant filtering and redundancy removal")
	# Filter phased VCF by gene region
	gene_filtered_vcfs = {}
	unphased_het_counts = {}
	for gene in config['genes_of_interest']:
		gff_gene_name = convert_gene_name_for_gff(gene)
		gff_file = os.path.join(config['gff_dir'], gff_gene_name + "_gene.gff3")
		gff_lines = [line.split("\t") for line in open(gff_file, "r").read().splitlines() if not line.startswith("#")]
		gff_record = gff_lines[0]
		chromosome, start, stop = gff_record[0], gff_record[3], gff_record[4]
		filter_region = f"{chromosome}:{start}-{stop}"

		# Create gene-specific output file paths in filtered_vcf directory
		gene_symbolic_vcf = os.path.join(config['filtered_vcf_dir'], f"{config['sample_ID']}_{gene}.symbolic.vcf.gz")
		gene_pass_vcf = os.path.join(config['filtered_vcf_dir'], f"{config['sample_ID']}_{gene}_PASS.vcf.gz")
		gene_fail_vcf = os.path.join(config['filtered_vcf_dir'], f"{config['sample_ID']}_{gene}_FAIL.vcf.gz")
		gene_sv_overlap_vcf = os.path.join(config['filtered_vcf_dir'], f"{config['sample_ID']}_{gene}_SV_OVERLAP.vcf.gz")
		gene_pass_unphased_vcf = os.path.join(config['filtered_vcf_dir'], f"{config['sample_ID']}_{gene}_PASS_UNPHASED.vcf.gz")
		gene_filtered_vcf = os.path.join(config['filtered_vcf_dir'], f"{config['sample_ID']}_{gene}_PASS_phased.vcf.gz")

		# Derive genotyper label from split caller config
		if config['snp_caller'] == config['indel_caller']:
			genotyper = config['snp_caller']
		else:
			genotyper = "hybrid"

		unphased_hets = filter_vcf_gene(
			input_vcf=input_vcf,
			gene=gene,
			filter_region=filter_region,
			symbolic_vcf=gene_symbolic_vcf,
			pass_vcf=gene_pass_vcf,
			fail_vcf=gene_fail_vcf,
			sv_overlap_vcf=gene_sv_overlap_vcf,
			pass_unphased=gene_pass_unphased_vcf,
			filtered_vcf=gene_filtered_vcf,
			platform=config['platform'],
			genotyper=genotyper,
			force_include_unphased=(gene in cds_rescued_genes)
		)
		
		gene_filtered_vcfs[gene] = gene_filtered_vcf
		if unphased_hets:
			unphased_het_counts[gene] = unphased_hets

	if unphased_het_counts:
		counts = ", ".join(f"{gene} ({count})" for gene, count in unphased_het_counts.items())
		print(f"Unphased PASS heterozygous variants: {counts}")
		print("Unphased variants written to: " + os.path.join(config['filtered_vcf_dir'],
												 f"{config['sample_ID']}_<gene>_PASS_UNPHASED.vcf.gz"))
	else:
		print("Unphased PASS heterozygous variants: none")
	print()

	stage("Gene sequence reconstruction")
	# Reset vcf2fasta_out_dir for sequential runs 
	if any(os.scandir(config['vcf2fasta_out_dir'])):
		shutil.rmtree(config['vcf2fasta_out_dir'])
		os.makedirs(config['vcf2fasta_out_dir'], exist_ok=True)

	for gene in config['genes_of_interest']:
		if gene in sufficient_coverage_genes and gene not in do_not_type_genes:
			gene_filtered_vcf = gene_filtered_vcfs.get(gene)
				
			gff_gene_name = convert_gene_name_for_gff(gene)
			run_vcf2fasta(
				input_vcf=gene_filtered_vcf,
				output_dir=os.path.join(config['vcf2fasta_out_dir'], gene),
				input_gff=os.path.join(config['gff_dir'], gff_gene_name + "_gene.gff3"),
				reference_genome=config['reference_genome'],
				gene=gene,
				feature="gene")

			run_vcf2fasta(
				input_vcf=gene_filtered_vcf,
				output_dir=os.path.join(config['vcf2fasta_out_dir'], gene),
				input_gff=os.path.join(config['gff_dir'], gff_gene_name + "_cds_sorted.gff3"),
				reference_genome=config['reference_genome'],
				gene=gene,
				feature="CDS")
	
	parse_fastas(
		sample_ID=config['sample_ID'],
		vcf2fasta_output_dir=config['vcf2fasta_out_dir'],
		outfile_gene=config['hla_gene_fasta'],
		outfile_CDS=config['hla_cds_fasta'],
		DNA_bases=config['DNA_bases'],
		stop_codons=config['stop_codons'],
		unphased_genes=unphased_genes,
		gene_dict=config['gene_dict'],
		CDS_dict=config['CDS_dict'],
		gff_dir=config['gff_dir'],
		cds_rescued_genes=cds_rescued_genes,
		ARS_dict=config.get('ARS_dict', None),
		CLASS_I_GENES=config.get('CLASS_I_GENES', None),
		gene_filtered_vcfs=gene_filtered_vcfs,
		reference_genome=config['reference_genome']
	)
	
	# DR/DQ read re-consensus context (PacBio only; gated by config flag).
	from .reconsensus_drdq import DRDQ_GENES
	reconsensus_ctx = None
	if getattr(hla_config, "reconsensus_drdq", False) and config['platform'] == "PACBIO":
		reconsensus_ctx = {
			"enabled": True,
			"mode": getattr(hla_config, "reconsensus_read_assignment", "hp_tag"),
			"bam": config['hg38_rmdup_chr6_haplotag_bam'],
			"gene_vcfs": {gene: gene_filtered_vcfs.get(gene) for gene in DRDQ_GENES},
			"gene_dict": config['gene_dict'],
			"gene_fasta": config['hla_gene_fasta'],
			"cds_fasta": config['hla_cds_fasta'],
		}

	# Step 3: HLA typing
	stage("IPD-IMGT/HLA database matching")
	print("Matching reconstructed sequences against the IPD-IMGT/HLA database...")
	classifications = classify_hla_alleles(
		reference_xml_file=config['IMGT_XML'],
		hla_fasta_dir=config['hla_fasta_dir'],
		sample_ID=config['sample_ID'],
		generate_query_ref_comp=True,
		reconsensus_ctx=reconsensus_ctx,
		output_dir=config['hla_typing_dir']
	)

	print()
	finish_stage()
	print()

	print_results(config)
	print()

	detail(f"HLA typing result files located in dir: {config['hla_typing_dir']}/")
	detail("")
	detail("g group results written to: g_group_output.csv")
	detail("three field results written to: 3_field_allele_output.csv")
	detail("four field results written to: allele_output.csv")
	detail("")
	detail("Results with ambiguities in the format of genotype list strings written to:")
	detail("g_group_output_full.csv")
	detail("3_field_allele_output_full.csv")
	detail("allele_output_full.csv")
	detail("")
	detail("Debugging files written to:")
	detail("g_group_assignment.log")
	detail("3_field_allele_assignment.log")
	detail("allele_assignment.log")
	detail("sample_ref_comp.csv")
	detail("\n")

	print_output_files(config, input_vcf)
	print()

	print("HLA allele resolution workflow completed")

	return {
		"classifications": classifications,
		"coverage_stats": coverage_stats,
		"phased_genes": fully_phased,
		"cds_rescued_genes": cds_rescued_genes,
	}
