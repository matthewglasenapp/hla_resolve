#!/usr/bin/env python3
"""
Extract RefCall indels from .dv.vcf.gz files that overlap HLA gene regions.
Prints sample, gene, position, REF, ALT, DP, AD, and VAF.
"""

import subprocess
import os
import sys

BASE_DIR = "/hb/scratch/mglasena/hla_resolve_results/pacbio"

# HLA full gene coordinates (from hla_genes.bed)
HLA_GENES = {
    "HLA-A":    (29941259, 29949572),
    "HLA-C":    (31268748, 31272130),
    "HLA-B":    (31353871, 31357681),
    "HLA-DRB1": (32577901, 32589848),
    "HLA-DQA1": (32632743, 32643685),
    "HLA-DQB1": (32659466, 32668383),
    "HLA-DPA1": (33064568, 33080775),
    "HLA-DPB1": (33075989, 33089696),
}

# Build a BED-style regions string for bcftools
REGIONS = ",".join(f"chr6:{s+1}-{e}" for s, e in HLA_GENES.values())

def get_gene(pos):
    """Return the HLA gene name for a given position, or None."""
    for gene, (start, end) in HLA_GENES.items():
        if start < pos <= end:
            return gene
    return None

def get_sample_dirs(base_dir):
    sample_dirs = []
    for entry in sorted(os.listdir(base_dir)):
        full_path = os.path.join(base_dir, entry)
        if os.path.isdir(full_path) and not entry.startswith("."):
            sample_dirs.append(entry)
    return sample_dirs

def extract_refcall_indels(base_dir):
    samples = get_sample_dirs(base_dir)

    print("\t".join(["SAMPLE", "GENE", "CHROM", "POS", "REF", "ALT", "DP", "AD", "VAF"]))

    for sample in samples:
        vcf_path = os.path.join(base_dir, sample, "genotype_calls", f"{sample}.dv.vcf.gz")

        if not os.path.exists(vcf_path):
            print(f"WARNING: {vcf_path} not found, skipping", file=sys.stderr)
            continue

        print(f"Processing {sample}...", file=sys.stderr)

        # Filter for indels with RefCall in FILTER, restricted to HLA gene regions
        cmd = (
            f"bcftools view -v indels -f RefCall -r {REGIONS} {vcf_path} | "
            f"bcftools query -f '%CHROM\t%POS\t%REF\t%ALT[\t%DP\t%AD\t%VAF]\n'"
        )

        result = subprocess.run(cmd, shell=True, capture_output=True, text=True)

        if result.returncode != 0:
            print(f"ERROR processing {sample}: {result.stderr.strip()}", file=sys.stderr)
            continue

        for line in result.stdout.strip().split("\n"):
            if not line:
                continue
            fields = line.split("\t")
            pos = int(fields[1])
            gene = get_gene(pos)
            if gene:
                print(f"{sample}\t{gene}\t{line}")

if __name__ == "__main__":
    extract_refcall_indels(BASE_DIR)
