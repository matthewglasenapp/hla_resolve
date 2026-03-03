#!/usr/bin/env python3
"""
Shift BED coordinates by +4959 bp for entries after HLA-DRB1 on chr6.

When the GRCh38 reference is augmented with DRB1 ALT scaffolds, coordinates
downstream of DRB1 (>= 32589742) are shifted by 4959 bp. This script applies
that offset to produce a shifted BED file for use with the augmented reference.

Input:  data/bed_files/parse_haploblocks/parse_haploblocks_bed.bed
Output: data/bed_files/parse_haploblocks/parse_haploblocks_bed.shifted.bed
"""
import os

SHIFT = 4959
DRB1_POS = 32589742

script_dir = os.path.dirname(__file__)
bed_dir = os.path.join(script_dir, "..", "..", "glasenapp_2026_manuscript", "data", "bed_files", "parse_haploblocks")

infile = os.path.join(bed_dir, "parse_haploblocks_bed.bed")
outfile = os.path.join(bed_dir, "parse_haploblocks_bed.shifted.bed")

with open(infile) as fin, open(outfile, "w") as fout:
    for line in fin:
        if line.strip() == "" or line.startswith("#"):
            fout.write(line)
            continue
        fields = line.strip().split("\t")
        chrom, start, end = fields[0], int(fields[1]), int(fields[2])
        if chrom == "chr6" and start >= DRB1_POS:
            start += SHIFT
            end += SHIFT
        fields[1] = str(start)
        fields[2] = str(end)
        fout.write("\t".join(fields) + "\n")

print(f"Shifted BED written to {outfile}")
