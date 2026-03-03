#!/usr/bin/env python3
"""
Shift repeat BED file coordinates downstream of the DRB1 insertion point.

Applies +4959 bp shift to chr6 entries with start >= 32,589,742 to account
for the augmented DRB1 region in augmented_hg38_with_long_drb1.fa.

Usage:
    python3 shift_repeat_beds.py <input_bed> <output_bed>

Example:
    python3 shift_repeat_beds.py \
        /path/to/hla_resolve/data/repeats_bed/human_GRCh38_no_alt_analysis_set.trf.bed \
        /hb/scratch/mglasena/drb1_experiment/human_GRCh38_no_alt_analysis_set.trf.shifted.bed

    python3 shift_repeat_beds.py \
        /path/to/hla_resolve/data/repeats_bed/chr6_polymorphic_repeats.hg38.bed \
        /hb/scratch/mglasena/drb1_experiment/chr6_polymorphic_repeats.hg38.shifted.bed
"""
import sys

SHIFT = 4959
DRB1_POS = 32589742  # DRB1 start codon — shift everything at or beyond this position

if len(sys.argv) != 3:
    print(f"Usage: {sys.argv[0]} <input_bed> <output_bed>")
    sys.exit(1)

input_bed = sys.argv[1]
output_bed = sys.argv[2]

shifted_count = 0
total_count = 0

with open(input_bed) as fin, open(output_bed, "w") as fout:
    for line in fin:
        if not line.strip() or line.startswith("#") or line.startswith("track"):
            fout.write(line)
            continue

        fields = line.rstrip("\n").split("\t")
        total_count += 1

        try:
            chrom = fields[0]
            start = int(fields[1])
            end = int(fields[2])

            if chrom == "chr6" and start >= DRB1_POS:
                fields[1] = str(start + SHIFT)
                fields[2] = str(end + SHIFT)
                shifted_count += 1

            fout.write("\t".join(fields) + "\n")
        except (ValueError, IndexError):
            fout.write(line)

print(f"Done: {shifted_count}/{total_count} entries shifted by +{SHIFT} bp")
print(f"Output: {output_bed}")
