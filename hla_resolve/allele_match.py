#!/usr/bin/env python3

import edlib
import parasail
import csv
import os
import subprocess
import re
from collections import defaultdict
from Bio import SeqIO

# ---------------------------
# Configuration
# ---------------------------

CAPTURE_FASTA = "/hb/scratch/mglasena/hla_resolve_results/pacbio/HLA_haplotypes_full.fa"
HPRC_FASTA_DIR = "/hb/groups/cornejo_lab/matt/hla_capture/input_data/hprc/"
PLATFORM = "pacbio"

GENES = {
    "HLA-A", "HLA-B", "HLA-C",
    "HLA-DPA1", "HLA-DPB1",
    "HLA-DQA1", "HLA-DQB1",
    "HLA-DRB1"
}

DISTANCE_ALERT_THRESHOLD = 2000
PRETTY_ALIGNMENT_WIDTH = 100

# Parasail scoring parameters
PARASAIL_MATRIX = parasail.matrix_create("ACGT", 1, -1)
PARASAIL_GAP_OPEN = 5
PARASAIL_GAP_EXTEND = 1

# ---------------------------
# Load HPRC reference haplotypes
# ---------------------------

HPRC_dict = defaultdict(lambda: defaultdict(tuple))

def create_HPRC_dict():
    find_cmd = f"find {HPRC_FASTA_DIR} -name '*.fa' -type f"
    fasta_files = subprocess.run(
        find_cmd, capture_output=True, text=True, shell=True
    ).stdout.strip().split("\n")

    for file in fasta_files:
        if not file:
            continue
        sample = os.path.basename(file).split(".")[0]

        for record in SeqIO.parse(file, "fasta"):
            gene = record.id.split("::")[0]
            seq = str(record.seq)
            HPRC_dict[gene][sample] += (seq,)

# ---------------------------
# Load capture haplotypes
# ---------------------------

def load_capture_haplotypes(fasta_path):
    capture = defaultdict(lambda: defaultdict(list))

    for record in SeqIO.parse(fasta_path, "fasta"):
        try:
            sample, gene, hap = record.id.split("_", 2)
        except ValueError:
            raise ValueError(f"Unexpected FASTA header: {record.id}")

        if gene not in GENES:
            continue

        capture[gene][sample].append(str(record.seq))

    for gene in capture:
        for sample, haps in capture[gene].items():
            if len(haps) != 2:
                raise ValueError(
                    f"{sample} {gene} has {len(haps)} capture haplotypes (expected 2)"
                )

    return capture

# ---------------------------
# Pretty alignment helper
# ---------------------------

def pretty_alignment(seq1, seq2, cigar, width=100):
    i = j = 0
    a1, a2, mid = [], [], []

    for length, op in re.findall(r"(\d+)([=XID])", cigar):
        length = int(length)

        if op == "=":
            for _ in range(length):
                a1.append(seq1[i])
                a2.append(seq2[j])
                mid.append("|")
                i += 1
                j += 1
        elif op == "X":
            for _ in range(length):
                a1.append(seq1[i])
                a2.append(seq2[j])
                mid.append("*")
                i += 1
                j += 1
        elif op == "I":
            for _ in range(length):
                a1.append(seq1[i])
                a2.append("-")
                mid.append(" ")
                i += 1
        elif op == "D":
            for _ in range(length):
                a1.append("-")
                a2.append(seq2[j])
                mid.append(" ")
                j += 1

    lines = []
    for k in range(0, len(a1), width):
        lines.append("".join(a1[k:k+width]))
        lines.append("".join(mid[k:k+width]))
        lines.append("".join(a2[k:k+width]))
        lines.append("")

    return "\n".join(lines)

# ---------------------------
# Alignment + distances
# ---------------------------

def compute_alignment(seq1, seq2, context=None):
    # Shorter sequence = query (fully aligned),
    # longer sequence = database (free begin/end gaps)
    if len(seq1) <= len(seq2):
        shorter, longer = seq1, seq2
    else:
        shorter, longer = seq2, seq1

    # --- edlib: true minimum edit distance + gap-compressed ---
    # HW mode: shorter is fully aligned within longer (database ends free)
    edlib_result = edlib.align(shorter, longer, mode="HW", task="path")
    edlib_raw = edlib_result["editDistance"]
    edlib_gc = edlib_raw
    for length, op in re.findall(r"(\d+)([ID])", edlib_result["cigar"]):
        length = int(length)
        if length > 1:
            edlib_gc -= (length - 1)

    # --- parasail: affine-gap NW on edlib-extracted region ---
    # edlib's HW mode reliably locates where the shorter sequence maps
    # within the longer. We extract that region and run parasail NW
    # (global alignment) on the extracted region for biologically
    # realistic affine-gap distances. This avoids NW force-aligning
    # flanking sequence when there are large length differences.
    start, end = edlib_result["locations"][0]
    db_region = longer[start:end + 1]

    parasail_result = parasail.nw_trace_striped_sat(
        shorter, db_region, PARASAIL_GAP_OPEN, PARASAIL_GAP_EXTEND, PARASAIL_MATRIX
    )
    cigar = parasail_result.cigar.decode.decode("utf-8")

    # Parasail raw = total mismatches + total indel bases
    # Parasail gap-compressed = mismatches + indel runs (each run counts as 1)
    parasail_raw = 0
    parasail_gc = 0
    for length, op in re.findall(r"(\d+)([=XID])", cigar):
        length = int(length)
        if op == "X":
            parasail_raw += length
            parasail_gc += length
        elif op in ("I", "D"):
            parasail_raw += length
            parasail_gc += 1

    aln_len = sum(int(l) for l, op in re.findall(r"(\d+)([=XID])", cigar))
    identity = 1 - (parasail_gc / aln_len) if aln_len > 0 else 0

    if parasail_gc > DISTANCE_ALERT_THRESHOLD:
        print("\n" + "=" * 80)
        print("LARGE GAP-COMPRESSED DISTANCE")
        if context:
            print("Context:", context)
        print("Edlib raw / GC:", edlib_raw, "/", edlib_gc)
        print("Parasail raw / GC:", parasail_raw, "/", parasail_gc)
        print("Alignment length:", aln_len)
        print("CIGAR:", cigar)
        print(pretty_alignment(shorter, db_region, cigar, PRETTY_ALIGNMENT_WIDTH))
        print("=" * 80 + "\n")

    return edlib_raw, parasail_raw, edlib_gc, parasail_gc, aln_len, identity

# ---------------------------
# Core comparison logic
# ---------------------------

def compare_haplotypes(sample, gene, capture_alleles):
    if gene not in HPRC_dict or sample not in HPRC_dict[gene]:
        return None

    ref = HPRC_dict[gene][sample]

    # Haploid reference
    if len(ref) == 1:
        a1 = compute_alignment(
            capture_alleles[0], ref[0],
            context=f"{sample} {gene} hap1 vs ref"
        )
        a2 = compute_alignment(
            capture_alleles[1], ref[0],
            context=f"{sample} {gene} hap2 vs ref"
        )
        return a1 if a1[3] <= a2[3] else a2

    # Diploid reference — pair by parasail gap-compressed distance (index 3)
    one_one = compute_alignment(
        capture_alleles[0], ref[0],
        context=f"{sample} {gene} hap1 vs ref1"
    )
    two_two = compute_alignment(
        capture_alleles[1], ref[1],
        context=f"{sample} {gene} hap2 vs ref2"
    )
    sum1 = one_one[3] + two_two[3]

    one_two = compute_alignment(
        capture_alleles[0], ref[1],
        context=f"{sample} {gene} hap1 vs ref2"
    )
    two_one = compute_alignment(
        capture_alleles[1], ref[0],
        context=f"{sample} {gene} hap2 vs ref1"
    )
    sum2 = one_two[3] + two_one[3]

    return (
        (*one_one, *two_two)
        if sum1 <= sum2
        else (*one_two, *two_one)
    )

# ---------------------------
# Main
# ---------------------------

def main():
    create_HPRC_dict()
    capture = load_capture_haplotypes(CAPTURE_FASTA)

    with open("allele_match.csv", "w", newline="") as f:
        writer = csv.writer(f)
        writer.writerow([
            "platform", "gene", "sample",
            "hap_1_edlib_raw", "hap_2_edlib_raw",
            "hap_1_parasail_raw", "hap_2_parasail_raw",
            "hap_1_edlib_gc", "hap_2_edlib_gc",
            "hap_1_parasail_gc", "hap_2_parasail_gc",
            "alignment_length_1", "alignment_length_2",
            "hap_1_identity", "hap_2_identity"
        ])

        for gene in capture:
            for sample, haps in capture[gene].items():
                result = compare_haplotypes(sample, gene, haps)
                if result is None:
                    continue

                if len(HPRC_dict[gene][sample]) == 1:
                    er1, pr1, egc1, pgc1, l1, i1 = result
                    er2 = pr2 = egc2 = pgc2 = l2 = i2 = "NA"
                else:
                    er1, pr1, egc1, pgc1, l1, i1, er2, pr2, egc2, pgc2, l2, i2 = result

                writer.writerow([
                    PLATFORM, gene, sample,
                    er1, er2,
                    pr1, pr2,
                    egc1, egc2,
                    pgc1, pgc2,
                    l1, l2,
                    i1, i2
                ])

if __name__ == "__main__":
    main()
