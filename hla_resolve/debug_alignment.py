#!/usr/bin/env python3
"""
Debug script: compare edlib and parasail CIGARs side-by-side.
Writes detailed alignment info to debug_alignment.log
"""

import edlib
import parasail
import os
import subprocess
import re
import sys
from collections import defaultdict
from Bio import SeqIO

# ---------------------------
# Configuration
# ---------------------------

CAPTURE_FASTA = "/hb/scratch/mglasena/hla_resolve_results/pacbio/HLA_haplotypes_full.fa"
HPRC_FASTA_DIR = "/hb/groups/cornejo_lab/matt/hla_capture/input_data/hprc/"

GENES = {
    "HLA-A", "HLA-B", "HLA-C",
    "HLA-DPA1", "HLA-DPB1",
    "HLA-DQA1", "HLA-DQB1",
    "HLA-DRB1"
}

PARASAIL_MATRIX = parasail.matrix_create("ACGT", 1, -1)
PARASAIL_GAP_OPEN = 5
PARASAIL_GAP_EXTEND = 1

LOG_FILE = "debug_alignment.log"

# ---------------------------
# Data loading (same as allele_match.py)
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

def load_capture_haplotypes(fasta_path):
    capture = defaultdict(lambda: defaultdict(list))
    for record in SeqIO.parse(fasta_path, "fasta"):
        try:
            sample, gene, hap = record.id.split("_", 2)
        except ValueError:
            continue
        if gene not in GENES:
            continue
        capture[gene][sample].append(str(record.seq))
    return capture

# ---------------------------
# Debug alignment
# ---------------------------

def debug_alignment(seq1, seq2, label, log):
    if len(seq1) <= len(seq2):
        shorter, longer = seq1, seq2
        shorter_label, longer_label = "seq1", "seq2"
    else:
        shorter, longer = seq2, seq1
        shorter_label, longer_label = "seq2", "seq1"

    log.write(f"\n{'='*80}\n")
    log.write(f"ALIGNMENT: {label}\n")
    log.write(f"  seq1 length: {len(seq1)}\n")
    log.write(f"  seq2 length: {len(seq2)}\n")
    log.write(f"  shorter ({shorter_label}): {len(shorter)}bp\n")
    log.write(f"  longer  ({longer_label}): {len(longer)}bp\n")
    log.write(f"  length diff: {len(longer) - len(shorter)}bp\n")
    log.write(f"\n")

    # --- edlib ---
    edlib_result = edlib.align(shorter, longer, mode="HW", task="path")
    edlib_raw = edlib_result["editDistance"]
    edlib_cigar = edlib_result["cigar"]
    start, end = edlib_result["locations"][0]

    edlib_gc = edlib_raw
    for length, op in re.findall(r"(\d+)([ID])", edlib_cigar):
        length = int(length)
        if length > 1:
            edlib_gc -= (length - 1)

    edlib_aln_len = sum(int(l) for l, op in re.findall(r"(\d+)([=XID])", edlib_cigar))

    log.write(f"  EDLIB (HW mode):\n")
    log.write(f"    raw edit distance: {edlib_raw}\n")
    log.write(f"    gc edit distance:  {edlib_gc}\n")
    log.write(f"    location in longer: [{start}, {end}] ({end - start + 1}bp)\n")
    log.write(f"    alignment length:  {edlib_aln_len}\n")
    log.write(f"    CIGAR length:      {len(edlib_cigar)} chars\n")
    if len(edlib_cigar) < 500:
        log.write(f"    CIGAR: {edlib_cigar}\n")
    else:
        log.write(f"    CIGAR (first 200): {edlib_cigar[:200]}...\n")
        log.write(f"    CIGAR (last 200):  ...{edlib_cigar[-200:]}\n")
    log.write(f"\n")

    # --- parasail NW (full sequences) ---
    parasail_result = parasail.nw_trace_striped_sat(
        shorter, longer, PARASAIL_GAP_OPEN, PARASAIL_GAP_EXTEND, PARASAIL_MATRIX
    )
    parasail_cigar_full = parasail_result.cigar.decode.decode("utf-8")
    parasail_score = parasail_result.score

    log.write(f"  PARASAIL NW (global, full sequences):\n")
    log.write(f"    score: {parasail_score}\n")
    log.write(f"    CIGAR length: {len(parasail_cigar_full)} chars\n")
    if len(parasail_cigar_full) < 500:
        log.write(f"    CIGAR: {parasail_cigar_full}\n")
    else:
        log.write(f"    CIGAR (first 200): {parasail_cigar_full[:200]}...\n")
        log.write(f"    CIGAR (last 200):  ...{parasail_cigar_full[-200:]}\n")

    # Analyze first/last operations
    ops = re.findall(r"(\d+)([=XID])", parasail_cigar_full)
    if ops:
        log.write(f"    first 5 ops: {ops[:5]}\n")
        log.write(f"    last 5 ops:  {ops[-5:]}\n")
    log.write(f"\n")

    # --- stripping ---
    leading = re.match(r"^(\d+[ID])+", parasail_cigar_full)
    trailing = re.search(r"(\d+[ID])+$", parasail_cigar_full)
    log.write(f"  STRIPPING:\n")
    log.write(f"    leading indels match: {leading.group() if leading else 'NONE'}\n")
    log.write(f"    trailing indels match: {trailing.group() if trailing else 'NONE'}\n")

    cigar_stripped = re.sub(r"^(\d+[ID])+", "", parasail_cigar_full)
    cigar_stripped = re.sub(r"(\d+[ID])+$", "", cigar_stripped)

    ops_stripped = re.findall(r"(\d+)([=XID])", cigar_stripped)
    if ops_stripped:
        log.write(f"    stripped first 5 ops: {ops_stripped[:5]}\n")
        log.write(f"    stripped last 5 ops:  {ops_stripped[-5:]}\n")
    log.write(f"\n")

    # --- parasail metrics after stripping ---
    parasail_raw = 0
    parasail_gc = 0
    for length, op in re.findall(r"(\d+)([=XID])", cigar_stripped):
        length = int(length)
        if op == "X":
            parasail_raw += length
            parasail_gc += length
        elif op in ("I", "D"):
            parasail_raw += length
            parasail_gc += 1

    parasail_aln_len = sum(int(l) for l, op in re.findall(r"(\d+)([=XID])", cigar_stripped))

    log.write(f"  PARASAIL METRICS (after stripping):\n")
    log.write(f"    parasail raw:  {parasail_raw}\n")
    log.write(f"    parasail gc:   {parasail_gc}\n")
    log.write(f"    alignment len: {parasail_aln_len}\n")
    log.write(f"\n")

    # --- parasail NW on edlib-extracted region ---
    db_region = longer[start:end + 1]
    parasail_extracted = parasail.nw_trace_striped_sat(
        shorter, db_region, PARASAIL_GAP_OPEN, PARASAIL_GAP_EXTEND, PARASAIL_MATRIX
    )
    cigar_extracted = parasail_extracted.cigar.decode.decode("utf-8")

    ext_raw = 0
    ext_gc = 0
    for length, op in re.findall(r"(\d+)([=XID])", cigar_extracted):
        length = int(length)
        if op == "X":
            ext_raw += length
            ext_gc += length
        elif op in ("I", "D"):
            ext_raw += length
            ext_gc += 1

    ext_aln_len = sum(int(l) for l, op in re.findall(r"(\d+)([=XID])", cigar_extracted))

    log.write(f"  PARASAIL NW (edlib-extracted region [{start}:{end+1}], {len(db_region)}bp):\n")
    log.write(f"    score: {parasail_extracted.score}\n")
    log.write(f"    parasail raw:  {ext_raw}\n")
    log.write(f"    parasail gc:   {ext_gc}\n")
    log.write(f"    alignment len: {ext_aln_len}\n")
    if len(cigar_extracted) < 500:
        log.write(f"    CIGAR: {cigar_extracted}\n")
    else:
        log.write(f"    CIGAR (first 200): {cigar_extracted[:200]}...\n")
        log.write(f"    CIGAR (last 200):  ...{cigar_extracted[-200:]}\n")
    log.write(f"\n")

    # --- summary ---
    log.write(f"  COMPARISON:\n")
    log.write(f"    {'metric':<25} {'edlib':>10} {'parasail_NW':>15} {'parasail_extract':>18}\n")
    log.write(f"    {'raw edit dist':<25} {edlib_raw:>10} {parasail_raw:>15} {ext_raw:>18}\n")
    log.write(f"    {'gc edit dist':<25} {edlib_gc:>10} {parasail_gc:>15} {ext_gc:>18}\n")
    log.write(f"    {'alignment length':<25} {edlib_aln_len:>10} {parasail_aln_len:>15} {ext_aln_len:>18}\n")
    log.write(f"{'='*80}\n")

# ---------------------------
# Main
# ---------------------------

def main():
    create_HPRC_dict()
    capture = load_capture_haplotypes(CAPTURE_FASTA)

    with open(LOG_FILE, "w") as log:
        log.write("DEBUG ALIGNMENT LOG\n")
        log.write(f"Parasail scoring: match=1, mismatch=-1, gap_open={PARASAIL_GAP_OPEN}, gap_extend={PARASAIL_GAP_EXTEND}\n\n")

        for gene in sorted(capture.keys()):
            for sample in sorted(capture[gene].keys()):
                haps = capture[gene][sample]
                if gene not in HPRC_dict or sample not in HPRC_dict[gene]:
                    continue

                ref_haps = HPRC_dict[gene][sample]

                # Just do first capture hap vs first ref hap for debug
                debug_alignment(
                    haps[0], ref_haps[0],
                    f"{sample} {gene} capture_hap1 vs ref_hap1",
                    log
                )

    print(f"Debug log written to {LOG_FILE}")

if __name__ == "__main__":
    main()
