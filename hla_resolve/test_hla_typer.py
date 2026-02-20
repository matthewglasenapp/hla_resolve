"""
test_hla_typer.py — Compare parasail vs edlib for pass 3 (4-field selection).

Skips passes 1/2 entirely. Takes the existing pipeline output to get 3-field
prefixes, then re-runs only pass 3 with both aligners against truth.

Usage:
    python test_hla_typer.py \
        --xml /Users/matt/Desktop/IPD_IMGT_XML/hla.xml \
        --full /Users/matt/Downloads/HLA_haplotypes_full.fa \
        --allele-output /Users/matt/Desktop/allele_output_drb.csv \
        --truth /Users/matt/Desktop/ihw_truth_full_edited.csv
"""

import sys
import os
import re
import argparse
import time

sys.path.insert(0, os.path.dirname(os.path.abspath(__file__)))

from Bio import SeqIO
import parasail
import pandas as pd

import edlib
import parasail

from hla_typer import (
    build_g_group_dict,
    produce_allele_seq_db,
    get_distance,
    assign_classification_to_sample_full_seq,
)


# ─────────────────────────────────────────────────────────
# Load existing pipeline output as pass 2 results
# ─────────────────────────────────────────────────────────

def load_pass2_from_output(allele_output_csv):
    """
    Read allele_output.csv (4-field calls from existing pipeline run).
    Return dict: {sample_name: 3-field prefix} for each sample/gene/hap.
    e.g. {"IHW09021_HLA-A_1": "HLA-A*30:01:01"}
    """
    df = pd.read_csv(allele_output_csv)
    pass2 = {}

    for _, row in df.iterrows():
        sample_id = row['sample']
        for col in df.columns:
            if col == 'sample':
                continue
            allele_4field = str(row[col]).strip()
            if allele_4field in ('nan', 'NA', 'None', '') or '*' not in allele_4field:
                continue

            # Build sample_name like "IHW09021_HLA-A_1"
            # col is like "HLA-A_1"
            sample_name = f"{sample_id}_{col}"

            # Strip 4th field to get 3-field prefix
            fields = allele_4field.split(":")
            if len(fields) >= 4:
                three_field = ":".join(fields[:3])
            else:
                three_field = allele_4field  # Already < 4 fields

            pass2[sample_name] = three_field

    return pass2


# ─────────────────────────────────────────────────────────
# Pass 3 candidate generation (shared)
# ─────────────────────────────────────────────────────────

def get_candidates(three_field_prefix, sequence_data):
    """Find all alleles in sequence_data matching the 3-field prefix."""
    matches = []
    for allele in sequence_data.keys():
        dist = get_distance(three_field_prefix, allele)
        if dist == 0:
            matches.append(allele)
    return matches


# ─────────────────────────────────────────────────────────
# Parasail pass 3
# ─────────────────────────────────────────────────────────

def pass_3_parasail(sequence_data, pass2_prefixes, full_samples,
                    gap_open=10, gap_extend=1, matrix=None, quiet=False):
    """Pass 3 using parasail alignment score (maximize)."""
    if not quiet:
        print(f"INFO: Pass 3 (parasail, go={gap_open}, ge={gap_extend})...", flush=True)
    start = time.time()

    if matrix is None:
        matrix = parasail.dnafull

    results = {}

    for sample_name, three_field in pass2_prefixes.items():
        if sample_name not in full_samples:
            continue

        candidates = get_candidates(three_field, sequence_data)
        if not candidates:
            results[sample_name] = (three_field, 0)
            continue

        candidate_db = produce_allele_seq_db(sequence_data, selected_alleles=candidates, exon_only=False)
        query = str(full_samples[sample_name])

        best_allele = None
        best_score = None

        for allele_name, ref_seq in candidate_db.items():
            ref = str(ref_seq)

            # sg_qx: shorter maps within longer, no penalty for reference overhangs
            if len(query) <= len(ref):
                result = parasail.sg_qx_trace_scan_sat(query, ref, gap_open, gap_extend, matrix)
            else:
                result = parasail.sg_qx_trace_scan_sat(ref, query, gap_open, gap_extend, matrix)

            if best_score is None or result.score > best_score:
                best_score = result.score
                best_allele = allele_name

        results[sample_name] = (best_allele, best_score)

    if not quiet:
        print(f"INFO: Pass 3 (parasail) done in {time.time()-start:.0f}s", flush=True)
    return results


# ─────────────────────────────────────────────────────────
# Hybrid edlib+parasail pass 3 (gap-compressed)
# ─────────────────────────────────────────────────────────

PARASAIL_MATRIX = parasail.matrix_create("ACGT", 1, -1)

def compute_parasail_gc(query, ref, gap_open=5, gap_extend=1):
    """Two-stage: edlib HW to locate region, parasail NW for affine-gap scoring."""
    if len(query) <= len(ref):
        shorter, longer = query, ref
    else:
        shorter, longer = ref, query

    # Stage 1: edlib HW to find mapping location
    edlib_result = edlib.align(shorter, longer, mode="HW", task="path")
    start, end = edlib_result["locations"][0]
    db_region = longer[start:end + 1]

    # Stage 2: parasail NW on extracted region
    parasail_result = parasail.nw_trace_striped_sat(
        shorter, db_region, gap_open, gap_extend, PARASAIL_MATRIX
    )
    cigar = parasail_result.cigar.decode.decode("utf-8")

    # Gap-compressed: mismatches count individually, each indel run = 1
    parasail_gc = 0
    for length, op in re.findall(r"(\d+)([=XID])", cigar):
        length = int(length)
        if op == "X":
            parasail_gc += length
        elif op in ("I", "D"):
            parasail_gc += 1

    aln_len = sum(int(l) for l, op in re.findall(r"(\d+)([=XID])", cigar))
    return parasail_gc, aln_len


def pass_3_edlib_parasail(sequence_data, pass2_prefixes, full_samples,
                          gap_open=5, gap_extend=1, quiet=False):
    """Pass 3 using edlib+parasail hybrid: edlib locates region, parasail NW scores it."""
    if not quiet:
        print(f"INFO: Pass 3 (edlib+parasail, go={gap_open}, ge={gap_extend})...", flush=True)
    start = time.time()

    results = {}

    for sample_name, three_field in pass2_prefixes.items():
        if sample_name not in full_samples:
            continue

        candidates = get_candidates(three_field, sequence_data)
        if not candidates:
            results[sample_name] = (three_field, 0, 0, 0, 0, [])
            continue

        candidate_db = produce_allele_seq_db(sequence_data, selected_alleles=candidates, exon_only=False)
        query = str(full_samples[sample_name])

        best_allele = None
        best_gc = None
        best_aln_len = 0
        same_dist = []

        for allele_name, ref_seq in candidate_db.items():
            ref = str(ref_seq)
            gc, aln_len = compute_parasail_gc(query, ref, gap_open, gap_extend)

            if best_gc is None or gc < best_gc:
                best_gc = gc
                best_allele = allele_name
                best_aln_len = aln_len
                same_dist = [allele_name]
            elif gc == best_gc:
                # Tiebreaker: prefer longer alignment
                if aln_len > best_aln_len:
                    best_allele = allele_name
                    best_aln_len = aln_len
                    same_dist = [allele_name]
                elif aln_len == best_aln_len:
                    same_dist.append(allele_name)

        identity = 1 - (best_gc / best_aln_len) if best_aln_len > 0 else 0
        results[sample_name] = (best_allele, best_gc, best_aln_len, identity, 0, same_dist)

    if not quiet:
        print(f"INFO: Pass 3 (edlib+parasail) done in {time.time()-start:.0f}s", flush=True)
    return results


# ─────────────────────────────────────────────────────────
# Edlib pass 3 (baseline)
# ─────────────────────────────────────────────────────────

def pass_3_edlib(sequence_data, pass2_prefixes, full_samples, metric="mismatch_identity", quiet=False):
    """Pass 3 using edlib + mismatch_identity (current production logic)."""
    if not quiet:
        print(f"INFO: Pass 3 (edlib, metric={metric})...", flush=True)
    start = time.time()

    results = {}

    for sample_name, three_field in pass2_prefixes.items():
        if sample_name not in full_samples:
            continue

        candidates = get_candidates(three_field, sequence_data)
        if not candidates:
            results[sample_name] = (three_field, 0, 0, 0, 0, [])
            continue

        candidate_db = produce_allele_seq_db(sequence_data, selected_alleles=candidates, exon_only=False)

        result = assign_classification_to_sample_full_seq(
            candidate_db, full_samples[sample_name], sample_name, eval_metric=metric
        )
        results[sample_name] = result

    if not quiet:
        print(f"INFO: Pass 3 (edlib) done in {time.time()-start:.0f}s", flush=True)
    return results


# ─────────────────────────────────────────────────────────
# Truth comparison
# ─────────────────────────────────────────────────────────

PROBLEMATIC_SAMPLES = ["IHW09117"]


def truncate_allele(allele, n):
    if not isinstance(allele, str) or '*' not in allele:
        return None
    head, tail = allele.split('*', 1)
    head = head.replace('HLA-', '')
    fields = tail.split(':')
    if len(fields) < n:
        return None
    padded = []
    for f in fields[:n]:
        m = re.match(r'(\d+)(.*)', f)
        if m:
            padded.append(f"{int(m.group(1)):02d}{m.group(2)}")
        else:
            padded.append(f)
    return f"{head}*{':'.join(padded)}"


def load_truth(path):
    truth_df = pd.read_csv(path)
    truth = {}
    for _, row in truth_df.iterrows():
        sample = row['sample']
        if sample in PROBLEMATIC_SAMPLES:
            continue
        for col in truth_df.columns:
            if col == 'sample':
                continue
            parts = col.rsplit('_', 1)
            gene_short = parts[0].replace('HLA-', '')
            hap = parts[1]
            val = str(row[col]).strip()
            if val in ('nan', 'NA', 'UNKNOWN', 'no call', 'None', ''):
                continue
            if '*' not in val:
                continue
            options = re.split(r'[/|]', val)
            options = [o.strip() for o in options if '*' in o]
            truth[(sample, gene_short, hap)] = options
    return truth


def score_concordance(results_dict, truth, n_fields=4):
    from collections import defaultdict, Counter

    called = {}
    for sample_name, result in results_dict.items():
        allele = result[0]
        if allele is None:
            continue
        parts = sample_name.split('_')
        sample_id = parts[0]
        gene = '_'.join(parts[1:-1]).replace('HLA-', '')
        hap = parts[-1]
        if sample_id in PROBLEMATIC_SAMPLES:
            continue
        called[(sample_id, gene, hap)] = allele

    sample_genes = defaultdict(dict)
    for (sample, gene, hap), allele in called.items():
        sample_genes[(sample, gene)][hap] = allele

    correct = 0
    total = 0
    gene_correct = Counter()
    gene_total = Counter()
    wrong = []

    for (sample, gene), hap_alleles in sample_genes.items():
        truth_alleles_all = []
        for hap in ['1', '2']:
            tkey = (sample, gene, hap)
            if tkey in truth:
                truth_alleles_all.append((hap, truth[tkey]))

        test_alleles = [hap_alleles[h] for h in sorted(hap_alleles.keys())]
        if not truth_alleles_all or not test_alleles:
            continue

        used_tests = set()
        for hap, truth_options in truth_alleles_all:
            truth_truncs = [truncate_allele(t, n_fields) for t in truth_options if truncate_allele(t, n_fields)]
            if not truth_truncs:
                continue

            total += 1
            gene_total[gene] += 1
            test_truncs = [truncate_allele(a, n_fields) for a in test_alleles if truncate_allele(a, n_fields)]

            all_truth_truncs = []
            for _, topts in truth_alleles_all:
                for t in topts:
                    tt = truncate_allele(t, n_fields)
                    if tt:
                        all_truth_truncs.append(tt)

            matched = False
            for i, test_allele in enumerate(test_alleles):
                test_trunc = truncate_allele(test_allele, n_fields)
                if not test_trunc:
                    continue
                if all_truth_truncs.count(test_trunc) >= 2 and test_truncs.count(test_trunc) >= 2:
                    if test_trunc in truth_truncs:
                        matched = True
                        break
                elif i not in used_tests and test_trunc in truth_truncs:
                    matched = True
                    used_tests.add(i)
                    break

            if matched:
                correct += 1
                gene_correct[gene] += 1
            else:
                called_trunc = truncate_allele(hap_alleles.get(hap, '?'), n_fields) or hap_alleles.get(hap, '?')
                wrong.append(f"  {sample} {gene}_{hap}: called={called_trunc}, truth={truth_options}")

    return correct, total, gene_correct, gene_total, wrong


def print_concordance(label, correct, total, gene_correct, gene_total, wrong):
    pct = correct / total * 100 if total > 0 else 0
    print(f"\n{'='*60}")
    print(f"  {label}")
    print(f"  4-field concordance: {correct}/{total} = {pct:.1f}%")
    print(f"{'='*60}")
    for gene in sorted(gene_total.keys()):
        c = gene_correct[gene]
        t = gene_total[gene]
        p = c / t * 100 if t > 0 else 0
        print(f"  {gene:<8}: {c}/{t} = {p:.1f}%")
    if wrong:
        print(f"\n  Wrong calls ({len(wrong)}):")
        for w in sorted(wrong):
            print(w)


# ─────────────────────────────────────────────────────────
# Main
# ─────────────────────────────────────────────────────────

def main():
    parser = argparse.ArgumentParser(description='Test parasail vs edlib for HLA pass 3')
    parser.add_argument('--xml', required=True, help='IMGT XML reference file')
    parser.add_argument('--full', required=True, help='Sample full-sequence FASTA')
    parser.add_argument('--allele-output', required=True, help='Existing pipeline allele_output.csv (provides 3-field prefixes)')
    parser.add_argument('--truth', required=True, help='Truth CSV file')
    parser.add_argument('--gap-open', type=int, default=10, help='Parasail gap open penalty (default: 10)')
    parser.add_argument('--gap-extend', type=int, default=1, help='Parasail gap extend penalty (default: 1)')
    parser.add_argument('--sweep', action='store_true', help='Sweep gap parameters to find optimal')
    args = parser.parse_args()

    # ── Load XML reference (needed for candidate sequences) ──
    print("INFO: Loading XML reference...", flush=True)
    g_group_dict, p_group_dict, sequence_data = build_g_group_dict(args.xml)
    print(f"INFO: {len(sequence_data)} reference alleles loaded", flush=True)

    # ── Load full-sequence samples ──
    print("INFO: Loading full-sequence FASTA...", flush=True)
    full_samples = {}
    for record in SeqIO.parse(args.full, "fasta"):
        full_samples[record.id] = record.seq
    print(f"INFO: {len(full_samples)} full samples loaded", flush=True)

    # ── Load 3-field prefixes from existing pipeline output ──
    print("INFO: Loading 3-field prefixes from allele output...", flush=True)
    pass2_prefixes = load_pass2_from_output(args.allele_output)
    print(f"INFO: {len(pass2_prefixes)} pass 2 entries loaded", flush=True)

    # ── Load truth ──
    truth = load_truth(args.truth)
    print(f"INFO: {len(truth)} truth entries loaded", flush=True)

    # ── Edlib baseline ──
    edlib_results = pass_3_edlib(sequence_data, pass2_prefixes, full_samples, metric="mismatch_identity")
    c, t, gc, gt, w = score_concordance(edlib_results, truth)
    print_concordance("EDLIB (mismatch_identity)", c, t, gc, gt, w)

    # ── Hybrid edlib+parasail ──
    hybrid_results = pass_3_edlib_parasail(sequence_data, pass2_prefixes, full_samples,
                                           gap_open=5, gap_extend=1)
    c, t, gc, gt, w = score_concordance(hybrid_results, truth)
    print_concordance("EDLIB+PARASAIL hybrid (go=5, ge=1, minimize parasail_gc)", c, t, gc, gt, w)

    # ── Parasail ──
    if args.sweep:
        print(f"\n{'='*60}")
        print("  PARAMETER SWEEP")
        print(f"{'='*60}")
        print(f"{'gap_open':>10} {'gap_extend':>10} {'correct':>8} {'total':>6} {'pct':>7}")
        print("-" * 45)

        best_params = None
        best_correct = 0

        for go in [3, 5, 6, 8, 10, 12, 15, 20]:
            for ge in [1, 2, 3]:
                pr = pass_3_parasail(sequence_data, pass2_prefixes, full_samples,
                                     gap_open=go, gap_extend=ge, quiet=True)
                c, t, gc, gt, w = score_concordance(pr, truth)
                pct = c / t * 100 if t > 0 else 0
                print(f"{go:>10} {ge:>10} {c:>8} {t:>6} {pct:>6.1f}%", flush=True)
                if c > best_correct:
                    best_correct = c
                    best_params = (go, ge)

        print(f"\nBest: gap_open={best_params[0]}, gap_extend={best_params[1]} -> {best_correct}/{t}")

        # Detailed results for best
        pr = pass_3_parasail(sequence_data, pass2_prefixes, full_samples,
                             gap_open=best_params[0], gap_extend=best_params[1])
        c, t, gc, gt, w = score_concordance(pr, truth)
        print_concordance(f"PARASAIL (best: go={best_params[0]}, ge={best_params[1]})", c, t, gc, gt, w)
    else:
        pr = pass_3_parasail(sequence_data, pass2_prefixes, full_samples,
                             gap_open=args.gap_open, gap_extend=args.gap_extend)
        c, t, gc, gt, w = score_concordance(pr, truth)
        print_concordance(f"PARASAIL (go={args.gap_open}, ge={args.gap_extend})", c, t, gc, gt, w)


if __name__ == "__main__":
    main()
