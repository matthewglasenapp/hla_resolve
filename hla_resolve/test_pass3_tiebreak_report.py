"""Pass 3 tiebreaking analysis across all HLA genes (IHW only, excl IHW09117).

Compares four pass 3 scoring chains (all use edlib, iteration order as final tiebreaker):
  1. MI → ML → iter  (production baseline: mismatch_identity → match_length → XML order)
  2. ML → iter       (match_length primary, XML order secondary)
  3. MI → iter       (mismatch_identity primary only, XML order secondary — no ML)
  4. ML → MI → iter  (match_length primary, mismatch_identity secondary, XML order tertiary)
"""
import sys, os, re, time
from collections import defaultdict
sys.path.insert(0, os.path.dirname(os.path.abspath(__file__)))

import pandas as pd
from hla_typer import (
    build_g_group_dict, get_g_group_exons, load_test_data,
    pass_1_classification, pass_2_classification,
    produce_allele_seq_db, get_distance, get_sampleid,
)

# ── Config ──
ALL_GENES = ["HLA-A", "HLA-B", "HLA-C", "HLA-DPA1",
             "HLA-DPB1", "HLA-DQA1", "HLA-DQB1", "HLA-DRB1"]
XML       = "/Users/matt/Desktop/IPD_IMGT_XML/hla.xml"
EXON_FASTA = "/Users/matt/Desktop/HLA_haplotypes.fa"
FULL_FASTA = "/Users/matt/Desktop/HLA_haplotypes_full.fa"
TRUTH_CSV  = "/Users/matt/Desktop/ihw_truth_full_edited.csv"
SKIP_SAMPLES = {"IHW09117"}


# ── Helpers ──
def truncate(allele, n):
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
        padded.append(f"{int(m.group(1)):02d}{m.group(2)}" if m else f)
    return f"{head}*{':'.join(padded)}"


def _score_all(query, candidate_db):
    """Compute (allele_name, match_len, mi) for every candidate in iteration order."""
    scores = []
    for allele_name, ref_seq in candidate_db.items():
        ref = str(ref_seq)
        distance, match_len, mismatch_len = get_distance(
            query, ref, get_length=True, gap_compressed=False)
        mi = match_len / (match_len + mismatch_len) if (match_len + mismatch_len) > 0 else 0
        scores.append((allele_name, match_len, mi))
    return scores  # preserved in XML iteration order


# ── Scoring functions (iteration order as final tiebreaker) ──

def chain_mi_ml_iter(query, candidate_db, sample_name):
    """MI → ML → iter  (production chain)"""
    scores = _score_all(query, candidate_db)
    best_mi = max(s[2] for s in scores)
    mi_group = [s for s in scores if s[2] == best_mi]      # preserves XML order
    best_ml = max(s[1] for s in mi_group)
    ml_group = [s for s in mi_group if s[1] == best_ml]    # preserves XML order
    selected = ml_group[0][0]                               # first = XML order winner
    tie_info = {"num_tied": len(ml_group)}
    return selected, tie_info


def chain_ml_iter(query, candidate_db, sample_name):
    """ML → iter  (match_length primary, XML order secondary)"""
    scores = _score_all(query, candidate_db)
    best_ml = max(s[1] for s in scores)
    ml_group = [s for s in scores if s[1] == best_ml]      # preserves XML order
    selected = ml_group[0][0]
    tie_info = {"num_tied": len(ml_group)}
    return selected, tie_info


def chain_mi_iter(query, candidate_db, sample_name):
    """MI → iter  (mismatch_identity primary only, XML order secondary — no ML)"""
    scores = _score_all(query, candidate_db)
    best_mi = max(s[2] for s in scores)
    mi_group = [s for s in scores if s[2] == best_mi]      # preserves XML order
    selected = mi_group[0][0]
    tie_info = {"num_tied": len(mi_group)}
    return selected, tie_info


def chain_ml_mi_iter(query, candidate_db, sample_name):
    """ML → MI → iter  (match_length primary, mismatch_identity secondary, XML order tertiary)"""
    scores = _score_all(query, candidate_db)
    best_ml = max(s[1] for s in scores)
    ml_group = [s for s in scores if s[1] == best_ml]      # preserves XML order
    best_mi = max(s[2] for s in ml_group)
    mi_group = [s for s in ml_group if s[2] == best_mi]    # preserves XML order
    selected = mi_group[0][0]
    tie_info = {"num_tied": len(mi_group)}
    return selected, tie_info


METHODS = [
    ("mi_ml_iter",  "MI→ML→iter (production)", chain_mi_ml_iter),
    ("ml_iter",     "ML→iter",                 chain_ml_iter),
    ("mi_iter",     "MI→iter",                 chain_mi_iter),
    ("ml_mi_iter",  "ML→MI→iter",              chain_ml_mi_iter),
]


# ── Pass 3: precompute candidates, then score with each method ──
def precompute_gene(gene, pass2, seq_data, full_samps):
    """Build candidate_db for each IHW sample once per gene."""
    gene_pass2 = {
        k: v for k, v in pass2.items()
        if gene in k and get_sampleid(k).startswith("IHW")
    }
    precomputed = {}
    for sample_name, classification in gene_pass2.items():
        classified_allele = classification[0]
        if classified_allele is None:
            continue
        fields = classified_allele.split(":")
        if len(fields) > 3:
            three_field = ":".join(fields[:-1])
        else:
            precomputed[sample_name] = {
                "type": "short_circuit", "allele": classified_allele,
                "candidates": [], "candidate_db": None, "query": None,
            }
            continue
        candidates = [a for a in seq_data.keys() if get_distance(three_field, a) == 0]
        if not candidates:
            precomputed[sample_name] = {
                "type": "no_candidates", "allele": classified_allele,
                "candidates": [], "candidate_db": None, "query": None,
            }
            continue
        candidate_db = produce_allele_seq_db(seq_data, selected_alleles=candidates, exon_only=False)
        if sample_name not in full_samps:
            precomputed[sample_name] = {
                "type": "no_full_seq", "allele": classified_allele,
                "candidates": candidates, "candidate_db": candidate_db, "query": None,
            }
            continue
        precomputed[sample_name] = {
            "type": "score", "allele": classified_allele,
            "candidates": candidates, "candidate_db": candidate_db,
            "query": full_samps[sample_name],
        }
    return precomputed


def score_gene(precomputed, score_fn):
    """Score all samples for a gene using one scoring chain."""
    results = {}
    for sample_name, pre in precomputed.items():
        if pre["type"] == "score":
            selected, tie_info = score_fn(pre["query"], pre["candidate_db"], sample_name)
            results[sample_name] = {"allele": selected, "tie_info": tie_info}
        else:
            results[sample_name] = {
                "allele": pre["allele"],
                "tie_info": {"num_tied": 1},
            }
    return results


# ── Cross-matching evaluation (matches alex_test_against_truth logic) ──
def evaluate_concordance(gene, results, truth_by_sample):
    n = 4
    correct, total = 0, 0
    wrong_details = []
    for sid, hap_truth in truth_by_sample.items():
        called_truth, called_test = [], []
        for hap_idx in ['1', '2']:
            if hap_idx in hap_truth:
                called_truth.append(hap_truth[hap_idx])
            key = f"{sid}_{gene}_{hap_idx}"
            if key in results:
                called_test.append(results[key]["allele"])
        if not called_truth or not called_test:
            continue
        test_truncs = [truncate(a, n) for a in called_test if truncate(a, n)]
        used_tests = set()
        for truth_options in called_truth:
            sample_truth = truth_options[0]
            if not isinstance(sample_truth, str) or '*' not in sample_truth:
                continue
            if len(sample_truth.split('*', 1)[1].split(':')) < n:
                continue
            truth_truncs = [truncate(t, n) for t in truth_options if truncate(t, n)]
            if not truth_truncs:
                continue
            total += 1
            matched = False
            for i, test_a in enumerate(called_test):
                test_trunc = truncate(test_a, n)
                if not test_trunc:
                    continue
                truth_all_truncs = []
                for t_opts in called_truth:
                    if t_opts:
                        for t in t_opts:
                            tr = truncate(t, n)
                            if tr:
                                truth_all_truncs.append(tr)
                if truth_all_truncs.count(test_trunc) == 2 and test_truncs.count(test_trunc) == 2:
                    if test_trunc in truth_truncs:
                        matched = True
                        break
                elif i not in used_tests and test_trunc in truth_truncs:
                    matched = True
                    used_tests.add(i)
                    break
            if matched:
                correct += 1
            else:
                wrong_details.append({
                    "sid": sid, "truth_options": truth_options,
                    "called_test": called_test,
                    "truth_truncs": truth_truncs,
                    "called_truncs": [truncate(a, n) for a in called_test],
                })
    return correct, total, wrong_details


# ════════════════════════════════════════════════════════════════════
#  MAIN
# ════════════════════════════════════════════════════════════════════
if __name__ == "__main__":
    print("=" * 70)
    print("  PASS 3 TIEBREAKING CHAIN COMPARISON — ALL GENES, IHW ONLY")
    print("  Chains: MI→ML→iter | ML→iter | MI→iter | ML→MI→iter")
    print("=" * 70, flush=True)

    # Load truth
    print("\nLoading truth data...", flush=True)
    truth_df = pd.read_csv(TRUTH_CSV)
    truth_all = {}
    for _, row in truth_df.iterrows():
        sid = row['sample']
        if sid in SKIP_SAMPLES:
            continue
        for col in truth_df.columns:
            if col in ('sample', 'source'):
                continue
            val = str(row[col]).strip()
            if val in ('nan', 'NA', '') or '*' not in val:
                continue
            options = [o.strip() for o in re.split(r'[/|]', val) if '*' in o]
            parts = col.rsplit('_', 1)
            if len(parts) != 2:
                continue
            gene_prefix, hap = parts
            truth_all.setdefault(gene_prefix, {}).setdefault(sid, {})[hap] = options
    for g in ALL_GENES:
        print(f"  {g}: {len(truth_all.get(g, {}))} samples with truth")

    # Load reference
    print("\nLoading XML (ignore_incomplete=True)...", flush=True)
    t0 = time.time()
    g_group_dict, _, sequence_data = build_g_group_dict(XML, ignore_incomplete=True)
    print(f"  {len(sequence_data)} alleles in {time.time()-t0:.0f}s")
    print("Building G-group exon sequences...", flush=True)
    g_group_common_seqs = get_g_group_exons(g_group_dict, sequence_data)

    # Load FASTAs
    print("Loading exon FASTA...", flush=True)
    exon_samples = load_test_data(EXON_FASTA)
    exon_samples = {k: v for k, v in exon_samples.items()
                    if get_sampleid(k) not in SKIP_SAMPLES}
    print(f"  {len(exon_samples)} samples")
    print("Loading full FASTA...", flush=True)
    full_samples = load_test_data(FULL_FASTA)
    full_samples = {k: v for k, v in full_samples.items()
                    if get_sampleid(k) not in SKIP_SAMPLES}
    print(f"  {len(full_samples)} samples")

    # Pass 1 + 2
    print("\nRunning pass 1...", flush=True)
    t0 = time.time()
    g_group_results = pass_1_classification(g_group_common_seqs, exon_samples, g_group_dict)
    print(f"  Done in {time.time()-t0:.0f}s")
    print("Running pass 2...", flush=True)
    t0 = time.time()
    pass2_results = pass_2_classification(
        sequence_data, g_group_dict, g_group_results, exon_samples, metric="edit_distance")
    print(f"  Done in {time.time()-t0:.0f}s")

    # Run all genes × methods
    all_concordance = {}
    grand_t0 = time.time()

    for gene in ALL_GENES:
        print(f"\n{'='*60}")
        print(f"  {gene}", flush=True)
        print(f"{'='*60}")

        truth_gene = truth_all.get(gene, {})
        all_concordance[gene] = {}

        # Precompute candidate_dbs once per gene
        print(f"  Precomputing candidates...", flush=True)
        t0 = time.time()
        precomputed = precompute_gene(gene, pass2_results, sequence_data, full_samples)
        n_scoreable = sum(1 for p in precomputed.values() if p["type"] == "score")
        print(f"    {len(precomputed)} samples, {n_scoreable} scoreable in {time.time()-t0:.0f}s")

        for mk, mn, score_fn in METHODS:
            print(f"  {mn}...", end="", flush=True)
            t0 = time.time()
            results = score_gene(precomputed, score_fn)
            elapsed = time.time() - t0
            c, t, wd = evaluate_concordance(gene, results, truth_gene)
            pct = c / t * 100 if t > 0 else 0
            print(f" {c}/{t} ({pct:.1f}%) in {elapsed:.0f}s", flush=True)
            all_concordance[gene][mk] = (c, t, wd)

    total_elapsed = time.time() - grand_t0
    print(f"\n  Total pass 3 time: {total_elapsed:.0f}s")

    # ════════════════════════════════════════════════════════════════
    #  REPORT
    # ════════════════════════════════════════════════════════════════
    print(f"\n\n{'#'*70}")
    print(f"  4-FIELD CONCORDANCE: GENE × CHAIN")
    print(f"{'#'*70}\n")

    # Header
    col_w = 22
    header = f"  {'Gene':<12}"
    for _, mn, _ in METHODS:
        header += f" {mn:<{col_w}}"
    print(header)
    print(f"  {'-'*12} " + " ".join(['-'*col_w]*len(METHODS)))

    for gene in ALL_GENES:
        row = f"  {gene:<12}"
        for mk, _, _ in METHODS:
            c, t, _ = all_concordance[gene][mk]
            pct = c / t * 100 if t > 0 else 0
            cell = f"{c}/{t} ({pct:.1f}%)"
            row += f" {cell:<{col_w}}"
        print(row)

    print(f"  {'-'*12} " + " ".join(['-'*col_w]*len(METHODS)))
    row = f"  {'TOTAL':<12}"
    totals = {}
    for mk, _, _ in METHODS:
        tc = sum(all_concordance[g][mk][0] for g in ALL_GENES)
        tt = sum(all_concordance[g][mk][1] for g in ALL_GENES)
        pct = tc / tt * 100 if tt > 0 else 0
        totals[mk] = (tc, tt)
        cell = f"{tc}/{tt} ({pct:.1f}%)"
        row += f" {cell:<{col_w}}"
    print(row)

    # Diff vs production
    print(f"\n  Difference vs MI→ML→iter (production):")
    base_c, base_t = totals["mi_ml_iter"]
    for mk, mn, _ in METHODS:
        if mk == "mi_ml_iter":
            continue
        tc, tt = totals[mk]
        diff = tc - base_c
        sign = "+" if diff >= 0 else ""
        print(f"    {mn:<22}: {sign}{diff:+d} ({tc}/{tt} vs {base_c}/{base_t})")

    # Per-gene diffs
    print(f"\n  Per-gene differences vs production:")
    print(f"  {'Gene':<12}", end="")
    for mk, mn, _ in METHODS:
        if mk != "mi_ml_iter":
            print(f" {mn:<22}", end="")
    print()
    for gene in ALL_GENES:
        base_c2, base_t2, _ = all_concordance[gene]["mi_ml_iter"]
        row = f"  {gene:<12}"
        any_diff = False
        for mk, _, _ in METHODS:
            if mk == "mi_ml_iter":
                continue
            c, t, _ = all_concordance[gene][mk]
            diff = c - base_c2
            sign = "+" if diff >= 0 else ""
            cell = f"{sign}{diff} ({c}/{t})"
            row += f" {cell:<22}"
            if diff != 0:
                any_diff = True
        print(row + (" *" if any_diff else ""))

    print(f"\n{'#'*70}")
    print(f"  END OF REPORT (total time: {time.time()-grand_t0:.0f}s)")
    print(f"{'#'*70}")
