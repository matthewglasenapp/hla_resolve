"""Analyze 4th-field concordance by truth allele type (:01 vs non-:01)
and tiebreaker resolution level from production hla_typer.py output.

Uses allele_output_trgt.csv (production) + ihw_truth_full_edited.csv (truth)
+ replicates pass 3 tiebreaking to determine resolution level per call.
"""
import sys, os, re
sys.path.insert(0, os.path.dirname(os.path.abspath(__file__)))

import pandas as pd
from hla_typer import (
    build_g_group_dict, produce_allele_seq_db, get_distance, get_sampleid,
)

# ── Config ──────────────────────────────────────────────────────────────
XML = "/Users/matt/Desktop/IPD_IMGT_XML/hla.xml"
FULL_FASTA = "/Users/matt/Desktop/HLA_haplotypes_full.fa"
TRUTH_CSV = "/Users/matt/Desktop/ihw_truth_full_edited.csv"
OUTPUT_CSV = "/Users/matt/Desktop/allele_output_trgt.csv"
SKIP_SAMPLES = {"IHW09117"}
GENES = ["HLA-A", "HLA-B", "HLA-C", "HLA-DPA1",
         "HLA-DPB1", "HLA-DQA1", "HLA-DQB1", "HLA-DRB1"]


def is_called(x):
    return isinstance(x, str) and "*" in x

def truncate(allele, n):
    if not is_called(allele):
        return None
    # Strip HLA- prefix if present
    a = allele.replace("HLA-", "")
    if "*" not in a:
        return None
    head, tail = a.split("*", 1)
    fields = tail.split(":")
    if len(fields) < n:
        return None
    padded = []
    for f in fields[:n]:
        m = re.match(r"(\d+)(.*)", f)
        padded.append(f"{int(m.group(1)):02d}{m.group(2)}" if m else f)
    return f"{head}*{':'.join(padded)}"

def get_4th_field(allele):
    """Extract raw 4th field from allele string."""
    if not is_called(allele):
        return None
    a = allele.replace("HLA-", "")
    if "*" not in a:
        return None
    _, tail = a.split("*", 1)
    fields = tail.split(":")
    if len(fields) < 4:
        return None
    # Return just the numeric part
    m = re.match(r"(\d+)", fields[3])
    return int(m.group(1)) if m else None

def split_options(val):
    if not isinstance(val, str):
        return []
    return [o.strip() for o in re.split(r"[/|]", val) if "*" in o]


# ── Part A: concordance by truth 4th-field value (CSV only) ────────────
print("=" * 70)
print("  PART A: 4-field concordance by truth 4th-field value")
print("=" * 70)

truth_df = pd.read_csv(TRUTH_CSV)
output_df = pd.read_csv(OUTPUT_CSV)

# Build truth dict: truth[sample][gene_hap] = [list of option alleles]
truth = {}
for _, row in truth_df.iterrows():
    sid = str(row["sample"]).strip()
    if sid in SKIP_SAMPLES:
        continue
    for col in truth_df.columns:
        if col in ("sample", "source"):
            continue
        val = str(row[col]).strip()
        if val in ("nan", "NA", "") or "*" not in val:
            continue
        truth.setdefault(sid, {})[col] = split_options(val)

# Build output dict: output[sample][gene_hap] = allele
output = {}
for _, row in output_df.iterrows():
    sid = str(row["sample"]).strip()
    if sid in SKIP_SAMPLES or not sid.startswith("IHW"):
        continue
    for col in output_df.columns:
        if col == "sample":
            continue
        val = str(row[col]).strip()
        if val in ("nan", "NA", ""):
            continue
        output.setdefault(sid, {})[col] = val

# Cross-matching concordance at 4-field, split by truth 4th field
results_01 = {"correct": 0, "wrong": 0, "details": []}
results_non01 = {"correct": 0, "wrong": 0, "details": []}
all_gene_results = {g: {"01_correct": 0, "01_wrong": 0, "non01_correct": 0, "non01_wrong": 0} for g in GENES}

for sid in sorted(set(truth.keys()) & set(output.keys())):
    for gene in GENES:
        gene_short = gene.replace("HLA-", "")
        hap_cols = [f"{gene}_{i}" for i in [1, 2]]

        # Get truth and test alleles for this gene
        truth_alleles = []
        for hc in hap_cols:
            if hc in truth[sid]:
                truth_alleles.append(truth[sid][hc])
            else:
                truth_alleles.append([])

        test_alleles = []
        for hc in hap_cols:
            if hc in output[sid]:
                test_alleles.append(output[sid][hc])
            else:
                test_alleles.append(None)

        # Cross-matching: try both assignments
        # For each truth haplotype, check if any option matches either test allele
        for hap_idx, truth_opts in enumerate(truth_alleles):
            if not truth_opts:
                continue

            # Check truth has 4 fields
            truth_4f_list = [truncate(t, 4) for t in truth_opts if truncate(t, 4)]
            if not truth_4f_list:
                continue

            # Get test 4-field truncations
            test_4f = [truncate(t, 4) for t in test_alleles if t and truncate(t, 4)]
            if not test_4f:
                continue

            # Determine truth 4th field value (use first option)
            truth_4th = get_4th_field(truth_opts[0])
            if truth_4th is None:
                continue

            is_01 = (truth_4th == 1)

            # Cross-match: does any truth option match any test allele?
            matched = any(t4 in test_4f for t4 in truth_4f_list)

            hap_label = f"{sid} {gene_short}_{hap_idx+1}"
            called = test_alleles[hap_idx] if test_alleles[hap_idx] else "N/A"
            called_4f = truncate(called, 4) if called != "N/A" else "N/A"

            if is_01:
                if matched:
                    results_01["correct"] += 1
                    all_gene_results[gene]["01_correct"] += 1
                else:
                    results_01["wrong"] += 1
                    results_01["details"].append(
                        f"  {hap_label}: called={called_4f}, truth={truth_4f_list}")
                    all_gene_results[gene]["01_wrong"] += 1
            else:
                if matched:
                    results_non01["correct"] += 1
                    all_gene_results[gene]["non01_correct"] += 1
                else:
                    results_non01["wrong"] += 1
                    results_non01["details"].append(
                        f"  {hap_label}: called={called_4f}, truth={truth_4f_list}")
                    all_gene_results[gene]["non01_wrong"] += 1

# Print Part A results
n01 = results_01["correct"] + results_01["wrong"]
nn01 = results_non01["correct"] + results_non01["wrong"]

print(f"\n  Truth 4th field = :01")
print(f"    Concordance: {results_01['correct']}/{n01} ({results_01['correct']/n01*100:.1f}%)" if n01 else "    No samples")
print(f"\n  Truth 4th field ≠ :01")
print(f"    Concordance: {results_non01['correct']}/{nn01} ({results_non01['correct']/nn01*100:.1f}%)" if nn01 else "    No samples")

print(f"\n  Per gene breakdown:")
print(f"  {'Gene':<12} {'Truth :01':>20} {'Truth ≠:01':>20}")
print(f"  {'-'*12} {'-'*20} {'-'*20}")
for gene in GENES:
    g = all_gene_results[gene]
    t01 = g["01_correct"] + g["01_wrong"]
    tn01 = g["non01_correct"] + g["non01_wrong"]
    s01 = f"{g['01_correct']}/{t01} ({g['01_correct']/t01*100:.1f}%)" if t01 else "N/A"
    sn01 = f"{g['non01_correct']}/{tn01} ({g['non01_correct']/tn01*100:.1f}%)" if tn01 else "N/A"
    print(f"  {gene:<12} {s01:>20} {sn01:>20}")

print(f"\n  Wrong calls where truth ≠ :01 ({results_non01['wrong']}):")
for d in results_non01["details"]:
    print(d)

print(f"\n  Wrong calls where truth = :01 ({results_01['wrong']}):")
for d in results_01["details"]:
    print(d)


# ── Part B: resolution level for each production call ──────────────────
# Truth-driven: iterate same 162 haplotypes as Part A, with cross-matching.
# For each truth haplotype, find the matching test allele (may be cross-matched
# to the other haplotype position), and determine its resolution level.
print(f"\n{'=' * 70}")
print(f"  PART B: Resolution level per production call")
print(f"  (replicate pass 3 tiebreaking chain on same candidates)")
print(f"{'=' * 70}")

from hla_typer import load_test_data

print("\nLoading sequence data...", flush=True)
_, _, sequence_data = build_g_group_dict(XML, ignore_incomplete=True)
full_samples = {k: v for k, v in load_test_data(FULL_FASTA).items()
                if get_sampleid(k) not in SKIP_SAMPLES}


def get_resolution_for_allele(sid, gene, hap_idx, called_allele, sequence_data, full_samples):
    """Determine MI/ML/iteration resolution level for a specific called allele."""
    gene_short = gene.replace("HLA-", "")
    clean = called_allele.replace("HLA-", "")
    fields = clean.split(":")
    if len(fields) < 4:
        return None  # can't determine resolution for <4 field calls

    three_field = ":".join(fields[:3])

    # Find sample key in full_samples
    sample_key = None
    for k in full_samples:
        if get_sampleid(k) == sid and gene in k and k.endswith(f"_{hap_idx}"):
            sample_key = k
            break
    if sample_key is None:
        return None

    # Find candidates matching 3-field prefix
    candidates = [a for a in sequence_data.keys()
                  if get_distance(three_field, a) == 0]
    if not candidates:
        return None

    candidate_db = produce_allele_seq_db(sequence_data,
                                          selected_alleles=candidates,
                                          exon_only=False)
    query = full_samples[sample_key]

    # Score all candidates
    scores = []
    for allele_name, ref_seq in candidate_db.items():
        if gene_short not in allele_name:
            continue
        ref = str(ref_seq)
        distance, match_len, mismatch_len = get_distance(
            query, ref, get_length=True, gap_compressed=False)
        mi = match_len / (match_len + mismatch_len) if (match_len + mismatch_len) > 0 else 0
        scores.append((allele_name, mi, match_len))

    if not scores:
        return None

    # Current production order: MI → ML → iteration
    best_mi = max(s[1] for s in scores)
    mi_tied = [s for s in scores if s[1] == best_mi]
    best_ml = max(s[2] for s in mi_tied)
    ml_tied = [s for s in mi_tied if s[2] == best_ml]

    if len(mi_tied) == 1:
        resolution = "MI"
    elif len(ml_tied) == 1:
        resolution = "match_length"
    else:
        resolution = "iteration_order"

    # Production winner: first in ml_tied (iteration order)
    prod_winner = truncate(ml_tied[0][0], 4)

    # Chain 1: MI → iteration order (no ML)
    mi_iter_winner = truncate(mi_tied[0][0], 4)

    # Chain 2: ML → iteration order (no MI)
    ml_best = max(s[2] for s in scores)
    ml_tied_all = [s for s in scores if s[2] == ml_best]
    ml_iter_winner = truncate(ml_tied_all[0][0], 4)

    return {
        "resolution": resolution,
        "num_candidates": len(candidates),
        "num_mi_tied": len(mi_tied),
        "num_ml_tied": len(ml_tied),
        "prod_winner": prod_winner,
        "mi_iter_winner": mi_iter_winner,
        "ml_iter_winner": ml_iter_winner,
    }


# Step 1: Compute alt winners for every (sample, gene, hap) query sequence.
# Key: (sid, gene, hap_1indexed) -> res_info dict (or None)
alt_winners = {}
skipped = []

for sid in sorted(set(truth.keys()) & set(output.keys())):
    for gene in GENES:
        gene_short = gene.replace("HLA-", "")
        for hap_1idx in [1, 2]:
            col = f"{gene}_{hap_1idx}"
            called = output.get(sid, {}).get(col)
            if not called or not is_called(called):
                continue
            clean = called.replace("HLA-", "")
            if len(clean.split(":")) < 4:
                # Production output < 4 fields — pass 3 wasn't applied or failed
                alt_winners[(sid, gene, hap_1idx)] = None
                continue
            res_info = get_resolution_for_allele(
                sid, gene, hap_1idx, called, sequence_data, full_samples)
            alt_winners[(sid, gene, hap_1idx)] = res_info

# Step 2: Truth-driven iteration (same 162 as Part A), with proper
# cross-matching for production, MI→iter, and ML→iter.
level_results = []

for sid in sorted(set(truth.keys()) & set(output.keys())):
    for gene in GENES:
        gene_short = gene.replace("HLA-", "")
        hap_cols = [f"{gene}_{i}" for i in [1, 2]]

        # Get truth alleles for this gene
        truth_alleles = []
        for hc in hap_cols:
            if hc in truth.get(sid, {}):
                truth_alleles.append(truth[sid][hc])
            else:
                truth_alleles.append([])

        # Get production test alleles
        test_alleles = []
        for hc in hap_cols:
            if hc in output.get(sid, {}):
                test_alleles.append(output[sid][hc])
            else:
                test_alleles.append(None)

        # Build 4-field lists for production and alt chains
        prod_4f = [truncate(t, 4) for t in test_alleles]  # [pos0, pos1], may be None
        mi_4f = [None, None]
        ml_4f = [None, None]
        res_by_hap = [None, None]
        for h in [0, 1]:
            ri = alt_winners.get((sid, gene, h + 1))
            res_by_hap[h] = ri
            if ri:
                mi_4f[h] = ri["mi_iter_winner"]
                ml_4f[h] = ri["ml_iter_winner"]

        # Skip if no production test alleles have 4 fields
        prod_4f_valid = [p for p in prod_4f if p]
        if not prod_4f_valid:
            continue

        # For each scoreable truth haplotype, do cross-matching
        for hap_idx, truth_opts in enumerate(truth_alleles):
            if not truth_opts:
                continue
            truth_4f_list = [truncate(t, 4) for t in truth_opts if truncate(t, 4)]
            if not truth_4f_list:
                continue
            truth_4th = get_4th_field(truth_opts[0])
            if truth_4th is None:
                continue

            # Production cross-match: truth matches any production test allele?
            prod_matched = any(t4 in prod_4f_valid for t4 in truth_4f_list)

            # MI→iter cross-match: truth matches any MI→iter allele?
            mi_4f_valid = [m for m in mi_4f if m]
            mi_matched = any(t4 in mi_4f_valid for t4 in truth_4f_list)

            # ML→iter cross-match: truth matches any ML→iter allele?
            ml_4f_valid = [m for m in ml_4f if m]
            ml_matched = any(t4 in ml_4f_valid for t4 in truth_4f_list)

            # Determine resolution level using same-pos or cross-matched query
            actual_hap = hap_idx
            same_pos_allele = test_alleles[actual_hap]
            same_pos_4f = truncate(same_pos_allele, 4) if same_pos_allele else None

            if same_pos_4f:
                res_info = res_by_hap[actual_hap]
            else:
                other_hap = 1 - actual_hap
                other_allele = test_alleles[other_hap]
                other_4f = truncate(other_allele, 4) if other_allele else None
                if other_4f:
                    res_info = res_by_hap[other_hap]
                else:
                    res_info = None

            if res_info is None:
                skipped.append(f"{sid} {gene_short}_{actual_hap+1}")
                level_results.append({
                    "sample": sid, "gene": gene, "hap": actual_hap + 1,
                    "truth_4f": truth_4f_list, "truth_4th_field": truth_4th,
                    "is_01": truth_4th == 1, "correct": prod_matched,
                    "resolution": "unknown",
                    "num_mi_tied": 0, "num_ml_tied": 0,
                    "mi_iter_correct": mi_matched,
                    "ml_iter_correct": ml_matched,
                })
                continue

            level_results.append({
                "sample": sid, "gene": gene, "hap": actual_hap + 1,
                "truth_4f": truth_4f_list, "truth_4th_field": truth_4th,
                "is_01": truth_4th == 1, "correct": prod_matched,
                "resolution": res_info["resolution"],
                "num_mi_tied": res_info["num_mi_tied"],
                "num_ml_tied": res_info["num_ml_tied"],
                "mi_iter_correct": mi_matched,
                "ml_iter_correct": ml_matched,
            })

if skipped:
    print(f"\n  WARNING: {len(skipped)} haplotypes could not determine resolution:")
    for s in skipped:
        print(f"    {s}")

# ── Part B summary ─────────────────────────────────────────────────────
total = len(level_results)
print(f"\n  {total} scoreable 4-field haplotypes\n")

# Resolution level counts
for res in ["MI", "match_length", "iteration_order", "unknown"]:
    n = sum(1 for r in level_results if r["resolution"] == res)
    if n == 0:
        continue
    correct = sum(1 for r in level_results if r["resolution"] == res and r["correct"])
    pct = n / total * 100 if total else 0
    cpct = correct / n * 100 if n else 0
    print(f"  {res:<20}: {n:3d}/{total} ({pct:5.1f}%)  concordance: {correct}/{n} ({cpct:.1f}%)")

# Resolution × truth type
print(f"\n  Resolution level × truth 4th-field type:")
print(f"  {'Resolution':<20} {'Truth :01':>25} {'Truth ≠:01':>25}")
print(f"  {'-'*20} {'-'*25} {'-'*25}")
for res in ["MI", "match_length", "iteration_order", "unknown"]:
    n_res = sum(1 for r in level_results if r["resolution"] == res)
    if n_res == 0:
        continue
    n01 = sum(1 for r in level_results if r["resolution"] == res and r["is_01"])
    c01 = sum(1 for r in level_results if r["resolution"] == res and r["is_01"] and r["correct"])
    s01 = f"{c01}/{n01} ({c01/n01*100:.1f}%)" if n01 else "N/A"
    nn = sum(1 for r in level_results if r["resolution"] == res and not r["is_01"])
    cn = sum(1 for r in level_results if r["resolution"] == res and not r["is_01"] and r["correct"])
    sn = f"{cn}/{nn} ({cn/nn*100:.1f}%)" if nn else "N/A"
    print(f"  {res:<20} {s01:>25} {sn:>25}")

# Per-gene resolution breakdown
print(f"\n  Per-gene resolution breakdown:")
for gene in GENES:
    gr = [r for r in level_results if r["gene"] == gene]
    if not gr:
        continue
    n = len(gr)
    mi = sum(1 for r in gr if r["resolution"] == "MI")
    ml = sum(1 for r in gr if r["resolution"] == "match_length")
    io = sum(1 for r in gr if r["resolution"] == "iteration_order")
    unk = sum(1 for r in gr if r["resolution"] == "unknown")
    extra = f" unk={unk}" if unk else ""
    print(f"  {gene:<12}: MI={mi:2d} ML={ml:2d} iter_order={io:2d}{extra} (total {n})")

# Detail: wrong calls with resolution info
print(f"\n  Wrong calls with resolution level:")
for r in sorted(level_results, key=lambda x: (x["gene"], x["sample"])):
    if r["correct"]:
        continue
    tag = ":01" if r["is_01"] else f":{r['truth_4th_field']:02d}"
    print(f"  {r['sample']} {r['gene'].replace('HLA-','')}_{r['hap']}: "
          f"res={r['resolution']:<20} "
          f"truth={r['truth_4f']} (truth_4th={tag}), "
          f"MI_tied={r['num_mi_tied']}, ML_tied={r['num_ml_tied']}")

# ── Part C: Compare two-step chains (with proper cross-matching) ─────
print(f"\n{'=' * 70}")
print(f"  PART C: Two-step tiebreaker comparison")
print(f"  Production (MI→ML→iter) vs MI→iter vs ML→iter")
print(f"  (cross-matching applied to all chains)")
print(f"{'=' * 70}")

total_c = len(level_results)

prod_c = sum(1 for r in level_results if r["correct"])
mi_c = sum(1 for r in level_results if r["mi_iter_correct"])
ml_c = sum(1 for r in level_results if r["ml_iter_correct"])

print(f"\n  Overall ({total_c} haplotypes):")
print(f"    Production (MI→ML→iter):  {prod_c}/{total_c} ({prod_c/total_c*100:.1f}%)")
print(f"    MI → iter order:          {mi_c}/{total_c} ({mi_c/total_c*100:.1f}%)  diff={mi_c-prod_c:+d}")
print(f"    ML → iter order:          {ml_c}/{total_c} ({ml_c/total_c*100:.1f}%)  diff={ml_c-prod_c:+d}")

# Breakdown by truth type
print(f"\n  By truth 4th-field type:")
print(f"  {'':25} {'Production':>15} {'MI→iter':>15} {'ML→iter':>15}")
print(f"  {'-'*25} {'-'*15} {'-'*15} {'-'*15}")
for label, filt in [("Truth = :01", lambda r: r["is_01"]), ("Truth ≠ :01", lambda r: not r["is_01"])]:
    sub = [r for r in level_results if filt(r)]
    n = len(sub)
    if n == 0:
        continue
    pc = sum(1 for r in sub if r["correct"])
    mic = sum(1 for r in sub if r["mi_iter_correct"])
    mlc = sum(1 for r in sub if r["ml_iter_correct"])
    print(f"  {label:<25} {pc:>6}/{n} ({pc/n*100:.1f}%) {mic:>6}/{n} ({mic/n*100:.1f}%) {mlc:>6}/{n} ({mlc/n*100:.1f}%)")

# Per-gene comparison
print(f"\n  Per-gene comparison:")
print(f"  {'Gene':<12} {'Production':>18} {'MI→iter':>18} {'ML→iter':>18}")
print(f"  {'-'*12} {'-'*18} {'-'*18} {'-'*18}")
for gene in GENES:
    gr = [r for r in level_results if r["gene"] == gene]
    if not gr:
        continue
    n = len(gr)
    pc = sum(1 for r in gr if r["correct"])
    mic = sum(1 for r in gr if r["mi_iter_correct"])
    mlc = sum(1 for r in gr if r["ml_iter_correct"])
    print(f"  {gene:<12} {pc:>7}/{n} ({pc/n*100:.1f}%) {mic:>7}/{n} ({mic/n*100:.1f}%) {mlc:>7}/{n} ({mlc/n*100:.1f}%)")

print(f"\n{'=' * 70}")
print(f"  DONE")
print(f"{'=' * 70}")
