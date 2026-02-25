"""Analyze 4th-field concordance by truth allele type (:01 vs non-:01)
and tiebreaker resolution level from production hla_typer.py output.

Uses allele_output_trgt.csv (production) + ihw_truth_full_edited.csv (truth)
+ replicates pass 3 tiebreaking to determine resolution level per call.
"""
import sys, os, re
sys.path.insert(0, os.path.dirname(os.path.abspath(__file__)))

import edlib
import parasail
import pandas as pd
from hla_typer import (
    build_g_group_dict, produce_allele_seq_db, get_distance, get_sampleid,
)

# ── Parasail config ──────────────────────────────────────────────────────
_PARASAIL_MATRIX     = parasail.matrix_create("ACGT", 1, -1)
_PARASAIL_GAP_OPEN   = 5
_PARASAIL_GAP_EXTEND = 1

def get_distance_hybrid(query, ref):
    """
    Hybrid edlib + parasail alignment.
      1. edlib HW (shorter fully within longer) locates the best window.
      2. parasail NW with affine gap runs on that extracted window,
         producing a biologically realistic CIGAR.
    Returns:
        p_mi        – match_len / (match_len + mismatch_len)  [higher = better]
        p_match_len – count of '=' bases                      [higher = better]
        p_gc        – mismatches + indel run count             [lower  = better]
        p_score     – parasail alignment score                 [higher = better]
    """
    query_s = str(query)
    ref_s   = str(ref)

    if len(query_s) <= len(ref_s):
        shorter, longer = query_s, ref_s
    else:
        shorter, longer = ref_s, query_s

    edlib_result = edlib.align(shorter, longer, mode="HW", task="path")
    start, end   = edlib_result["locations"][0]
    region       = longer[start:end + 1]

    ps = parasail.nw_trace_striped_sat(
        shorter, region,
        _PARASAIL_GAP_OPEN, _PARASAIL_GAP_EXTEND, _PARASAIL_MATRIX
    )
    cigar = ps.cigar.decode.decode("utf-8")
    score = ps.score

    match_len = mismatch_len = p_gc = 0
    for length, op in re.findall(r"(\d+)([=XID])", cigar):
        length = int(length)
        if   op == "=": match_len     += length
        elif op == "X": mismatch_len  += length; p_gc += length
        elif op in ("I", "D"): p_gc  += 1        # each indel run = 1

    p_mi = match_len / (match_len + mismatch_len) if (match_len + mismatch_len) > 0 else 0
    return p_mi, match_len, p_gc, score

# ── Config ──────────────────────────────────────────────────────────────
XML = "/Users/matt/Desktop/IPD_IMGT_XML/hla.xml"
FULL_FASTA = "/Users/matt/Desktop/HLA_haplotypes_full.fa"
TRUTH_CSV = "/Users/matt/Desktop/ihw_truth_full_edited.csv"
OUTPUT_CSV = "/Users/matt/Desktop/allele_output_latest.csv"
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


def strict_cross_match(truth_4f_list, test_4f_indexed, test_4f_valid,
                       all_truth_first_4f, used_test_indices):
    """Cross-match one truth haplotype against test alleles with used-test tracking.

    Matches alex_test_against_truth.py logic: each test allele index can only be
    consumed once, UNLESS both truth haplotypes and both test alleles share the
    same 4-field value (homozygous exception).

    Args:
        truth_4f_list: truncated 4-field options for this truth haplotype
        test_4f_indexed: list of 4-field values per test position (may contain None)
        test_4f_valid: list of non-None 4-field test values
        all_truth_first_4f: first-option 4-field value per truth haplotype (for homo check)
        used_test_indices: set of already-consumed test indices (modified in-place)
    Returns:
        True if matched.
    """
    for i, test_val in enumerate(test_4f_indexed):
        if test_val is None:
            continue
        # Homozygous exception: both truth haps AND both test alleles have this value
        if (all_truth_first_4f.count(test_val) >= 2
                and test_4f_valid.count(test_val) >= 2
                and test_val in truth_4f_list):
            return True
        elif i not in used_test_indices and test_val in truth_4f_list:
            used_test_indices.add(i)
            return True
    return False


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

        # Pre-compute test 4-field values with indices for strict tracking
        test_4f_indexed = [truncate(ta, 4) if ta else None for ta in test_alleles]
        test_4f_valid = [t for t in test_4f_indexed if t]

        # Pre-compute all truth first-option 4f values for homozygous check
        all_truth_first_4f = []
        for ta_opts in truth_alleles:
            if ta_opts:
                all_truth_first_4f.append(truncate(ta_opts[0], 4))
            else:
                all_truth_first_4f.append(None)

        # Cross-matching with proper used-test tracking
        used_test_indices = set()
        for hap_idx, truth_opts in enumerate(truth_alleles):
            if not truth_opts:
                continue

            # Check truth has 4 fields
            truth_4f_list = [truncate(t, 4) for t in truth_opts if truncate(t, 4)]
            if not truth_4f_list:
                continue

            if not test_4f_valid:
                continue

            # Determine truth 4th field value (use first option)
            truth_4th = get_4th_field(truth_opts[0])
            if truth_4th is None:
                continue

            is_01 = (truth_4th == 1)

            # Cross-match with used-test index tracking
            matched = strict_cross_match(truth_4f_list, test_4f_indexed,
                                         test_4f_valid, all_truth_first_4f,
                                         used_test_indices)

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
                if get_sampleid(k) not in SKIP_SAMPLES
                and get_sampleid(k).startswith("IHW")}


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

    # Score all candidates — edlib metrics + parasail hybrid metrics
    # Tuple: (allele_name, e_mi, e_ml, e_raw, e_gc, p_mi, p_ml, p_gc, p_score)
    scores = []
    for allele_name, ref_seq in candidate_db.items():
        if gene_short not in allele_name:
            continue
        ref = str(ref_seq)
        _, _, _, e_raw, e_gc, e_mi, e_ml = get_distance(query, ref, get_detailed=True)
        p_mi, p_ml, p_gc, p_score = get_distance_hybrid(query, ref)
        scores.append((allele_name, e_mi, e_ml, e_raw, e_gc, p_mi, p_ml, p_gc, p_score))

    if not scores:
        return None

    # ── edlib-based chains ──────────────────────────────────────────────

    # Production: edlib MI → ML → iter
    best_mi = max(s[1] for s in scores)
    mi_tied = [s for s in scores if s[1] == best_mi]
    best_ml = max(s[2] for s in mi_tied)
    ml_tied = [s for s in mi_tied if s[2] == best_ml]

    resolution = "MI" if len(mi_tied) == 1 else ("match_length" if len(ml_tied) == 1 else "iteration_order")
    prod_winner    = truncate(ml_tied[0][0], 4)

    # Chain 1: edlib MI → iter
    mi_iter_winner = truncate(mi_tied[0][0], 4)

    # Chain 2: edlib ML → iter
    ml_best        = max(s[2] for s in scores)
    ml_tied_all    = [s for s in scores if s[2] == ml_best]
    ml_iter_winner = truncate(ml_tied_all[0][0], 4)

    # Chain 3: edlib raw edit → iter (minimize)
    best_raw       = min(s[3] for s in scores)
    raw_iter_winner = truncate(next(s for s in scores if s[3] == best_raw)[0], 4)

    # Chain 4: edlib gc edit → iter (minimize)
    best_egc       = min(s[4] for s in scores)
    gc_iter_winner  = truncate(next(s for s in scores if s[4] == best_egc)[0], 4)

    # ── parasail-based chains ───────────────────────────────────────────

    # Chain 5: parasail MI → iter
    best_pmi       = max(s[5] for s in scores)
    pmi_tied       = [s for s in scores if s[5] == best_pmi]
    pmi_iter_winner = truncate(pmi_tied[0][0], 4)

    # Chain 6: parasail MI → parasail gc → iter
    best_pgc_in_pmi = min(s[7] for s in pmi_tied)
    pmi_pgc_tied    = [s for s in pmi_tied if s[7] == best_pgc_in_pmi]
    pmi_pgc_winner  = truncate(pmi_pgc_tied[0][0], 4)

    # Chain 7: parasail MI → parasail ML → iter  (direct analog of production)
    best_pml_in_pmi = max(s[6] for s in pmi_tied)
    pmi_pml_tied    = [s for s in pmi_tied if s[6] == best_pml_in_pmi]
    pmi_pml_winner  = truncate(pmi_pml_tied[0][0], 4)

    # Chain 8: parasail score → iter (alignment score, higher = better)
    best_pscore     = max(s[8] for s in scores)
    pscore_winner   = truncate(next(s for s in scores if s[8] == best_pscore)[0], 4)

    # Chain 9: edlib MI → parasail score → iter  (edlib primary, parasail tiebreaker)
    best_ps_in_mi   = max(s[8] for s in mi_tied)
    mi_pscore_tied  = [s for s in mi_tied if s[8] == best_ps_in_mi]
    mi_pscore_winner = truncate(mi_pscore_tied[0][0], 4)

    # Chain 10: edlib MI → edlib ML → parasail score → iter  (production + parasail final tiebreak)
    best_ps_in_ml   = max(s[8] for s in ml_tied)
    ml_pscore_tied  = [s for s in ml_tied if s[8] == best_ps_in_ml]
    prod_pscore_winner = truncate(ml_pscore_tied[0][0], 4)

    return {
        "resolution": resolution,
        "num_candidates": len(candidates),
        "num_mi_tied": len(mi_tied),
        "num_ml_tied": len(ml_tied),
        "prod_winner":       prod_winner,
        "mi_iter_winner":    mi_iter_winner,
        "ml_iter_winner":    ml_iter_winner,
        "raw_iter_winner":   raw_iter_winner,
        "gc_iter_winner":    gc_iter_winner,
        "pmi_iter_winner":   pmi_iter_winner,
        "pmi_pgc_winner":    pmi_pgc_winner,
        "pmi_pml_winner":    pmi_pml_winner,
        "pscore_winner":     pscore_winner,
        "mi_pscore_winner":  mi_pscore_winner,
        "prod_pscore_winner": prod_pscore_winner,
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
        mi_4f        = [None, None]
        ml_4f        = [None, None]
        raw_4f       = [None, None]
        gc_4f        = [None, None]
        pmi_4f       = [None, None]
        pmi_pgc_4f   = [None, None]
        pmi_pml_4f   = [None, None]
        pscore_4f    = [None, None]
        mi_pscore_4f = [None, None]
        prod_pscore_4f = [None, None]
        res_by_hap = [None, None]
        for h in [0, 1]:
            ri = alt_winners.get((sid, gene, h + 1))
            res_by_hap[h] = ri
            if ri:
                mi_4f[h]          = ri["mi_iter_winner"]
                ml_4f[h]          = ri["ml_iter_winner"]
                raw_4f[h]         = ri["raw_iter_winner"]
                gc_4f[h]          = ri["gc_iter_winner"]
                pmi_4f[h]         = ri["pmi_iter_winner"]
                pmi_pgc_4f[h]     = ri["pmi_pgc_winner"]
                pmi_pml_4f[h]     = ri["pmi_pml_winner"]
                pscore_4f[h]      = ri["pscore_winner"]
                mi_pscore_4f[h]   = ri["mi_pscore_winner"]
                prod_pscore_4f[h] = ri["prod_pscore_winner"]

        # Skip if no production test alleles have 4 fields
        prod_4f_valid = [p for p in prod_4f if p]
        if not prod_4f_valid:
            continue

        # Build per-chain test-allele lists and pre-compute for strict matching
        chain_test_4f = {
            "correct":              prod_4f,
            "mi_iter_correct":      mi_4f,
            "ml_iter_correct":      ml_4f,
            "raw_iter_correct":     raw_4f,
            "gc_iter_correct":      gc_4f,
            "pmi_iter_correct":     pmi_4f,
            "pmi_pgc_correct":      pmi_pgc_4f,
            "pmi_pml_correct":      pmi_pml_4f,
            "pscore_correct":       pscore_4f,
            "mi_pscore_correct":    mi_pscore_4f,
            "prod_pscore_correct":  prod_pscore_4f,
        }
        chain_used = {k: set() for k in chain_test_4f}

        # Pre-compute truth first-option 4f for homozygous check
        all_truth_first_4f = []
        for ta_opts in truth_alleles:
            if ta_opts:
                all_truth_first_4f.append(truncate(ta_opts[0], 4))
            else:
                all_truth_first_4f.append(None)

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

            # Strict cross-match each chain with proper used-test tracking
            chain_matched = {}
            for chain_key, chain_4f in chain_test_4f.items():
                chain_4f_valid = [p for p in chain_4f if p]
                chain_matched[chain_key] = strict_cross_match(
                    truth_4f_list, chain_4f, chain_4f_valid,
                    all_truth_first_4f, chain_used[chain_key])

            prod_matched = chain_matched["correct"]

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

            _common = {k: chain_matched[k] for k in chain_matched if k != "correct"}

            if res_info is None:
                skipped.append(f"{sid} {gene_short}_{actual_hap+1}")
                level_results.append({
                    "sample": sid, "gene": gene, "hap": actual_hap + 1,
                    "truth_4f": truth_4f_list, "truth_4th_field": truth_4th,
                    "is_01": truth_4th == 1, "correct": prod_matched,
                    "resolution": "unknown",
                    "num_mi_tied": 0, "num_ml_tied": 0,
                    **_common,
                })
                continue

            level_results.append({
                "sample": sid, "gene": gene, "hap": actual_hap + 1,
                "truth_4f": truth_4f_list, "truth_4th_field": truth_4th,
                "is_01": truth_4th == 1, "correct": prod_matched,
                "resolution": res_info["resolution"],
                "num_mi_tied": res_info["num_mi_tied"],
                "num_ml_tied": res_info["num_ml_tied"],
                **_common,
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

# ── Part C: Compare metric chains (with proper cross-matching) ─────────
print(f"\n{'=' * 70}")
print(f"  PART C: Primary metric comparison  (all chains, cross-matching applied)")
print(f"{'=' * 70}")

total_c = len(level_results)

# Collect counts for all chains
chains = [
    ("Production (eMI→eML→iter)",   "correct"),
    ("eMI → iter",                  "mi_iter_correct"),
    ("eML → iter",                  "ml_iter_correct"),
    ("e_raw → iter",                "raw_iter_correct"),
    ("e_gc → iter",                 "gc_iter_correct"),
    ("pMI → iter",                  "pmi_iter_correct"),
    ("pMI → pGC → iter",            "pmi_pgc_correct"),
    ("pMI → pML → iter",            "pmi_pml_correct"),
    ("p_score → iter",              "pscore_correct"),
    ("eMI → p_score → iter",        "mi_pscore_correct"),
    ("eMI→eML → p_score → iter",    "prod_pscore_correct"),
]

def count(key, rows=None):
    rows = rows or level_results
    return sum(1 for r in rows if r[key])

print(f"\n  Overall ({total_c} haplotypes):")
print(f"  {'Chain':<30} {'Correct':>10}  {'%':>6}  {'vs Prod':>8}")
print(f"  {'-'*30} {'-'*10}  {'-'*6}  {'-'*8}")
prod_n = count("correct")
for label, key in chains:
    n = count(key)
    diff = n - prod_n
    marker = " ◄" if n > prod_n else ("" if n == prod_n else "")
    print(f"  {label:<30} {n:>5}/{total_c:<5}  {n/total_c*100:>5.1f}%  {diff:>+7d}{marker}")

# By truth type
print(f"\n  By truth 4th-field type:")
w = 16
header = f"  {'Chain':<30}"
for lbl in ("Truth=:01", "Truth≠:01"):
    header += f" {lbl:>{w}}"
print(header)
print(f"  {'-'*30}" + f" {'-'*w}" * 2)
for label, key in chains:
    row = f"  {label:<30}"
    for filt_key, filt in [("is_01", True), ("is_01", False)]:
        sub = [r for r in level_results if r["is_01"] == filt]
        n   = len(sub)
        c   = count(key, sub)
        row += f" {c}/{n} ({c/n*100:.1f}%)".rjust(w)
    print(row)

# Per-gene
print(f"\n  Per-gene comparison:")
chain_keys = [k for _, k in chains]
chain_labels_short = ["Prod","eMI","eML","eRaw","eGC","pMI","pMI+pGC","pMI+pML","pScr","eMI+pScr","Prod+pScr"]
gene_header = f"  {'Gene':<12}" + "".join(f" {l:>11}" for l in chain_labels_short)
print(gene_header)
print(f"  {'-'*12}" + f" {'-'*11}" * len(chain_labels_short))
for gene in GENES:
    gr = [r for r in level_results if r["gene"] == gene]
    if not gr:
        continue
    n = len(gr)
    row = f"  {gene:<12}"
    for key in chain_keys:
        c = count(key, gr)
        row += f" {c}/{n}({c/n*100:.0f}%)".rjust(11)
    print(row)

print(f"\n{'=' * 70}")
print(f"  DONE")
print(f"{'=' * 70}")
