"""
Detect ties in pass 2 and evaluate whether a match-length tiebreaker
would improve 3-field concordance vs the current alphabetical tiebreaker.

Excludes IHW09117 (contaminated sample).
Cross-references against truth data.
"""
import sys, os, re
import pandas as pd
sys.path.insert(0, os.path.dirname(os.path.abspath(__file__)))

from hla_typer_claude import (
    build_g_group_dict,
    get_g_group_exons,
    load_test_data,
    pass_1_classification,
    produce_allele_seq_db,
    assign_classification_to_sample_full_seq,
    generate_allele_dict,
    get_gene,
    get_sampleid,
    get_distance,
)

XML = "/Users/matt/Desktop/IPD_IMGT_XML/hla.xml"
SAMPLES = "/Users/matt/Desktop/HLA_haplotypes.fa"
TRUTH = "/Users/matt/Desktop/ihw_truth_full_edited.csv"
EXCLUDE = {"IHW09117"}

# ── Helpers ──────────────────────────────────────────────────────────

def truncate(allele, n_fields):
    """Truncate an allele name to n fields, with zero-padded first field."""
    if not isinstance(allele, str) or "*" not in allele:
        return None
    head, tail = allele.split("*", 1)
    fields = tail.split(":")
    if len(fields) < n_fields:
        return None
    padded = []
    for f in fields[:n_fields]:
        m = re.match(r"(\d+)(.*)", f)
        if m:
            padded.append(f"{int(m.group(1)):02d}{m.group(2)}")
        else:
            padded.append(f)
    return f"{head}*{':'.join(padded)}"

def normalize_allele(a):
    """Strip HLA- prefix and expression suffixes for comparison."""
    if not isinstance(a, str):
        return None
    a = a.replace("HLA-", "")
    # Strip trailing expression characters (N, L, S, C, A, Q)
    a = re.sub(r"[NLSCAQ]$", "", a)
    return a

def truth_matches_at_n(called_allele, truth_options, n):
    """Check if called_allele matches any truth option at n-field level."""
    called_trunc = truncate(normalize_allele(called_allele), n)
    if called_trunc is None:
        return None  # can't evaluate
    for opt in truth_options:
        truth_trunc = truncate(normalize_allele(opt), n)
        if truth_trunc is not None and called_trunc == truth_trunc:
            return True
    return False

# ── Load data ────────────────────────────────────────────────────────

print("Loading XML reference...")
g_group_dict, _, sequence_data = build_g_group_dict(XML, ignore_incomplete=True)
print(f"  {len(sequence_data)} alleles loaded")

print("Extracting G-group common sequences...")
g_group_common_seqs = get_g_group_exons(g_group_dict, sequence_data)

print("Loading samples...")
samples = load_test_data(SAMPLES)
print(f"  {len(samples)} samples loaded")

# Filter out excluded samples
samples = {k: v for k, v in samples.items() if get_sampleid(k) not in EXCLUDE}
print(f"  {len(samples)} samples after excluding {EXCLUDE}")

print("Loading truth data...")
truth_df = pd.read_csv(TRUTH)
# Build truth lookup: {sample_id: {gene_col: [options]}}
truth_lookup = {}
for _, row in truth_df.iterrows():
    sid = row["sample"]
    if sid in EXCLUDE:
        continue
    truth_lookup[sid] = {}
    for col in truth_df.columns:
        if col in ("sample", "source"):
            continue
        val = row[col]
        if pd.isna(val):
            continue
        # Split on | or / for multiple truth options
        options = [x.strip() for x in re.split(r"[|/]", str(val))]
        truth_lookup[sid][col] = options

# ── Run pass 1 ───────────────────────────────────────────────────────

print("\nRunning pass 1...")
g_group_results = pass_1_classification(g_group_common_seqs, samples, g_group_dict)

failed_pass1 = {n: r for n, r in g_group_results.items() if r[0] is None}
perfect_pass1 = {n: r for n, r in g_group_results.items() if r[0] is not None}
print(f"Pass 1: {len(perfect_pass1)} perfect, {len(failed_pass1)} failed")

# ── Pass 2 with detailed tie analysis ────────────────────────────────

g_group_to_allele = generate_allele_dict(g_group_dict)
all_allele_seq_db = produce_allele_seq_db(sequence_data, exon_only=True)
METRIC = "edit_distance"

def get_truth_for_sample(sample_name):
    """Get truth options for a sample, trying both allele indices."""
    sid = get_sampleid(sample_name)
    if sid not in truth_lookup:
        return None

    # Extract gene column from sample name: e.g. "IHW09049_HLA-DQA1_1" -> "HLA-DQA1_1"
    hla_idx = sample_name.find("HLA")
    if hla_idx < 0:
        return None
    gene_col = sample_name[hla_idx:]

    return truth_lookup[sid].get(gene_col)

def run_pass2_for_sample(sample_name, classification):
    """Run pass 2 for a single sample, return (result, search_type, allele_seq_db)."""
    g_group = classification[0]

    if g_group is not None and ';' in g_group and classification[1] == 0 and classification[2] == 1:
        return None, "shortcircuit", None

    search_all = g_group not in g_group_to_allele
    if not search_all:
        alleles = g_group_to_allele[g_group]
        db = produce_allele_seq_db(sequence_data, selected_alleles=alleles, exon_only=True)
        search_type = "g_group_restricted"
    else:
        db = all_allele_seq_db
        search_type = "full_search"

    result = assign_classification_to_sample_full_seq(
        db, samples[sample_name], sample_name, eval_metric=METRIC
    )
    return result, search_type, db

def pick_by_match_length(equidistant, db, sample_seq, sample_name):
    """Among tied alleles, pick the one with longest match length."""
    best_allele = None
    best_match_len = -1
    for allele in equidistant:
        if allele not in db:
            continue
        _, match_len = get_distance(sample_seq, db[allele], get_length=True, gap_compressed=False)
        if match_len > best_match_len:
            best_match_len = match_len
            best_allele = allele
    return best_allele, best_match_len

# ── Analyze all samples ──────────────────────────────────────────────

print(f"\n{'='*70}")
print(f"  ANALYZING TIES: alphabetical vs match-length tiebreaker")
print(f"{'='*70}\n")

results = []
current = 0
total = len(g_group_results)

for sample_name, classification in g_group_results.items():
    current += 1
    if current % 50 == 0:
        print(f"  Processing {current}/{total}...")

    result, search_type, db = run_pass2_for_sample(sample_name, classification)
    if result is None:
        continue

    equidistant = result[-1]
    truth_options = get_truth_for_sample(sample_name)

    # Current winner: alphabetical
    alpha_winner = sorted(equidistant)[0] if len(equidistant) > 1 else result[0]

    # Match-length winner
    if len(equidistant) > 1 and db is not None:
        ml_winner, ml_len = pick_by_match_length(equidistant, db, samples[sample_name], sample_name)
    else:
        ml_winner = result[0]
        ml_len = result[2]

    # Check how many unique 3-field alleles are in the tie set
    unique_3field = set()
    for a in equidistant:
        t = truncate(normalize_allele(a), 3)
        if t:
            unique_3field.add(t)

    results.append({
        "sample": sample_name,
        "search_type": search_type,
        "edit_distance": result[1],
        "num_tied": len(equidistant),
        "num_unique_3field": len(unique_3field),
        "unique_3field_alleles": sorted(unique_3field),
        "alpha_winner": alpha_winner,
        "ml_winner": ml_winner,
        "ml_match_len": ml_len,
        "truth_options": truth_options,
        "has_truth": truth_options is not None,
    })

# ── Report: ties crossing 3-field boundaries ─────────────────────────

cross_3field = [r for r in results if r["num_unique_3field"] > 1]

print(f"\n{'='*70}")
print(f"  TIES CROSSING 3-FIELD BOUNDARIES (excluding IHW09117)")
print(f"{'='*70}")
print(f"  Total: {len(cross_3field)} samples\n")

for r in cross_3field:
    print(f"  {r['sample']} ({r['search_type']})")
    print(f"    Edit distance: {r['edit_distance']}, Tied: {r['num_tied']}, Unique 3-field groups: {r['num_unique_3field']}")
    print(f"    3-field groups: {r['unique_3field_alleles']}")
    print(f"    Alpha winner:  {r['alpha_winner']}")
    print(f"    ML winner:     {r['ml_winner']} (match_len={r['ml_match_len']})")
    if r["truth_options"]:
        # Check both at 3-field
        alpha_correct = truth_matches_at_n(r["alpha_winner"], r["truth_options"], 3)
        ml_correct = truth_matches_at_n(r["ml_winner"], r["truth_options"], 3)
        print(f"    Truth: {r['truth_options']}")
        print(f"    Alpha correct at 3-field: {alpha_correct}")
        print(f"    ML correct at 3-field:    {ml_correct}")
    else:
        print(f"    (No truth data for this gene)")
    print()

# ── Report: 3-field concordance comparison ────────────────────────────

print(f"\n{'='*70}")
print(f"  3-FIELD CONCORDANCE: alphabetical vs match-length tiebreaker")
print(f"  (only samples with truth data AND ties)")
print(f"{'='*70}")

tied_with_truth = [r for r in results if r["num_tied"] > 1 and r["has_truth"]]
alpha_correct_3 = 0
ml_correct_3 = 0
evaluable = 0
disagree_cases = []

for r in tied_with_truth:
    ac = truth_matches_at_n(r["alpha_winner"], r["truth_options"], 3)
    mc = truth_matches_at_n(r["ml_winner"], r["truth_options"], 3)
    if ac is None or mc is None:
        continue
    evaluable += 1
    if ac:
        alpha_correct_3 += 1
    if mc:
        ml_correct_3 += 1
    if ac != mc:
        disagree_cases.append((r, ac, mc))

print(f"\n  Evaluable tied samples with truth: {evaluable}")
print(f"  Alphabetical correct at 3-field:  {alpha_correct_3} ({alpha_correct_3/evaluable*100:.1f}%)" if evaluable else "")
print(f"  Match-length correct at 3-field:  {ml_correct_3} ({ml_correct_3/evaluable*100:.1f}%)" if evaluable else "")

if disagree_cases:
    print(f"\n  Cases where tiebreakers DISAGREE ({len(disagree_cases)}):")
    for r, ac, mc in disagree_cases:
        print(f"\n    {r['sample']} (edit_dist={r['edit_distance']}, {r['num_tied']} tied)")
        print(f"      Truth:         {r['truth_options']}")
        print(f"      Alpha: {r['alpha_winner']} -> {'CORRECT' if ac else 'WRONG'}")
        print(f"      ML:    {r['ml_winner']} (len={r['ml_match_len']}) -> {'CORRECT' if mc else 'WRONG'}")
else:
    print(f"\n  No disagreements between tiebreakers.")

# ── Overall summary ──────────────────────────────────────────────────

all_tied = [r for r in results if r["num_tied"] > 1]
full_search_tied = [r for r in all_tied if r["search_type"] == "full_search"]
group_tied = [r for r in all_tied if r["search_type"] == "g_group_restricted"]

print(f"\n{'='*70}")
print(f"  SUMMARY")
print(f"{'='*70}")
print(f"  Total samples (excl IHW09117):     {len(results)}")
print(f"  Samples with ties:                 {len(all_tied)}")
print(f"    - in full search:                {len(full_search_tied)}")
print(f"    - in G-group restricted:         {len(group_tied)}")
print(f"  Ties crossing 3-field boundary:    {len(cross_3field)}")
print(f"  Tiebreaker disagreements vs truth: {len(disagree_cases)}")
