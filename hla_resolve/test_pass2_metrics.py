"""
test_pass2_metrics.py — Compare 5 different pass 2 metrics for 3-field concordance.

Runs pass 1 once, then pass 2 five times (one per metric), scoring each
at all field levels using the concordance logic from alex_test_against_truth.py.

Usage:
    python test_pass2_metrics.py \
        --xml /Users/matt/Desktop/IPD_IMGT_XML/hla.xml \
        --samples /Users/matt/Desktop/HLA_haplotypes.fa \
        --truth /Users/matt/Desktop/ihw_truth_full_edited.csv
"""

import sys
import os
import re
import argparse
import pandas as pd

sys.path.insert(0, os.path.dirname(os.path.abspath(__file__)))

# Import from hla_typer_claude (bug-fixed get_distance)
from hla_typer_claude import (
    build_g_group_dict,
    get_g_group_exons,
    load_test_data,
    pass_1_classification,
    pass_2_classification,
    produce_allele_seq_db,
    get_distance,
    get_gene,
    get_sampleid,
    generate_allele_dict,
)


# =============================================================================
# Concordance logic (from alex_test_against_truth.py)
# =============================================================================

PROBLEMATIC_SAMPLES = ["IHW09117"]
PROBLEMATIC_GENES = []

NO_CALL_PAT = re.compile(r"(?:^|\b)(?:no\s*call|na|none)$", re.IGNORECASE)


def is_called(x):
    return isinstance(x, str) and "*" in x and not NO_CALL_PAT.search(x)


def split_truth_options(a):
    if not isinstance(a, str):
        return []
    parts = re.split(r"[|/]", a)
    return [p.strip() for p in parts]


def truncate_to_fields(allele, n):
    if not is_called(allele):
        return None
    head, tail = allele.split("*", 1)
    fields = tail.split(":")
    if len(fields) < n:
        return None
    padded = []
    for f in fields[:n]:
        match = re.match(r"(\d+)(.*)", f)
        if match:
            padded.append(f"{int(match.group(1)):02d}{match.group(2)}")
        else:
            padded.append(f)
    return f"{head}*{':'.join(padded)}"


def get_percentages(truth_df, class_df, genes, genes_numbered, source="IHW"):
    """Score concordance at all field levels (1-4). Returns (concordance, fields) dicts."""
    # Filter by tool/platform if columns exist
    if "tool" in class_df.columns:
        class_df = class_df[class_df["tool"] == "hla_resolve"]
    if "platform" in class_df.columns:
        class_df = class_df[class_df["platform"] == "pacbio"]

    # Fill NaN with "no call"
    for col in genes:
        for suffix in ["_1", "_2"]:
            colname = col + suffix
            if colname in truth_df.columns:
                null_entries = truth_df[colname].isna()
                truth_df.loc[null_entries, colname] = "no call"
            if colname in class_df.columns:
                null_entries = class_df[colname].isna()
                class_df.loc[null_entries, colname] = "no call"

    # Remove HLA- prefix from classifier allele values if present
    for col in genes_numbered:
        if col not in class_df.columns:
            continue
        letter = col.split("_")[0]
        mask = class_df[col].str.contains(f"HLA-{letter}\\*", regex=True, na=True)
        class_df.loc[~mask, col] = class_df.loc[~mask, col].astype(str).str.removeprefix('HLA-')

    # Common samples
    common_samples = set(truth_df["sample"]) & set(class_df["sample"])
    for ps in set(PROBLEMATIC_SAMPLES) & common_samples:
        common_samples.remove(ps)

    common_truth = truth_df[truth_df["sample"].isin(common_samples)]
    common_class = class_df[class_df["sample"].isin(common_samples)]

    # Select truth source if column exists
    if "source" in common_truth.columns:
        src_common_truth = common_truth[common_truth["source"] == source]
    else:
        src_common_truth = common_truth

    concordance = {f"{n}-field": {g: 0 for g in genes} for n in range(1, 5)}
    fields = {f"{n}-field": {g: 0 for g in genes} for n in range(1, 5)}

    for sample in common_samples:
        truth = src_common_truth[src_common_truth["sample"] == sample].iloc[0]
        classification = common_class[common_class["sample"] == sample].iloc[0]

        for gene in genes:
            gene_name = gene.split("HLA-")[1] if "HLA-" in gene else gene
            truth_alleles = [
                truth.get(f"{gene}_{i}", "NA") for i in [1, 2]
                if f"{sample} {gene_name}_{i}" not in PROBLEMATIC_GENES
            ]
            test_alleles = [
                classification.get(f"{gene}_{i}", "NA") for i in [1, 2]
                if f"{sample} {gene_name}_{i}" not in PROBLEMATIC_GENES
            ]
            truth_alleles = [a for a in truth_alleles if is_called(a)]
            test_alleles = [a for a in test_alleles if is_called(a)]
            if not truth_alleles or not test_alleles:
                continue

            for n in range(1, 5):
                used_tests = set()

                for truth_allele in truth_alleles:
                    if not is_called(truth_allele):
                        continue
                    if len(truth_allele.split("*", 1)[1].split(":")) < n:
                        continue

                    truth_opts = split_truth_options(truth_allele) or [truth_allele]
                    truth_truncs = [truncate_to_fields(a, n) for a in truth_opts if truncate_to_fields(a, n)]
                    if not truth_truncs:
                        continue

                    test_truncs = [truncate_to_fields(a, n) for a in test_alleles if truncate_to_fields(a, n)]

                    fields[f"{n}-field"][gene] += 1

                    matched = False
                    for i, test_allele in enumerate(test_alleles):
                        test_trunc = truncate_to_fields(test_allele, n)
                        if not test_trunc:
                            continue

                        truth_options = [truncate_to_fields(x, n) for x in truth_alleles]
                        if truth_options.count(test_trunc) == 2 and test_truncs.count(test_trunc) == 2:
                            if test_trunc in truth_truncs:
                                matched = True
                                break
                        elif i not in used_tests and test_trunc in truth_truncs:
                            matched = True
                            used_tests.add(i)
                            break

                    if matched:
                        concordance[f"{n}-field"][gene] += 1

    return concordance, fields


# =============================================================================
# Raw edit distance classifier (not natively supported)
# =============================================================================

def assign_classification_raw_edit(full_sequence, sequence, full_sample_name):
    """Classify using raw (non-gap-compressed) edit distance, minimizing."""
    best = (None, None, None, None, None)
    same_dist = []
    hla_gene = get_gene(full_sample_name)

    for class_name, seq_data in full_sequence.items():
        if hla_gene not in class_name:
            continue

        raw_dist, match_len, mismatches = get_distance(
            sequence, seq_data, get_length=True, gap_compressed=False, get_mismatches=True
        )

        if best[1] is None or raw_dist < best[1]:
            best = (class_name, raw_dist, match_len, 0, 0)
            same_dist = [class_name]
        elif raw_dist == best[1]:
            same_dist.append(class_name)

    return (*best, same_dist)


def pass_2_raw_edit(sequence_data, allele_to_g_groups, results_dict, samples):
    """Run pass 2 using raw (non-gap-compressed) edit distance."""
    g_group_to_allele = generate_allele_dict(allele_to_g_groups)
    all_allele_sequence_db = produce_allele_seq_db(sequence_data, exon_only=True)

    results = {}
    current = 0
    for sample_name, classification in results_dict.items():
        current += 1
        if current % 20 == 0:
            print(f"INFO: Processing {current}/{len(results_dict)} ({(current/len(results_dict))*100:.2f}%)")

        g_group = classification[0]

        if g_group is not None and ';' in g_group and \
           classification[1] == 0 and classification[2] == 1:
            results[sample_name] = (g_group.split(";")[1], 0, 0, 0, 0, [g_group.split(";")[1]])
            continue

        search_all = g_group not in g_group_to_allele
        if not search_all:
            alleles = g_group_to_allele[g_group]
            allele_sequence_db = produce_allele_seq_db(sequence_data, selected_alleles=alleles, exon_only=True)
        else:
            allele_sequence_db = all_allele_sequence_db

        result = assign_classification_raw_edit(allele_sequence_db, samples[sample_name], sample_name)

        # Equidistant tiebreaker: pick alphabetically lowest
        equidistant = result[-1]
        if len(equidistant) > 1:
            equidistant = sorted(equidistant)
            result = (equidistant[0], *result[1:])

        results[sample_name] = result

    return results


# =============================================================================
# Convert results dict to DataFrame for concordance scoring
# =============================================================================

def results_to_dataframe(results):
    """
    Convert pass 2 results dict {sample_name: (allele, ...)} to a DataFrame
    with columns [sample, HLA-A_1, HLA-A_2, ...] matching truth format.
    """
    data = {}
    for sample_name, result in results.items():
        allele_name = result[0] if isinstance(result, tuple) else result
        if allele_name is None:
            continue

        sample_id = get_sampleid(sample_name)
        # Extract gene column: "HLA-DQB1_1" from "HG002_HLA-DQB1_1"
        hla_idx = sample_name.find("HLA")
        gene_col = sample_name[hla_idx:]

        if sample_id not in data:
            data[sample_id] = {}
        data[sample_id][gene_col] = allele_name

    rows = []
    for sample_id in sorted(data.keys()):
        row = {"sample": sample_id}
        row.update(data[sample_id])
        rows.append(row)

    return pd.DataFrame(rows)


# =============================================================================
# Main
# =============================================================================

METRICS = [
    ("match_length",             "match_length", False),
    ("raw_edit_distance",        None,           True),
    ("gap_compressed_edit_dist", "edit_distance", False),
    ("sequence_identity",        "identity",     False),
    ("mismatch_identity",        "mismatch",     False),
]


def main():
    parser = argparse.ArgumentParser(description='Compare pass 2 metrics for 3-field concordance')
    parser.add_argument('--xml', required=True, help='IMGT XML reference file')
    parser.add_argument('--samples', required=True, help='CDS FASTA (concatenated exons)')
    parser.add_argument('--truth', required=True, help='Truth CSV file')
    args = parser.parse_args()

    # ── Load shared data (done once) ──
    print("INFO: Loading XML reference...")
    g_group_dict, p_group_dict, sequence_data = build_g_group_dict(args.xml)
    print(f"INFO: {len(sequence_data)} alleles loaded")

    print("INFO: Extracting G-group common sequences...")
    g_group_common_seqs = get_g_group_exons(g_group_dict, sequence_data)
    if g_group_common_seqs is None:
        print("ERR: G-group exon extraction failed")
        sys.exit(1)

    print("INFO: Loading samples...")
    samples = load_test_data(args.samples)
    print(f"INFO: {len(samples)} samples loaded")

    print("INFO: Loading truth data...")
    truth_df = pd.read_csv(args.truth)

    # Derive gene lists from truth file
    genes_numbered = [c for c in truth_df.columns if c != "sample" and c != "source"]
    genes = sorted(list(set(c.rsplit("_", 1)[0] for c in genes_numbered)))

    # ── Run pass 1 once ──
    print("\nINFO: Running pass 1 (G-group assignment)...")
    g_group_results = pass_1_classification(g_group_common_seqs, samples, g_group_dict)

    # ── Run pass 2 for each metric and score ──
    all_concordances = {}

    for label, metric_str, is_custom in METRICS:
        print(f"\n{'='*60}")
        print(f"  Running pass 2 with metric: {label}")
        print(f"{'='*60}")

        if is_custom:
            results = pass_2_raw_edit(sequence_data, g_group_dict, g_group_results, samples)
        else:
            results = pass_2_classification(
                sequence_data, g_group_dict, g_group_results, samples,
                metric=metric_str
            )

        # Convert to DataFrame
        class_df = results_to_dataframe(results)

        # Score concordance
        concordance, field_counts = get_percentages(
            truth_df.copy(), class_df, genes, genes_numbered
        )
        all_concordances[label] = (concordance, field_counts)

    # ── Print results ──
    print(f"\n\n{'='*70}")
    print(f"  RESULTS: 3-FIELD CONCORDANCE COMPARISON")
    print(f"{'='*70}")

    # Summary table
    print(f"\n{'Metric':<30} {'Correct':>8} {'Total':>6} {'Pct':>7}")
    print("-" * 55)
    for label, (concordance, field_counts) in all_concordances.items():
        correct = sum(concordance["3-field"].values())
        total = sum(field_counts["3-field"].values())
        pct = correct / total * 100 if total > 0 else 0
        print(f"{label:<30} {correct:>8} {total:>6} {pct:>6.1f}%")

    # Per-gene breakdown for each metric
    for label, (concordance, field_counts) in all_concordances.items():
        correct = sum(concordance["3-field"].values())
        total = sum(field_counts["3-field"].values())
        pct = correct / total * 100 if total > 0 else 0
        print(f"\n{'='*60}")
        print(f"  {label}")
        print(f"  3-field concordance: {correct}/{total} = {pct:.1f}%")
        print(f"{'='*60}")
        for gene in sorted(field_counts["3-field"].keys()):
            c = concordance["3-field"][gene]
            t = field_counts["3-field"][gene]
            p = c / t * 100 if t > 0 else 0
            print(f"  {gene:<12}: {c}/{t} = {p:.1f}%")

    # Also show all field levels in a compact table
    print(f"\n\n{'='*70}")
    print(f"  ALL FIELD LEVELS")
    print(f"{'='*70}")
    for n in range(1, 5):
        field_key = f"{n}-field"
        print(f"\n  {field_key}:")
        print(f"  {'Metric':<30} {'Correct':>8} {'Total':>6} {'Pct':>7}")
        print(f"  {'-'*53}")
        for label, (concordance, field_counts) in all_concordances.items():
            correct = sum(concordance[field_key].values())
            total = sum(field_counts[field_key].values())
            pct = correct / total * 100 if total > 0 else 0
            print(f"  {label:<30} {correct:>8} {total:>6} {pct:>6.1f}%")


if __name__ == "__main__":
    main()
