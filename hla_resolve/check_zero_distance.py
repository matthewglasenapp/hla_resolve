"""
check_zero_distance.py — Count how many pass 2 assignments have zero edit distance.

Usage:
    python check_zero_distance.py \
        --xml /Users/matt/Desktop/IPD_IMGT_XML/hla.xml \
        --samples /Users/matt/Desktop/HLA_haplotypes.fa
"""

import sys, os, argparse
sys.path.insert(0, os.path.dirname(os.path.abspath(__file__)))

from hla_typer import (
    build_g_group_dict, get_g_group_exons, load_test_data,
    pass_1_classification, produce_allele_seq_db, get_distance,
    get_gene, get_sampleid, generate_allele_dict,
)


def pass_2_with_distances(sequence_data, allele_to_g_groups, results_dict, samples):
    """Run pass 2 with raw edit distance, returning per-sample distances."""
    g_group_to_allele = generate_allele_dict(allele_to_g_groups)
    all_allele_sequence_db = produce_allele_seq_db(sequence_data, exon_only=True)

    results = {}
    for sample_name, classification in results_dict.items():
        g_group = classification[0]

        # Perfect G-group match — distance is 0 by definition
        if g_group is not None and ';' in g_group and \
           classification[1] == 0 and classification[2] == 1:
            results[sample_name] = (g_group.split(";")[1], 0)
            continue

        search_all = g_group not in g_group_to_allele
        if not search_all:
            alleles = g_group_to_allele[g_group]
            allele_sequence_db = produce_allele_seq_db(sequence_data, selected_alleles=alleles, exon_only=True)
        else:
            allele_sequence_db = all_allele_sequence_db

        hla_gene = get_gene(sample_name)
        best_allele = None
        best_dist = None

        for class_name, seq_data in allele_sequence_db.items():
            if hla_gene not in class_name:
                continue
            raw_dist, match_len, mismatches = get_distance(
                samples[sample_name], seq_data, get_length=True, gap_compressed=False
            )
            if best_dist is None or raw_dist < best_dist:
                best_allele = class_name
                best_dist = raw_dist

        results[sample_name] = (best_allele, best_dist)

    return results


def main():
    parser = argparse.ArgumentParser()
    parser.add_argument('--xml', required=True)
    parser.add_argument('--samples', required=True)
    parser.add_argument('--ignore-incomplete', action='store_true', default=True)
    args = parser.parse_args()

    print("Loading XML reference...")
    g_group_dict, _, sequence_data = build_g_group_dict(args.xml, ignore_incomplete=args.ignore_incomplete)
    print(f"{len(sequence_data)} alleles loaded")

    print("Extracting G-group exons...")
    g_group_common_seqs = get_g_group_exons(g_group_dict, sequence_data)

    print("Loading samples...")
    samples = load_test_data(args.samples)
    # Exclude problematic samples
    samples = {k: v for k, v in samples.items() if "IHW09117" not in k}
    print(f"{len(samples)} samples loaded (excluding IHW09117)")

    print("\nRunning pass 1...")
    g_group_results = pass_1_classification(g_group_common_seqs, samples, g_group_dict)

    print("Running pass 2 (raw edit distance)...")
    results = pass_2_with_distances(sequence_data, g_group_dict, g_group_results, samples)

    # Count and report
    zero = []
    nonzero = []
    for sample_name, (allele, dist) in sorted(results.items()):
        gene = get_gene(sample_name)
        sample_id = get_sampleid(sample_name)
        if dist == 0:
            zero.append((sample_id, gene, sample_name, allele))
        else:
            nonzero.append((sample_id, gene, sample_name, allele, dist))

    print(f"\n{'='*60}")
    print(f"  ZERO EDIT DISTANCE: {len(zero)}/{len(results)}")
    print(f"  NON-ZERO:           {len(nonzero)}/{len(results)}")
    print(f"{'='*60}")

    # Per-gene summary
    from collections import Counter
    gene_zero = Counter(g for _, g, _, _ in zero)
    gene_nonzero = Counter(g for _, g, _, _, _ in nonzero)
    all_genes = sorted(set(list(gene_zero.keys()) + list(gene_nonzero.keys())))

    print(f"\n{'Gene':<12} {'Zero':>6} {'Non-zero':>10} {'Total':>7} {'% Zero':>8}")
    print("-" * 47)
    for gene in all_genes:
        z = gene_zero[gene]
        nz = gene_nonzero[gene]
        t = z + nz
        pct = z / t * 100 if t > 0 else 0
        print(f"{gene:<12} {z:>6} {nz:>10} {t:>7} {pct:>7.1f}%")

    total_z = len(zero)
    total_nz = len(nonzero)
    total = total_z + total_nz
    print("-" * 47)
    print(f"{'TOTAL':<12} {total_z:>6} {total_nz:>10} {total:>7} {total_z/total*100:>7.1f}%")

    # List non-zero cases
    if nonzero:
        print(f"\nNon-zero edit distance assignments:")
        print(f"{'Sample':<12} {'Gene':<12} {'Allele':<30} {'Dist':>5}")
        print("-" * 62)
        for sample_id, gene, sample_name, allele, dist in sorted(nonzero):
            print(f"{sample_id:<12} {gene:<12} {allele:<30} {dist:>5}")


if __name__ == "__main__":
    main()
