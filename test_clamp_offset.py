#!/usr/bin/env python3
"""
Test script to validate indel-aware clamping logic for DPA1 (minus strand).

Compares old method (reference-based indices) vs new method (indel-aware)
for subsetting the vcf2fasta gene output.

Usage:
    python test_clamp_offset.py
"""

import pysam
import glob

# ── Paths (CLQC1901_19QC-1 DPA1) ──
VCF_PATH = "/hb/groups/cornejo_lab/matt/ihiws_2026/data/hla_resolve_results/CLQC1901_19QC-1/filtered_vcf/CLQC1901_19QC-1_HLA-DPA1_PASS_phased.vcf.gz"
RAW_FASTA_DIR = "/hb/groups/cornejo_lab/matt/ihiws_2026/data/hla_resolve_results/CLQC1901_19QC-1/vcf2fasta_out/HLA-DPA1_gene/"
REFERENCE = "/hb/groups/cornejo_lab/matt/hla_resolve/hla_resolve/data/reference/GCA_000001405.15_GRCh38_no_alt_analysis_set.fna"
GENE_COORDS_FILE = "/hb/groups/cornejo_lab/matt/hla_resolve/hla_resolve/data/hla_gff/hla_dpa1_gene_coords.txt"

# ── DPA1 parameters from the pipeline output ──
GENE_START = 33064569   # lowest genomic coordinate (GFF start)
GENE_STOP = 33080775    # highest genomic coordinate (GFF stop)
CLAMPED_START = 33065644  # low end of clamped interval
CLAMPED_STOP = 33080775   # high end of clamped interval


def read_fasta_alleles(fasta_dir):
    """Read both haplotype sequences from vcf2fasta gene output."""
    fas_files = glob.glob(fasta_dir + "/*.fas")
    if not fas_files:
        raise FileNotFoundError(f"No .fas files in {fasta_dir}")
    alleles = []
    with open(fas_files[0]) as f:
        content = f.read().split(">")
    for entry in content[1:]:
        seq = "".join(entry.split("\n")[1:]).replace("-", "").strip()
        alleles.append(seq)
    return alleles[0], alleles[1]


def revcomp(seq):
    """Reverse complement a DNA sequence."""
    comp = str.maketrans("ACGTacgt", "TGCAtgca")
    return seq.translate(comp)[::-1]


def compute_offset_at_position(vcf_path, haplotype, target_pos):
    """
    Cumulative indel offset for the given haplotype from all variants
    with pos < target_pos. This is in forward-strand coordinate space.
    """
    offset = 0
    vcf = pysam.VariantFile(vcf_path)
    for rec in vcf:
        if rec.pos >= target_pos:
            break
        sample = list(rec.samples.values())[0]
        gt = sample['GT']
        if gt is None or None in gt:
            continue
        allele_idx = gt[haplotype]
        if allele_idx == 0:
            continue
        ref = rec.ref
        alt = rec.alleles[allele_idx]
        offset += len(alt) - len(ref)
    vcf.close()
    return offset


def main():
    gene_coords = [int(x) for x in open(GENE_COORDS_FILE).read().splitlines()]
    is_minus_strand = gene_coords[0] > gene_coords[-1]

    allele_1, allele_2 = read_fasta_alleles(RAW_FASTA_DIR)
    total_len = len(allele_1)  # both haplotypes same length for this case

    print(f"Gene coords: {len(gene_coords)} positions, range {gene_coords[0]}-{gene_coords[-1]}")
    print(f"Minus strand: {is_minus_strand}")
    print(f"Raw vcf2fasta allele_1: {len(allele_1)} bp")
    print(f"Raw vcf2fasta allele_2: {len(allele_2)} bp")
    print(f"Reference gene length:  {len(gene_coords)} bp")

    # ══════════════════════════════════════════════════════
    # Old method (reference-based indices from gene_coords)
    # ══════════════════════════════════════════════════════
    ref_idx_start = gene_coords.index(CLAMPED_START)
    ref_idx_stop = gene_coords.index(CLAMPED_STOP)
    old_fasta_start, old_fasta_stop = sorted((ref_idx_start, ref_idx_stop))
    old_allele_1 = allele_1[old_fasta_start:old_fasta_stop + 1]
    old_allele_2 = allele_2[old_fasta_start:old_fasta_stop + 1]

    print(f"\n{'='*60}")
    print(f"OLD METHOD (reference-based indices)")
    print(f"{'='*60}")
    print(f"gene_coords.index({CLAMPED_START}) = {ref_idx_start}")
    print(f"gene_coords.index({CLAMPED_STOP}) = {ref_idx_stop}")
    print(f"fasta_start={old_fasta_start}, fasta_stop={old_fasta_stop}")
    print(f"Allele 1: {len(old_allele_1)} bp")
    print(f"  starts with: {old_allele_1[:20]}")
    print(f"  ends with:   {old_allele_1[-20:]}")
    print(f"Allele 2: {len(old_allele_2)} bp")
    print(f"  starts with: {old_allele_2[:20]}")
    print(f"  ends with:   {old_allele_2[-20:]}")

    # ══════════════════════════════════════════════════════
    # New method (indel-aware, strand-aware)
    # ══════════════════════════════════════════════════════
    #
    # vcf2fasta for minus strand genes:
    #   1. Fetches forward-strand reference for [gene_start, gene_stop]
    #   2. Applies variants (forward-strand coordinates)
    #   3. Reverse complements the result
    #
    # So in the intermediate (pre-RC) forward-strand sequence:
    #   fwd_index(pos) = (pos - GENE_START) + offset_before(pos)
    #
    # After RC, FASTA index for a given forward-strand position:
    #   fasta_index(pos) = total_len - 1 - fwd_index(pos)
    #
    # We want genomic interval [CLAMPED_START, CLAMPED_STOP].
    # In forward strand:
    #   fwd_low  = fwd_index(CLAMPED_START)   (lowest genomic coord in interval)
    #   fwd_high = fwd_index(CLAMPED_STOP)    (highest genomic coord in interval)
    # After RC:
    #   fasta_start = total_len - 1 - fwd_high  (high coord → start of RC'd FASTA)
    #   fasta_stop  = total_len - 1 - fwd_low   (low coord → end of RC'd FASTA)

    print(f"\n{'='*60}")
    print(f"NEW METHOD (indel-aware, strand-aware)")
    print(f"{'='*60}")

    for hap_idx, allele, label in [(0, allele_1, "Allele 1"), (1, allele_2, "Allele 2")]:
        # Offset from indels before each boundary
        offset_low = compute_offset_at_position(VCF_PATH, hap_idx, CLAMPED_START)
        offset_high = compute_offset_at_position(VCF_PATH, hap_idx, CLAMPED_STOP + 1)

        # Forward-strand indices
        fwd_low = (CLAMPED_START - GENE_START) + offset_low
        fwd_high = (CLAMPED_STOP - GENE_START) + offset_high

        print(f"\n{label} (haplotype {hap_idx}):")
        print(f"  offset_before({CLAMPED_START}) = {offset_low}")
        print(f"  offset_before({CLAMPED_STOP + 1}) = {offset_high}")
        print(f"  fwd_low (intermediate idx for {CLAMPED_START}) = {fwd_low}")
        print(f"  fwd_high (intermediate idx for {CLAMPED_STOP}) = {fwd_high}")

        if is_minus_strand:
            fasta_start = total_len - 1 - fwd_high
            fasta_stop = total_len - 1 - fwd_low
            print(f"  Minus strand transform:")
            print(f"    fasta_start = {total_len} - 1 - {fwd_high} = {fasta_start}")
            print(f"    fasta_stop  = {total_len} - 1 - {fwd_low} = {fasta_stop}")
        else:
            fasta_start = fwd_low
            fasta_stop = fwd_high
            print(f"  Plus strand (direct mapping):")
            print(f"    fasta_start = {fasta_start}")
            print(f"    fasta_stop  = {fasta_stop}")

        clamped = allele[fasta_start:fasta_stop + 1]
        print(f"  Clamped: {len(clamped)} bp")
        print(f"    starts with: {clamped[:20]}")
        print(f"    ends with:   {clamped[-20:]}")

    # ══════════════════════════════════════════════════════
    # Validation against reference
    # ══════════════════════════════════════════════════════
    print(f"\n{'='*60}")
    print(f"VALIDATION")
    print(f"{'='*60}")

    ref = pysam.FastaFile(REFERENCE)

    # For minus strand: FASTA starts at CLAMPED_STOP (high coord, RC'd)
    # The first bases of the FASTA should match RC of reference at the high end
    ref_at_high = ref.fetch("chr6", CLAMPED_STOP - 10, CLAMPED_STOP)  # 0-based
    rc_ref_at_high = revcomp(ref_at_high)
    print(f"\nReference at gene high end (chr6:{CLAMPED_STOP-9}-{CLAMPED_STOP}): {ref_at_high}")
    print(f"RC of that: {rc_ref_at_high}")
    print(f"Old allele_1 starts with: {old_allele_1[:10]}")

    # FASTA ends at CLAMPED_START (low coord, RC'd)
    # Last bases should match RC of reference at clamp start
    ref_at_low = ref.fetch("chr6", CLAMPED_START - 1, CLAMPED_START + 9)  # 0-based
    rc_ref_at_low = revcomp(ref_at_low)
    print(f"\nReference at clamp start (chr6:{CLAMPED_START}-{CLAMPED_START+9}): {ref_at_low}")
    print(f"RC of that: {rc_ref_at_low}")
    print(f"Old allele_1 ends with:   {old_allele_1[-10:]}")

    # Also check what's at the het position we're trying to exclude
    ref_at_het = ref.fetch("chr6", 33065643 - 1, 33065643 + 9)
    rc_ref_at_het = revcomp(ref_at_het)
    print(f"\nReference at het pos (chr6:33065643-33065652): {ref_at_het}")
    print(f"RC of that: {rc_ref_at_het}")

    ref.close()

    # ══════════════════════════════════════════════════════
    # Summary
    # ══════════════════════════════════════════════════════
    print(f"\n{'='*60}")
    print(f"SUMMARY")
    print(f"{'='*60}")
    print(f"Old method: FASTA[{old_fasta_start}:{old_fasta_stop + 1}] = {len(old_allele_1)} bp")

    # Recompute new for allele_1 for comparison
    offset_low = compute_offset_at_position(VCF_PATH, 0, CLAMPED_START)
    offset_high = compute_offset_at_position(VCF_PATH, 0, CLAMPED_STOP + 1)
    fwd_low = (CLAMPED_START - GENE_START) + offset_low
    fwd_high = (CLAMPED_STOP - GENE_START) + offset_high
    new_fasta_start = total_len - 1 - fwd_high
    new_fasta_stop = total_len - 1 - fwd_low
    new_len = new_fasta_stop - new_fasta_start + 1

    print(f"New method: FASTA[{new_fasta_start}:{new_fasta_stop + 1}] = {new_len} bp")
    print(f"Difference: {len(old_allele_1) - new_len} bp")
    print(f"Old stop index: {old_fasta_stop}, New stop index: {new_fasta_stop}")
    if old_fasta_stop != new_fasta_stop:
        print(f"The old method includes {old_fasta_stop - new_fasta_stop} extra base(s) at the stop end")
        print(f"These bases are from the genomic region BELOW {CLAMPED_START} (outside safe interval)")


if __name__ == "__main__":
    main()
