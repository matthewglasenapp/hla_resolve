#!/usr/bin/env python3
"""
Test script to validate indel-aware clamping for CLQC1903 HLA-DQA1 (plus strand).

DQA1 is partially phased with haploblock 32640816-32655563.
Gene is 32628179-32647062 (plus strand).
Effective clamp: 32640816-32647062 (trims ~12.6kb from 5' end).
39 indels before the clamp point.

Usage:
    python test_clamp_offset_dqa1.py
"""

import pysam
import glob

# ── Paths (CLQC1903_19QC-1 DQA1) ──
VCF_PATH = "/hb/groups/cornejo_lab/matt/ihiws_2026/data/hla_resolve_results/CLQC1903_19QC-1/filtered_vcf/CLQC1903_19QC-1_HLA-DQA1_PASS_phased.vcf.gz"
RAW_FASTA_DIR = "/hb/groups/cornejo_lab/matt/ihiws_2026/data/hla_resolve_results/CLQC1903_19QC-1/vcf2fasta_out/HLA-DQA1_gene/"
REFERENCE = "/hb/groups/cornejo_lab/matt/hla_resolve/hla_resolve/data/reference/GCA_000001405.15_GRCh38_no_alt_analysis_set.fna"
GENE_COORDS_FILE = "/hb/groups/cornejo_lab/matt/hla_resolve/hla_resolve/data/hla_gff/hla_dqa1_gene_coords.txt"

# ── DQA1 parameters ──
GENE_START = 32628179   # GFF start
GENE_STOP = 32647062    # GFF stop
# Haploblock is 32640816-32655563, clamped to gene boundary
CLAMPED_START = 32640816  # haploblock start
CLAMPED_STOP = 32647062   # gene stop (haploblock extends beyond)


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
    comp = str.maketrans("ACGTacgt", "TGCAtgca")
    return seq.translate(comp)[::-1]


def compute_offset_at_position(vcf_path, haplotype, target_pos):
    """
    Cumulative indel offset for the given haplotype from all variants
    with pos < target_pos.
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

    print(f"Gene: HLA-DQA1 (CLQC1903_19QC-1)")
    print(f"Gene coords: {len(gene_coords)} positions, range {gene_coords[0]}-{gene_coords[-1]}")
    print(f"Strand: {'minus' if is_minus_strand else 'plus'}")
    print(f"Raw vcf2fasta allele_1: {len(allele_1)} bp")
    print(f"Raw vcf2fasta allele_2: {len(allele_2)} bp")
    print(f"Reference gene length:  {len(gene_coords)} bp")
    print(f"Clamp interval: {CLAMPED_START}-{CLAMPED_STOP}")
    print(f"Expected clamp size (reference): {CLAMPED_STOP - CLAMPED_START + 1} bp")

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
    print(f"  starts with: {old_allele_1[:30]}")
    print(f"  ends with:   {old_allele_1[-30:]}")
    print(f"Allele 2: {len(old_allele_2)} bp")
    print(f"  starts with: {old_allele_2[:30]}")
    print(f"  ends with:   {old_allele_2[-30:]}")

    # ══════════════════════════════════════════════════════
    # New method (indel-aware, strand-aware)
    # ══════════════════════════════════════════════════════
    print(f"\n{'='*60}")
    print(f"NEW METHOD (indel-aware, strand-aware)")
    print(f"{'='*60}")

    for hap_idx, allele, label in [(0, allele_1, "Allele 1"), (1, allele_2, "Allele 2")]:
        total_len = len(allele)

        offset_start = compute_offset_at_position(VCF_PATH, hap_idx, CLAMPED_START)
        offset_stop = compute_offset_at_position(VCF_PATH, hap_idx, CLAMPED_STOP + 1)

        if is_minus_strand:
            fwd_low = (CLAMPED_START - GENE_START) + offset_start
            fwd_high = (CLAMPED_STOP - GENE_START) + offset_stop
            fasta_start = total_len - 1 - fwd_high
            fasta_stop = total_len - 1 - fwd_low
        else:
            # Plus strand: direct mapping
            fasta_start = (CLAMPED_START - GENE_START) + offset_start
            fasta_stop = (CLAMPED_STOP - GENE_START) + offset_stop

        clamped = allele[fasta_start:fasta_stop + 1]

        print(f"\n{label} (haplotype {hap_idx}):")
        print(f"  offset_before({CLAMPED_START}) = {offset_start}")
        print(f"  offset_before({CLAMPED_STOP + 1}) = {offset_stop}")
        print(f"  fasta_start = {fasta_start}")
        print(f"  fasta_stop  = {fasta_stop}")
        print(f"  Clamped: {len(clamped)} bp")
        print(f"    starts with: {clamped[:30]}")
        print(f"    ends with:   {clamped[-30:]}")

    # ══════════════════════════════════════════════════════
    # Validation against reference
    # ══════════════════════════════════════════════════════
    print(f"\n{'='*60}")
    print(f"VALIDATION")
    print(f"{'='*60}")

    ref = pysam.FastaFile(REFERENCE)

    # For plus strand: FASTA start should correspond to CLAMPED_START
    ref_at_start = ref.fetch("chr6", CLAMPED_START - 1, CLAMPED_START + 19)
    print(f"\nReference at clamp start (chr6:{CLAMPED_START}-{CLAMPED_START+19}): {ref_at_start}")
    print(f"Old allele_1 starts with:  {old_allele_1[:20]}")

    # New method allele_1
    offset_s = compute_offset_at_position(VCF_PATH, 0, CLAMPED_START)
    offset_e = compute_offset_at_position(VCF_PATH, 0, CLAMPED_STOP + 1)
    new_start = (CLAMPED_START - GENE_START) + offset_s
    new_stop = (CLAMPED_STOP - GENE_START) + offset_e
    new_allele_1 = allele_1[new_start:new_stop + 1]
    print(f"New allele_1 starts with:  {new_allele_1[:20]}")

    # End of clamp
    ref_at_stop = ref.fetch("chr6", CLAMPED_STOP - 19, CLAMPED_STOP)
    print(f"\nReference at clamp stop (chr6:{CLAMPED_STOP-19}-{CLAMPED_STOP}): {ref_at_stop}")
    print(f"Old allele_1 ends with:    {old_allele_1[-20:]}")
    print(f"New allele_1 ends with:    {new_allele_1[-20:]}")

    ref.close()

    # ══════════════════════════════════════════════════════
    # Summary
    # ══════════════════════════════════════════════════════
    print(f"\n{'='*60}")
    print(f"SUMMARY")
    print(f"{'='*60}")
    print(f"Old method: FASTA[{old_fasta_start}:{old_fasta_stop + 1}] = {len(old_allele_1)} bp")
    print(f"New method: FASTA[{new_start}:{new_stop + 1}] = {len(new_allele_1)} bp")
    diff = new_start - old_fasta_start
    print(f"Start index shift: {diff} bp (old={old_fasta_start}, new={new_start})")
    print(f"Length difference: {len(new_allele_1) - len(old_allele_1)} bp")
    if diff != 0:
        print(f"\nThe old method starts {abs(diff)} bp {'early' if diff > 0 else 'late'} in the FASTA")
        print(f"For plus strand, this means the old method includes {abs(diff)} bp")
        print(f"from OUTSIDE the haploblock (genomic positions before {CLAMPED_START})")


if __name__ == "__main__":
    main()
