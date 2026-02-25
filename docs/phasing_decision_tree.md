# Phasing Decision Tree

## Overview

This document describes the decision tree used to classify genes by phasing status and determine how allele sequences are reconstructed. The logic is implemented across `investigate_haploblocks_methods.py` (classification) and `reconstruct_fasta_methods.py` (VCF filtering and FASTA extraction).

---

## Decision Tree

1. **Gene has ≤1 het genotype**
   - Classify as `fully_phased`
   - Whitelist the unphased genotype, run vcf2fasta
   - Match with full output

2. **Gene fully spanned by single extended haploblock**
   - Classify as `fully_phased`
   - Remove unphased genotypes, run vcf2fasta
   - Match with full output

3. **Gene NOT fully spanned** — enter rescue mode

   a. **ARS fully spanned by a single extended haploblock**
      - Find largest ARS-spanning block
      - Remove unphased genotypes, run vcf2fasta
      - Subset output to haploblock coordinates

   b. **ARS NOT spanned** — CDS-aware rescue

      Count all QC-pass het sites (phased + unphased, QUAL>=10 / GQ>=20) in CDS and ARS CDS regions. These come from the existing `gene_het_sites`, which is already filtered by `parse_haploblocks()`.

      i. **Total CDS hets <= 1**
         - Whitelist all unphased hets in `filter_vcf_gene`
         - Run vcf2fasta
         - CDS: write full concatenated exon output to `_CDS.fasta`
         - Gene:
           - If CDS hets = 1: anchor on the CDS het position. In the sorted list of all het positions in the gene, find the het immediately before and immediately after the CDS het. The interval is (prev_het + 1) to (next_het - 1), using gene boundaries as fallback if no prev/next het exists. If this interval fully contains ARS, write to `_gene.fasta`. If not, do not write.
           - If CDS hets = 0: all hets are intronic. Find the largest ARS-overlapping interval containing <=1 het:
             1. Sort all het positions in the gene: h1, h2, ..., hn
             2. For each het h_i, compute interval_i = [h_(i-1) + 1, h_(i+1) - 1] (gene boundaries as fallback) — this is the maximal region containing only h_i
             3. Filter to intervals overlapping ARS
             4. Pick the largest; clamp gene output to that interval

             If no interval overlaps ARS, do not write gene record. (With <=1 het, phase is irrelevant — vcf2fasta random assignment produces both haplotypes correctly.)

      ii. **Total CDS hets > 1**

          A. **ARS CDS hets <= 1**
             - Whitelist all hets in `filter_vcf_gene`
             - Run vcf2fasta
             - CDS: extract ARS CDS sequence only, write to `_CDS.fasta`
             - Gene:
               - Class II: extend ARS outward to the nearest het (any het) in each direction; write that interval to `_gene.fasta`. (No intron concern since Class II ARS is a single exon.)
               - Class I: check if any het (phased or unphased) exists in intron 2 (gap between the two ARS CDS exons).
                 - No intron 2 het: extend ARS outward to the nearest het in each direction; write that interval to `_gene.fasta`
                 - Het in intron 2: do not write gene record

          B. **ARS CDS hets > 1**
             - Classify as `do_not_type_genes` (no typing performed)
