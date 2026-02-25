# Phasing Decision Tree

## Overview

This document describes the decision tree used to classify genes by phasing status and determine how allele sequences are reconstructed. The logic is implemented across `investigate_haploblocks_methods.py` (classification) and `reconstruct_fasta_methods.py` (VCF filtering and FASTA extraction).

---

## Decision Tree

1. **Gene has ≤1 het genotype**
   - Classify as `fully_phased`
   - Whitelist the unphased genotype, run vcf2fasta
   - Match to the HLA database with the full vcf2fasta output

2. **Gene fully spanned by single extended haploblock**
   - Classify as `fully_phased`
   - Remove unphased genotypes, run vcf2fasta
   - Match to the HLA database with the full vcf2fasta output

3. **Gene NOT fully spanned** — enter rescue mode

   a. **ARS fully spanned by a single extended haploblock**
      - Find the largest ARS-spanning block
      - Remove unphased genotypes, run vcf2fasta
      - Subset the vcf2fasta output to the haploblock coordinates.
      - Match to the HLA database with the subsetted sequence. 

   b. **ARS NOT spanned** — CDS-aware rescue

      Count all QC-pass heterozygous genotypes in the CDS and ARS CDS regions. These come from the existing `gene_het_sites`, which is already filtered by `parse_haploblocks()`.

      i. **Total CDS hets <= 1**
         - Whitelist all unphased heterozygous genotypes in `filter_vcf_gene`
         - Run vcf2fasta
         - CDS: write full concatenated exon output (vcf2fasta --feat CDS) to `_CDS.fasta`
         - Gene (vcf2fasta --feat gene):
           - If there is 1 CDS heterozygous genotype, anchor on the CDS heterozygous genotype position. In the sorted list of all heterozygous genotype positions in the gene, find the heterozygous genotype immediately before and immediately after the CDS heterozygous genotype. The interval is (prev_het + 1) to (next_het - 1), using gene boundaries as fallback if no prev/next het exists. If this interval fully contains ARS, write to `_gene.fasta`. If not, do not write a record to `_gene.fasta` (no 4th field assignment will be made).
           - If there are zero CDS heterozygous genotypes, all heterozygous genotypes are intronic. Find the largest ARS-overlapping interval containing <=1 heterozygous genotype:
             1. Sort all het positions in the gene: h1, h2, ..., hn
             2. For each het h_i, compute interval_i = [h_(i-1) + 1, h_(i+1) - 1] (gene boundaries as fallback) — this is the maximal region containing only h_i
             3. Filter to intervals overlapping ARS
             4. Pick the largest; clamp gene output to that interval

             If no interval overlaps ARS, do not write gene record. (With <=1 het, phase is irrelevant — vcf2fasta random assignment produces both haplotypes correctly.)

      ii. **Total CDS heterozygous genotypes > 1**

          A. **ARS CDS hets <= 1**
             - Whitelist all heterozygous genotypes in `filter_vcf_gene`
             - Run vcf2fasta
             - CDS: extract ARS CDS sequence only, write to `_CDS.fasta`
             - Gene:
               - Class II: extend ARS (exon 2) outward to the nearest heterozygous genotypes in each direction; write that interval to `_gene.fasta`.
               - Class I: check if any heterozygous genotype exists in intron 2 (gap between the two ARS CDS exons).
                 - No intron 2 heterozygous genotype: extend ARS outward to the nearest heterozygous genotype in each direction; write that interval to `_gene.fasta`
                 - Heterozygous genotype in intron 2: do not write gene record

          B. **ARS CDS hets > 1**
             - Classify as `do_not_type_genes` (no HLA typing performed)
