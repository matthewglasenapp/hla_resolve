# HLA-Resolve Technical Reference

Detailed documentation on the algorithms, decision logic, and tools used by HLA-Resolve.

---

## 1. Preparing Phased VCF for vcf2fasta

Following phasing with HiPhase or longphase, the SNV and SV VCFs are merged and normalized with bcftools. The following steps are taken to prepare the VCF records for use with vcf2fasta.
1. Filter the joint phased VCF by gene. The gene coordinates are defined by the GRCh38 GFF3 records.

2. Send all symbolic variants to SAMPLE_GENE.symbolic.vcf.gz. These variant types (e.g., TRID, BND, DUP, INV) are not currently compatible with vcf2fasta. These variants will not participate in quality filtering or be used downstream, but are kept for record.

3. Suppress large indels overlapping PASS pbsv structural variants (haplotype-aware). This step addresses the known gray area in variant size and representation where small-variant genotypers and structural-variant genotypers may both emit overlapping calls for the same underlying insertion/deletion, leading to redundant or conflicting genotypes in the merged VCF. During per-gene post-phasing filtering, the pipeline identifies PASS pbsv structural variants and evaluates overlap with indels from bcftools or DeepVariant. Indels from small-variant genotypers are suppressed in rare cases only if they are ≥30 bp, overlap a pbsv structural variant on the reference, and both variants are assigned to the same haplotype (i.e., the indel must be phased). SNPs and small indels are never suppressed. This prevents redundant or conflicting variant calls within regions affected by structural variation. This logic is applied only to pbsv SV genotypes that are non-symbolic, FILTER=PASS, and phased or homozygous alternate. Suppressed indels are written to SAMPLE_GENE_sv_overlap.vcf.gz. These records are excluded from downstream processing (including vcf2fasta).

4. Apply hard filters to all non-symbolic variants and send QC-pass variants to SAMPLE_GENE_PASS.vcf.gz and QC-fail variants to SAMPLE_GENE_FAIL.vcf.gz. The hard filters for SNVs are QUAL ≥ 10, GQ ≥ 20, and DP ≥ 2. The hard filters for bcftools indels are GQ ≥ 10 and DP ≥ 2, with GQ allowed to be missing. No QUAL threshold is applied to indels. Pbsv non-symbolic SV genotypes (e.g., insertion, deletion) are filtered separately, scanning only for FILTER=PASS. If using DeepVariant as the small variant genotyper, genotypes are filtered separately using only FILTER=PASS, with no additional QUAL, GQ, or DP thresholds applied.

5. Count heterozygous genotypes to address edge-case of a single heterozygous genotype that is unphased. If there is only one heterozygous genotype and it is unphased, generate a whitelist expression that allows it to be retained in the QC-pass phased genotype VCF.

6. Filter the QC-pass genotypes (SAMPLE_GENE_PASS.vcf.gz) to remove unphased heterozygous genotypes. If there is only one heterozygous genotype and it is unphased, it will be retained due to the whitelist rule. The QC-pass genotypes with unphased heterozygous genotypes removed are written to SAMPLE_GENE_PASS_phased.vcf.gz. Unphased heterozygous genotypes are written to SAMPLE_GENE_PASS_UNPHASED.vcf.gz.

7. Both symbolic genotypes and unphased QC-pass heterozygous genotypes that overlap an HLA gene will be reported. These records cannot be incorporated into the vcf2fasta haplotype reconstruction, but may confer important information.

---

## 2. Notes on the Use of vcf2fasta

HLA-Resolve uses a forked version of vcf2fasta, originally written by Santiago Sanchez-Ramirez (https://github.com/santiagosnchez/vcf2fasta, MIT license). The fork is maintained at https://github.com/matthewglasenapp/vcf2fasta and is installed as a pip dependency via `environment.yml`. The fork includes several bug fixes described below.

### CDS ordering and GFF preprocessing

vcf2fasta concatenates CDS records in the order they appear in the GFF file, then applies reverse complement for minus strand genes. To ensure correct ordering, `supplementary_scripts/sort_cds.py` extracts CDS records from the raw GRCh38 GFF3 annotations and sorts them by genomic coordinate (low to high). The pre-sorted GFF3 files are included in `hla_resolve/data/hla_gff`. In that directory, `hla_<gene>.gff3` is the unaltered annotation from GRCh38.p14, and `hla_<gene>_cds_sorted.gff3` is the CDS-only sorted version used as input for vcf2fasta. `hla_<gene>_gene.gff3` contains the gene-level record used for full-gene reconstruction.

### Forced phased mode

vcf2fasta checks the first genotype of the overall VCF to determine phase status. If the first variant is unphased, it treats the whole VCF as unphased and outputs one diploid sequence with IUPAC codes representing heterozygous bases. The fork hardcodes vcf2fasta to operate in "phased" mode. The HLA-Resolve pipeline ensures this is appropriate, as it first ensures each gene is completely spanned by a haplotype block, and features special treatment for genes that have incomplete phasing.

```diff
- phased = v2f_helper.getPhased(vcf)
- if not phased:
-     print('No phased genotypes found on first variant. Treating as "unphased"')
- else:
-     print('Phased genotypes found on first variant. Treating as "phased"')
+ # Patch: forcibly treat all input as phased
+ phased = True
+ print('Treating as "phased"')
```

When vcf2fasta treats a VCF as "phased," it treats unphased heterozygous genotypes as phased and arbitrarily assigns each allele to a haplotype. This means that unphased heterozygous genotypes must be removed from the phased VCF prior to running vcf2fasta. However, if there is only one heterozygous genotype in the single-gene VCF and it remains unphased, it should not be filtered before vcf2fasta. Because we care only about phased haplotypes within genes, it doesn't matter which haplotype gets which allele in this scenario. The two reconstructed haplotypes will be correct. However, the standard pipeline removes unphased heterozygous genotypes from the VCF before vcf2fasta, because when operating in "phased" mode, if there is a mix of phased and unphased heterozygous genotypes, vcf2fasta randomly assigns alleles from the unphased heterozygous genotypes to haplotypes, resulting in inevitable switch errors. The edge case of a single unphased heterozygous genotype is addressed in `reconstruct_fasta_methods.py`.

### Bug fixes in the fork

1. **`getPloidy()` and `getPhased()` robustness**: Originally, both functions relied on the first variant record to infer ploidy and phasing status. If that first record contained missing genotypes or was of an unexpected format, the program exited with an error. The updated versions iterate through the VCF until they encounter a valid genotype.

    **`getPloidy()` (`v2f/functions.py`):**

    ```diff
    - def getPloidy(vcf):
    -     var = [ y for x,y in next(vcf.fetch()).samples.items() ]
    -     p = [ len(v.get('GT')) for v in var if v.get('GT')[0] is not None ]
    -     return p[0]
    + def getPloidy(vcf):
    +     for rec in vcf.fetch():
    +         for sample in rec.samples.values():
    +             gt = sample.get('GT')
    +             if gt and all(g is not None for g in gt):
    +                 return len(gt)
    +     # Default to diploid if no variants found (empty VCF)
    +     return 2
    ```

    **`getPhased()` (`v2f/functions.py`):**

    ```diff
    - def getPhased(vcf):
    -     var = [ y for x,y in next(vcf.fetch()).samples.items() ]
    -     p = any([ not v.phased for v in var ])
    -     return not p
    + def getPhased(vcf):
    +     for rec in vcf.fetch():
    +         for sample in rec.samples.values():
    +             gt = sample.get('GT')
    +             if gt and all(g is not None for g in gt):
    +                 return sample.phased
    +     raise ValueError("No phased genotypes found to determine phasing.")
    ```

2. **Right-to-left variant application**: The original indel-handling code applied variants from low to high genomic coordinates using a running offset to adjust positions after insertions and deletions. This can break when a SNP and an indel overlap, causing the indel to be skipped. The fix removes the shifting offset and applies variants from right to left (highest to lowest coordinate).

    **`getSequences()` (`v2f/functions.py`):**

    ```diff
    -        posadd = 0
    -        for vrec in vcf.fetch(chrom, start, end):
    -            pos = vrec.pos - start - 1 + posadd
    -            alleles, max_len = getAlleles(vrec, ploidy, phased, args.addref)
    -            ref_len = len(vrec.ref)
    -            for sample in alleles.keys():
    -                tmpseqs[sample] = UpdateSeq(alleles, sample, pos, ref_len, tmpseqs[sample])
    -            posadd += (max_len - ref_len)
    +        vrec_list = list(vcf.fetch(chrom, start, end))
    +        vrec_list.sort(key=lambda r: r.pos, reverse=True)
    +
    +        for vrec in vrec_list:
    +            pos = vrec.pos - start - 1
    +            alleles, max_len = getAlleles(vrec, ploidy, phased, args.addref)
    +            ref_len = len(vrec.ref)
    +            for sample in alleles.keys():
    +                tmpseqs[sample] = UpdateSeq(
    +                    alleles,
    +                    sample,
    +                    pos,
    +                    ref_len,
    +                    tmpseqs[sample]
    +                )
    ```

3. **Empty VCF handling**: The original code attempted to get sample names from the first variant record, which fails with a `StopIteration` error when the VCF has no variants. The fix retrieves sample names from the VCF header, which is always present.

    **`cli.py`:**

    ```diff
    - samples = v2f_helper.get_samples(vcf)
    + # Get a list of samples from header (works even if VCF has no variants)
    + samples = list(vcf.header.samples)
    ```

4. **Default ploidy for empty VCFs**: `getPloidy()` returns a default ploidy of 2 (diploid) when the VCF contains no variants, rather than raising an error. This allows vcf2fasta to output reference sequences for genes with no variants. (See `getPloidy()` diff above.)

### Additional vcf2fasta notes

1. Multiple transcripts exist for all HLA genes in the hg38 GFF3. If you don't select one, vcf2fasta will make FASTA files for all annotated transcripts.

2. vcf2fasta applies all variants it sees when generating the FASTA. It does not check the FILTER status of records. Therefore, it is important to filter the VCF prior to running vcf2fasta.

3. Before running vcf2fasta in phased mode, you need to ensure that the genes you're making FASTAs for are fully contained within a single phase set. If there are two or more phase sets contained within your gene, it could cause a switch error.

---

## 3. Phasing Decision Tree

### Overview

This section describes the decision tree used to classify genes by phasing status and determine how allele sequences are reconstructed. The logic is implemented across `investigate_haploblocks_methods.py` (classification) and `reconstruct_fasta_methods.py` (VCF filtering and FASTA extraction).

### Decision Tree

1. **Gene has ≤1 het genotype**
   - Classify as `fully_phased`
   - Whitelist the unphased genotype, run vcf2fasta
   - Match to the HLA database with the full vcf2fasta output

2. **Gene fully spanned by single extended haploblock**
   - Classify as `fully_phased`
   - Remove unphased genotypes, run vcf2fasta
   - Match to the HLA database with the full vcf2fasta output

3. **Gene NOT fully spanned** — enter rescue mode

   - **a. ARS fully spanned by a single extended haploblock**
     - Find the largest ARS-spanning block
     - Remove unphased genotypes, run vcf2fasta
     - Subset the vcf2fasta output to the haploblock coordinates.
     - Match to the HLA database with the subsetted sequence.

     Note: Because vcf2fasta applies indels that change the sequence length relative to the reference, subsetting the output by genomic coordinates requires computing per-haplotype indel offsets from the filtered VCF. For minus-strand genes, the vcf2fasta output is reverse-complemented, so the coordinate transformation must account for strand orientation. This applies to all gene-level FASTA subsetting (sections 3a and 3b).

   - **b. ARS NOT spanned** — CDS-aware rescue

     Count all QC-pass heterozygous genotypes in the CDS and ARS CDS regions. These come from the existing `gene_het_sites`, which is already filtered by `parse_haploblocks()`.

     - **i. Total CDS hets <= 1**
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

     - **ii. Total CDS heterozygous genotypes > 1**

       - **A. ARS CDS hets <= 1**
         - Whitelist all heterozygous genotypes in `filter_vcf_gene`
         - Run vcf2fasta
         - CDS: extract ARS CDS sequence only, write to `_CDS.fasta`
         - Gene:
           - Class II: extend ARS (exon 2) outward to the nearest heterozygous genotypes in each direction; write that interval to `_gene.fasta`.
           - Class I: check if any heterozygous genotype exists in intron 2 (gap between the two ARS CDS exons).
             - No intron 2 heterozygous genotype: extend ARS outward to the nearest heterozygous genotype in each direction; write that interval to `_gene.fasta`
             - Heterozygous genotype in intron 2: do not write gene record

       - **B. ARS CDS hets > 1**
         - Classify as `do_not_type_genes` (no HLA typing performed)

---

## 4. Allele Classification Methodology (`hla_typer.py`)

This section describes the allele classification algorithm implemented in `hla_typer.py`, which performs HLA allele assignment after haplotype reconstruction from phased FASTA sequences.

The typing workflow operates downstream of variant calling, phasing, and FASTA reconstruction, and is responsible for assigning HLA alleles at G-group, three-field, and four-field resolution.

---

### 4.1 Reference Data

Reference data are derived from the [IPD-IMGT/HLA](https://www.ebi.ac.uk/ipd/imgt/hla/) database.

The XML reference file (`hla.xml`) is parsed to extract:
- Allele nucleotide sequences
- Exon, intron, and UTR annotations
- Peptide-binding domain definitions
- G-group membership metadata

This metadata is used to:
- Reconstruct reference exon and full-gene sequences
- Generate allele lookup tables by G group
- Validate that all alleles within a G group share identical peptide-binding domain sequences

---

### 4.2 Inputs to the Typing Algorithm

The typing script consumes three primary inputs, all generated by upstream stages of the HLA-Resolve pipeline:

1. **IMGT/HLA reference database**
   - XML format (`hla.xml`)
   - Contains all reference allele sequences and annotations

2. **Concatenated exon FASTA**
   - Required input
   - Contains reconstructed, phased exon sequences for each sample allele
   - Used for G-group assignment and three-field typing

3. **Full-gene FASTA (optional)**
   - Includes introns and UTRs
   - Enables refinement to four-field resolution
   - If not provided, four-field typing is skipped

The pipeline will run successfully with exon-only FASTA input, but accurate four-field assignments require full-gene sequence data.

---

### 4.3 Typing Outputs

The typing program produces both CSV result files and detailed log files documenting classification decisions.

#### Primary Output Files

- **`g_group_output.csv`**
  G-group classification for each allele.
  Unassigned entries are left blank.
  Alleles are formatted using standard HLA nomenclature (without the `HLA-` prefix).

- **`3_field_allele_output.csv`**
  Three-field–accurate allele assignments for each sample.
  Includes a placeholder fourth field for compatibility with reference-matching tools.
  These fourth-field values do **not** necessarily reflect true four-field assignments.

- **`allele_output.csv`**
  Final four-field–accurate allele assignments (when full-gene sequence data are provided).
  Uses standard HLA nomenclature including the `HLA-` prefix.

#### Log Files

- **`g_group_assignment.log`**
  Records whether G-group assignment succeeded for each allele, including edit distance values.

- **`3_field_allele_assignment.log`**
  Logs the assignment metric used, metric values, and any equidistant candidate alleles.

- **`allele_assignment.log`**
  Documents fourth-field refinement logic, including:
  - Whether refinement was attempted
  - Alleles considered during classification
  - Final assignment metrics

#### Auxiliary Data

- **`database_data.json`**
  Serialized reference data for debugging and testing, including:
  - Feature-level sequence data
  - G-group membership tables

---

### 4.4 Hierarchical Classification Strategy

The typing algorithm employs a **three-pass hierarchical approach**:

#### Pass 1: G-group Assignment
- Compares the query peptide-binding domain sequence to stored G-group peptide-binding domain sequences
- Uses gap-compressed edit distance
- A G group is assigned only if edit distance equals zero
- If no G group is assigned, the algorithm proceeds with an unrestricted search

#### Pass 2: Three-field Assignment
- Compares concatenated exon sequences of the query allele against reconstructed exon sequences of candidate reference alleles
- Search space is restricted to the assigned G group when available
- Produces accurate three-field allele assignments

#### Pass 3: Four-field Refinement (optional)
- Performed only if full-gene FASTA input is provided
- Ignores the previously assigned fourth field and searches using a wildcard (e.g., `C*02:02:02:XX`)
- Compares full-gene sequences including introns and UTRs
- Uses sequence identity rather than raw edit distance

---

### 4.5 Assignment Metrics

#### Gap-compressed Edit Distance
During passes 1 and 2, similarity is evaluated using a modified Levenshtein edit distance computed with [edlib](https://github.com/Martinsos/edlib):

- Calculated in infix mode (free end gaps)
- Continuous insertions or deletions of any length are penalized as a single event
- This reflects a parsimonious mutation model for indels

#### Sequence Identity (Pass 3)
For four-field refinement, sequence identity is used: identity = 1 - (edit_distance / match_length)


This avoids bias toward incomplete reference alleles when comparing noncoding regions, where higher mutation rates are expected.

---

### 4.6 Reference Sequence Reconstruction

Reference exon and full-gene sequences are reconstructed dynamically using annotation metadata from `hla.xml`.

This includes:
- Ordered exon concatenation
- Intron and UTR integration for full-gene sequences
- Peptide-binding domain extraction

Sanity checks ensure consistency between G-group definitions and peptide-binding domain sequences.
