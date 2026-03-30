<br/>

<p align="center">
  <img src="images/hla_resolve.png" alt="HLA-RESOLVE Logo" width="340"/><br/>
  <b>HLA Typing from PacBio Reads</b>
</p>

HLA-Resolve is a command-line tool for high-resolution (four-field) HLA typing from PacBio HiFi sequencing reads. It reconstructs phased, full-gene haplotypes for the eight classical HLA loci (HLA-A, -B, -C, -DPA1, -DPB1, -DQA1, -DQB1, -DRB1) and queries the [IPD-IMGT/HLA database](https://www.ebi.ac.uk/ipd/imgt/hla/) to return star allele calls. The tool is compatible with whole-genome, whole-exome, and targeted sequencing data.

**Authors:** [Matthew Glasenapp](https://github.com/matthewglasenapp), [Alex Symons](https://github.com/FlyingFish800), [Omar Cornejo](https://github.com/oeco28)

**⚠️ Note:** HLA-Resolve is intended for high-coverage PacBio HiFi reads. ONT support is still in development. The HLA-Resolve manuscript is under peer review.

## Requirements

- **Linux (x86_64)** — Several dependencies (pbmarkdup, hiphase, trgt, pbsv, pbmm2) are distributed as precompiled Linux binaries via Bioconda and are not available for macOS.
- **Conda** and **pip** — Used to install all dependencies (see [Installation](#installation)).

## Overview

#### Input
A raw, single-sample (demultiplexed) PacBio sequencing file in FASTQ or unmapped BAM format (compressed or uncompressed). The tool is compatible with WGS, WES, and targeted sequencing schemes.

#### Output(s)

**Primary Results**

HLA star-allele calls for the following genes:

```
HLA-A, HLA-B, HLA-C, HLA-DPA1, HLA-DPB1, HLA-DQA1, HLA-DQB1, HLA-DRB1
```

**Intermediate Files**
- Haplotagged, mapped BAMs for chromosome 6 (for visualization in genome browsers such as IGV)
- Phased VCFs (chromosome 6 and individual gene)
- Reconstructed haplotype nucleotide sequences for each HLA gene in FASTA format

#### Runtime and Required Resources
Runtime depends heavily on input file size and available compute resources. Targeted MHC capture data typically completes in **<30 minutes** using **6 CPUs and 20 GB RAM**. Runtime increases for high-coverage WGS or WES datasets, as all reads must be mapped to the human reference genome prior to restricting downstream analysis to the MHC region on chromosome 6.

## Installation
```
git clone https://github.com/matthewglasenapp/hla_resolve
conda env create -f hla_resolve/environment.yml
conda activate hla_resolve
pip install -e hla_resolve
hla_resolve --help
```
The first time ``hla_resolve`` is executed, it will automatically download the following required files:

| File | Source |
|------|--------|
| GRCh38 reference genome | NCBI |
| picard.jar | Broad Institute |
| LongPhase binary | GitHub |
| hla.xml (IPD-IMGT/HLA database) | IMGTHLA |
| Clair3 Singularity image | Docker Hub |
| DeepVariant Singularity image | Docker Hub |

**Note:** These downloads are large. Ensure sufficient disk space is available in the install directory before the first run.


## Quick Start

```
usage: hla_resolve [-h] --input_file INPUT_FILE --sample_name SAMPLE_NAME --platform {pacbio,ont} --scheme
              {WGS,WES,targeted} --output_dir OUTPUT_DIR [--trim_adapters] [--adapter_file ADAPTER_FILE]
              [--threads THREADS] [--read_group_string READ_GROUP_STRING] [--clean-up]

Run HLA-Resolve

optional arguments:
  -h, --help            show this help message and exit
  --input_file INPUT_FILE
                        Path to the raw sequencing reads file (default: None)
  --sample_name SAMPLE_NAME
                        Override the parsed sample name (default: None)
  --platform {pacbio,ont}
                        Specify sequencing platform (pacbio, ont) (default: None)
  --scheme {WGS,WES,targeted}
                        Sequencing scheme (default: None)
  --output_dir OUTPUT_DIR
                        Output Directory (default: None)
  --trim_adapters       Enable adapter trimming before processing (default: False)
  --adapter_file ADAPTER_FILE
                        Path to a file with custom adapter sequences (FASTA/FASTQ). If not provided, default
                        adapters will be used. (default: None)
  --threads THREADS     Number of threads to use (default: 6)
  --read_group_string READ_GROUP_STRING
                        Override the parsed read group string (default: None)
  --clean-up            Remove intermediate files (default: False)

Example: hla_resolve --input_file reads.bam --sample_name HG002 --platform pacbio --scheme targeted --output_dir out --threads 10

```

## Demo

Input data: PacBio Revio HiFi targeted sequencing reads from HG002 (Ashkenazi Son), a sample from the GIAB and HPRC benchmarks.

```
hla_resolve \
--input_file hla_resolve/demo/HG002.hifi_reads.fastq.gz \
--sample_name HG002 \
--platform pacbio \
--scheme targeted \
--output_dir test \
--trim_adapters \
--adapter_file hla_resolve/demo/adapters.fasta \
--threads 6
```

The command will print the final star allele calls to STDOUT, along with important logging information, including coverage depth metrics, heterozygous genotypes that could not be phased, and the paths of intermediate files (e.g., BAM, VCF).

Intermediate files will be written to the following dirctories. The user can specify the ```--clean-up``` option if they do not want intermediate files, such as mapped BAM, phased genotypes (VCFs), or fasta haplotype nucleotide sequences for the HLA genes.


| Directory                | Description                                                                 |
|---------------------------|-----------------------------------------------------------------------------|
| `fastq_raw/`              | Raw fastq. Converted from BAM format if input is BAM. Copied from raw file if input is fastq |
| `fastq_trimmed/`          | Fastq reads with adapters/barcodes trimmed, if specified by user. If no trimming is specified, will be a copy of the reads in `fastq_raw/` |
| `mapped_bam/`             | Contains BAM files from reference genome alignments                        |
| `genotype_calls/`         | Contains the raw small variant genotype calls (`.vcf.gz`) from the user-specified genotyping tool |
| `structural_variant_vcf/` | Contains the SV genotype calls from either Sniffles (ONT) or pbsv (PacBio) |
| `pbtrgt_vcf/`             | Contains the tandem repeat genotypes from TRGT (PacBio-only)               |
| `phased_vcf/`             | Contains phased genotype calls from joint phasing of small variants, structural variants, and tandem repeat genotypes |
| `mosdepth/`               | Contains coverage depth output files from mosdepth for the HLA genes        |
| `haploblocks/`            | Contains a list of fully-phased MHC genes                                  |
| `filtered_vcf/`           | Contains the final, filtered VCF of variants to be applied during fasta haplotype reconstruction |
| `vcf2fasta_out/`          | Contains the raw sequence output from vcf2fasta                            |
| `hla_fasta_haplotypes/`   | Contains fasta files of full gene and CDS sequences for each HLA gene       |
| `hla_typing_results/`     | Contains the final results of HLA typing                                   |

## Workflow and Dependencies

HLA-Resolve takes raw PacBio HiFi reads (FASTQ or uBAM) as input and executes the following steps to produce four-field HLA star allele calls.

#### 1. Adapter Trimming
Adapter and barcode sequences are removed from raw reads using [cutadapt](https://doi.org/10.14806/ej.17.1.200) (when adapter sequences are provided) or [fastplong](https://doi.org/10.1002/imt2.107) (auto-detection mode).

#### 2. PCR Duplicate Removal
PCR duplicates are identified and removed from the trimmed reads using [pbmarkdup](https://github.com/PacificBiosciences/pbmarkdup).

#### 3. Reference Genome Alignment
Deduplicated reads are aligned to a modified GRCh38 reference genome (no-alt analysis set) using [minimap2](https://doi.org/10.1093/bioinformatics/bty191). The modified reference includes an additional scaffold containing the HLA-Y/HLA-OLI insertion to prevent mismapping of HLA-Y reads to HLA-A.

#### 4. HLA-DRB Paralog Filtering
A separate alignment step maps reads against a multi-allele HLA-DRB reference (containing HLA-DRB1, -DRB3, and -DRB4 sequences from the IPD-IMGT/HLA database) to identify and remove HLA-DRB3/DRB4 reads that would otherwise mismap to HLA-DRB1.

#### 5. Read Filtering
Aligned reads are filtered to retain only primary alignments on chromosome 6, with HLA-DRB3/DRB4 reads excluded.

#### 6. Small Variant Calling
SNVs are called with [bcftools](https://doi.org/10.1093/bioinformatics/btr509) and indels are called with [DeepVariant](https://doi.org/10.1038/nbt.4235). DeepVariant RefCall genotypes with sufficient read support are rescued and reclassified as heterozygous or homozygous ALT based on variant allele frequency.

#### 7. Structural Variant Calling
Structural variants are called from the aligned reads using [pbsv](https://github.com/PacificBiosciences/pbsv).

#### 8. Tandem Repeat Genotyping
Tandem repeats within the MHC region are genotyped using [TRGT](https://doi.org/10.1038/s41587-023-02057-3).

#### 9. Joint Phasing
Small variants, structural variants, and tandem repeat genotypes are jointly phased with [HiPhase](https://doi.org/10.1093/bioinformatics/btae042), producing haplotagged BAMs, phased VCFs, and haplotype block coordinates.

#### 10. Coverage Assessment
Per-gene coverage depth and breadth are calculated with [mosdepth](https://doi.org/10.1093/bioinformatics/btx699). Genes failing minimum coverage thresholds are excluded from HLA typing.

#### 11. Haploblock Evaluation
Phased haplotype blocks are evaluated to determine whether each HLA gene is fully spanned by a single phase set. Genes with internal phasing breaks enter a rescue pipeline that attempts to recover coding sequence or antigen recognition sequence (ARS) haplotypes.

#### 12. Variant Filtering and Redundancy Removal
Phased genotypes are filtered by gene to remove redundant calls from overlapping variant callers (e.g., DeepVariant indels overlapping pbsv structural variants, or non-TRGT variants within TRGT tandem repeat regions). Symbolic and complex structural variant types (BND, INV, DUP) are excluded.

#### 13. Haplotype Reconstruction
Phased, filtered genotypes are applied to GRCh38 gene models using [vcf2fasta](https://github.com/sanchez-ramirez/vcf2fasta) to reconstruct full-gene and coding-sequence haplotype FASTA files for each HLA gene.

#### 14. IPD-IMGT/HLA Database Matching
Reconstructed haplotypes are compared against alleles in the [IPD-IMGT/HLA database](https://doi.org/10.1093/nar/gkac1011) with a three-pass hierarchical classification algorithm using [edlib](https://doi.org/10.1093/bioinformatics/btw753):

   1. **G-group assignment** — The antigen recognition sequence (ARS exons) is matched to G-group reference sequences by edit distance. An exact match is required.
   2. **Three-field allele assignment** — The full concatenated exon sequence is compared against alleles within the assigned G group, ranked by edit distance.
   3. **Four-field refinement** — The full-gene haplotype (including introns and UTRs) is compared against candidate alleles, ranked by mismatch identity (the proportion of matching bases at 1:1-aligned positions), which avoids penalizing insertions and deletions from unreliable intronic reconstruction. Ties are broken by match length, then by lowest fourth-field value.

#### Note
For WGS and WES input, the pipeline skips adapter trimming (step 1) and pre-alignment duplicate removal (step 2), and uses [pbmm2](https://github.com/PacificBiosciences/pbmm2) instead of minimap2 for reference genome alignment.

## Planned Features (In Development)

1. HLA typing at P-group resolution
2. HLA typing for additional HLA Class I protein-coding genes and pseudogenes
   (HLA-E, HLA-F, HLA-G; HLA-H, HLA-J, HLA-K, HLA-L, HLA-S, HLA-V, HLA-W)
3. HLA typing for additional HLA Class II protein-coding genes
   (HLA-DRB3, HLA-DRB4, HLA-DRB5)
4. CYP21A2 star-allele calling

#### Technical Reference
For detailed documentation on the algorithms, decision logic, and tools used internally by HLA-Resolve, see the
[Technical Reference](https://github.com/matthewglasenapp/hla_resolve/blob/main/docs/technical_reference.md).

## Citation

If you use HLA-Resolve, please cite:

> Glasenapp, M.R., Yee, M.-C., Symons, A.E., Garcia, O.A. & Cornejo, O.E. HLA-Resolve: High-Resolution HLA Haplotyping Using Long-Read Hybrid Capture. *medRxiv* (2026). DOI: [pending]

## License

HLA-Resolve is released under the [UC Santa Cruz Noncommercial License](LICENSE.txt).
