<br/>

<p align="center">
  <picture>
    <source media="(prefers-color-scheme: dark)" srcset="images/hla_resolve.png">
    <img src="images/hla_resolve_light.png" alt="HLA-RESOLVE Logo" width="340"/>
  </picture>
  <br/>
  <b>HLA Typing from PacBio Reads</b>
</p>

HLA-Resolve is a command-line tool for high-resolution HLA typing from high-coverage PacBio HiFi sequencing reads. It reconstructs phased, full-gene sequences for the eight classical HLA loci (HLA-A, -B, -C, -DPA1, -DPB1, -DQA1, -DQB1, -DRB1) and queries the [IPD-IMGT/HLA database](https://www.ebi.ac.uk/ipd/imgt/hla/) to assign HLA allele calls.

HLA-Resolve was designed for and fully validated on PacBio hybrid-capture libraries (read N50 ~4 kb). It should also work with PacBio whole-genome (WGS), whole-exome (WES), and amplicon data. WGS support has been validated on high-coverage PacBio HiFi libraries from the GIAB and HPRC benchmarks (see [Validated WGS Libraries](#validated-wgs-libraries)).

**Authors:** [Matthew Glasenapp](https://github.com/matthewglasenapp), [Alex Symons](https://github.com/FlyingFish800), [Omar Cornejo](https://github.com/oeco28)

**⚠️ Note:** HLA-Resolve is intended for high-coverage PacBio HiFi reads. ONT support is still in development. The HLA-Resolve [manuscript](https://doi.org/10.64898/2026.03.27.26349549) is under peer review.

> **Disclaimer:** HLA-Resolve is pre-release software in active development. It is intended for research use only and not for use in diagnostic procedures.

## Table of Contents

- [Requirements](#requirements)
- [Overview](#overview)
- [Installation](#installation)
- [Updating](#updating)
- [Quick Start](#quick-start)
- [Demo](#demo)
- [Validated WGS Libraries](#validated-wgs-libraries)
- [Workflow and Dependencies](#workflow-and-dependencies)
- [Planned Features (In Development)](#planned-features-in-development)
- [Citation](#citation)
- [License](#license)

## Requirements

- **Linux (x86_64)** — Several dependencies (pbmarkdup, hiphase, trgt, pbsv) are distributed as precompiled Linux binaries via Bioconda and are not available for macOS.
- **Conda** and **pip** — Used to install all dependencies (see [Installation](#installation)).

## Overview

#### Input
A raw, single-sample (demultiplexed) PacBio sequencing file in FASTQ or unmapped BAM format (compressed or uncompressed). The tool is compatible with WGS, WES, hybrid-capture, and amplicon sequencing schemes.

#### Output(s)

**Primary Results**

HLA allele calls for the following genes:

```
HLA-A, HLA-B, HLA-C, HLA-DPA1, HLA-DPB1, HLA-DQA1, HLA-DQB1, HLA-DRB1
```

**Intermediate Files**
- Haplotagged, mapped BAMs for chromosome 6 (for visualization in genome browsers such as IGV)
- Phased VCFs (chromosome 6 and individual gene)
- Reconstructed haplotype nucleotide sequences for each HLA gene in FASTA format

#### Runtime and Required Resources
Runtime depends heavily on input file size and available compute resources. Targeted HLA capture data typically completes in **<30 minutes** using **6 CPUs and 20 GB RAM**. Runtime increases for high-coverage WGS or WES datasets, as all reads must be mapped to the human reference genome prior to restricting downstream analysis to the HLA region on chromosome 6.

Reference genome alignment is the rate-limiting step and is multithreaded, so increasing the thread count with `--threads` (default **6**) provides the largest runtime reduction, particularly for high-coverage WGS or WES inputs.

## Installation
```text
git clone https://github.com/matthewglasenapp/hla_resolve
cd hla_resolve        # the repository directory created by the clone above
conda env create -f environment.yml
conda activate hla_resolve
pip install -e .
hla_resolve setup
```
``hla_resolve setup`` downloads and builds every required dependency once, up front:

| File | Source |
|------|--------|
| GRCh38 reference genome | NCBI |
| Picard | Broad Institute |
| LongPhase binary | GitHub |
| hla.xml ([IPD-IMGT/HLA database](https://github.com/ANHIG/IMGTHLA)) | IMGTHLA |
| Clair3 Singularity image | Docker Hub |
| DeepVariant Singularity image | Docker Hub |

**Note:** These downloads are large. Ensure sufficient disk space is available in the install directory before the first run.

## Updating
Please ensure you are running the latest version. To update an existing installation to the latest version, run ``update.sh`` from the root of your cloned ``hla_resolve`` repository:
```text
chmod a+x update.sh
bash update.sh
```


## Quick Start

```
usage: hla_resolve [-h] [--version] --input_file INPUT_FILE --sample_name SAMPLE_NAME --platform {pacbio,ont} --scheme {WGS,WES,hybrid_capture,amplicon} --output_dir OUTPUT_DIR
                   [--trim_adapters] [--adapter_file ADAPTER_FILE] [--threads THREADS] [--read_group_string READ_GROUP_STRING] [--clean-up] [--clair3_model CLAIR3_MODEL]

Run HLA-Resolve

options:
  -h, --help            show this help message and exit
  --version             show program's version number and exit
  --input_file INPUT_FILE
                        Path to the raw sequencing reads file (default: None)
  --sample_name SAMPLE_NAME
                        Override the parsed sample name (default: None)
  --platform {pacbio,ont}
                        Specify sequencing platform (pacbio, ont) (default: None)
  --scheme {WGS,WES,hybrid_capture,amplicon}
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
  --clair3_model CLAIR3_MODEL
                        Clair3 model name (bundled in SIF). Defaults to r1041_e82_400bps_sup_v500 for ONT
                        and hifi_revio for PacBio. (default: None)

Example: hla_resolve --input_file reads.bam --sample_name HG002 --platform pacbio --scheme hybrid_capture --output_dir out --threads 10

```

## Demo

Input data: PacBio Revio HiFi hybrid-capture sequencing reads from HG002 (Ashkenazi Son), a sample from the GIAB and HPRC benchmarks. Run from the repository root:

```text
hla_resolve \
  --input_file demo/HG002.hifi_reads.fastq.gz \
  --sample_name HG002 \
  --platform pacbio \
  --scheme hybrid_capture \
  --output_dir test \
  --trim_adapters \
  --adapter_file demo/adapters.fasta \
  --threads 6
```

The command will print the final HLA allele calls to STDOUT, along with important logging information, including coverage depth metrics, heterozygous genotypes that could not be phased, and the paths of intermediate files (e.g., BAM, VCF).

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
| `haploblocks/`            | Contains a list of fully-phased HLA genes                                  |
| `filtered_vcf/`           | Contains the final, filtered VCF of variants to be applied during fasta haplotype reconstruction |
| `vcf2fasta_out/`          | Contains the raw sequence output from vcf2fasta                            |
| `hla_fasta_haplotypes/`   | Contains fasta files of full gene and CDS sequences for each HLA gene       |
| `hla_typing_results/`     | Contains the final results of HLA typing                                   |

## Validated WGS Libraries

HLA-Resolve has produced high-quality calls for the following whole-genome (WGS) PacBio HiFi libraries:

| Sample | Source | Instrument | HLA&nbsp;Coverage | Concordance | File |
|--------|--------|------------|--------------|-------------|------|
| HG002 | GIAB | Revio | ~30× | 1–3&nbsp;field:&nbsp;100%<br>4-field:&nbsp;7/8 | `https://ftp-trace.ncbi.nlm.nih.gov/ReferenceSamples/giab/data/AshkenazimTrio/HG002_NA24385_son/PacBio_HiFi-Revio_20231031/HG002_PacBio-Revio_m84039_231005_222902_s1.hifi_reads.bam` |
| HG002 | GIAB | Revio | ~30× | 1–3&nbsp;field:&nbsp;100%<br>4-field:&nbsp;100% | `https://ftp-trace.ncbi.nlm.nih.gov/ReferenceSamples/giab/data/AshkenazimTrio/HG002_NA24385_son/PacBio_HiFi-Revio_20231031/HG002_PacBio-Revio_m84039_230928_213653_s3.hifi_reads.bam` |
| HG002 | HPRC | Revio | ~33× | 1–3&nbsp;field:&nbsp;100%<br>4-field:&nbsp;100% | `s3://human-pangenomics/working/HPRC_PLUS/HG002/raw_data/PacBio_HiFi/wMods/m84011_220902_175841_s1.hifi_reads.bam` |
| HG01258 | HPRC | Revio | ~19× | 1–3&nbsp;field:&nbsp;100%<br>4-field:&nbsp;7/8 | `s3://human-pangenomics/working/HPRC/HG01258/raw_data/PacBio_HiFi/m84046_231202_090949_s3.hifi_reads.bc2054.bam` |
| HG03579 | HPRC | Sequel II | ~15× | 1–3&nbsp;field:&nbsp;100%<br>4-field:&nbsp;7/8 | `s3://human-pangenomics/working/HPRC/HG03579/raw_data/PacBio_HiFi/m64043_200516_230634.ccs.bam` |

HPRC libraries can be downloaded without credentials using the AWS CLI:

`aws s3 cp --no-sign-request <s3-path> .`

**Note:** Concordance was evaluated against ground-truth HLA annotations provided by Lai et al. 2023 ([DOI: 10.1016/j.csbj.2024.03.030](https://doi.org/10.1016/j.csbj.2024.03.030); [Supplementary File 6](docs/Lai_Supplementary-6.xlsx)). All 4-field discordances are single-field miscalls in the fourth field; 1–3 field concordance is 100% across all libraries.

The `hla_resolve` command used to analyze these libraries was:

```bash
hla_resolve --input_file <INPUT_uBAM> --sample_name <sample_name> --platform pacbio --scheme WGS --output_dir <output_dir> --threads <threads>
```

## Workflow and Dependencies

HLA-Resolve takes raw PacBio HiFi reads (FASTQ or uBAM) as input and executes the following steps to produce four-field HLA allele assignments.

#### 1. Adapter Trimming
Adapter and barcode sequences are removed from raw reads using [cutadapt](https://doi.org/10.14806/ej.17.1.200) (when adapter sequences are provided) or [fastplong](https://doi.org/10.1002/imt2.107) (auto-detection mode).

#### 2. PCR Duplicate Removal
PCR duplicates are identified and removed from the trimmed reads using [pbmarkdup](https://github.com/PacificBiosciences/pbmarkdup).

#### 3. Reference Genome Alignment
Deduplicated reads are aligned to a modified GRCh38 reference genome (no-alt analysis set) using [rammap](https://doi.org/10.64898/2026.05.26.726289) (Wang and Li, 2026), a memory-safe Rust reimplementation of [minimap2](https://doi.org/10.1093/bioinformatics/bty191) that produces identical alignments. The modified reference includes an additional scaffold containing the HLA-Y/HLA-OLI insertion to prevent mismapping of HLA-Y reads to HLA-A. 

#### 4. HLA-DRB Paralog Filtering
A separate alignment step maps reads against a multi-allele HLA-DRB reference (`DRB_reference.fa`) continaing 13 HLA-DRB1 alleles, 3 HLA-DRB3 alleles, and 1 HLA-DRB4 allele from IPD-IMGT/HLA, as well as the HLA-DRB5, -DRB6, and -DRB9 sequences from GRCh38. Reads with primary alignments to anything other than an HLA-DRB1 are flagged for removal. 

#### 5. Read Filtering
Aligned reads are filtered to retain only primary alignments on chromosome 6. Reads with primary alignments to HLA-DRB1 paralogs (identified in step 4) are removed. Additionally, reads with primary alignments upstream of HLA-DRB1 (chr6:32439878-32572902), the region containing HLA-DRA, -DRB9, -DRB5, and -DRB6 are removed to prevent spurious SV calls by pbsv based on split alignments across the DRB paralogs.

#### 6. Small Variant Calling
SNVs are called with [bcftools](https://doi.org/10.1093/bioinformatics/btr509) and indels are called with [DeepVariant](https://doi.org/10.1038/nbt.4235). DeepVariant RefCall genotypes with sufficient read support are rescued and reclassified as heterozygous or homozygous ALT based on variant allele frequency.

#### 7. Structural Variant Calling
Structural variants are called from the aligned reads using [pbsv](https://github.com/PacificBiosciences/pbsv).

#### 8. Tandem Repeat Genotyping
Tandem repeats within the HLA region are genotyped using [TRGT](https://doi.org/10.1038/s41587-023-02057-3).

#### 9. Joint Phasing
Small variants, structural variants, and tandem repeat genotypes are jointly phased with [HiPhase](https://doi.org/10.1093/bioinformatics/btae042), producing haplotagged BAMs, phased VCFs, and haplotype block coordinates.

#### 10. Coverage Assessment
Per-gene coverage depth and breadth are calculated with [mosdepth](https://doi.org/10.1093/bioinformatics/btx699). Genes failing minimum coverage thresholds are excluded from HLA typing.

#### 11. Haploblock Evaluation
Phased haplotype blocks are evaluated to determine whether each HLA gene is fully spanned by a single phase set. Genes with internal phasing breaks enter a rescue pipeline that attempts to recover coding sequence or antigen recognition sequence (ARS) haplotypes.

#### 12. Variant Filtering and Redundancy Removal
Phased genotypes are filtered by gene to remove redundant calls from overlapping variant callers (e.g., DeepVariant indels overlapping pbsv structural variants, or non-TRGT variants within TRGT tandem repeat regions). Symbolic and complex structural variant types (BND, INV, DUP) are excluded.

#### 13. Haplotype Reconstruction
Phased, filtered genotypes are applied to GRCh38 gene models using [vcf2fasta](https://github.com/santiagosnchez/vcf2fasta) to reconstruct full-gene and coding-sequence haplotype FASTA files for each HLA gene.

#### 14. IPD-IMGT/HLA Database Matching
Reconstructed haplotypes are compared against alleles in the [IPD-IMGT/HLA database](https://doi.org/10.1093/nar/gkac1011) with a three-pass hierarchical classification algorithm using [edlib](https://doi.org/10.1093/bioinformatics/btw753):

   1. **G-group assignment** — The antigen recognition sequence (ARS exons) is matched to G-group reference sequences by edit distance. An exact match is required.
   2. **Three-field allele assignment** — The full concatenated exon sequence is compared against alleles within the assigned G group, ranked by edit distance.
   3. **Four-field refinement** — The full-gene haplotype (including introns and UTRs) is compared against candidate alleles, ranked by mismatch identity (the proportion of matching bases at 1:1-aligned positions), which avoids penalizing insertions and deletions from unreliable intronic reconstruction. Ties are broken by match length, then by lowest fourth-field value.
   4. **DR/DQ re-consensus refinement** — For HLA-DQA1, HLA-DQB1, and HLA-DRB1, the fourth-field call is re-derived by re-mapping the underlying reads to the best-guess allele, rebuilding a consensus, and re-matching within the same three-field group.

#### Note
For WGS and WES input, the pipeline skips adapter trimming (step 1) and pre-alignment duplicate removal (step 2).

## Planned Features (In Development)

1. HLA typing at P-group resolution
2. HLA typing for additional HLA Class I protein-coding genes and pseudogenes
   (HLA-E, HLA-F, HLA-G; HLA-H, HLA-J, HLA-K, HLA-L, HLA-S, HLA-V, HLA-W)
3. HLA typing for additional HLA Class II protein-coding genes
   (HLA-DRB3, HLA-DRB4, HLA-DRB5)

#### Technical Reference
For detailed documentation on the algorithms, decision logic, and tools used internally by HLA-Resolve, see the
[Technical Reference](https://github.com/matthewglasenapp/hla_resolve/blob/main/docs/technical_reference.md).

## Citation

If you use HLA-Resolve, please cite:

> Glasenapp, M.R., Yee, M.-C., Symons, A.E., Cornejo, O.E. & Garcia, O.A. HLA-Resolve: High-Resolution HLA Haplotyping Using Long-Read Hybrid Capture. *medRxiv* (2026). https://doi.org/10.64898/2026.03.27.26349549

## License

HLA-Resolve is released under the [UC Santa Cruz Noncommercial License](LICENSE.txt).
