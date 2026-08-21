<br/>

<p align="center">
  <picture>
    <source media="(prefers-color-scheme: dark)" srcset="images/hla_resolve.png">
    <img src="images/hla_resolve_light.png" alt="HLA-RESOLVE Logo" width="340"/>
  </picture>
  <br/>
  <b>HLA Typing from PacBio Reads</b>
</p>

<p align="center">
  <img src="https://img.shields.io/github/v/release/matthewglasenapp/hla_resolve?label=version&color=blue" alt="Version">
  <img src="https://img.shields.io/badge/platform-linux--64-lightgrey" alt="Platform">
  <a href="LICENSE.txt"><img src="https://img.shields.io/badge/license-UCSC%20Noncommercial-green" alt="License"></a>
  <a href="https://doi.org/10.64898/2026.03.27.26349549"><img src="https://img.shields.io/badge/medRxiv-10.64898-b31b1b" alt="Preprint"></a>
</p>

<br/>

<p align="center">
  <a href="#installation"><picture><source media="(prefers-color-scheme: dark)" srcset="https://img.shields.io/badge/Installation-1F2937?style=for-the-badge&logo=anaconda&logoColor=white"><img src="https://img.shields.io/badge/Installation-E5E7EB?style=for-the-badge&logo=anaconda&logoColor=1F2937" alt="Installation"></picture></a>
  &nbsp;&nbsp;&nbsp;&nbsp;
  <a href="#quick-start-and-demo"><picture><source media="(prefers-color-scheme: dark)" srcset="https://img.shields.io/badge/Quick_Start-1F2937?style=for-the-badge&logo=gnubash&logoColor=white"><img src="https://img.shields.io/badge/Quick_Start-E5E7EB?style=for-the-badge&logo=gnubash&logoColor=1F2937" alt="Quick Start"></picture></a>
  &nbsp;&nbsp;&nbsp;&nbsp;
  <a href="#technical-reference"><picture><source media="(prefers-color-scheme: dark)" srcset="https://img.shields.io/badge/Technical_Reference-1F2937?style=for-the-badge&logo=readthedocs&logoColor=white"><img src="https://img.shields.io/badge/Technical_Reference-E5E7EB?style=for-the-badge&logo=readthedocs&logoColor=1F2937" alt="Technical Reference"></picture></a>
  &nbsp;&nbsp;&nbsp;&nbsp;
  <a href="#benchmarks"><picture><source media="(prefers-color-scheme: dark)" srcset="https://img.shields.io/badge/Benchmark-1F2937?style=for-the-badge&logo=databricks&logoColor=white"><img src="https://img.shields.io/badge/Benchmark-E5E7EB?style=for-the-badge&logo=databricks&logoColor=1F2937" alt="Benchmark"></picture></a>
  &nbsp;&nbsp;&nbsp;&nbsp;
  <a href="#citation"><picture><source media="(prefers-color-scheme: dark)" srcset="https://img.shields.io/badge/Citation-1F2937?style=for-the-badge&logo=googlescholar&logoColor=white"><img src="https://img.shields.io/badge/Citation-E5E7EB?style=for-the-badge&logo=googlescholar&logoColor=1F2937" alt="Citation"></picture></a>
</p>

<h4 align="center">
  <a href="https://github.com/matthewglasenapp">Matthew Glasenapp</a> &nbsp;·&nbsp;
  <a href="https://github.com/FlyingFish800">Alex Symons</a> &nbsp;·&nbsp;
  <a href="https://github.com/oeco28">Omar Cornejo</a>
</h4>

<br/>

## Introduction

HLA-Resolve is a command-line tool for high-resolution HLA typing from high-coverage PacBio sequencing reads. It reconstructs phased, coding and full-gene sequences for the eight classical HLA loci (HLA-A, -B, -C, -DPA1, -DPB1, -DQA1, -DQB1, -DRB1) and queries the [IPD-IMGT/HLA database](https://www.ebi.ac.uk/ipd/imgt/hla/) to assign HLA allele calls at four different levels of resolution (G group, P group, three-field, and four-field)

HLA-Resolve was designed for and fully validated on PacBio hybrid-capture libraries (read N50 ~4 kb). WGS support has been validated on PacBio whole-genome sequencing reads from the GIAB and HPRC benchmarks (see [Benchmarks](#benchmarks)). HLA-Resolve has not been tested on amplicon sequencing data yet. 

> [!IMPORTANT]
> 1. HLA-Resolve is pre-release software in active development. It is intended for high-coverage PacBio reads. A gene is typed only if its peptide-binding domain reaches at least 8× mean coverage depth.
> 2. ONT support is still in development, and `--platform ont` is rejected at runtime until it lands.
> 3. The software is for research use only and not for use in diagnostic procedures. The HLA-Resolve [manuscript](https://doi.org/10.64898/2026.03.27.26349549) is under peer review.

<br/>

<details>
<summary><b>Table of Contents</b></summary>

- [Introduction](#introduction)
- [Tool Overview](#tool-overview)
  - [Input](#input)
  - [Outputs](#outputs)
  - [Runtime and Required Resources](#runtime-and-required-resources)
- [Installation](#installation)
- [Quick Start and Demo](#quick-start-and-demo)
  - [Minimal Command](#minimal-command)
  - [Demo](#demo)
- [Technical Reference](#technical-reference)
- [Benchmarks](#benchmarks)
  - [Hybrid Capture](#hybrid-capture)
  - [Whole Genome Sequencing](#whole-genome-sequencing)
- [Planned Features (In Development)](#planned-features-in-development)
- [Citation](#citation)
- [Support](#support)
- [License](#license)

</details>

---

## Tool Overview

### Input

<table>
<tr><td><b>Reads</b></td><td>A raw, single-sample (demultiplexed) PacBio sequencing file</td></tr>
<tr><td><b>Format</b></td><td>FASTQ or unmapped BAM, compressed or uncompressed</td></tr>
</table>

### Outputs

#### Primary Results

HLA allele calls for the following eight genes:

<table>
<tr><td><b>Class I</b></td><td>HLA-A, HLA-B, HLA-C</td></tr>
<tr><td><b>Class II</b></td><td>HLA-DPA1, HLA-DPB1, HLA-DQA1, HLA-DQB1, HLA-DRB1</td></tr>
</table>

**Example Output**

<table>
<tr>
<th>sample</th>
<th>HLA-A_1</th>
<th>HLA-A_2</th>
<th>HLA-B_1</th>
<th>HLA-B_2</th>
<th>...</th>
</tr>
<tr>
<td>HG002</td>
<td>HLA-A&#42;01:01:01:01</td>
<td>HLA-A&#42;26:01:01:01</td>
<td>HLA-B&#42;38:01:01:01</td>
<td>HLA-B&#42;35:08:01:01</td>
<td>...</td>
</tr>
</table>

The calls are written to `<output_dir>/<sample>/hla_typing_results/`, one row per sample.

> [!NOTE]
> Allele order within each gene is arbitrary and is not consistent between genes.

Four output files hold the same calls at different levels of resolution.

| Resolution | Best Guess | Ambiguities Reported |
|------------|------------|----------------------|
| P group | `p_group_output.csv` | `p_group_output_full.csv` |
| G group | `g_group_output.csv` | `g_group_output_full.csv` |
| Three field | `3_field_allele_output.csv` | `3_field_allele_output_full.csv` |
| Four field | `allele_output.csv` | `allele_output_full.csv` |

The Best Guess file forces a single best guess for every allele. The Ambiguities Reported file gives the calls as genotype list strings when multiple candidate alleles cannot be distinguished at the sequence level.

P and G groups are read from the reconstructed sequence rather than looked up from the allele call. A P group is identity of the protein encoded by the peptide-binding domain, exon 2 for class II and exons 2 and 3 for class I, excluding the codons split across an exon border. A G group is identity of the same region at the nucleotide level. Where neither is observed, the three-field call is reported instead, so a name ending in P or G is a group and anything else is an allele.

#### Supporting Files

- A haplotagged BAM of the reads mapped to the MHC region of chromosome 6
- Phased VCFs, both chromosome 6 and per gene
- Nucleotide sequences (FASTA) for each HLA gene

### Runtime and Required Resources

Runtime depends heavily on input file size and available compute resources. Target capture data typically completes in **<15 minutes** using **6 CPUs and 20 GB RAM**. Runtime increases for high-coverage WGS datasets, as all reads must be mapped to the human reference genome prior to restricting downstream analysis to the HLA region on chromosome 6.

Reference genome alignment is the rate-limiting step and is multithreaded, so increasing the thread count with `--threads` provides the largest runtime reduction, especially for high-coverage WGS inputs. The default number of threads is **6**, or the number of CPUs available if fewer, read from the Slurm allocation where there is one. Memory requirements rise with the thread count, so raise the job's memory alongside `--threads`. `--threads` sets the thread count for the tools HLA-Resolve calls directly. DeepVariant also parallelizes internally, so on a cluster its CPU use is bounded by the job allocation and not by this flag.

---

## Installation

**Requirements**

- **Linux (x86_64)** — Several dependencies (pbmarkdup, hiphase, trgt, pbsv) are distributed as precompiled Linux binaries via Bioconda and are not available for macOS.
- **Conda** and **pip** — Used to install all dependencies.

**Install**

```bash
git clone https://github.com/matthewglasenapp/hla_resolve
cd hla_resolve        # the repository directory created by the clone above
conda env create -f environment.yml
conda activate hla_resolve
pip install -e .
hla_resolve setup
```
`hla_resolve setup` downloads and builds every required dependency once, up front:

| File | Version | Source |
|------|---------|--------|
| GRCh38 reference genome (no-alt analysis set) | GCA_000001405.15 | NCBI |
| rammap binary | v1.0.0 | GitHub |
| hla.xml ([IPD-IMGT/HLA database](https://github.com/ANHIG/IMGTHLA)) | 3.64.0 | IMGTHLA |
| DeepVariant Singularity image | 1.6.1 | Docker Hub |

> [!NOTE]
> Setup needs about **9 GB** free in the install directory while it runs, and leaves about **6 GB** in place once it finishes. The Singularity image pull also needs temporary space, so set `SINGULARITY_TMPDIR` if `/tmp` is small on your nodes.

**Updating an existing installation**

Please ensure you are running the latest version. To update, run `update.sh` from the root of your cloned `hla_resolve` repository:

```bash
chmod a+x update.sh
bash update.sh
```

---

## Quick Start and Demo

### Minimal Command

```bash
hla_resolve \
  --input_file <reads.fastq.gz | reads.bam> \
  --sample_name <sample_name> \
  --platform pacbio \
  --scheme <WGS | WES | hybrid_capture | amplicon> \
  --output_dir <output_dir>
```

### Demo

The repository includes a demo dataset of PacBio HLA hybrid capture sequencing reads from HG002 (Ashkenazi Son), a sample from the GIAB and HPRC benchmarks. Run this from the repository root:

<table>
<tr>
<th align="left">Command</th>
<th align="left">Example Output</th>
</tr>
<tr>
<td valign="top">

```bash
hla_resolve \
  --input_file demo/HG002.hifi_reads.fastq.gz \
  --sample_name HG002 \
  --platform pacbio \
  --scheme hybrid_capture \
  --output_dir demo_out \
  --trim_adapters \
  --adapter_file demo/adapters.fasta \
  --threads 6
```

</td>
<td valign="top">

```text
Sample: HG002

Four-field resolution
gene       _1                 _2
HLA-A      A*01:01:01:01      A*26:01:01:01
HLA-B      B*38:01:01:01      B*35:08:01:01
HLA-C      C*04:01:01:06      C*12:03:01:01
HLA-DPA1   DPA1*01:03:01:02   DPA1*01:03:01:04
HLA-DPB1   DPB1*04:01:01:01   DPB1*04:01:01:03
HLA-DQA1   DQA1*03:01:01:01   DQA1*01:05:01:01
HLA-DQB1   DQB1*05:01:01:05   DQB1*03:02:01:01
HLA-DRB1   DRB1*04:02:01      DRB1*10:01:01:03

Note: Allele order within each gene is
arbitrary and is not consistent
between genes.
Note: Ambiguous calls are written as
genotype list strings to the
ambiguities files below.

Finished HG002 in 9m 18s (status: ok)
```

</td>
</tr>
</table>

The command prints the HLA allele calls at all four resolutions to STDOUT, together with the paths of the main result files and logging information such as coverage depth metrics. Genes that could not be reconstructed are shown as `not_typed`. The paths of intermediate files go to the log file only. Pass `--verbose` to see them on screen as well.

The same HLA allele calls are written to `demo_out/HG002/hla_typing_results/`, the primary results directory for this run. It holds the eight result files described in [Primary Results](#primary-results), with `allele_output.csv` giving the four-field calls shown above.

The reconstructed sequences that produced those calls are in `demo_out/HG002/hla_fasta_haplotypes/`. `HG002_HLA_haplotypes_gene.fasta` holds the full gene sequences and `HG002_HLA_haplotypes_CDS.fasta` holds the coding sequences (CDS), two records per gene, one for each haplotype.

<details>
<summary><b>Full command-line options</b></summary>

```
usage: hla_resolve [-h] [--version] --input_file INPUT_FILE --sample_name
                   SAMPLE_NAME --platform {pacbio,ont} --scheme
                   {WGS,WES,hybrid_capture,amplicon} --output_dir OUTPUT_DIR
                   [--trim_adapters] [--adapter_file ADAPTER_FILE]
                   [--threads THREADS] [--read_group_string READ_GROUP_STRING]
                   [--clean_up] [--keep_all_intermediates] [--keep_full_bam]
                   [--clair3_model CLAIR3_MODEL] [--verbose] [--quiet]

Run HLA-Resolve

optional arguments:
  -h, --help            show this help message and exit
  --version             show program's version number and exit
  --input_file INPUT_FILE
                        Path to the raw sequencing reads file (default: None)
  --sample_name SAMPLE_NAME
                        Name for this sample. Used for output filenames and the
                        read group (default: None)
  --platform {pacbio,ont}
                        Sequencing platform. Only pacbio is supported; ont is
                        not yet available (default: None)
  --scheme {WGS,WES,hybrid_capture,amplicon}
                        Sequencing scheme (default: None)
  --output_dir OUTPUT_DIR
                        Output directory. Results are written to
                        <output_dir>/<sample_name>/ (default: None)
  --trim_adapters       Enable adapter trimming before processing (default:
                        False)
  --adapter_file ADAPTER_FILE
                        Path to a file with custom adapter sequences
                        (FASTA/FASTQ). If not provided, fastplong auto-
                        detection will be used. (default: None)
  --threads THREADS     Number of threads to use, lowered to the CPU count when
                        fewer are available (default: 6)
  --read_group_string READ_GROUP_STRING
                        Override the parsed read group string (default: None)
  --clean_up            Keep only the HLA typing results and the run log.
                        Everything else HLA-Resolve wrote, including the
                        haplotagged BAM, the VCFs, and the FASTA haplotypes, is
                        removed at the end of the run. Files it did not write
                        are reported and left in place (default: False)
  --keep_all_intermediates
                        Retain all intermediate files. Only use for debugging
                        (default: False)
  --keep_full_bam       Keep the genome-wide mapped BAM. It is deleted by
                        default once reads are filtered to the MHC (default:
                        False)
  --clair3_model CLAIR3_MODEL
                        Clair3 model name (bundled in SIF). Defaults to
                        r1041_e82_400bps_sup_v500 for ONT and hifi_revio for
                        PacBio. (default: None)
  --verbose             Print intermediate file paths and detailed per-variant
                        diagnostic output (overlap suppression, RefCall rescue,
                        unphased het records, CDS sanity check) (default:
                        False)
  --quiet               Print only stage headers, warnings, the final results
                        tables, and the output file paths. The full log is
                        still written to the log file (default: False)

One-time setup (downloads references, binaries, and images):
  hla_resolve setup

Example run:
  hla_resolve --input_file reads.bam --sample_name HG002 --platform pacbio --scheme hybrid_capture --output_dir out

HLA-Resolve is pre-release software intended for research use only
and not for use in diagnostic procedures.
```

</details>

<details>
<summary><b>Output Directories</b></summary>

A run keeps the directories below. Working files that a later stage supersedes are removed as the run goes. `--keep_all_intermediates` retains them.

| Directory | Description |
|-----------|-------------|
| `fastq_raw/`            | Read quality reports from fastplong. Written for a fastq input |
| `mapped_bam/`           | The haplotagged BAM of reads mapped to the MHC region, and the list of HLA-DRB paralog reads removed before variant calling |
| `phased_vcf/`           | The merged VCF of small variants, structural variants, and tandem repeat genotypes following joint phasing, including the phase block report |
| `filtered_vcf/`         | The single-gene phased VCFs used in gene sequence reconstruction, and any unphased heterozygous variants that could not be applied |
| `hla_fasta_haplotypes/` | Fasta files of full gene and CDS sequences for each HLA gene |
| `hla_typing_results/`   | The final results of HLA typing |

</details>

---

## Technical Reference

The [Technical Reference](https://github.com/matthewglasenapp/hla_resolve/blob/main/docs/technical_reference.md) gives the full workflow, the dependencies used at each step, and detailed documentation on the algorithms and decision logic used internally by HLA-Resolve.

---

## Benchmarks

### Hybrid Capture

HLA-Resolve was run on PacBio hybrid capture reads for 31 samples, at a mean coverage of 365× across the eight classical HLA genes. Every gene in every sample was typed (496/496 alleles). Concordance was measured for the 27 samples that carry a reference typing, 15 from the International Histocompatibility Working Group (IHWG) and 12 from the Human Pangenome Reference Consortium (HPRC). The HPRC reference typings are from [Lai et al. 2024](https://doi.org/10.1016/j.csbj.2024.03.030).

| Resolution | IHWG | HPRC | Combined |
|------------|------|------|----------|
| One field | 98.7% (227/230) | 100% (192/192) | 99.3% (419/422) |
| Two field | 98.6% (219/222) | 100% (192/192) | 99.3% (411/414) |
| Three field | 99.0% (196/198) | 100% (190/190) | 99.5% (386/388) |
| Four field | 83.3% (135/162) | 96.8% (180/186) | 90.5% (315/348) |

**[Browse the full benchmark →](docs/capture_validation.md)**

- Concordance by gene and field of resolution
- Mean HLA coverage depth by sample
- The SRA accession for each sample

Raw reads are available under BioProject [PRJNA1417783](https://www.ncbi.nlm.nih.gov/bioproject/PRJNA1417783).

### Whole Genome Sequencing

HLA-Resolve was run on whole-genome PacBio sequencing reads for 39 HPRC samples, at a mean coverage of 35.2× across the eight classical HLA genes. Concordance with the reference typings of Lai et al. was **100% at one- through three-field resolution** (610/610 alleles) and **92.8% at four-field resolution** (555/598), with a call rate of 99.4% (620/624).

**[Browse the full benchmark →](docs/wgs_validation.md)**

- Per-sample concordance by field of resolution
- Mean HLA coverage depth by sample
- The raw PacBio input file(s) used for each sample

All inputs are publicly available, so the benchmark can be reproduced end to end.

---

## Planned Features (In Development)

| | Feature |
|---|---------|
| 1 | HLA typing at P-group resolution |
| 2 | HLA typing for additional HLA Class I protein-coding genes and pseudogenes<br>(HLA-E, HLA-F, HLA-G; HLA-H, HLA-J, HLA-K, HLA-L, HLA-S, HLA-V, HLA-W) |
| 3 | HLA typing for additional HLA Class II protein-coding genes<br>(HLA-DRB3, HLA-DRB4, HLA-DRB5) |

---

## Citation

If you use HLA-Resolve, please cite:

> Glasenapp, M.R., Yee, M.-C., Symons, A.E., Cornejo, O.E. & Garcia, O.A. HLA-Resolve: High-Resolution HLA Haplotyping Using Long-Read Hybrid Capture. *medRxiv* (2026). https://doi.org/10.64898/2026.03.27.26349549

---

## Support

Questions, bug reports, and feature requests are welcome. Please [open an issue](https://github.com/matthewglasenapp/hla_resolve/issues) on GitHub.

---

## License

HLA-Resolve is released under the [UC Santa Cruz Noncommercial License](LICENSE.txt).
