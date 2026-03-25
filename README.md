<br/>

<p align="center">
  <img src="images/hla_resolve.png" alt="HLA-RESOLVE Logo" width="340"/><br/>
  <b>HLA Typing from PacBio Reads</b>
</p>

**Authors:** [Matthew Glasenapp](https://github.com/matthewglasenapp), [Alex Symons](https://github.com/FlyingFish800), [Omar Cornejo](https://github.com/oeco28)

**⚠️ Note:** HLA-Resolve is intended for high-coverage PacBio HiFi reads. ONT support is still in development. The HLA-Resolve manuscript is under peer review. 

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

#### Planned Features (In Development)

1. HLA typing at P-group resolution
2. HLA typing for additional HLA Class I protein-coding genes and pseudogenes  
   (HLA-E, HLA-F, HLA-G; HLA-H, HLA-J, HLA-K, HLA-L, HLA-S, HLA-V, HLA-W)
3. HLA typing for additional HLA Class II protein-coding genes  
   (HLA-DRB3, HLA-DRB4, HLA-DRB5)
4. CYP21A2 star-allele calling

# Installation
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
| hla.xml (IPD-IMGT/HLA database) | ANHIG/IMGTHLA |
| Clair3 Singularity image | Docker Hub |
| DeepVariant Singularity image | Docker Hub |

**Note:** These downloads are large. Ensure sufficient disk space is available in the install directory before the first run.

Note: HLA-Resolve depends on several PacBio command-line utilities that are distributed only as precompiled Linux binaries through the Bioconda channel (pbmarkdup, pbtk, hiphase, trgt). Because these packages are not built or supported for macOS (osx-64 / osx-arm64), HLA-Resolve cannot be fully installed or executed natively on a MacBook or other macOS device.

# Quick Start

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

# Demo

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

#### Technical Reference
For detailed documentation on the algorithms, decision logic, and tools used internally by HLA-Resolve, see the
[Technical Reference](https://github.com/matthewglasenapp/hla_resolve/blob/main/docs/technical_reference.md).
