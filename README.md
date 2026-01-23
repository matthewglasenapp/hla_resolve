<br/>

<p align="center">
  <img src="images/hla_resolve.png" alt="HLA-RESOLVE Logo" width="340"/><br/>
  <b>HLA Typing from Long Reads</b>
</p>

**Authors:** [Matthew Glasenapp](https://github.com/matthewglasenapp), [Alex Symons](https://github.com/FlyingFish800), [Omar Cornejo](https://github.com/oeco28)

Input: A raw, single-sample (demultiplexed) long-read sequencing read file in FASTQ or unmapped BAM format. The program automatically detects the input file format and can handle both compressed and uncompressed input. The sequencing platforms currently supported are PacBio and ONT. The tool is compatible with WGS, WES, and targeted sequencing schemes.  

Output: HLA star allele calls based on the latest IPD-IMGT/HLA database. Star allele calls are provided for the following genes:  

```
HLA-A, HLA-B, HLA-C, HLA-DPA1, HLA-DPB1, HLA-DQA1, HLA-DQB1, HLA-DRB1 
```
Haplotagged, mapped BAM files for chomosome 6 are provided for visualization with genome browsers such as IGV. Phased VCFs for chomrosome 6 are provided. Reconstructed, haploid (phased) nucleotide sequences are provided for each gene in fasta format. 

Runtime: Depends heavily on the size of the raw sequence reads file and CPU allocation. MHC target capture data should run in < 30min with 6CPU and 20GB RAM. Runtime will take longer with high-coverage WGS and WES, as all reads must be mapped to the human reference genome before restricting downstream analysis to the MHC region of chromosome 6. 

### Planned Features (In Development)

1. HLA typing at P-group resolution
2. HLA typing for additional HLA Class I protein-coding genes and pseudogenes  
   (HLA-E, HLA-F, HLA-G; HLA-H, HLA-J, HLA-K, HLA-L, HLA-S, HLA-V, HLA-W)
3. HLA typing for additional HLA Class II protein-coding genes  
   (HLA-DRB3, HLA-DRB4, HLA-DRB5)
4. CYP21A2 star-allele calling
5. Option to use DeepVariant or Clair3 for genotyping


# Installation
```
git clone https://github.com/matthewglasenapp/hla_resolve
conda env create -f hla_resolve/environment.yml
conda activate hla_resolve
pip install -e hla_resolve
hla_resolve --help
```
The first trime ``hla_resolve`` is executed, it will download the required reference genome from NCBI, picard.jar from the Broad Institute, and hla.xml from the IPD-IMGT/HLA database. 

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

Input data: PacBio Revio HiFi targeted sequencing reads from IHW09200, a Southeast Asian sample from Thailand that was part of the 4th Asia-Oceania Histocompatibility Workshop (4AOHW) cell line panel.

```
hla_resolve \
--input_file hla_resolve/demo/IHW09200.hifi_reads.fastq.gz \
--sample_name IHW09200 \
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
| `structural_variant_vcf/` | Contains the SV genotype calls from either Sniffles (ONT) or sawfish (PacBio) |
| `pbtrgt_vcf/`             | Contains the tandem repeat genotypes from TRGT (PacBio-only)               |
| `phased_vcf/`             | Contains phased genotype calls from joint phasing of small variants, structural variants, and tandem repeat genotypes (PacBio only) |
| `mosdepth/`               | Contains coverage depth output files from mosdepth for the HLA genes        |
| `haploblocks/`            | Contains a list of fully-phased MHC genes                                  |
| `filtered_vcf/`           | Contains the final, filtered VCF of variants to be applied during fasta haplotype reconstruction |
| `vcf2fasta_out/`          | Contains the raw sequence output from vcf2fasta                            |
| `hla_fasta_haplotypes/`   | Contains fasta files of full gene and CDS sequences for each HLA gene       |
| `hla_typing_results/`     | Contains the final results of HLA typing                                   |
