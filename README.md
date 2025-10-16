# HLA Resolve
HLA typing from raw long-read sequencing data (FASTQ or unmapped BAM)

**Authors:** [Matthew Glasenapp](https://github.com/matthewglasenapp), [Alex Symons](https://github.com/FlyingFish800), [Omar Cornejo](https://github.com/oeco28)

Input: A single-sample, demultiplexed, raw long-read sequencing read file. The sequencing platforms currently supported are PacBio and ONT. The tool is compatible with WGS, WES, and targeted sequencing. 

Output(s): HLA star allele calls based on the latest IPD-IMGT/HLA database. Star allele calls are provided for the following genes:
```
HLA-A, HLA-B, HLA-C, HLA-DPA1, HLA-DPB1, HLA-DQA1, HLA-DQB1, HLA-DRB1 
```
We are working to include the other HLA-DRB paralogs in a future release. 

Haplotagged, mapped BAM files for chomosome 6 are provided for visualization with genome browsers such as IGV. Phased VCFs for chomrosome 6 are provided. Reconstructed, haploid (phased) nucleotide sequences are provided for each gene in fasta format. 

Runtime: Depends heavily on the size of the raw sequence reads file and CPU allocation. MHC target capture data should run in < 30min with 6CPU and 20GB RAM. Runtime will take longer with high-coverage WGS and WES, as all reads must be mapped to the human reference genome before restricting downstream analysis to the MHC region of chromosome 6. 

```
usage: cli.py [-h] --input_file INPUT_FILE --sample_name SAMPLE_NAME --platform {pacbio,ont} --scheme
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

Example: python3 hla_resolve/cli.py --input_file reads.bam --sample_name HG002 --platform pacbio --scheme targeted --output_dir out --threads 10

```

Demo example:

Input data: PacBio Revio HiFi targeted sequencing reads from IHW09200, a Southeast Asian sample from Thailand that was part of the 4th Asia-Oceania Histocompatibility Workshop (4AOHW) cell line panel.

```
python3 hla_resolve/cli.py \
--input_file demo/IHW09200.hifi_reads.fastq.gz \
--sample_name IHW09200 \
--platform pacbio \
--scheme targeted \
--output_dir test \
--trim_adapters \
--adapter_file demo/adapters.fasta \
--threads 6
```

The command will print the final star allele calls to STDOUT, along with important logging information, including coverage depth metrics, heterozygous genotypes that could not be phased, and the paths of intermediate files (e.g., BAM, VCF).

Intermediate files will be written to the following dirctory. The user can specify the ```--clean-up``` option if they do not want intermediate files, such as mapped BAM, phased genotypes (VCFs), or fasta haplotype nucleotide sequences for the HLA genes. 


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
