# WGS Validation

HLA-Resolve v0.9.8 (IPD-IMGT/HLA release 3.64.0) was run on public whole-genome PacBio reads for **143 sequencing libraries from 129 samples**. All libraries were sequenced by the Human Pangenome Reference Consortium (HPRC; cohorts HPRC, HPRC_PLUS) except the two Genome in a Bottle (GIAB) HG002 libraries. Every library was typed independently. Libraries sequenced across multiple SMRT cells were combined to a target genome-wide coverage depth of 30×. Samples with two libraries appear twice as independent measurements.

Allele calls were compared to two reference sets:

- **Lai et al. 2024** ([doi.org/10.1016/j.csbj.2024.03.030](https://doi.org/10.1016/j.csbj.2024.03.030); [Supplementary File 6](Lai_Supplementary-6.xlsx))
  - 42 HPRC release 1 samples, 49 libraries
  - HLA-A, -B, -C, -DPA1, -DPB1, -DQA1, -DQB1, -DRB1
  - Resolution: up through four-field (full) resolution
- **1000 Genomes Project HLA panel** ([20181129_HLA_types_full_1000_Genomes_Project_panel.txt](https://ftp.1000genomes.ebi.ac.uk/vol1/ftp/data_collections/HLA_types/20181129_HLA_types_full_1000_Genomes_Project_panel.txt))
  - 87 samples, 94 libraries
  - HLA-A, -B, -C, -DQB1, -DRB1
  - Resolution: two-field resolution
  - Exome-based (IPD-IMGT/HLA 3.28)

The following command was used for each sample:

```bash
hla_resolve --input_file <LIBRARY_uBAM> --sample_name <SAMPLE> \
            --platform pacbio --scheme WGS \
            --output_dir <OUTPUT_DIR> --threads <THREADS>
```

## Scoring

Concordance is reported among alleles called. An allele was evaluated at a given field resolution only if the reference specified it to that resolution. A call resolved to fewer fields than the reference is considered discordant at that level of resolution. By default, genes do not type if their antigen recognition site coverage depth sits below 8x, so call rate was measured as the proportion of total possible allele calls emitted by HLA-Resolve.

## Lai et al. reference (49 libraries, 42 samples)

| Metric | Value |
|---|---:|
| Total possible allele calls | 784 |
| Alleles called | 784 |
| Call rate | 100.0% |

| Resolution | Concordance | Concordant alleles | Alleles evaluated |
|---|---:|---:|---:|
| 1-field | 99.9% | 782 | 783 |
| 2-field | 99.7% | 781 | 783 |
| 3-field | 99.5% | 768 | 772 |
| 4-field | 92.7% | 701 | 756 |

Alleles evaluated is lower than alleles called because 1 uncertain reference allele was excluded (HG01358 DRB1*04:92#), and not all reference alleles have third and fourth field values.

## 1000 Genomes panel (94 libraries, 87 samples)

| Metric | Value |
|---|---:|
| Total possible allele calls | 940 |
| Alleles called | 940 |
| Call rate | 100.0% |

| Resolution | Concordance | Concordant alleles | Alleles evaluated |
|---|---:|---:|---:|
| 1-field | 99.9% | 911 | 912 |
| 2-field | 99.0% | 903 | 912 |

Alleles evaluated is lower than alleles called because the reference has missing HLA-DQB1 calls for 14 samples (28 alleles).

## All libraries

Both reference sets pooled, 143 libraries. The number of alleles evaluated falls with resolution because the 1000 Genomes panel stops at two fields and covers five genes, and because a reference allele is evaluated only at the fields it specifies.

| Metric | Value |
|---|---:|
| Total possible allele calls | 1724 |
| Alleles called | 1724 |
| Call rate | 100.0% |

| Resolution | Concordance | Concordant alleles | Alleles evaluated |
|---|---:|---:|---:|
| 1-field | 99.9% | 1693 | 1695 |
| 2-field | 99.4% | 1684 | 1695 |
| 3-field | 99.5% | 768 | 772 |
| 4-field | 92.7% | 701 | 756 |

## Results by library

One row per sequencing library, grouped by sample.

- The 3- and 4-field columns are blank where the reference is two-field.
- HLA coverage is the mean depth across the eight classical HLA genes from the HLA-Resolve run log (all 143 libraries: mean 33.2×, median 34.0×, range 19.6–46.8×).

| Sample | Library | Cohort | Instrument | Reference | HLA&nbsp;coverage | Genes&nbsp;typed | Alleles&nbsp;called | 1-field | 2-field | 3-field | 4-field |
|---|---|---|---|---|---:|---|---|---|---|---|---|
| HG00097 | PG00097.HFSS | HPRC | Revio | 1kGP | 41.3× | 8/8 | 16/16 | 10/10 | 10/10 |  |  |
| HG00099 | PG00099_1.HFSS | HPRC | Revio + Sequel II | 1kGP | 37.4× | 8/8 | 16/16 | 10/10 | 10/10 |  |  |
| HG00099 | PG00099_2.HFSS | HPRC | Sequel II | 1kGP | 31.2× | 8/8 | 16/16 | 10/10 | 10/10 |  |  |
| HG00126 | PG00126.HFSS | HPRC | Revio | 1kGP | 27.2× | 8/8 | 16/16 | 10/10 | 10/10 |  |  |
| HG00128 | PG00128.HFSS | HPRC | Revio | 1kGP | 33.5× | 8/8 | 16/16 | 10/10 | 10/10 |  |  |
| HG00128 | PG00128.HFSS2 | HPRC | Revio | 1kGP | 25.6× | 8/8 | 16/16 | 10/10 | 10/10 |  |  |
| HG00133 | PG00133.HFSS | HPRC | Revio | 1kGP | 26.6× | 8/8 | 16/16 | 10/10 | 10/10 |  |  |
| HG00140 | HG00140_lib1 | HPRC | Sequel II | 1kGP | 33.4× | 8/8 | 16/16 | 8/8 | 8/8 |  |  |
| HG00146 | HG00146_lib1 | HPRC | Revio | 1kGP | 37.5× | 8/8 | 16/16 | 10/10 | 10/10 |  |  |
| HG00232 | HG00232_lib1 | HPRC | Revio | 1kGP | 34.4× | 8/8 | 16/16 | 8/8 | 8/8 |  |  |
| HG00235 | HG00235_PB1 | HPRC | Revio | 1kGP | 31.0× | 8/8 | 16/16 | 10/10 | 10/10 |  |  |
| HG00253 | HG00253_PB1 | HPRC | Revio | 1kGP | 35.5× | 8/8 | 16/16 | 8/8 | 8/8 |  |  |
| HG00272 | PG00272.HFSS | HPRC | Revio | 1kGP | 42.3× | 8/8 | 16/16 | 10/10 | 10/10 |  |  |
| HG00280 | PG00280.HFSS | HPRC | Sequel II | 1kGP | 26.1× | 8/8 | 16/16 | 8/8 | 8/8 |  |  |
| HG00290 | PG00290.HFSS | HPRC | Revio | 1kGP | 33.7× | 8/8 | 16/16 | 10/10 | 10/10 |  |  |
| HG00320 | HG00320_lib1 | HPRC | Revio | 1kGP | 30.3× | 8/8 | 16/16 | 8/8 | 8/8 |  |  |
| HG00321 | HG00321_lib1 | HPRC | Revio | 1kGP | 33.1× | 8/8 | 16/16 | 8/8 | 8/8 |  |  |
| HG00323 | HG00323_lib1 | HPRC | Sequel II | 1kGP | 32.7× | 8/8 | 16/16 | 8/8 | 8/8 |  |  |
| HG00344 | PG00344.HFSS2 | HPRC | Revio | 1kGP | 20.7× | 8/8 | 16/16 | 8/8 | 8/8 |  |  |
| HG00350 | PG00350.HFSS | HPRC | Revio | 1kGP | 27.8× | 8/8 | 16/16 | 10/10 | 10/10 |  |  |
| HG00438 | HG00438_lib1 | HPRC | Sequel II | Lai | 29.1× | 8/8 | 16/16 | 16/16 | 16/16 | 16/16 | 15/16 |
| HG00621 | HG00621_lib1 | HPRC | Sequel II | Lai | 39.4× | 8/8 | 16/16 | 16/16 | 16/16 | 16/16 | 15/16 |
| HG00673 | HG00673_lib1 | HPRC | Sequel II | Lai | 40.0× | 8/8 | 16/16 | 16/16 | 16/16 | 16/16 | 14/16 |
| HG00733 | HG00733.HFSS | HPRC_PLUS | Sequel II | Lai | 39.7× | 8/8 | 16/16 | 16/16 | 16/16 | 16/16 | 14/16 |
| HG00733 | HG00733:untagged | HPRC_PLUS | Sequel II | Lai | 37.1× | 8/8 | 16/16 | 16/16 | 16/16 | 16/16 | 15/16 |
| HG00735 | HG00735_lib1 | HPRC | Sequel II | Lai | 28.6× | 8/8 | 16/16 | 16/16 | 16/16 | 16/16 | 15/16 |
| HG00741 | HG00741_lib1 | HPRC | Sequel II | Lai | 35.9× | 8/8 | 16/16 | 16/16 | 16/16 | 16/16 | 14/16 |
| HG01071 | HG01071_lib1 | HPRC | Sequel II | Lai | 27.6× | 8/8 | 16/16 | 16/16 | 16/16 | 16/16 | 15/16 |
| HG01106 | HG01106_lib1 | HPRC | Sequel II | Lai | 32.4× | 8/8 | 16/16 | 16/16 | 16/16 | 16/16 | 15/16 |
| HG01109 | HG01109:untagged | HPRC_PLUS | Sequel II | Lai | 38.1× | 8/8 | 16/16 | 16/16 | 16/16 | 16/16 | 15/16 |
| HG01167 | HG01167_PB1 | HPRC | Revio | 1kGP | 33.1× | 8/8 | 16/16 | 10/10 | 10/10 |  |  |
| HG01175 | HG01175_lib1 | HPRC | Sequel II | Lai | 36.0× | 8/8 | 16/16 | 16/16 | 16/16 | 16/16 | 15/16 |
| HG01243 | HG01243:untagged | HPRC_PLUS | Sequel II | Lai | 35.6× | 8/8 | 16/16 | 16/16 | 16/16 | 16/16 | 15/16 |
| HG01258 | HG01258.HiFiEx_f1 | HPRC | Sequel II | Lai | 41.6× | 8/8 | 16/16 | 16/16 | 16/16 | 15/15 | 14/15 |
| HG01358 | HG01358_HiFiEx_f1 | HPRC | Sequel II | Lai | 37.0× | 8/8 | 16/16 | 15/15 | 15/15 | 15/15 | 14/14 |
| HG01361 | HG01361.HFSS3 | HPRC | Revio | Lai | 24.0× | 8/8 | 16/16 | 15/16 | 14/16 | 14/16 | 12/15 |
| HG01361 | HG01361.HiFiEx_f2 | HPRC | Sequel II | Lai | 33.2× | 8/8 | 16/16 | 16/16 | 16/16 | 16/16 | 13/15 |
| HG01530 | HG01530_PB1 | HPRC | Sequel II | 1kGP | 30.3× | 8/8 | 16/16 | 10/10 | 10/10 |  |  |
| HG01784 | HG01784_PB1 | HPRC | Sequel II | 1kGP | 33.7× | 8/8 | 16/16 | 10/10 | 10/10 |  |  |
| HG01784 | HG01784_PB2 | HPRC | Revio | 1kGP | 20.7× | 8/8 | 16/16 | 10/10 | 10/10 |  |  |
| HG01786 | HG01786_lib1 | HPRC | Revio | 1kGP | 37.9× | 8/8 | 16/16 | 10/10 | 10/10 |  |  |
| HG01891 | HG01891.HFSS | HPRC | Revio | Lai | 24.5× | 8/8 | 16/16 | 16/16 | 16/16 | 16/16 | 15/16 |
| HG01891 | HG01891.HiFiEx_f2 | HPRC | Sequel II | Lai | 38.8× | 8/8 | 16/16 | 16/16 | 16/16 | 16/16 | 15/16 |
| HG01928 | HG01928_lib1 | HPRC | Sequel II | Lai | 29.1× | 8/8 | 16/16 | 16/16 | 16/16 | 13/15 | 12/15 |
| HG01952 | HG01952_lib1 | HPRC | Sequel II | Lai | 29.6× | 8/8 | 16/16 | 16/16 | 16/16 | 16/16 | 13/14 |
| HG01978 | HG01978_lib1 | HPRC | Sequel II | Lai | 33.1× | 8/8 | 16/16 | 16/16 | 16/16 | 16/16 | 14/15 |
| HG02040 | HG02040.HFSS2 | HPRC | Sequel II | 1kGP | 25.8× | 8/8 | 16/16 | 10/10 | 9/10 |  |  |
| HG02055 | HG02055_Revio_validated | HPRC_PLUS | Revio | Lai | 34.0× | 8/8 | 16/16 | 16/16 | 16/16 | 16/16 | 15/16 |
| HG02080 | HG02080_ELF2_480c6e | HPRC_PLUS | Sequel II | Lai | 36.1× | 8/8 | 16/16 | 16/16 | 16/16 | 16/16 | 16/16 |
| HG02109 | HG02109_Revio_validated | HPRC_PLUS | Revio | Lai | 35.5× | 8/8 | 16/16 | 16/16 | 16/16 | 14/14 | 14/14 |
| HG02145 | HG02145_Revio_validated | HPRC_PLUS | Revio | Lai | 33.1× | 8/8 | 16/16 | 16/16 | 16/16 | 16/16 | 13/16 |
| HG02148 | HG02148_lib1 | HPRC | Sequel II | Lai | 27.6× | 8/8 | 16/16 | 16/16 | 16/16 | 16/16 | 15/16 |
| HG02155 | HG02155_lib1 | HPRC | Sequel II | 1kGP | 34.2× | 8/8 | 16/16 | 10/10 | 10/10 |  |  |
| HG02165 | PG02165.HFSS2 | HPRC | Sequel II | 1kGP | 30.4× | 8/8 | 16/16 | 10/10 | 10/10 |  |  |
| HG02178 | HG02178_PB1 | HPRC | Revio | 1kGP | 34.1× | 8/8 | 16/16 | 10/10 | 10/10 |  |  |
| HG02257 | HG02257.HFSS | HPRC | Revio | Lai | 22.4× | 8/8 | 16/16 | 16/16 | 16/16 | 16/16 | 14/15 |
| HG02257 | HG02257.HiFiEx_f2 | HPRC | Sequel II | Lai | 37.0× | 8/8 | 16/16 | 16/16 | 16/16 | 16/16 | 14/15 |
| HG02391 | HG02391_PB1 | HPRC | Sequel II | 1kGP | 30.0× | 8/8 | 16/16 | 10/10 | 10/10 |  |  |
| HG02392 | HG02392_PB1 | HPRC | Revio | 1kGP | 32.1× | 8/8 | 16/16 | 10/10 | 10/10 |  |  |
| HG02572 | HG02572.HFSS3 | HPRC | Revio | Lai | 24.3× | 8/8 | 16/16 | 16/16 | 16/16 | 15/15 | 12/14 |
| HG02572 | HG02572.HiFiEx_f2 | HPRC | Sequel II | Lai | 27.6× | 8/8 | 16/16 | 16/16 | 16/16 | 15/15 | 12/14 |
| HG02583 | HG02583_PB1 | HPRC | Revio | 1kGP | 33.5× | 8/8 | 16/16 | 10/10 | 10/10 |  |  |
| HG02622 | HG02622_lib1 | HPRC | Sequel II | Lai | 32.8× | 8/8 | 16/16 | 16/16 | 16/16 | 16/16 | 13/16 |
| HG02630 | HG02630_lib1 | HPRC | Sequel II | Lai | 37.2× | 8/8 | 16/16 | 16/16 | 16/16 | 16/16 | 15/16 |
| HG02717 | HG02717_lib1 | HPRC | Sequel II | Lai | 33.5× | 8/8 | 16/16 | 16/16 | 16/16 | 15/15 | 11/15 |
| HG02723 | HG02723_Revio_validated | HPRC_PLUS | Revio | Lai | 34.4× | 8/8 | 16/16 | 16/16 | 16/16 | 15/15 | 14/15 |
| HG02723 | HG02723:untagged | HPRC_PLUS | Sequel II | Lai | 43.4× | 8/8 | 16/16 | 16/16 | 16/16 | 15/15 | 14/15 |
| HG02818 | HG02818:untagged | HPRC_PLUS | Sequel II | Lai | 37.1× | 8/8 | 16/16 | 16/16 | 16/16 | 16/16 | 15/16 |
| HG02886 | HG02886_lib1 | HPRC | Sequel II | Lai | 35.2× | 8/8 | 16/16 | 16/16 | 16/16 | 16/16 | 14/16 |
| HG02922 | PG02922.HFSS2 | HPRC | Revio + Sequel II | 1kGP | 33.4× | 8/8 | 16/16 | 10/10 | 10/10 |  |  |
| HG02965 | HG02965_lib1 | HPRC | Sequel II | 1kGP | 34.5× | 8/8 | 16/16 | 10/10 | 10/10 |  |  |
| HG02976 | HG02976_lib1 | HPRC | Sequel II | 1kGP | 27.6× | 8/8 | 16/16 | 10/10 | 9/10 |  |  |
| HG03098 | HG03098_Fraction2_Fraction3_480cnp | HPRC_PLUS | Sequel II | Lai | 35.0× | 8/8 | 16/16 | 16/16 | 16/16 | 15/15 | 15/15 |
| HG03130 | HG03130.HFSS | HPRC | Sequel II | 1kGP | 39.6× | 8/8 | 16/16 | 10/10 | 10/10 |  |  |
| HG03139 | HG03139.HFSS | HPRC | Revio + Sequel II | 1kGP | 33.6× | 8/8 | 16/16 | 10/10 | 10/10 |  |  |
| HG03195 | HG03195_lib1 | HPRC | Sequel II | 1kGP | 32.1× | 8/8 | 16/16 | 10/10 | 9/10 |  |  |
| HG03209 | HG03209.HFSS | HPRC | Sequel II | 1kGP | 34.6× | 8/8 | 16/16 | 10/10 | 10/10 |  |  |
| HG03225 | HG03225_lib1 | HPRC | Sequel II | 1kGP | 36.5× | 8/8 | 16/16 | 10/10 | 10/10 |  |  |
| HG03270 | HG03270_PB1 | HPRC | Revio | 1kGP | 37.5× | 8/8 | 16/16 | 10/10 | 10/10 |  |  |
| HG03369 | HG03369_PB1 | HPRC | Revio | 1kGP | 32.9× | 8/8 | 16/16 | 10/10 | 10/10 |  |  |
| HG03453 | HG03453_lib1 | HPRC | Sequel II | Lai | 38.7× | 8/8 | 16/16 | 16/16 | 16/16 | 16/16 | 14/15 |
| HG03470 | HG03470_PB1 | HPRC | Revio + Sequel II | 1kGP | 35.1× | 8/8 | 16/16 | 10/10 | 10/10 |  |  |
| HG03486 | HG03486:untagged | HPRC_PLUS | Sequel II | Lai | 39.5× | 8/8 | 16/16 | 16/16 | 16/16 | 16/16 | 16/16 |
| HG03492 | HG03492_ELF2_ELF3_480c6f | HPRC_PLUS | Sequel II | Lai | 37.0× | 8/8 | 16/16 | 16/16 | 16/16 | 16/16 | 15/15 |
| HG03516 | HG03516.HFSS | HPRC | Revio | Lai | 24.4× | 8/8 | 16/16 | 16/16 | 16/16 | 16/16 | 16/16 |
| HG03516 | HG03516_HiFiEx_mix | HPRC | Sequel II | Lai | 39.7× | 8/8 | 16/16 | 16/16 | 16/16 | 16/16 | 15/16 |
| HG03521 | PG03521.HFSS2 | HPRC | Revio | 1kGP | 36.6× | 8/8 | 16/16 | 10/10 | 10/10 |  |  |
| HG03540 | HG03540_lib1 | HPRC | Sequel II | Lai | 41.7× | 8/8 | 16/16 | 16/16 | 16/16 | 15/15 | 14/14 |
| HG03579 | HG03579_lib1 | HPRC | Sequel II | Lai | 39.8× | 8/8 | 16/16 | 16/16 | 16/16 | 16/16 | 16/16 |
| HG03583 | HG03583_PB1 | HPRC | Revio | 1kGP | 30.1× | 8/8 | 16/16 | 10/10 | 10/10 |  |  |
| HG03742 | HG03742_PB1 | HPRC | Sequel II | 1kGP | 26.5× | 8/8 | 16/16 | 10/10 | 9/10 |  |  |
| HG03784 | HG03784_lib1 | HPRC | Revio | 1kGP | 35.7× | 8/8 | 16/16 | 10/10 | 10/10 |  |  |
| HG03874 | HG03874_PB1 | HPRC | Revio | 1kGP | 27.0× | 8/8 | 16/16 | 10/10 | 10/10 |  |  |
| NA18505 | NA18505_PB1 | HPRC | Revio | 1kGP | 36.2× | 8/8 | 16/16 | 10/10 | 10/10 |  |  |
| NA18508 | PG18508.HFSS2 | HPRC | Revio | 1kGP | 20.6× | 8/8 | 16/16 | 10/10 | 10/10 |  |  |
| NA18522 | NA18522.HFSS | HPRC | Revio + Sequel II | 1kGP | 34.0× | 8/8 | 16/16 | 10/10 | 10/10 |  |  |
| NA18565 | NA18565_lib1 | HPRC | Revio | 1kGP | 34.6× | 8/8 | 16/16 | 10/10 | 10/10 |  |  |
| NA18570 | NA18570_PB1 | HPRC | Sequel II | 1kGP | 28.8× | 8/8 | 16/16 | 8/8 | 8/8 |  |  |
| NA18608 | NA18608_lib1 | HPRC | Revio | 1kGP | 44.0× | 8/8 | 16/16 | 8/8 | 8/8 |  |  |
| NA18620 | NA18620_PB1 | HPRC | Revio | 1kGP | 36.6× | 8/8 | 16/16 | 8/8 | 7/8 |  |  |
| NA18747 | PG18747_1.HFSS | HPRC | Revio + Sequel II | 1kGP | 38.0× | 8/8 | 16/16 | 10/10 | 10/10 |  |  |
| NA18747 | PG18747_2.HFSS | HPRC | Sequel II | 1kGP | 20.3× | 8/8 | 16/16 | 10/10 | 10/10 |  |  |
| NA18879 | PG18879.HFSS | HPRC | Revio | 1kGP | 32.7× | 8/8 | 16/16 | 10/10 | 10/10 |  |  |
| NA18906 | NA18906:untagged | HPRC_PLUS | Sequel II | Lai | 46.8× | 8/8 | 16/16 | 16/16 | 16/16 | 16/16 | 15/15 |
| NA18952 | NA18952_lib1 | HPRC | Revio | 1kGP | 35.9× | 8/8 | 16/16 | 8/8 | 8/8 |  |  |
| NA18971 | NA18971_lib1 | HPRC | Sequel II | 1kGP | 35.7× | 8/8 | 16/16 | 10/10 | 10/10 |  |  |
| NA18974 | NA18974_lib1 | HPRC | Revio | 1kGP | 40.1× | 8/8 | 16/16 | 10/10 | 10/10 |  |  |
| NA18976 | NA18976_lib1 | HPRC | Revio | 1kGP | 38.1× | 8/8 | 16/16 | 10/10 | 10/10 |  |  |
| NA18983 | PG18983.HFSS2 | HPRC | Sequel II | 1kGP | 25.0× | 8/8 | 16/16 | 10/10 | 10/10 |  |  |
| NA19036 | NA19036_PB1 | HPRC | Revio | 1kGP | 34.9× | 8/8 | 16/16 | 10/10 | 10/10 |  |  |
| NA19043 | PG19043.HFSS | HPRC | Revio + Sequel II | 1kGP | 38.1× | 8/8 | 16/16 | 10/10 | 10/10 |  |  |
| NA19087 | NA19087_PB1 | HPRC | Sequel II | 1kGP | 25.1× | 8/8 | 16/16 | 8/8 | 8/8 |  |  |
| NA19159 | NA19159_PB1 | HPRC | Revio + Sequel II | 1kGP | 34.7× | 8/8 | 16/16 | 10/10 | 9/10 |  |  |
| NA19185 | NA19185_PB1 | HPRC | Sequel II | 1kGP | 29.5× | 8/8 | 16/16 | 10/10 | 9/10 |  |  |
| NA19240 | NA19240:untagged | HPRC_PLUS | Sequel II | Lai | 41.8× | 8/8 | 16/16 | 16/16 | 16/16 | 16/16 | 15/16 |
| NA19338 | NA19338_PB1 | HPRC | Sequel II | 1kGP | 35.2× | 8/8 | 16/16 | 10/10 | 10/10 |  |  |
| NA19338 | NA19338_PB2 | HPRC | Revio | 1kGP | 19.6× | 8/8 | 16/16 | 9/10 | 9/10 |  |  |
| NA19391 | NA19391_PB1 | HPRC | Sequel II | 1kGP | 31.3× | 8/8 | 16/16 | 10/10 | 10/10 |  |  |
| NA19391 | NA19391_PB3 | HPRC | Revio | 1kGP | 22.0× | 8/8 | 16/16 | 10/10 | 10/10 |  |  |
| NA19443 | PG19443.HFSS | HPRC | Revio | 1kGP | 34.0× | 8/8 | 16/16 | 10/10 | 10/10 |  |  |
| NA19468 | NA19468_PB1 | HPRC | Sequel II | 1kGP | 30.9× | 8/8 | 16/16 | 10/10 | 10/10 |  |  |
| NA19682 | PG19682.HFSS | HPRC | Revio | 1kGP | 30.6× | 8/8 | 16/16 | 10/10 | 10/10 |  |  |
| NA19700 | PG19700.HFSS | HPRC | Revio | 1kGP | 32.2× | 8/8 | 16/16 | 10/10 | 10/10 |  |  |
| NA19776 | PG19776.HFSS | HPRC | Revio | 1kGP | 21.6× | 8/8 | 16/16 | 10/10 | 10/10 |  |  |
| NA19835 | NA19835_PB1 | HPRC | Revio | 1kGP | 31.0× | 8/8 | 16/16 | 10/10 | 10/10 |  |  |
| NA19909 | NA19909_Fraction2_Fraction3 | HPRC | Revio | 1kGP | 35.0× | 8/8 | 16/16 | 10/10 | 10/10 |  |  |
| NA20129 | NA20129:untagged | HPRC_PLUS | Sequel II | Lai | 38.2× | 8/8 | 16/16 | 16/16 | 16/16 | 16/16 | 15/15 |
| NA20282 | NA20282_lib1 | HPRC | Revio | 1kGP | 35.8× | 8/8 | 16/16 | 10/10 | 10/10 |  |  |
| NA20346 | NA20346_lib1 | HPRC | Revio | 1kGP | 27.5× | 8/8 | 16/16 | 10/10 | 10/10 |  |  |
| NA20752 | PG20752.HFSS | HPRC | Sequel II | 1kGP | 34.0× | 8/8 | 16/16 | 8/8 | 8/8 |  |  |
| NA20799 | NA20799_PB1 | HPRC | Sequel II | 1kGP | 28.4× | 8/8 | 16/16 | 10/10 | 9/10 |  |  |
| NA20799 | NA20799_PB2 | HPRC | Revio | 1kGP | 24.8× | 8/8 | 16/16 | 10/10 | 10/10 |  |  |
| NA20805 | NA20805_lib1 | HPRC | Sequel II | 1kGP | 40.4× | 8/8 | 16/16 | 10/10 | 10/10 |  |  |
| NA20809 | NA20809_lib1 | HPRC | Revio | 1kGP | 32.0× | 8/8 | 16/16 | 10/10 | 10/10 |  |  |
| NA20850 | PG20850.HFSS | HPRC | Revio | 1kGP | 23.1× | 8/8 | 16/16 | 10/10 | 10/10 |  |  |
| NA20870 | NA20870_PB1 | HPRC | Revio | 1kGP | 38.8× | 8/8 | 16/16 | 10/10 | 10/10 |  |  |
| NA20905 | PG20905.HFSS | HPRC | Revio + Sequel II | 1kGP | 36.4× | 8/8 | 16/16 | 10/10 | 10/10 |  |  |
| NA21093 | PG21093.HFSS | HPRC | Revio | 1kGP | 36.3× | 8/8 | 16/16 | 10/10 | 10/10 |  |  |
| NA21102 | NA21102_lib1 | HPRC | Revio | 1kGP | 31.6× | 8/8 | 16/16 | 10/10 | 10/10 |  |  |
| NA21106 | PG21106.HFSS | HPRC | Revio | 1kGP | 42.4× | 8/8 | 16/16 | 10/10 | 10/10 |  |  |
| NA21110 | NA21110_PB1 | HPRC | Revio | 1kGP | 34.8× | 8/8 | 16/16 | 10/10 | 10/10 |  |  |
| NA21144 | NA21144_lib1 | HPRC | Revio | 1kGP | 40.0× | 8/8 | 16/16 | 10/10 | 10/10 |  |  |
| NA21309 | NA21309:untagged | HPRC_PLUS | Sequel II | Lai | 37.2× | 8/8 | 16/16 | 16/16 | 16/16 | 16/16 | 15/15 |

## Discordant alleles

Gene × library combinations discordant at 1 or 2 or 3 fields, with the shallowest discordant field, the gene's mean ARS depth, the calls and the reference alleles.

| Sample | Library | Gene | Field | ARS&nbsp;depth | Called | Reference | Note |
|---|---|---|---|---:|---|---|---|
| HG01361 | HG01361.HFSS3 | DRB1 | 1 | 18.0× | DRB1*15:01:32 / DRB1*11:34 | DRB1*11:02:01:02 / DRB1*07:01:01:01 | ARS depth under 20x |
| NA19338 | NA19338_PB2 | DRB1 | 1 | 16.0× | DRB1*15:03:01:03 / DRB1*15:03:01:03 | DRB1*13:02 / DRB1*15:03 | ARS depth under 20x |
| HG02040 | HG02040.HFSS2 | A | 2 | 24.0× | A*02:03:01:01 / A*29:01:01:01 | A*02:148/A*02:281/A*02:370/A*02:427/A*02:544/A*02:595/A*02:634 / A*29:01 |  |
| HG02976 | HG02976_lib1 | DRB1 | 2 | 23.0× | DRB1*15:03:01:03 / DRB1*11:01:02:03 | DRB1*11:10 / DRB1*15:03 |  |
| HG03195 | HG03195_lib1 | B | 2 | 18.4× | B*35:598 / B*15:03:01:02 | B*15:03 / B*35:01 | ARS depth under 20x |
| HG03742 | HG03742_PB1 | B | 2 | 30.0× | B*52:01:01:09 / B*37:110 | B*37:01 / B*52:01 |  |
| NA18620 | NA18620_PB1 | C | 2 | 47.0× | C*04:82:01 / C*07:02:01:15 | C*04:01 / C*07:02 |  |
| NA19159 | NA19159_PB1 | DRB1 | 2 | 25.9× | DRB1*07:01:01:01 / DRB1*13:01:01:04 | DRB1*13:177 / DRB1*07:01 | panel allele flagged exome-only (*) |
| NA19185 | NA19185_PB1 | C | 2 | 27.6× | C*17:01:01:02 / C*16:01:01:01 | C*16:01 / C*17:03 |  |
| NA20799 | NA20799_PB1 | DRB1 | 2 | 31.3× | DRB1*11:11:01 / DRB1*01:02:01:01 | DRB1*01:02 / DRB1*11:01 |  |
| HG01928 | HG01928_lib1 | A | 3 | 24.0× | A*02:01:52 / A*02:01:52 | A*02:01:01:01 / A*02:01:01:01 |  |

## Input files

Each library was typed from one unaligned BAM made by concatenating the files listed. HPRC files are public on the AWS Open Data registry and need no credential. Build a path by substituting the Cohort, Sample and Files values from the table into:

```
s3://human-pangenomics/working/<COHORT>/<SAMPLE>/raw_data/PacBio_HiFi/<FILE>
```

For example, the first file listed for HG00097 (cohort HPRC) is at:

```
s3://human-pangenomics/working/HPRC/HG00097/raw_data/PacBio_HiFi/m84046_230716_051716_s3.hifi_reads.bc2086.bam
```

Files are public and need no credentials:

```bash
aws s3 cp --no-sign-request \
  s3://human-pangenomics/working/HPRC/HG00097/raw_data/PacBio_HiFi/m84046_230716_051716_s3.hifi_reads.bc2086.bam .
```

Some files sit one folder deeper (`primrose/`, `wMods/`). The manifest below gives the exact key.

| Sample | Library | Cohort | Files |
|---|---|---|---|
| HG00097 | PG00097.HFSS | HPRC | m84046_230716_051716_s3.hifi_reads.bc2086.bam<br>m84046_230716_054822_s4.hifi_reads.bc2086.bam |
| HG00099 | PG00099_1.HFSS | HPRC | m54329U_220825_174247-bc2012.5mc.hifi_reads.bam<br>m84046_230517_211630_s4.hifi_reads.bc2012.bam |
| HG00099 | PG00099_2.HFSS | HPRC | m54329U_220827_143814-bc2050.5mc.hifi_reads.bam<br>m54329U_220829_095708-bc2050.5mc.hifi_reads.bam |
| HG00126 | PG00126.HFSS | HPRC | m84046_230602_203207_s4.hifi_reads.bc2043.bam<br>m84046_230617_043254_s3.hifi_reads.bc2043.bam<br>m84046_230617_050400_s4.hifi_reads.bc2043.bam |
| HG00128 | PG00128.HFSS | HPRC | m84046_230602_203207_s4.hifi_reads.bc2044.bam<br>m84046_230617_043254_s3.hifi_reads.bc2044.bam<br>m84046_230617_050400_s4.hifi_reads.bc2044.bam |
| HG00128 | PG00128.HFSS2 | HPRC | m84046_230630_233157_s3.bc2069--bc2069.bam |
| HG00133 | PG00133.HFSS | HPRC | m84046_230602_203207_s4.hifi_reads.bc2045.bam<br>m84046_230617_043254_s3.hifi_reads.bc2045.bam<br>m84046_230617_050400_s4.hifi_reads.bc2045.bam |
| HG00140 | HG00140_lib1 | HPRC | m64043_220728_173215-bc1018.5mc.hifi_reads.bam<br>m64136_220715_182717-bc1018.5mc.hifi_reads.bam<br>m64136_220717_152248-bc1018.5mc.hifi_reads.bam<br>m64136_220719_122056-bc1018.5mc.hifi_reads.bam |
| HG00146 | HG00146_lib1 | HPRC | m84081_230616_182824_s3.hifi_reads.bc2013.bam |
| HG00232 | HG00232_lib1 | HPRC | m84081_230623_202140_s1.hifi_reads.bc2014.bam<br>m84081_230714_201817_s1.hifi_reads.bc2014.bam<br>m84081_230714_205519_s2.hifi_reads.bc2014.bam |
| HG00235 | HG00235_PB1 | HPRC | m84091_230721_141900_s2.hifi_reads.bc1008.bam |
| HG00253 | HG00253_PB1 | HPRC | m84091_230721_134835_s1.hifi_reads.bc1003.bam |
| HG00272 | PG00272.HFSS | HPRC | m84046_230716_051716_s3.hifi_reads.bc2087.bam<br>m84046_230716_054822_s4.hifi_reads.bc2087.bam |
| HG00280 | PG00280.HFSS | HPRC | m54329U_220901_221341-bc2051.5mc.hifi_reads.bam<br>m64076_220831_191646-bc2051.5mc.hifi_reads.bam |
| HG00290 | PG00290.HFSS | HPRC | m84046_230623_203534_s4.hifi_reads.bc2052.bam<br>m84046_230627_201017_s1.hifi_reads.bc2052.bam |
| HG00320 | HG00320_lib1 | HPRC | m84081_230610_202430_s3.hifi_reads.bc2011.bam |
| HG00321 | HG00321_lib1 | HPRC | m84081_230623_212309_s3.hifi_reads.bc2015.bam<br>m84081_230630_190803_s1.hifi_reads.bc2015.bam |
| HG00323 | HG00323_lib1 | HPRC | m64043_220728_173215-bc1016.5mc.hifi_reads.bam<br>m64136_220708_205611-bc1016.5mc.hifi_reads.bam<br>m64136_220710_175208-bc1016.5mc.hifi_reads.bam<br>m64136_220712_144943-bc1016.5mc.hifi_reads.bam |
| HG00344 | PG00344.HFSS2 | HPRC | m84046_230715_064007_s2.hifi_reads.bc2083.bam |
| HG00350 | PG00350.HFSS | HPRC | m84046_230701_000218_s4.hifi_reads.bc2048.bam<br>m84046_230712_231732_s2.hifi_reads.bc2048.bam |
| HG00438 | HG00438_lib1 | HPRC | m64043_200710_174426.ccs.bam<br>m64043_200711_235708.ccs.bam<br>m64043_200713_062240.ccs.bam |
| HG00621 | HG00621_lib1 | HPRC | m64136_200710_174522.ccs.bam<br>m64136_200711_235843.ccs.bam<br>m64136_200713_062514.ccs.bam<br>m64136_200714_125149.ccs.bam |
| HG00673 | HG00673_lib1 | HPRC | m64043_200716_182902.ccs.bam<br>m64043_200718_004213.ccs.bam<br>m64043_200719_070806.ccs.bam<br>m64043_200720_133355.ccs.bam |
| HG00733 | HG00733.HFSS | HPRC_PLUS | m64076_211214_012715.hifi_reads.bam<br>m64076_211215_225159.hifi_reads.bam<br>m64076_211217_081943.hifi_reads.bam |
| HG00733 | HG00733:untagged | HPRC_PLUS | HG00733_20190925_EEE_m54329U_190607_185248.ccs.bam<br>HG00733_20190925_EEE_m54329U_190615_010947.ccs.bam<br>HG00733_20190925_EEE_m54329U_190617_231905.ccs.bam<br>HG00733_20190925_EEE_m54329U_190619_052546.ccs.bam<br>HG00733_20190925_EEE_m54329U_190629_180018.ccs.bam<br>HG00733_20190925_EEE_m54329U_190701_222759.ccs.bam<br>HG00733_20190925_EEE_m54329U_190827_173812.ccs.bam |
| HG00735 | HG00735_lib1 | HPRC | m64043_200702_173033.ccs.bam<br>m64043_200703_234328.ccs.bam<br>m64043_200705_060840.ccs.bam |
| HG00741 | HG00741_lib1 | HPRC | m64136_200625_174949.ccs.bam<br>m64136_200627_000247.ccs.bam<br>m64136_200628_062837.ccs.bam<br>m64136_200629_125431.ccs.bam |
| HG01071 | HG01071_lib1 | HPRC | m64136_200702_173125.ccs.bam<br>m64136_200703_234438.ccs.bam<br>m64136_200705_061033.ccs.bam<br>m64136_200706_123635.ccs.bam |
| HG01106 | HG01106_lib1 | HPRC | m64043_200625_174853.ccs.bam<br>m64043_200627_000137.ccs.bam<br>m64043_200628_062711.ccs.bam |
| HG01109 | HG01109:untagged | HPRC_PLUS | m64043_200827_191459.ccs.bam<br>m64043_200829_012836.ccs.bam<br>m64043_200830_075523.ccs.bam |
| HG01167 | HG01167_PB1 | HPRC | m84091_230725_164143_s3.hifi_reads.bc1012.bam |
| HG01175 | HG01175_lib1 | HPRC | m64043_200618_201934.ccs.bam<br>m64043_200620_173220.ccs.bam<br>m64043_200621_234442.ccs.bam<br>m64043_200623_060946.ccs.bam |
| HG01243 | HG01243:untagged | HPRC_PLUS | m64136_200827_191603.ccs.bam<br>m64136_200829_012933.ccs.bam<br>m64136_200830_075556.ccs.bam |
| HG01258 | HG01258.HiFiEx_f1 | HPRC | wMods/m64076_200306_185917.ccs_reprocessed.jasmine.5mc.hifi_reads.bam<br>wMods/m64076_200308_194406.ccs_reprocessed.jasmine.5mc.hifi_reads.bam<br>wMods/m64076_200310_015720.ccs_reprocessed.jasmine.5mc.hifi_reads.bam<br>wMods/m64076_200311_082315.ccs_reprocessed.jasmine.5mc.hifi_reads.bam |
| HG01358 | HG01358_HiFiEx_f1 | HPRC | m64076_200201_051547.ccs.bam<br>m64076_200203_181219.ccs.bam<br>m64076_200206_215943.ccs.bam<br>m64076_200208_041234.ccs.bam |
| HG01361 | HG01361.HFSS3 | HPRC | m84046_231202_071034_s1.hifi_reads.bc2052.bam |
| HG01361 | HG01361.HiFiEx_f2 | HPRC | wMods/m54329U_200306_185930.ccs_reprocessed.jasmine.5mc.hifi_reads.bam<br>wMods/m54329U_200308_194417.ccs_reprocessed.jasmine.5mc.hifi_reads.bam<br>wMods/m54329U_200310_015838.ccs_reprocessed.jasmine.5mc.hifi_reads.bam |
| HG01530 | HG01530_PB1 | HPRC | primrose/m64055e_220915_170505.5mc.hifi_reads.bam<br>primrose/m64334e_220920_180418.5mc.hifi_reads.bam |
| HG01784 | HG01784_PB1 | HPRC | primrose/m64055e_220920_175231.5mc.hifi_reads.bam<br>primrose/m64055e_220922_130045.5mc.hifi_reads.bam<br>primrose/m64330e_220917_122517.5mc.hifi_reads.bam |
| HG01784 | HG01784_PB2 | HPRC | m84091_231120_182829_s1.hifi_reads.bc1019.bam |
| HG01786 | HG01786_lib1 | HPRC | m84081_230616_172650_s1.hifi_reads.bc2012.bam |
| HG01891 | HG01891.HFSS | HPRC | m84046_230712_224626_s1.hifi_reads.bc2066.bam |
| HG01891 | HG01891.HiFiEx_f2 | HPRC | m54329U_200124_193652.ccs.bam<br>m54329U_200127_180554.ccs.bam<br>m54329U_200129_001928.ccs.bam<br>m54329U_200130_064539.ccs.bam |
| HG01928 | HG01928_lib1 | HPRC | m64043_200612_200936.ccs.bam<br>m64043_200614_191756.ccs.bam<br>m64043_200616_013031.ccs.bam<br>m64043_200617_075530.ccs.bam |
| HG01952 | HG01952_lib1 | HPRC | m64136_200612_201033.ccs.bam<br>m64136_200614_192134.ccs.bam<br>m64136_200616_013426.ccs.bam |
| HG01978 | HG01978_lib1 | HPRC | m64136_200530_164818.ccs.bam<br>m64136_200602_153012.ccs.bam<br>m64136_200603_214308.ccs.bam<br>m64136_200605_040848.ccs.bam |
| HG02040 | HG02040.HFSS2 | HPRC | m64076_220908_204946-bc2056.5mc.hifi_reads.bam<br>m64076_220910_174812-bc2056.5mc.hifi_reads.bam |
| HG02055 | HG02055_Revio_validated | HPRC_PLUS | m84039_230303_005138_s2.hifi_reads.bc2010.bam<br>m84039_230311_024907_s1.hifi_reads.bc2010.bam<br>m84039_230314_210032_s1.hifi_reads.bc2010.bam<br>m84039_230316_185945_s1.hifi_reads.bc2010.bam |
| HG02080 | HG02080_ELF2_480c6e | HPRC_PLUS | m64043_200904_190723.ccs.bam<br>m64043_200906_012211.ccs.bam<br>m64043_200907_074948.ccs.bam |
| HG02109 | HG02109_Revio_validated | HPRC_PLUS | m84036_230317_175945_s2.hifi_reads.bc2013.bam<br>m84039_230303_012244_s3.hifi_reads.bc2013.bam<br>m84039_230314_213047_s2.hifi_reads.bc2013.bam<br>m84039_230316_193003_s2.hifi_reads.bc2013.bam |
| HG02145 | HG02145_Revio_validated | HPRC_PLUS | m84036_230317_175945_s2.hifi_reads.bc2014.bam<br>m84039_230303_012244_s3.hifi_reads.bc2014.bam<br>m84039_230314_213047_s2.hifi_reads.bc2014.bam |
| HG02148 | HG02148_lib1 | HPRC | m64136_200618_202033.ccs.bam<br>m64136_200620_173618.ccs.bam<br>m64136_200621_234916.ccs.bam |
| HG02155 | HG02155_lib1 | HPRC | m64043_220728_173215-bc1012.5mc.hifi_reads.bam<br>m64136_220701_181202-bc1012.5mc.hifi_reads.bam<br>m64136_220703_150713-bc1012.5mc.hifi_reads.bam<br>m64136_220705_120453-bc1012.5mc.hifi_reads.bam |
| HG02165 | PG02165.HFSS2 | HPRC | m64076_220902_143723-bc2054.5mc.hifi_reads.bam<br>m64076_220904_113359-bc2054.5mc.hifi_reads.bam |
| HG02178 | HG02178_PB1 | HPRC | m84091_230629_191221_s4.hifi_reads.bc1008.bam<br>m84091_230710_180033_s3.hifi_reads.bc1008.bam |
| HG02257 | HG02257.HFSS | HPRC | m84046_230628_182559_s4.hifi_reads.bc2064.bam |
| HG02257 | HG02257.HiFiEx_f2 | HPRC | m64076_200125_231256.ccs.bam<br>m64076_200127_180545.ccs.bam<br>m64076_200129_001835.ccs.bam<br>m64076_200130_064345.ccs.bam |
| HG02391 | HG02391_PB1 | HPRC | primrose/m54306Ue_220821_023435.5mc.hifi_reads.bam<br>primrose/m64055e_220821_033925.5mc.hifi_reads.bam |
| HG02392 | HG02392_PB1 | HPRC | m84091_230629_181009_s2.hifi_reads.bc1002.bam<br>m84091_230710_180033_s3.hifi_reads.bc1002.bam |
| HG02572 | HG02572.HFSS3 | HPRC | m84046_230711_233802_s3.hifi_reads.default.bam |
| HG02572 | HG02572.HiFiEx_f2 | HPRC | m54329U_200319_002813.ccs.bam<br>m64076_200313_161705.ccs.bam<br>m64076_200317_041201.ccs.bam<br>m64076_200318_103811.ccs.bam |
| HG02583 | HG02583_PB1 | HPRC | m84091_230707_183012_s4.hifi_reads.bc1016.bam |
| HG02622 | HG02622_lib1 | HPRC | m64043_200530_164723.ccs.bam<br>m64043_200601_191521.ccs.bam<br>m64043_200603_012738.ccs.bam |
| HG02630 | HG02630_lib1 | HPRC | m64043_200501_162248.ccs.bam<br>m64043_200502_223511.ccs.bam<br>m64043_200504_050026.ccs.bam |
| HG02717 | HG02717_lib1 | HPRC | m64043_200403_163826.ccs.bam<br>m64043_200405_180950.ccs.bam<br>m64043_200407_002219.ccs.bam |
| HG02723 | HG02723_Revio_validated | HPRC_PLUS | m84036_230317_175945_s2.hifi_reads.bc2012.bam<br>m84039_230303_012244_s3.hifi_reads.bc2012.bam<br>m84039_230314_213047_s2.hifi_reads.bc2012.bam<br>m84039_230316_193003_s2.hifi_reads.bc2012.bam |
| HG02723 | HG02723:untagged | HPRC_PLUS | m64043_191221_024136.ccs.bam<br>m64043_191223_180311.ccs.bam<br>m64043_191225_001554.ccs.bam<br>m64043_191226_064057.ccs.bam |
| HG02818 | HG02818:untagged | HPRC_PLUS | m64043_200206_173947.ccs.bam<br>m64043_200207_235213.ccs.bam<br>m64043_200209_061852.ccs.bam<br>m64043_200314_004623.ccs.bam<br>m64043_200315_071057.ccs.bam<br>m64043_200316_214923.ccs.bam<br>m64043_200318_040100.ccs.bam |
| HG02886 | HG02886_lib1 | HPRC | m64043_200410_214826.ccs.bam<br>m64043_200412_040054.ccs.bam<br>m64043_200413_102554.ccs.bam |
| HG02922 | PG02922.HFSS2 | HPRC | m54329U_220903_190900-bc2053.5mc.hifi_reads.bam<br>m54329U_220905_144015-bc2053.5mc.hifi_reads.bam<br>m84046_230617_040148_s2.hifi_reads.bc2053.bam |
| HG02965 | HG02965_lib1 | HPRC | m64043_220715_182700-bc1017.5mc.hifi_reads.bam<br>m64043_220717_152237-bc1017.5mc.hifi_reads.bam<br>m64043_220719_122051-bc1017.5mc.hifi_reads.bam |
| HG02976 | HG02976_lib1 | HPRC | m64043_220617_195135-bc1009.5mc.hifi_reads.bam<br>m64043_220619_164635-bc1009.5mc.hifi_reads.bam<br>m64043_220621_134351-bc1009.5mc.hifi_reads.bam |
| HG03098 | HG03098_Fraction2_Fraction3_480cnp | HPRC_PLUS | m64043_201128_031055.ccs.bam<br>m64043_201203_004011.ccs.bam<br>m64043_201123_083343.ccs.bam |
| HG03130 | HG03130.HFSS | HPRC | m64076_220513_215716-bc2018.5mc.hifi_reads.bam<br>m64076_220516_221911-bc2018.5mc.hifi_reads.bam<br>m64076_220518_191414-bc2018.5mc.hifi_reads.bam |
| HG03139 | HG03139.HFSS | HPRC | m54329U_220519_225150-bc2019.5mc.hifi_reads.bam<br>m54329U_220523_155616-bc2019.5mc.hifi_reads.bam<br>m84046_230602_220440_s3.hifi_reads.bc2019.bam |
| HG03195 | HG03195_lib1 | HPRC | m64043_220726_203720-bc1010.5mc.hifi_reads.bam<br>m64136_220617_195203-bc1010.5mc.hifi_reads.bam<br>m64136_220619_164654-bc1010.5mc.hifi_reads.bam<br>m64136_220621_134558-bc1010.5mc.hifi_reads.bam |
| HG03209 | HG03209.HFSS | HPRC | m64076_220520_234906-bc2020.5mc.hifi_reads.bam<br>m64076_220524_161808-bc2020.5mc.hifi_reads.bam<br>m64076_220526_115049-bc2020.5mc.hifi_reads.bam |
| HG03225 | HG03225_lib1 | HPRC | m64043_220516_150157-bc1008.5mc.hifi_reads.bam<br>m64043_220518_120044-bc1008.5mc.hifi_reads.bam<br>m64043_220726_203720-bc1008.5mc.hifi_reads.bam<br>m64136_220518_004520-bc1008.5mc.hifi_reads.bam |
| HG03270 | HG03270_PB1 | HPRC | m84091_230707_175906_s3.hifi_reads.bc1015.bam |
| HG03369 | HG03369_PB1 | HPRC | m84091_230721_145006_s3.hifi_reads.bc1009.bam |
| HG03453 | HG03453_lib1 | HPRC | m64043_200508_172634.ccs.bam<br>m64043_200509_233929.ccs.bam<br>m64043_200511_060458.ccs.bam |
| HG03470 | HG03470_PB1 | HPRC | primrose/m54306Ue_220817_100412.5mc.hifi_reads.bam<br>primrose/m54306Ue_220819_070258.5mc.hifi_reads.bam<br>m84091_230705_174256_s2.hifi_reads.bc1009.bam |
| HG03486 | HG03486:untagged | HPRC_PLUS | m64043_200424_162541.ccs.bam<br>m64043_200425_223840.ccs.bam<br>m64043_200428_155222.ccs.bam<br>m64043_200429_220517.ccs.bam |
| HG03492 | HG03492_ELF2_ELF3_480c6f | HPRC_PLUS | m64136_200904_190830.ccs.bam<br>m64136_200906_012331.ccs.bam<br>m64136_200907_075143.ccs.bam |
| HG03516 | HG03516.HFSS | HPRC | m84046_231202_110908_s4.hifi_reads.bc2056.bam |
| HG03516 | HG03516_HiFiEx_mix | HPRC | m54329U_200610_234222.ccs.bam<br>m54329U_200612_200443.ccs.bam<br>m54329U_200614_021746.ccs.bam<br>m54329U_200615_084313.ccs.bam |
| HG03521 | PG03521.HFSS2 | HPRC | m84046_230715_060901_s1.hifi_reads.bc2084.bam<br>m84046_230715_064007_s2.hifi_reads.bc2084.bam |
| HG03540 | HG03540_lib1 | HPRC | m64043_200521_171703.ccs.bam<br>m64043_200522_232930.ccs.bam<br>m64043_200524_055430.ccs.bam |
| HG03579 | HG03579_lib1 | HPRC | m64043_200515_165406.ccs.bam<br>m64043_200516_230634.ccs.bam<br>m64043_200518_053124.ccs.bam |
| HG03583 | HG03583_PB1 | HPRC | m84091_230712_170232_s2.hifi_reads.bc1020.bam<br>m84091_230817_153352_s4.hifi_reads.bc1020.bam |
| HG03742 | HG03742_PB1 | HPRC | primrose/m64055e_220917_121508.5mc.hifi_reads.bam<br>primrose/m64055e_220924_080740.5mc.hifi_reads.bam |
| HG03784 | HG03784_lib1 | HPRC | m84081_230624_213331_s3.hifi_reads.bc2017.bam |
| HG03874 | HG03874_PB1 | HPRC | m84091_230712_163213_s1.hifi_reads.bc1012.bam<br>m84091_230731_130312_s1.hifi_reads.bc1012.bam |
| NA18505 | NA18505_PB1 | HPRC | m84091_230718_154039_s1.hifi_reads.bc1021.bam<br>m84091_230724_152950_s1.hifi_reads.bc1021.bam |
| NA18508 | PG18508.HFSS2 | HPRC | m84046_230513_215634_s1.hifi_reads.bc2001.bam |
| NA18522 | NA18522.HFSS | HPRC | m54329U_220604_002013-bc2021.5mc.hifi_reads.bam<br>m54329U_220605_194957-bc2021.5mc.hifi_reads.bam<br>m84046_230531_214232_s2.hifi_reads.bc2021.bam |
| NA18565 | NA18565_lib1 | HPRC | m84081_230714_201817_s1.hifi_reads.bc2020.bam<br>m84081_230728_191731_s1.hifi_reads.bc2020.bam |
| NA18570 | NA18570_PB1 | HPRC | primrose/m64330e_220821_125227.5mc.hifi_reads.bam<br>primrose/m64334e_220815_142753.5mc.hifi_reads.bam |
| NA18608 | NA18608_lib1 | HPRC | m84081_230601_210338_s2.hifi_reads.bc2007.bam |
| NA18620 | NA18620_PB1 | HPRC | m84091_230703_172326_s3.hifi_reads.bc1009.bam |
| NA18747 | PG18747_1.HFSS | HPRC | m54329U_220823_191353-bc2011.5mc.hifi_reads.bam<br>m64076_221001_041132-bc2011.5mc.hifi_reads.bam<br>m84046_230531_214232_s2.hifi_reads.bc2011.bam |
| NA18747 | PG18747_2.HFSS | HPRC | m64076_220826_143529-bc2049.5mc.hifi_reads.bam<br>m64076_220828_113336-bc2049.5mc.hifi_reads.bam |
| NA18879 | PG18879.HFSS | HPRC | m84046_230727_225904_s1.hifi_reads.bc2089.bam<br>m84046_230728_195023_s2.hifi_reads.bc2089.bam<br>m84046_230728_202047_s3.hifi_reads.bc2089.bam |
| NA18906 | NA18906:untagged | HPRC_PLUS | m64136_200521_171936.ccs.bam<br>m64136_200523_195722.ccs.bam<br>m64136_200525_021027.ccs.bam<br>m64136_200526_083627.ccs.bam |
| NA18952 | NA18952_lib1 | HPRC | m84081_230609_201402_s3.hifi_reads.bc2009.bam |
| NA18971 | NA18971_lib1 | HPRC | m64043_220701_181144-bc1011.5mc.hifi_reads.bam<br>m64043_220703_150718-bc1011.5mc.hifi_reads.bam<br>m64043_220705_120954-bc1011.5mc.hifi_reads.bam |
| NA18974 | NA18974_lib1 | HPRC | m84081_230601_213444_s3.hifi_reads.bc2008.bam |
| NA18976 | NA18976_lib1 | HPRC | m84081_230610_192218_s1.hifi_reads.bc2010.bam |
| NA18983 | PG18983.HFSS2 | HPRC | m64457e_220902_165459-bc2055.5mc.hifi_reads.bam<br>m64457e_220904_135104-bc2055.5mc.hifi_reads.bam |
| NA19036 | NA19036_PB1 | HPRC | m84091_230718_161058_s2.hifi_reads.bc1022.bam |
| NA19043 | PG19043.HFSS | HPRC | m54329U_220916_162216-bc2061.5mc.hifi_reads.bam<br>m54329U_220918_114340-bc2061.5mc.hifi_reads.bam<br>m84046_230519_231055_s2.hifi_reads.bc2061.bam |
| NA19087 | NA19087_PB1 | HPRC | primrose/m54306Ue_220902_172744.5mc.hifi_reads.bam<br>primrose/m64055e_220911_090142.5mc.hifi_reads.bam |
| NA19159 | NA19159_PB1 | HPRC | primrose/m54306Ue_220907_175549.5mc.hifi_reads.bam<br>primrose/m64330e_220902_175209.5mc.hifi_reads.bam<br>m84091_230817_150246_s3.hifi_reads.bc1015.bam |
| NA19185 | NA19185_PB1 | HPRC | primrose/m64055e_220902_174345.5mc.hifi_reads.bam<br>primrose/m64055e_220907_174608.5mc.hifi_reads.bam |
| NA19240 | NA19240:untagged | HPRC_PLUS | m64043_200128_181438.ccs.bam<br>m64043_200130_175216.ccs.bam<br>m64043_200201_000449.ccs.bam<br>m64043_200202_062937.ccs.bam |
| NA19338 | NA19338_PB1 | HPRC | primrose/m54306Ue_220722_175158.5mc.hifi_reads.bam<br>primrose/m64055e_220729_170145.5mc.hifi_reads.bam<br>primrose/m64055e_220731_124800.5mc.hifi_reads.bam |
| NA19338 | NA19338_PB2 | HPRC | m84091_230705_174256_s2.hifi_reads.bc1002.bam<br>m84091_230705_184508_s4.hifi_reads.bc1002.bam<br>m84091_231120_175800_s4.hifi_reads.bc1002.bam |
| NA19391 | NA19391_PB1 | HPRC | primrose/m54306Ue_220815_144205.5mc.hifi_reads.bam<br>primrose/m64330e_220819_173115.5mc.hifi_reads.bam<br>primrose/m64334e_220819_174259.5mc.hifi_reads.bam |
| NA19391 | NA19391_PB3 | HPRC | m84091_231120_175800_s4.hifi_reads.bc1015.bam |
| NA19443 | PG19443.HFSS | HPRC | m84046_230623_203534_s4.hifi_reads.bc2053.bam<br>m84046_230627_201017_s1.hifi_reads.bc2053.bam |
| NA19468 | NA19468_PB1 | HPRC | primrose/m54306Ue_220927_182256.demultiplex.bc1008--bc1008.5mc.hifi_reads.bam<br>primrose/m64055e_220817_094430.5mc.hifi_reads.bam<br>primrose/m64055e_220819_064153.5mc.hifi_reads.bam |
| NA19682 | PG19682.HFSS | HPRC | m84046_230428_231421_s3.hifi_reads.bc2090.bam<br>m84046_230504_195327_s1.hifi_reads.bc2090.bam |
| NA19700 | PG19700.HFSS | HPRC | m84046_230428_224315_s2.hifi_reads.bc2091.bam<br>m84046_230502_203455_s3.hifi_reads.bc2091.bam |
| NA19776 | PG19776.HFSS | HPRC | m84046_230513_215634_s1.hifi_reads.bc2014.bam |
| NA19835 | NA19835_PB1 | HPRC | m84091_230718_164204_s3.hifi_reads.bc1001.bam<br>m84091_230731_133332_s2.hifi_reads.bc1001.bam |
| NA19909 | NA19909_Fraction2_Fraction3 | HPRC | m84081_230523_180945_s1.hifi_reads.bc2005.bam |
| NA20129 | NA20129:untagged | HPRC_PLUS | m64043_191227_185626.ccs.bam<br>m64043_191229_010753.ccs.bam<br>m64043_191230_073311.ccs.bam<br>m64043_200111_140530.ccs.bam<br>m64043_200114_192155.ccs.bam |
| NA20282 | NA20282_lib1 | HPRC | m84081_230516_191559_s1.hifi_reads.bc2004.bam |
| NA20346 | NA20346_lib1 | HPRC | m84081_230523_184006_s2.hifi_reads.bc2006.bam |
| NA20752 | PG20752.HFSS | HPRC | m64457e_220827_152226-bc2052.5mc.hifi_reads.bam<br>m64457e_220829_104034-bc2052.5mc.hifi_reads.bam |
| NA20799 | NA20799_PB1 | HPRC | primrose/m54306Ue_220920_181216.5mc.hifi_reads.bam<br>primrose/m54306Ue_220922_132103.5mc.hifi_reads.bam<br>primrose/m54306Ue_220927_182256.demultiplex.bc1018--bc1018.5mc.hifi_reads.bam |
| NA20799 | NA20799_PB2 | HPRC | m84091_231208_204928_s2.hifi_reads.bc1018.bam<br>m84091_231211_215746_s3.hifi_reads.bc1018.bam |
| NA20805 | NA20805_lib1 | HPRC | m64136_220726_203708-bc1020.5mc.hifi_reads.bam<br>m64136_220728_160902-bc1020.5mc.hifi_reads.bam<br>m64136_220730_113921-bc1020.5mc.hifi_reads.bam |
| NA20809 | NA20809_lib1 | HPRC | m84081_230715_203436_s1.hifi_reads.bc2022.bam<br>m84081_230728_191731_s1.hifi_reads.bc2022.bam |
| NA20850 | PG20850.HFSS | HPRC | m84046_230721_175304_s3.hifi_reads.bc2092.bam |
| NA20870 | NA20870_PB1 | HPRC | m84091_230703_175432_s4.hifi_reads.bc1010.bam |
| NA20905 | PG20905.HFSS | HPRC | m64076_220914_214248-bc2060.5mc.hifi_reads.bam<br>m64076_220916_183810-bc2060.5mc.hifi_reads.bam<br>m84046_230522_234602_s2.hifi_reads.bc2060.bam |
| NA21093 | PG21093.HFSS | HPRC | m84046_230428_231421_s3.hifi_reads.bc2094.bam<br>m84046_230504_195327_s1.hifi_reads.bc2094.bam |
| NA21102 | NA21102_lib1 | HPRC | m84081_230715_213648_s3.hifi_reads.bc2023.bam<br>m84081_230728_191731_s1.hifi_reads.bc2023.bam |
| NA21106 | PG21106.HFSS | HPRC | m84046_230501_151821_s4.hifi_reads.bc2089.bam<br>m84046_230504_205453_s3.hifi_reads.bc2089.bam |
| NA21110 | NA21110_PB1 | HPRC | m84091_230710_165908_s1.hifi_reads.bc1011.bam<br>m84091_230719_165228_s3.hifi_reads.bc1011.bam |
| NA21144 | NA21144_lib1 | HPRC | m84081_230629_184915_s1.hifi_reads.bc2018.bam |
| NA21309 | NA21309:untagged | HPRC_PLUS | m64043_191210_201113.ccs.bam<br>m64043_191213_191857.ccs.bam<br>m64043_191215_014401.ccs.bam<br>m64043_191219_192900.ccs.bam |

