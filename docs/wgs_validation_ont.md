# HLA-Resolve concordance on whole-genome ONT libraries

HLA-Resolve v0.12.0 (IPD-IMGT/HLA release 3.64.0) was run on public whole-genome ONT PromethION R10.4.1 reads for **72 sequencing libraries from 72 samples**. All libraries were sequenced by the Human Pangenome Reference Consortium (HPRC; cohorts HPP, HPRC, HPRC_PLUS). Every library was typed independently. Each sample had 43 to 99× of R10.4.1 reads. Flow cells were left out to bring each library to about 30× (21 to 39×), so that accuracy reflects typical WGS depth.

Allele calls were compared to two reference sets:

- **Lai et al. 2024** ([doi.org/10.1016/j.csbj.2024.03.030](https://doi.org/10.1016/j.csbj.2024.03.030); [Supplementary File 6](Lai_Supplementary-6.xlsx))
  - 6 HPRC release 1 samples, 6 libraries
  - HLA-A, -B, -C, -DPA1, -DPB1, -DQA1, -DQB1, -DRB1
  - Resolution: up through four-field (full) resolution
- **1000 Genomes Project HLA panel** ([20181129_HLA_types_full_1000_Genomes_Project_panel.txt](https://ftp.1000genomes.ebi.ac.uk/vol1/ftp/data_collections/HLA_types/20181129_HLA_types_full_1000_Genomes_Project_panel.txt))
  - 66 samples, 66 libraries
  - HLA-A, -B, -C, -DQB1, -DRB1
  - Resolution: two-field resolution
  - Exome-based (IPD-IMGT/HLA 3.28)

The following command was used for each sample:

```bash
hla_resolve --input_file <LIBRARY_uBAM> --sample_name <SAMPLE> \
            --platform ont --scheme WGS \
            --output_dir <OUTPUT_DIR> --threads <THREADS>
```

Libraries basecalled with Dorado sup v4.1.0 were run with `--clair3_model r1041_e82_400bps_sup_v410`, and libraries basecalled with v4.3.0 with `r1041_e82_400bps_sup_v500`, the default.

## Coverage depth

Mean depth across the antigen recognition site of each gene, over all 72 libraries.

| Gene | Mean | Median | Range |
|------|------|--------|-------|
| HLA-A | 30.2× | 29.4× | 15.9 to 50.6× |
| HLA-B | 31.0× | 30.5× | 16.8 to 53.8× |
| HLA-C | 30.9× | 30.9× | 12.8 to 54.6× |
| HLA-DPA1 | 31.5× | 31.6× | 17.6 to 50.4× |
| HLA-DPB1 | 31.5× | 30.6× | 16.0 to 50.1× |
| HLA-DQA1 | 29.0× | 27.8× | 15.5 to 47.6× |
| HLA-DQB1 | 30.5× | 29.8× | 13.7 to 51.7× |
| HLA-DRB1 | 23.7× | 23.8× | 11.7 to 48.1× |

## Scoring

Concordance is reported among alleles called. An allele was evaluated at a given field resolution only if the reference specified it to that resolution. A call resolved to fewer fields than the reference is considered discordant at that level of resolution. By default, genes do not type if their antigen recognition site coverage depth sits below 8x, so call rate was measured as the proportion of total possible allele calls emitted by HLA-Resolve.

## Lai et al. reference (6 libraries, 6 samples)

| Metric | Value |
|---|---:|
| Total possible allele calls | 96 |
| Alleles called | 96 |
| Call rate | 100.0% |

| Resolution | Concordance | Concordant alleles | Alleles evaluated |
|---|---:|---:|---:|
| 1-field | 100.0% | 96 | 96 |
| 2-field | 100.0% | 96 | 96 |
| 3-field | 100.0% | 96 | 96 |
| 4-field | 84.9% | 79 | 93 |

Alleles evaluated is lower than alleles called because not all reference alleles have third and fourth field values.

## 1000 Genomes panel (66 libraries, 66 samples)

| Metric | Value |
|---|---:|
| Total possible allele calls | 1056 |
| Alleles called | 1056 |
| Call rate | 100.0% |

| Resolution | Concordance | Concordant alleles | Alleles evaluated |
|---|---:|---:|---:|
| 1-field | 100.0% | 632 | 632 |
| 2-field | 99.4% | 628 | 632 |

Alleles evaluated is lower than alleles called because the panel covers five of the eight genes, and the reference has missing HLA-DQB1 calls for 14 samples (28 alleles).

## All libraries

Both reference sets pooled, 72 libraries. The number of alleles evaluated falls with resolution because the 1000 Genomes panel stops at two fields and covers five genes, and because a reference allele is evaluated only at the fields it specifies.

| Metric | Value |
|---|---:|
| Total possible allele calls | 1152 |
| Alleles called | 1152 |
| Call rate | 100.0% |

| Resolution | Concordance | Concordant alleles | Alleles evaluated |
|---|---:|---:|---:|
| 1-field | 100.0% | 728 | 728 |
| 2-field | 99.5% | 724 | 728 |
| 3-field | 100.0% | 96 | 96 |
| 4-field | 84.9% | 79 | 93 |

## Results by library

One row per sequencing library, grouped by sample.

- The 3- and 4-field columns are blank where the reference is two-field.
- HLA coverage is the mean depth across the eight classical HLA genes from the HLA-Resolve run log (all 72 libraries: mean 29.5×, median 28.0×, range 19.0–42.4×).

| Sample | Library | Cohort | Instrument | Reference | HLA&nbsp;coverage | Genes&nbsp;typed | Alleles&nbsp;called | 1-field | 2-field | 3-field | 4-field |
|---|---|---|---|---|---:|---|---|---|---|---|---|
| HG00097 | 10_31_23_R1041_HPRC_HG00097_1 | HPRC | PromethION | 1kGP | 37.3× | 8/8 | 16/16 | 10/10 | 10/10 |  |  |
| HG00126 | 08_29_23_R1041_HPRC_HG00126_1 | HPRC | PromethION | 1kGP | 26.1× | 8/8 | 16/16 | 10/10 | 10/10 |  |  |
| HG00128 | 11_7_23_R1041_HPRC_HG00128_1 | HPRC | PromethION | 1kGP | 37.8× | 8/8 | 16/16 | 10/10 | 10/10 |  |  |
| HG00133 | 08_29_23_R1041_HPRC_HG00133_1 | HPRC | PromethION | 1kGP | 19.6× | 8/8 | 16/16 | 10/10 | 10/10 |  |  |
| HG00146 | 01_23_24_R1041_HPRC_HG00146_1 | HPRC | PromethION | 1kGP | 33.9× | 8/8 | 16/16 | 10/10 | 10/10 |  |  |
| HG00232 | 01_23_24_R1041_HPRC_HG00232_1 | HPRC | PromethION | 1kGP | 39.4× | 8/8 | 16/16 | 8/8 | 8/8 |  |  |
| HG00235 | 01_16_24_R1041_HPRC_HG00235_1 | HPRC | PromethION | 1kGP | 35.1× | 8/8 | 16/16 | 10/10 | 10/10 |  |  |
| HG00253 | 01_16_24_R1041_HPRC_HG00253_1 | HPRC | PromethION | 1kGP | 27.2× | 8/8 | 16/16 | 8/8 | 8/8 |  |  |
| HG00272 | 08_08_23_R1041_HPRC_HG00272_1 | HPRC | PromethION | 1kGP | 22.7× | 8/8 | 16/16 | 10/10 | 10/10 |  |  |
| HG00290 | 01_23_24_R1041_HPRC_HG00290_1 | HPRC | PromethION | 1kGP | 29.3× | 8/8 | 16/16 | 10/10 | 10/10 |  |  |
| HG00320 | 12_12_23_R1041_HPRC_HG00320_1 | HPRC | PromethION | 1kGP | 22.3× | 8/8 | 16/16 | 8/8 | 8/8 |  |  |
| HG00321 | 01_16_24_R1041_HPRC_HG00321_1 | HPRC | PromethION | 1kGP | 25.9× | 8/8 | 16/16 | 8/8 | 7/8 |  |  |
| HG00329 | 10_24_23_R1041_HPRC_HG00329_1 | HPRC | PromethION | 1kGP | 25.8× | 8/8 | 16/16 | 8/8 | 8/8 |  |  |
| HG00344 | 10_31_23_R1041_HPRC_HG00344_1 | HPRC | PromethION | 1kGP | 36.5× | 8/8 | 16/16 | 8/8 | 8/8 |  |  |
| HG00350 | 10_24_23_R1041_HPRC_HG00350_1 | HPRC | PromethION | 1kGP | 22.3× | 8/8 | 16/16 | 10/10 | 10/10 |  |  |
| HG01167 | 08_29_23_R1041_HPRC_HG01167_1 | HPRC | PromethION | 1kGP | 23.2× | 8/8 | 16/16 | 10/10 | 10/10 |  |  |
| HG01786 | 11_28_23_R1041_HPRC_HG01786_1 | HPRC | PromethION | 1kGP | 25.5× | 8/8 | 16/16 | 10/10 | 10/10 |  |  |
| HG02178 | 12_12_23_R1041_HPRC_HG02178_1 | HPRC | PromethION | 1kGP | 39.4× | 8/8 | 16/16 | 10/10 | 10/10 |  |  |
| HG02392 | 10_31_23_R1041_HPRC_HG02392_1 | HPRC | PromethION | 1kGP | 29.3× | 8/8 | 16/16 | 10/10 | 10/10 |  |  |
| HG02583 | 01_16_24_R1041_HPRC_HG02583_1 | HPRC | PromethION | 1kGP | 27.9× | 8/8 | 16/16 | 10/10 | 10/10 |  |  |
| HG02818 | 02_22_23_R1041_HG02818_1_dorado0.7.2_sup4.1.0_5mCG_5hmCG | HPRC_PLUS | PromethION | Lai | 35.2× | 8/8 | 16/16 | 16/16 | 16/16 | 16/16 | 12/16 |
| HG03270 | 10_24_23_R1041_HPRC_HG03270_2 | HPRC | PromethION | 1kGP | 21.9× | 8/8 | 16/16 | 10/10 | 10/10 |  |  |
| HG03369 | 10_10_23_R1041_HPRC_HG03369_1 | HPRC | PromethION | 1kGP | 35.3× | 8/8 | 16/16 | 10/10 | 10/10 |  |  |
| HG03486 | 02_15_23_R1041_HG03486_1_dorado0.7.2_sup4.1.0_5mCG_5hmCG | HPRC_PLUS | PromethION | Lai | 27.6× | 8/8 | 16/16 | 16/16 | 16/16 | 16/16 | 15/16 |
| HG03521 | 10_24_23_R1041_HPRC_HG03521_1 | HPRC | PromethION | 1kGP | 29.1× | 8/8 | 16/16 | 10/10 | 10/10 |  |  |
| HG03583 | 10_31_23_R1041_HPRC_HG03583_1 | HPRC | PromethION | 1kGP | 22.1× | 8/8 | 16/16 | 10/10 | 10/10 |  |  |
| HG03784 | 10_17_23_R1041_HPRC_HG03784_1 | HPRC | PromethION | 1kGP | 39.1× | 8/8 | 16/16 | 10/10 | 10/10 |  |  |
| HG03874 | 01_23_24_R1041_HPRC_HG03874_1 | HPRC | PromethION | 1kGP | 28.0× | 8/8 | 16/16 | 10/10 | 10/10 |  |  |
| NA18505 | 10_24_23_R1041_HPRC_NA18505_1 | HPRC | PromethION | 1kGP | 27.7× | 8/8 | 16/16 | 10/10 | 10/10 |  |  |
| NA18508 | 09_12_23_R1041_HPRC_NA18508_1 | HPRC | PromethION | 1kGP | 26.9× | 8/8 | 16/16 | 10/10 | 9/10 |  |  |
| NA18565 | 10_10_23_R1041_HPRC_NA18565_5 | HPRC | PromethION | 1kGP | 22.8× | 8/8 | 16/16 | 10/10 | 10/10 |  |  |
| NA18608 | 10_31_23_R1041_HPRC_NA18608_1 | HPRC | PromethION | 1kGP | 28.9× | 8/8 | 16/16 | 8/8 | 7/8 |  |  |
| NA18620 | 01_16_24_R1041_HPRC_NA18620_1 | HPRC | PromethION | 1kGP | 26.4× | 8/8 | 16/16 | 8/8 | 7/8 |  |  |
| NA18879 | 08_08_23_R1041_HPRC_NA18879_1 | HPRC | PromethION | 1kGP | 40.8× | 8/8 | 16/16 | 10/10 | 10/10 |  |  |
| NA18906 | 02_15_23_R1041_NA18906_2_dorado0.7.2_sup4.1.0_5mCG_5hmCG | HPRC_PLUS | PromethION | Lai | 19.5× | 8/8 | 16/16 | 16/16 | 16/16 | 16/16 | 12/15 |
| NA18940 | 11_28_23_R1041_AMED-Collaboration-HPRC_GM18940_1 | HPP | PromethION | 1kGP | 19.2× | 8/8 | 16/16 | 8/8 | 8/8 |  |  |
| NA18943 | 11_28_23_R1041_AMED-Collaboration-HPRC_GM18943_3 | HPP | PromethION | 1kGP | 37.2× | 8/8 | 16/16 | 8/8 | 8/8 |  |  |
| NA18944 | 11_28_23_R1041_AMED-Collaboration-HPRC_GM18944_1 | HPP | PromethION | 1kGP | 30.9× | 8/8 | 16/16 | 8/8 | 8/8 |  |  |
| NA18945 | 12_05_23_R1041_AMED-Collaboration-HPRC_GM18945_1 | HPP | PromethION | 1kGP | 33.7× | 8/8 | 16/16 | 8/8 | 8/8 |  |  |
| NA18948 | 01_09_24_R1041_AMED-Collaboration-HPRC_GM18948_2 | HPP | PromethION | 1kGP | 27.9× | 8/8 | 16/16 | 8/8 | 8/8 |  |  |
| NA18952 | 11_28_23_R1041_HPRC_NA18952_1 | HPRC | PromethION | 1kGP | 30.0× | 8/8 | 16/16 | 8/8 | 8/8 |  |  |
| NA18959 | 01_09_24_R1041_AMED-Collaboration-HPRC_GM18959_1 | HPP | PromethION | 1kGP | 19.0× | 8/8 | 16/16 | 10/10 | 10/10 |  |  |
| NA18960 | 01_09_24_R1041_AMED-Collaboration-HPRC_GM18960_1 | HPP | PromethION | 1kGP | 26.4× | 8/8 | 16/16 | 10/10 | 10/10 |  |  |
| NA18967 | 01_09_24_R1041_AMED-Collaboration-HPRC_GM18967_1 | HPP | PromethION | 1kGP | 31.0× | 8/8 | 16/16 | 10/10 | 10/10 |  |  |
| NA18970 | 01_09_24_R1041_AMED-Collaboration-HPRC_GM18970_1 | HPP | PromethION | 1kGP | 27.8× | 8/8 | 16/16 | 10/10 | 10/10 |  |  |
| NA18974 | 11_28_23_R1041_HPRC_NA18974_1 | HPRC | PromethION | 1kGP | 41.6× | 8/8 | 16/16 | 10/10 | 10/10 |  |  |
| NA18976 | 12_12_23_R1041_HPRC_NA18976_1 | HPRC | PromethION | 1kGP | 27.1× | 8/8 | 16/16 | 10/10 | 10/10 |  |  |
| NA18982 | 01_09_24_R1041_AMED-Collaboration-HPRC_GM18982_1 | HPP | PromethION | 1kGP | 32.9× | 8/8 | 16/16 | 10/10 | 10/10 |  |  |
| NA19036 | 01_23_24_R1041_HPRC_NA19036_1 | HPRC | PromethION | 1kGP | 42.4× | 8/8 | 16/16 | 10/10 | 10/10 |  |  |
| NA19240 | 03_14_23_R1041_NA19240_1_dorado0.7.2_sup4.1.0_5mCG_5hmCG | HPRC_PLUS | PromethION | Lai | 24.8× | 8/8 | 16/16 | 16/16 | 16/16 | 16/16 | 14/16 |
| NA19443 | 01_23_24_R1041_HPRC_NA19443_1 | HPRC | PromethION | 1kGP | 41.9× | 8/8 | 16/16 | 10/10 | 10/10 |  |  |
| NA19682 | 11_28_23_R1041_HPRC_NA19682_1 | HPRC | PromethION | 1kGP | 32.5× | 8/8 | 16/16 | 10/10 | 10/10 |  |  |
| NA19700 | 01_23_24_R1041_HPRC_NA19700_1 | HPRC | PromethION | 1kGP | 28.6× | 8/8 | 16/16 | 10/10 | 10/10 |  |  |
| NA19776 | 10_24_23_R1041_HPRC_NA19776_1 | HPRC | PromethION | 1kGP | 32.9× | 8/8 | 16/16 | 10/10 | 10/10 |  |  |
| NA19835 | 10_31_23_R1041_HPRC_NA19835_1 | HPRC | PromethION | 1kGP | 34.2× | 8/8 | 16/16 | 10/10 | 10/10 |  |  |
| NA19909 | 11_14_23_R1041_HPRC_NA19909_1 | HPRC | PromethION | 1kGP | 27.1× | 8/8 | 16/16 | 10/10 | 10/10 |  |  |
| NA20129 | 03_08_23_R1041_NA20129_1_dorado0.7.2_sup4.1.0_5mCG_5hmCG | HPRC_PLUS | PromethION | Lai | 27.0× | 8/8 | 16/16 | 16/16 | 16/16 | 16/16 | 14/15 |
| NA20282 | 12_12_23_R1041_HPRC_NA20282_1 | HPRC | PromethION | 1kGP | 28.0× | 8/8 | 16/16 | 10/10 | 10/10 |  |  |
| NA20346 | 11_28_23_R1041_HPRC_NA20346_2 | HPRC | PromethION | 1kGP | 29.2× | 8/8 | 16/16 | 10/10 | 10/10 |  |  |
| NA20503 | 12_12_23_R1041_HPRC_NA20503_1 | HPRC | PromethION | 1kGP | 22.4× | 8/8 | 16/16 | 10/10 | 10/10 |  |  |
| NA20762 | 01_23_24_R1041_HPRC_NA20762_1 | HPRC | PromethION | 1kGP | 32.7× | 8/8 | 16/16 | 10/10 | 10/10 |  |  |
| NA20806 | 12_12_23_R1041_HPRC_NA20806_1 | HPRC | PromethION | 1kGP | 31.2× | 8/8 | 16/16 | 10/10 | 10/10 |  |  |
| NA20809 | 08_08_23_R1041_HPRC_NA20809_1 | HPRC | PromethION | 1kGP | 25.5× | 8/8 | 16/16 | 10/10 | 10/10 |  |  |
| NA20827 | 01_23_24_R1041_HPRC_NA20827_1 | HPRC | PromethION | 1kGP | 35.1× | 8/8 | 16/16 | 10/10 | 10/10 |  |  |
| NA20850 | 08_08_23_R1041_HPRC_NA20850_1 | HPRC | PromethION | 1kGP | 20.4× | 8/8 | 16/16 | 10/10 | 10/10 |  |  |
| NA20870 | 11_14_23_R1041_HPRC_NA20870_1 | HPRC | PromethION | 1kGP | 27.9× | 8/8 | 16/16 | 10/10 | 10/10 |  |  |
| NA21093 | 01_23_24_R1041_HPRC_NA21093_1 | HPRC | PromethION | 1kGP | 27.4× | 8/8 | 16/16 | 10/10 | 10/10 |  |  |
| NA21102 | 01_23_24_R1041_HPRC_NA21102_1 | HPRC | PromethION | 1kGP | 34.1× | 8/8 | 16/16 | 10/10 | 10/10 |  |  |
| NA21106 | 10_24_23_R1041_HPRC_NA21106_1 | HPRC | PromethION | 1kGP | 29.0× | 8/8 | 16/16 | 10/10 | 10/10 |  |  |
| NA21110 | 01_16_24_R1041_HPRC_NA21110_1 | HPRC | PromethION | 1kGP | 25.9× | 8/8 | 16/16 | 10/10 | 10/10 |  |  |
| NA21144 | 01_23_24_R1041_HPRC_NA21144_1 | HPRC | PromethION | 1kGP | 27.2× | 8/8 | 16/16 | 10/10 | 10/10 |  |  |
| NA21309 | 02_22_23_R1041_NA21309_1_dorado0.7.2_sup4.1.0_5mCG_5hmCG | HPRC_PLUS | PromethION | Lai | 31.2× | 8/8 | 16/16 | 16/16 | 16/16 | 16/16 | 12/15 |

## Input files

Each library was typed from one unaligned BAM made by concatenating the files listed. HPRC files are public on the AWS Open Data registry and need no credential. Build a path by substituting the Cohort, Sample and Files values from the table into:

```
s3://human-pangenomics/working/<COHORT>/<SAMPLE>/raw_data/nanopore/<FILE>
```

For example, the first file listed for HG00097 (cohort HPRC) is at:

```
s3://human-pangenomics/working/HPRC/HG00097/raw_data/nanopore/dorado0.6.0_sup4.3.0_5mCG_5hmCG/10_31_23_R1041_HPRC_HG00097_1_dorado0.6.0_sup4.3.0_5mCG_5hmCG.bam
```

Files are public and need no credentials:

```bash
aws s3 cp --no-sign-request \
  s3://human-pangenomics/working/HPRC/HG00097/raw_data/nanopore/dorado0.6.0_sup4.3.0_5mCG_5hmCG/10_31_23_R1041_HPRC_HG00097_1_dorado0.6.0_sup4.3.0_5mCG_5hmCG.bam .
```

Each file sits in a folder named after its basecaller version, which the Files column includes.

| Sample | Library | Cohort | Files |
|---|---|---|---|
| HG00097 | 10_31_23_R1041_HPRC_HG00097_1 | HPRC | dorado0.6.0_sup4.3.0_5mCG_5hmCG/10_31_23_R1041_HPRC_HG00097_1_dorado0.6.0_sup4.3.0_5mCG_5hmCG.bam<br>dorado0.6.0_sup4.3.0_5mCG_5hmCG/10_31_23_R1041_HPRC_HG00097_2_dorado0.6.0_sup4.3.0_5mCG_5hmCG.bam |
| HG00126 | 08_29_23_R1041_HPRC_HG00126_1 | HPRC | dorado0.6.0_sup4.3.0_5mCG_5hmCG/08_29_23_R1041_HPRC_HG00126_1_dorado0.6.0_sup4.3.0_5mCG_5hmCG.bam |
| HG00128 | 11_7_23_R1041_HPRC_HG00128_1 | HPRC | dorado0.6.0_sup4.3.0_5mCG_5hmCG/11_7_23_R1041_HPRC_HG00128_1_dorado0.6.0_sup4.3.0_5mCG_5hmCG.bam<br>dorado0.6.0_sup4.3.0_5mCG_5hmCG/11_7_23_R1041_HPRC_HG00128_2_dorado0.6.0_sup4.3.0_5mCG_5hmCG.bam |
| HG00133 | 08_29_23_R1041_HPRC_HG00133_1 | HPRC | dorado0.6.0_sup4.3.0_5mCG_5hmCG/08_29_23_R1041_HPRC_HG00133_1_dorado0.6.0_sup4.3.0_5mCG_5hmCG.bam |
| HG00146 | 01_23_24_R1041_HPRC_HG00146_1 | HPRC | dorado0.6.0_sup4.3.0_5mCG_5hmCG/01_23_24_R1041_HPRC_HG00146_1_dorado0.6.0_sup4.3.0_5mCG_5hmCG.bam<br>dorado0.6.0_sup4.3.0_5mCG_5hmCG/01_23_24_R1041_HPRC_HG00146_2_dorado0.6.0_sup4.3.0_5mCG_5hmCG.bam<br>dorado0.6.0_sup4.3.0_5mCG_5hmCG/01_23_24_R1041_HPRC_HG00146_3_dorado0.6.0_sup4.3.0_5mCG_5hmCG.bam |
| HG00232 | 01_23_24_R1041_HPRC_HG00232_1 | HPRC | dorado0.6.0_sup4.3.0_5mCG_5hmCG/01_23_24_R1041_HPRC_HG00232_1_dorado0.6.0_sup4.3.0_5mCG_5hmCG.bam<br>dorado0.6.0_sup4.3.0_5mCG_5hmCG/01_23_24_R1041_HPRC_HG00232_2_dorado0.6.0_sup4.3.0_5mCG_5hmCG.bam |
| HG00235 | 01_16_24_R1041_HPRC_HG00235_1 | HPRC | dorado0.6.0_sup4.3.0_5mCG_5hmCG/01_16_24_R1041_HPRC_HG00235_1_dorado0.6.0_sup4.3.0_5mCG_5hmCG.bam<br>dorado0.6.0_sup4.3.0_5mCG_5hmCG/01_16_24_R1041_HPRC_HG00235_3_dorado0.6.0_sup4.3.0_5mCG_5hmCG.bam |
| HG00253 | 01_16_24_R1041_HPRC_HG00253_1 | HPRC | dorado0.6.0_sup4.3.0_5mCG_5hmCG/01_16_24_R1041_HPRC_HG00253_1_dorado0.6.0_sup4.3.0_5mCG_5hmCG.bam<br>dorado0.6.0_sup4.3.0_5mCG_5hmCG/10_31_23_R1041_HPRC_HG00253_1_dorado0.6.0_sup4.3.0_5mCG_5hmCG.bam<br>dorado0.6.0_sup4.3.0_5mCG_5hmCG/10_31_23_R1041_HPRC_HG00253_2_dorado0.6.0_sup4.3.0_5mCG_5hmCG.bam<br>dorado0.6.0_sup4.3.0_5mCG_5hmCG/10_31_23_R1041_HPRC_HG00253_3_dorado0.6.0_sup4.3.0_5mCG_5hmCG.bam |
| HG00272 | 08_08_23_R1041_HPRC_HG00272_1 | HPRC | dorado0.6.0_sup4.3.0_5mCG_5hmCG/08_08_23_R1041_HPRC_HG00272_1_dorado0.6.0_sup4.3.0_5mCG_5hmCG.bam |
| HG00290 | 01_23_24_R1041_HPRC_HG00290_1 | HPRC | dorado0.6.0_sup4.3.0_5mCG_5hmCG/01_23_24_R1041_HPRC_HG00290_1_dorado0.6.0_sup4.3.0_5mCG_5hmCG.bam |
| HG00320 | 12_12_23_R1041_HPRC_HG00320_1 | HPRC | dorado0.6.0_sup4.3.0_5mCG_5hmCG/12_12_23_R1041_HPRC_HG00320_1_dorado0.6.0_sup4.3.0_5mCG_5hmCG.bam |
| HG00321 | 01_16_24_R1041_HPRC_HG00321_1 | HPRC | dorado0.6.0_sup4.3.0_5mCG_5hmCG/01_16_24_R1041_HPRC_HG00321_1_dorado0.6.0_sup4.3.0_5mCG_5hmCG.bam<br>dorado0.6.0_sup4.3.0_5mCG_5hmCG/11_7_23_R1041_HPRC_HG00321_1_dorado0.6.0_sup4.3.0_5mCG_5hmCG.bam<br>dorado0.6.0_sup4.3.0_5mCG_5hmCG/11_7_23_R1041_HPRC_HG00321_2_dorado0.6.0_sup4.3.0_5mCG_5hmCG.bam<br>dorado0.6.0_sup4.3.0_5mCG_5hmCG/11_7_23_R1041_HPRC_HG00321_3_dorado0.6.0_sup4.3.0_5mCG_5hmCG.bam |
| HG00329 | 10_24_23_R1041_HPRC_HG00329_1 | HPRC | dorado0.6.0_sup4.3.0_5mCG_5hmCG/10_24_23_R1041_HPRC_HG00329_1_dorado0.6.0_sup4.3.0_5mCG_5hmCG.bam |
| HG00344 | 10_31_23_R1041_HPRC_HG00344_1 | HPRC | dorado0.6.0_sup4.3.0_5mCG_5hmCG/10_31_23_R1041_HPRC_HG00344_1_dorado0.6.0_sup4.3.0_5mCG_5hmCG.bam<br>dorado0.6.0_sup4.3.0_5mCG_5hmCG/10_31_23_R1041_HPRC_HG00344_2_dorado0.6.0_sup4.3.0_5mCG_5hmCG.bam |
| HG00350 | 10_24_23_R1041_HPRC_HG00350_1 | HPRC | dorado0.6.0_sup4.3.0_5mCG_5hmCG/10_24_23_R1041_HPRC_HG00350_1_dorado0.6.0_sup4.3.0_5mCG_5hmCG.bam |
| HG01167 | 08_29_23_R1041_HPRC_HG01167_1 | HPRC | dorado0.6.0_sup4.3.0_5mCG_5hmCG/08_29_23_R1041_HPRC_HG01167_1_dorado0.6.0_sup4.3.0_5mCG_5hmCG.bam |
| HG01786 | 11_28_23_R1041_HPRC_HG01786_1 | HPRC | dorado0.6.0_sup4.3.0_5mCG_5hmCG/11_28_23_R1041_HPRC_HG01786_1_dorado0.6.0_sup4.3.0_5mCG_5hmCG.bam |
| HG02178 | 12_12_23_R1041_HPRC_HG02178_1 | HPRC | dorado0.6.0_sup4.3.0_5mCG_5hmCG/12_12_23_R1041_HPRC_HG02178_1_dorado0.6.0_sup4.3.0_5mCG_5hmCG.bam<br>dorado0.6.0_sup4.3.0_5mCG_5hmCG/12_12_23_R1041_HPRC_HG02178_2_dorado0.6.0_sup4.3.0_5mCG_5hmCG.bam |
| HG02392 | 10_31_23_R1041_HPRC_HG02392_1 | HPRC | dorado0.6.0_sup4.3.0_5mCG_5hmCG/10_31_23_R1041_HPRC_HG02392_1_dorado0.6.0_sup4.3.0_5mCG_5hmCG.bam |
| HG02583 | 01_16_24_R1041_HPRC_HG02583_1 | HPRC | dorado0.6.0_sup4.3.0_5mCG_5hmCG/01_16_24_R1041_HPRC_HG02583_1_dorado0.6.0_sup4.3.0_5mCG_5hmCG.bam<br>dorado0.6.0_sup4.3.0_5mCG_5hmCG/11_7_23_R1041_HPRC_HG02583_1_dorado0.6.0_sup4.3.0_5mCG_5hmCG.bam<br>dorado0.6.0_sup4.3.0_5mCG_5hmCG/11_7_23_R1041_HPRC_HG02583_2_dorado0.6.0_sup4.3.0_5mCG_5hmCG.bam<br>dorado0.6.0_sup4.3.0_5mCG_5hmCG/11_7_23_R1041_HPRC_HG02583_3_dorado0.6.0_sup4.3.0_5mCG_5hmCG.bam |
| HG02818 | 02_22_23_R1041_HG02818_1_dorado0.7.2_sup4.1.0_5mCG_5hmCG | HPRC_PLUS | dorado0.7.2_sup4.1.0_5mCG_5hmCG/02_22_23_R1041_HG02818_1_dorado0.7.2_sup4.1.0_5mCG_5hmCG.bam<br>dorado0.7.2_sup4.1.0_5mCG_5hmCG/02_22_23_R1041_HG02818_2_dorado0.7.2_sup4.1.0_5mCG_5hmCG.bam |
| HG03270 | 10_24_23_R1041_HPRC_HG03270_2 | HPRC | dorado0.6.0_sup4.3.0_5mCG_5hmCG/10_24_23_R1041_HPRC_HG03270_2_dorado0.6.0_sup4.3.0_5mCG_5hmCG.bam |
| HG03369 | 10_10_23_R1041_HPRC_HG03369_1 | HPRC | dorado0.6.0_sup4.3.0_5mCG_5hmCG/10_10_23_R1041_HPRC_HG03369_1_dorado0.6.0_sup4.3.0_5mCG_5hmCG.bam<br>dorado0.6.0_sup4.3.0_5mCG_5hmCG/10_10_23_R1041_HPRC_HG03369_2_dorado0.6.0_sup4.3.0_5mCG_5hmCG.bam |
| HG03486 | 02_15_23_R1041_HG03486_1_dorado0.7.2_sup4.1.0_5mCG_5hmCG | HPRC_PLUS | dorado0.7.2_sup4.1.0_5mCG_5hmCG/02_15_23_R1041_HG03486_1_dorado0.7.2_sup4.1.0_5mCG_5hmCG.bam<br>dorado0.7.2_sup4.1.0_5mCG_5hmCG/02_15_23_R1041_HG03486_2_dorado0.7.2_sup4.1.0_5mCG_5hmCG.bam |
| HG03521 | 10_24_23_R1041_HPRC_HG03521_1 | HPRC | dorado0.6.0_sup4.3.0_5mCG_5hmCG/10_24_23_R1041_HPRC_HG03521_1_dorado0.6.0_sup4.3.0_5mCG_5hmCG.bam |
| HG03583 | 10_31_23_R1041_HPRC_HG03583_1 | HPRC | dorado0.6.0_sup4.3.0_5mCG_5hmCG/10_31_23_R1041_HPRC_HG03583_1_dorado0.6.0_sup4.3.0_5mCG_5hmCG.bam<br>dorado0.6.0_sup4.3.0_5mCG_5hmCG/10_31_23_R1041_HPRC_HG03583_2_dorado0.6.0_sup4.3.0_5mCG_5hmCG.bam |
| HG03784 | 10_17_23_R1041_HPRC_HG03784_1 | HPRC | dorado0.6.0_sup4.3.0_5mCG_5hmCG/10_17_23_R1041_HPRC_HG03784_1_dorado0.6.0_sup4.3.0_5mCG_5hmCG.bam<br>dorado0.6.0_sup4.3.0_5mCG_5hmCG/10_17_23_R1041_HPRC_HG03784_3_dorado0.6.0_sup4.3.0_5mCG_5hmCG.bam |
| HG03874 | 01_23_24_R1041_HPRC_HG03874_1 | HPRC | dorado0.6.0_sup4.3.0_5mCG_5hmCG/01_23_24_R1041_HPRC_HG03874_1_dorado0.6.0_sup4.3.0_5mCG_5hmCG.bam<br>dorado0.6.0_sup4.3.0_5mCG_5hmCG/11_14_23_R1041_HPRC_HG03874_1_dorado0.6.0_sup4.3.0_5mCG_5hmCG.bam |
| NA18505 | 10_24_23_R1041_HPRC_NA18505_1 | HPRC | dorado0.6.0_sup4.3.0_5mCG_5hmCG/10_24_23_R1041_HPRC_NA18505_1_dorado0.6.0_sup4.3.0_5mCG_5hmCG.bam<br>dorado0.6.0_sup4.3.0_5mCG_5hmCG/10_24_23_R1041_HPRC_NA18505_4_dorado0.6.0_sup4.3.0_5mCG_5hmCG.bam |
| NA18508 | 09_12_23_R1041_HPRC_NA18508_1 | HPRC | dorado0.6.0_sup4.3.0_5mCG_5hmCG/09_12_23_R1041_HPRC_NA18508_1_dorado0.6.0_sup4.3.0_5mCG_5hmCG.bam |
| NA18565 | 10_10_23_R1041_HPRC_NA18565_5 | HPRC | dorado0.6.0_sup4.3.0_5mCG_5hmCG/10_10_23_R1041_HPRC_NA18565_5_dorado0.6.0_sup4.3.0_5mCG_5hmCG.bam |
| NA18608 | 10_31_23_R1041_HPRC_NA18608_1 | HPRC | dorado0.6.0_sup4.3.0_5mCG_5hmCG/10_31_23_R1041_HPRC_NA18608_1_dorado0.6.0_sup4.3.0_5mCG_5hmCG.bam |
| NA18620 | 01_16_24_R1041_HPRC_NA18620_1 | HPRC | dorado0.6.0_sup4.3.0_5mCG_5hmCG/01_16_24_R1041_HPRC_NA18620_1_dorado0.6.0_sup4.3.0_5mCG_5hmCG.bam<br>dorado0.6.0_sup4.3.0_5mCG_5hmCG/11_14_23_R1041_HPRC_NA18620_1_dorado0.6.0_sup4.3.0_5mCG_5hmCG.bam<br>dorado0.6.0_sup4.3.0_5mCG_5hmCG/11_14_23_R1041_HPRC_NA18620_2_dorado0.6.0_sup4.3.0_5mCG_5hmCG.bam<br>dorado0.6.0_sup4.3.0_5mCG_5hmCG/11_14_23_R1041_HPRC_NA18620_3_dorado0.6.0_sup4.3.0_5mCG_5hmCG.bam |
| NA18879 | 08_08_23_R1041_HPRC_NA18879_1 | HPRC | dorado0.6.0_sup4.3.0_5mCG_5hmCG/08_08_23_R1041_HPRC_NA18879_1_dorado0.6.0_sup4.3.0_5mCG_5hmCG.bam<br>dorado0.6.0_sup4.3.0_5mCG_5hmCG/08_08_23_R1041_HPRC_NA18879_2_dorado0.6.0_sup4.3.0_5mCG_5hmCG.bam |
| NA18906 | 02_15_23_R1041_NA18906_2_dorado0.7.2_sup4.1.0_5mCG_5hmCG | HPRC_PLUS | dorado0.7.2_sup4.1.0_5mCG_5hmCG/02_15_23_R1041_NA18906_2_dorado0.7.2_sup4.1.0_5mCG_5hmCG.bam |
| NA18940 | 11_28_23_R1041_AMED-Collaboration-HPRC_GM18940_1 | HPP | dorado0.5.3_sup4.3.0_5mCG_5hmCG/11_28_23_R1041_AMED-Collaboration-HPRC_NA18940_1_dorado0.5.3_sup4.3.0_5mCG_5hmCG.bam |
| NA18943 | 11_28_23_R1041_AMED-Collaboration-HPRC_GM18943_3 | HPP | dorado0.5.3_sup4.3.0_5mCG_5hmCG/11_28_23_R1041_AMED-Collaboration-HPRC_NA18943_3_dorado0.5.3_sup4.3.0_5mCG_5hmCG.bam<br>dorado0.5.3_sup4.3.0_5mCG_5hmCG/12_12_23_R1041_AMED-Collaboration-HPRC_NA18943_2_dorado0.5.3_sup4.3.0_5mCG_5hmCG.bam |
| NA18944 | 11_28_23_R1041_AMED-Collaboration-HPRC_GM18944_1 | HPP | dorado0.5.3_sup4.3.0_5mCG_5hmCG/11_28_23_R1041_AMED-Collaboration-HPRC_NA18944_1_dorado0.5.3_sup4.3.0_5mCG_5hmCG.bam |
| NA18945 | 12_05_23_R1041_AMED-Collaboration-HPRC_GM18945_1 | HPP | dorado0.5.3_sup4.3.0_5mCG_5hmCG/12_05_23_R1041_AMED-Collaboration-HPRC_NA18945_1_dorado0.5.3_sup4.3.0_5mCG_5hmCG.bam<br>dorado0.5.3_sup4.3.0_5mCG_5hmCG/12_05_23_R1041_AMED-Collaboration-HPRC_NA18945_2_dorado0.5.3_sup4.3.0_5mCG_5hmCG.bam<br>dorado0.5.3_sup4.3.0_5mCG_5hmCG/12_05_23_R1041_AMED-Collaboration-HPRC_NA18945_3_dorado0.5.3_sup4.3.0_5mCG_5hmCG.bam<br>dorado0.5.3_sup4.3.0_5mCG_5hmCG/12_05_23_R1041_AMED-Collaboration-HPRC_NA18945_4_dorado0.5.3_sup4.3.0_5mCG_5hmCG.bam<br>dorado0.5.3_sup4.3.0_5mCG_5hmCG/12_12_23_R1041_AMED-Collaboration-HPRC_NA18945_1_dorado0.5.3_sup4.3.0_5mCG_5hmCG.bam |
| NA18948 | 01_09_24_R1041_AMED-Collaboration-HPRC_GM18948_2 | HPP | dorado0.5.3_sup4.3.0_5mCG_5hmCG/01_09_24_R1041_AMED-Collaboration-HPRC_NA18948_2_dorado0.5.3_sup4.3.0_5mCG_5hmCG.bam |
| NA18952 | 11_28_23_R1041_HPRC_NA18952_1 | HPRC | dorado0.6.0_sup4.3.0_5mCG_5hmCG/11_28_23_R1041_HPRC_NA18952_1_dorado0.6.0_sup4.3.0_5mCG_5hmCG.bam |
| NA18959 | 01_09_24_R1041_AMED-Collaboration-HPRC_GM18959_1 | HPP | dorado0.5.3_sup4.3.0_5mCG_5hmCG/01_09_24_R1041_AMED-Collaboration-HPRC_NA18959_1_dorado0.5.3_sup4.3.0_5mCG_5hmCG.bam |
| NA18960 | 01_09_24_R1041_AMED-Collaboration-HPRC_GM18960_1 | HPP | dorado0.5.3_sup4.3.0_5mCG_5hmCG/01_09_24_R1041_AMED-Collaboration-HPRC_NA18960_1_dorado0.5.3_sup4.3.0_5mCG_5hmCG.bam |
| NA18967 | 01_09_24_R1041_AMED-Collaboration-HPRC_GM18967_1 | HPP | dorado0.5.3_sup4.3.0_5mCG_5hmCG/01_09_24_R1041_AMED-Collaboration-HPRC_NA18967_1_dorado0.5.3_sup4.3.0_5mCG_5hmCG.bam<br>dorado0.5.3_sup4.3.0_5mCG_5hmCG/01_09_24_R1041_AMED-Collaboration-HPRC_NA18967_2_dorado0.5.3_sup4.3.0_5mCG_5hmCG.bam |
| NA18970 | 01_09_24_R1041_AMED-Collaboration-HPRC_GM18970_1 | HPP | dorado0.5.3_sup4.3.0_5mCG_5hmCG/01_09_24_R1041_AMED-Collaboration-HPRC_NA18970_1_dorado0.5.3_sup4.3.0_5mCG_5hmCG.bam |
| NA18974 | 11_28_23_R1041_HPRC_NA18974_1 | HPRC | dorado0.6.0_sup4.3.0_5mCG_5hmCG/11_28_23_R1041_HPRC_NA18974_1_dorado0.6.0_sup4.3.0_5mCG_5hmCG.bam<br>dorado0.6.0_sup4.3.0_5mCG_5hmCG/11_28_23_R1041_HPRC_NA18974_3_dorado0.6.0_sup4.3.0_5mCG_5hmCG.bam |
| NA18976 | 12_12_23_R1041_HPRC_NA18976_1 | HPRC | dorado0.6.0_sup4.3.0_5mCG_5hmCG/12_12_23_R1041_HPRC_NA18976_1_dorado0.6.0_sup4.3.0_5mCG_5hmCG.bam |
| NA18982 | 01_09_24_R1041_AMED-Collaboration-HPRC_GM18982_1 | HPP | dorado0.5.3_sup4.3.0_5mCG_5hmCG/01_09_24_R1041_AMED-Collaboration-HPRC_NA18982_1_dorado0.5.3_sup4.3.0_5mCG_5hmCG.bam<br>dorado0.5.3_sup4.3.0_5mCG_5hmCG/01_09_24_R1041_AMED-Collaboration-HPRC_NA18982_2_dorado0.5.3_sup4.3.0_5mCG_5hmCG.bam |
| NA19036 | 01_23_24_R1041_HPRC_NA19036_1 | HPRC | dorado0.6.0_sup4.3.0_5mCG_5hmCG/01_23_24_R1041_HPRC_NA19036_1_dorado0.6.0_sup4.3.0_5mCG_5hmCG.bam<br>dorado0.6.0_sup4.3.0_5mCG_5hmCG/01_23_24_R1041_HPRC_NA19036_2_dorado0.6.0_sup4.3.0_5mCG_5hmCG.bam<br>dorado0.6.0_sup4.3.0_5mCG_5hmCG/01_23_24_R1041_HPRC_NA19036_3_dorado0.6.0_sup4.3.0_5mCG_5hmCG.bam |
| NA19240 | 03_14_23_R1041_NA19240_1_dorado0.7.2_sup4.1.0_5mCG_5hmCG | HPRC_PLUS | dorado0.7.2_sup4.1.0_5mCG_5hmCG/03_14_23_R1041_NA19240_1_dorado0.7.2_sup4.1.0_5mCG_5hmCG.bam |
| NA19443 | 01_23_24_R1041_HPRC_NA19443_1 | HPRC | dorado0.6.0_sup4.3.0_5mCG_5hmCG/01_23_24_R1041_HPRC_NA19443_1_dorado0.6.0_sup4.3.0_5mCG_5hmCG.bam<br>dorado0.6.0_sup4.3.0_5mCG_5hmCG/01_23_24_R1041_HPRC_NA19443_2_dorado0.6.0_sup4.3.0_5mCG_5hmCG.bam |
| NA19682 | 11_28_23_R1041_HPRC_NA19682_1 | HPRC | dorado0.6.0_sup4.3.0_5mCG_5hmCG/11_28_23_R1041_HPRC_NA19682_1_dorado0.6.0_sup4.3.0_5mCG_5hmCG.bam |
| NA19700 | 01_23_24_R1041_HPRC_NA19700_1 | HPRC | dorado0.6.0_sup4.3.0_5mCG_5hmCG/01_23_24_R1041_HPRC_NA19700_1_dorado0.6.0_sup4.3.0_5mCG_5hmCG.bam<br>dorado0.6.0_sup4.3.0_5mCG_5hmCG/01_23_24_R1041_HPRC_NA19700_2_dorado0.6.0_sup4.3.0_5mCG_5hmCG.bam |
| NA19776 | 10_24_23_R1041_HPRC_NA19776_1 | HPRC | dorado0.6.0_sup4.3.0_5mCG_5hmCG/10_24_23_R1041_HPRC_NA19776_1_dorado0.6.0_sup4.3.0_5mCG_5hmCG.bam |
| NA19835 | 10_31_23_R1041_HPRC_NA19835_1 | HPRC | dorado0.6.0_sup4.3.0_5mCG_5hmCG/10_31_23_R1041_HPRC_NA19835_1_dorado0.6.0_sup4.3.0_5mCG_5hmCG.bam |
| NA19909 | 11_14_23_R1041_HPRC_NA19909_1 | HPRC | dorado0.6.0_sup4.3.0_5mCG_5hmCG/11_14_23_R1041_HPRC_NA19909_1_dorado0.6.0_sup4.3.0_5mCG_5hmCG.bam |
| NA20129 | 03_08_23_R1041_NA20129_1_dorado0.7.2_sup4.1.0_5mCG_5hmCG | HPRC_PLUS | dorado0.7.2_sup4.1.0_5mCG_5hmCG/03_08_23_R1041_NA20129_1_dorado0.7.2_sup4.1.0_5mCG_5hmCG.bam |
| NA20282 | 12_12_23_R1041_HPRC_NA20282_1 | HPRC | dorado0.6.0_sup4.3.0_5mCG_5hmCG/12_12_23_R1041_HPRC_NA20282_1_dorado0.6.0_sup4.3.0_5mCG_5hmCG.bam |
| NA20346 | 11_28_23_R1041_HPRC_NA20346_2 | HPRC | dorado0.6.0_sup4.3.0_5mCG_5hmCG/11_28_23_R1041_HPRC_NA20346_2_dorado0.6.0_sup4.3.0_5mCG_5hmCG.bam |
| NA20503 | 12_12_23_R1041_HPRC_NA20503_1 | HPRC | dorado0.6.0_sup4.3.0_5mCG_5hmCG/12_12_23_R1041_HPRC_NA20503_1_dorado0.6.0_sup4.3.0_5mCG_5hmCG.bam |
| NA20762 | 01_23_24_R1041_HPRC_NA20762_1 | HPRC | dorado0.6.0_sup4.3.0_5mCG_5hmCG/01_23_24_R1041_HPRC_NA20762_1_dorado0.6.0_sup4.3.0_5mCG_5hmCG.bam<br>dorado0.6.0_sup4.3.0_5mCG_5hmCG/10_31_23_R1041_HPRC_NA20762_1_dorado0.6.0_sup4.3.0_5mCG_5hmCG.bam |
| NA20806 | 12_12_23_R1041_HPRC_NA20806_1 | HPRC | dorado0.6.0_sup4.3.0_5mCG_5hmCG/12_12_23_R1041_HPRC_NA20806_1_dorado0.6.0_sup4.3.0_5mCG_5hmCG.bam |
| NA20809 | 08_08_23_R1041_HPRC_NA20809_1 | HPRC | dorado0.6.0_sup4.3.0_5mCG_5hmCG/08_08_23_R1041_HPRC_NA20809_1_dorado0.6.0_sup4.3.0_5mCG_5hmCG.bam |
| NA20827 | 01_23_24_R1041_HPRC_NA20827_1 | HPRC | dorado0.6.0_sup4.3.0_5mCG_5hmCG/01_23_24_R1041_HPRC_NA20827_1_dorado0.6.0_sup4.3.0_5mCG_5hmCG.bam<br>dorado0.6.0_sup4.3.0_5mCG_5hmCG/01_23_24_R1041_HPRC_NA20827_2_dorado0.6.0_sup4.3.0_5mCG_5hmCG.bam |
| NA20850 | 08_08_23_R1041_HPRC_NA20850_1 | HPRC | dorado0.6.0_sup4.3.0_5mCG_5hmCG/08_08_23_R1041_HPRC_NA20850_1_dorado0.6.0_sup4.3.0_5mCG_5hmCG.bam |
| NA20870 | 11_14_23_R1041_HPRC_NA20870_1 | HPRC | dorado0.6.0_sup4.3.0_5mCG_5hmCG/11_14_23_R1041_HPRC_NA20870_1_dorado0.6.0_sup4.3.0_5mCG_5hmCG.bam |
| NA21093 | 01_23_24_R1041_HPRC_NA21093_1 | HPRC | dorado0.6.0_sup4.3.0_5mCG_5hmCG/01_23_24_R1041_HPRC_NA21093_1_dorado0.6.0_sup4.3.0_5mCG_5hmCG.bam |
| NA21102 | 01_23_24_R1041_HPRC_NA21102_1 | HPRC | dorado0.6.0_sup4.3.0_5mCG_5hmCG/01_23_24_R1041_HPRC_NA21102_1_dorado0.6.0_sup4.3.0_5mCG_5hmCG.bam<br>dorado0.6.0_sup4.3.0_5mCG_5hmCG/01_23_24_R1041_HPRC_NA21102_2_dorado0.6.0_sup4.3.0_5mCG_5hmCG.bam |
| NA21106 | 10_24_23_R1041_HPRC_NA21106_1 | HPRC | dorado0.6.0_sup4.3.0_5mCG_5hmCG/10_24_23_R1041_HPRC_NA21106_1_dorado0.6.0_sup4.3.0_5mCG_5hmCG.bam<br>dorado0.6.0_sup4.3.0_5mCG_5hmCG/10_24_23_R1041_HPRC_NA21106_4_dorado0.6.0_sup4.3.0_5mCG_5hmCG.bam |
| NA21110 | 01_16_24_R1041_HPRC_NA21110_1 | HPRC | dorado0.6.0_sup4.3.0_5mCG_5hmCG/01_16_24_R1041_HPRC_NA21110_1_dorado0.6.0_sup4.3.0_5mCG_5hmCG.bam |
| NA21144 | 01_23_24_R1041_HPRC_NA21144_1 | HPRC | dorado0.6.0_sup4.3.0_5mCG_5hmCG/01_23_24_R1041_HPRC_NA21144_1_dorado0.6.0_sup4.3.0_5mCG_5hmCG.bam<br>dorado0.6.0_sup4.3.0_5mCG_5hmCG/12_05_23_R1041_HPRC_NA21144_1_dorado0.6.0_sup4.3.0_5mCG_5hmCG.bam |
| NA21309 | 02_22_23_R1041_NA21309_1_dorado0.7.2_sup4.1.0_5mCG_5hmCG | HPRC_PLUS | dorado0.7.2_sup4.1.0_5mCG_5hmCG/02_22_23_R1041_NA21309_1_dorado0.7.2_sup4.1.0_5mCG_5hmCG.bam |

