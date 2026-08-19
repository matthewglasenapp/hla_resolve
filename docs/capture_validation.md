# Hybrid Capture Validation

HLA-Resolve was run on PacBio hybrid capture reads for 31 samples, using IPD-IMGT/HLA release v3.64.0. Every gene in every sample was typed, giving 496 of 496 possible alleles. Mean coverage across the eight classical HLA genes was 365×, ranging from 43× to 822×.

HLA allele concordance was measured for the 27 samples that carry a reference typing. Fifteen come from the International HLA and Immunogenetics Workshop (IHWG) and twelve from the Human Pangenome Reference Consortium (HPRC), scored against the reference typings of Lai et al. 2024 ([DOI: 10.1016/j.csbj.2024.03.030](https://doi.org/10.1016/j.csbj.2024.03.030); [Supplementary File 6](Lai_Supplementary-6.xlsx)). The two truth sets share no samples and are reported separately. Four samples have no reference typing and are typed but not scored.

The following command was used for each sample:

```bash
hla_resolve --input_file <FASTQ> --sample_name <SAMPLE> \
            --platform pacbio --scheme hybrid_capture \
            --output_dir <OUTPUT_DIR> --threads <THREADS>
```

Concordance is reported among alleles called. An allele was evaluated at a given field resolution only if the reference specified it to that resolution.

## Summary

| Resolution | IHWG | HPRC | Combined |
|------------|------|------|----------|
| One field | 227/230 (98.7%) | 192/192 (100%) | 419/422 (99.3%) |
| Two field | 219/222 (98.6%) | 192/192 (100%) | 411/414 (99.3%) |
| Three field | 196/198 (99.0%) | 190/190 (100%) | 386/388 (99.5%) |
| Four field | 135/162 (83.3%) | 180/186 (96.8%) | 315/348 (90.5%) |

## Sequencing reads

Raw reads are available from the NCBI Sequence Read Archive under BioProject [PRJNA1417783](https://www.ncbi.nlm.nih.gov/bioproject/PRJNA1417783). The run accessions below are the PacBio HiFi runs for each sample.

### Samples

| Sample | Reference typing | SRA run | HLA coverage | Alleles called |
|---|---|---|---|---|
| HG002 | HPRC | SRR37076010 | 822× | 16/16 |
| HG003 | none | SRR37076009 | 540× | 16/16 |
| HG004 | none | SRR37075998 | 372× | 16/16 |
| HG005 | HPRC | SRR37075987 | 327× | 16/16 |
| HG01106 | HPRC | SRR37075982 | 629× | 16/16 |
| HG01258 | HPRC | SRR37075981 | 732× | 16/16 |
| HG01928 | HPRC | SRR37075980 | 450× | 16/16 |
| HG02055 | HPRC | SRR37075979 | 436× | 16/16 |
| HG02630 | HPRC | SRR37076008 | 420× | 16/16 |
| HG03492 | HPRC | SRR37076007 | 369× | 16/16 |
| HG03579 | HPRC | SRR37076006 | 179× | 16/16 |
| IHW09021 | IHWG | SRR40178158 | 182× | 16/16 |
| IHW09049 | IHWG | SRR40178157 | 327× | 16/16 |
| IHW09071 | IHWG | SRR40178146 | 111× | 16/16 |
| IHW09118 | IHWG | SRR40178132 | 73× | 16/16 |
| IHW09122 | IHWG | SRR40178131 | 633× | 16/16 |
| IHW09125 | IHWG | SRR40178130 | 344× | 16/16 |
| IHW09175 | IHWG | SRR40178129 | 162× | 16/16 |
| IHW09198 | IHWG | SRR40178128 | 303× | 16/16 |
| IHW09200 | IHWG | SRR40178127 | 326× | 16/16 |
| IHW09224 | IHWG | SRR40178156 | 323× | 16/16 |
| IHW09245 | IHWG | SRR40178155 | 43× | 16/16 |
| IHW09251 | IHWG | SRR40178154 | 181× | 16/16 |
| IHW09359 | IHWG | SRR40178153 | 447× | 16/16 |
| IHW09364 | IHWG | SRR40178152 | 78× | 16/16 |
| IHW09409 | IHWG | SRR40178151 | 213× | 16/16 |
| NA19240 | HPRC | SRR37076005 | 203× | 16/16 |
| NA20129 | HPRC | SRR37076004 | 653× | 16/16 |
| NA21309 | HPRC | SRR37076003 | 567× | 16/16 |
| NA24694 | none | SRR37075984 | 615× | 16/16 |
| NA24695 | none | SRR37075983 | 248× | 16/16 |

### Concordance by gene, IHWG

| Gene | 1 field | 2 field | 3 field | 4 field |
|---|---|---|---|---|
| HLA-A | 30/30 (100%) | 28/28 (100%) | 28/28 (100%) | 24/24 (100%) |
| HLA-B | 30/30 (100%) | 28/28 (100%) | 26/26 (100%) | 16/19 (84%) |
| HLA-C | 30/30 (100%) | 28/28 (100%) | 27/27 (100%) | 20/24 (83%) |
| HLA-DPA1 | 26/26 (100%) | 26/26 (100%) | 25/25 (100%) | 18/23 (78%) |
| HLA-DPB1 | 27/28 (96%) | 27/28 (96%) | 21/21 (100%) | 15/16 (94%) |
| HLA-DQA1 | 27/28 (96%) | 27/28 (96%) | 23/24 (96%) | 11/20 (55%) |
| HLA-DQB1 | 27/28 (96%) | 27/28 (96%) | 25/26 (96%) | 18/20 (90%) |
| HLA-DRB1 | 30/30 (100%) | 28/28 (100%) | 21/21 (100%) | 13/16 (81%) |
| **Total** | **227/230 (98.7%)** | **219/222 (98.6%)** | **196/198 (99.0%)** | **135/162 (83.3%)** |

### Concordance by gene, HPRC

| Gene | 1 field | 2 field | 3 field | 4 field |
|---|---|---|---|---|
| HLA-A | 24/24 (100%) | 24/24 (100%) | 23/23 (100%) | 23/23 (100%) |
| HLA-B | 24/24 (100%) | 24/24 (100%) | 23/23 (100%) | 21/23 (91%) |
| HLA-C | 24/24 (100%) | 24/24 (100%) | 24/24 (100%) | 23/24 (96%) |
| HLA-DPA1 | 24/24 (100%) | 24/24 (100%) | 24/24 (100%) | 24/24 (100%) |
| HLA-DPB1 | 24/24 (100%) | 24/24 (100%) | 24/24 (100%) | 21/23 (91%) |
| HLA-DQA1 | 24/24 (100%) | 24/24 (100%) | 24/24 (100%) | 24/24 (100%) |
| HLA-DQB1 | 24/24 (100%) | 24/24 (100%) | 24/24 (100%) | 24/24 (100%) |
| HLA-DRB1 | 24/24 (100%) | 24/24 (100%) | 24/24 (100%) | 20/21 (95%) |
| **Total** | **192/192 (100.0%)** | **192/192 (100.0%)** | **190/190 (100.0%)** | **180/186 (96.8%)** |
