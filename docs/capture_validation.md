# Hybrid Capture Validation

HLA-Resolve was run on PacBio hybrid capture reads for 31 samples, using IPD-IMGT/HLA release v3.64.0. Every gene in every sample was typed, giving 496 of 496 possible alleles. Mean coverage across the eight classical HLA genes was 365×, ranging from 43× to 822×.

HLA allele concordance was measured for the 27 samples that carry a reference typing. Fifteen come from the International Histocompatibility Working Group (IHWG) and twelve from the Human Pangenome Reference Consortium (HPRC). The HPRC reference typings are from Lai et al. 2024 ([DOI: 10.1016/j.csbj.2024.03.030](https://doi.org/10.1016/j.csbj.2024.03.030); [Supplementary File 6](Lai_Supplementary-6.xlsx)). The two truth sets share no samples. Four samples have no reference typing and were typed but not scored.

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
| One field | 98.7% (227/230) | 100% (192/192) | 99.3% (419/422) |
| Two field | 98.6% (219/222) | 100% (192/192) | 99.3% (411/414) |
| Three field | 99.0% (196/198) | 100% (190/190) | 99.5% (386/388) |
| Four field | 83.3% (135/162) | 96.8% (180/186) | 90.5% (315/348) |

## Concordance by gene

Both truth sets pooled.

| Gene | One field | Two field | Three field | Four field |
|------|------|------|------|------|
| HLA-A | 100% (54/54) | 100% (52/52) | 100% (51/51) | 100% (47/47) |
| HLA-B | 100% (54/54) | 100% (52/52) | 100% (49/49) | 88.1% (37/42) |
| HLA-C | 100% (54/54) | 100% (52/52) | 100% (51/51) | 89.6% (43/48) |
| HLA-DPA1 | 100% (50/50) | 100% (50/50) | 100% (49/49) | 89.4% (42/47) |
| HLA-DPB1 | 98.1% (51/52) | 98.1% (51/52) | 100% (45/45) | 92.3% (36/39) |
| HLA-DQA1 | 98.1% (51/52) | 98.1% (51/52) | 97.9% (47/48) | 79.5% (35/44) |
| HLA-DQB1 | 98.1% (51/52) | 98.1% (51/52) | 98.0% (49/50) | 95.5% (42/44) |
| HLA-DRB1 | 100% (54/54) | 100% (52/52) | 100% (45/45) | 89.2% (33/37) |
| **Total** | **99.3% (419/422)** | **99.3% (411/414)** | **99.5% (386/388)** | **90.5% (315/348)** |

## Sequencing reads

Raw reads are available from the NCBI Sequence Read Archive under BioProject [PRJNA1417783](https://www.ncbi.nlm.nih.gov/bioproject/PRJNA1417783). The accessions below are the PacBio HiFi runs for each sample.

| Sample | SRA Accession | HLA Coverage |
|--------|---------------|--------------|
| HG002 | SRR37076010 | 822× |
| HG003 | SRR37076009 | 540× |
| HG004 | SRR37075998 | 372× |
| HG005 | SRR37075987 | 327× |
| HG01106 | SRR37075982 | 629× |
| HG01258 | SRR37075981 | 732× |
| HG01928 | SRR37075980 | 450× |
| HG02055 | SRR37075979 | 436× |
| HG02630 | SRR37076008 | 420× |
| HG03492 | SRR37076007 | 369× |
| HG03579 | SRR37076006 | 179× |
| IHW09021 | SRR40178158 | 182× |
| IHW09049 | SRR40178157 | 327× |
| IHW09071 | SRR40178146 | 111× |
| IHW09118 | SRR40178132 | 73× |
| IHW09122 | SRR40178131 | 633× |
| IHW09125 | SRR40178130 | 344× |
| IHW09175 | SRR40178129 | 162× |
| IHW09198 | SRR40178128 | 303× |
| IHW09200 | SRR40178127 | 326× |
| IHW09224 | SRR40178156 | 323× |
| IHW09245 | SRR40178155 | 43× |
| IHW09251 | SRR40178154 | 181× |
| IHW09359 | SRR40178153 | 447× |
| IHW09364 | SRR40178152 | 78× |
| IHW09409 | SRR40178151 | 213× |
| NA19240 | SRR37076005 | 203× |
| NA20129 | SRR37076004 | 653× |
| NA21309 | SRR37076003 | 567× |
| NA24694 | SRR37075984 | 615× |
| NA24695 | SRR37075983 | 248× |
