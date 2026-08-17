# WGS Validation

HLA-Resolve v0.9.0 was run on whole-genome PacBio raw sequencing reads for 39 samples from the Human Pangenome Reference Consortium (HPRC) and two additional HG002 runs from the Genome in a Bottle Consortium (GIAB), using IPD-IMGT/HLA release v3.64.0. HG002 appears in both sample sets. The GIAB runs are reported separately and are not included in the HPRC totals.

HLA allele concordance was measured against the reference typings provided by Lai et al. 2024 ([DOI: 10.1016/j.csbj.2024.03.030](https://doi.org/10.1016/j.csbj.2024.03.030); [Supplementary File 6](Lai_Supplementary-6.xlsx)). Where a sample was sequenced across multiple SMRT cells, reads were merged into a single unaligned BAM to provide adequate coverage of the HLA genes.

The following command was used for each sample:

```bash
hla_resolve --input_file <MERGED_uBAM> --sample_name <SAMPLE> \
            --platform pacbio --scheme WGS \
            --output_dir <OUTPUT_DIR> --threads <THREADS>
```


Concordance is reported among alleles called. An allele was evaluated at a given field resolution only if the reference specified it to that resolution, which is why the denominator differs between resolutions.


## Genome in a Bottle

Two PacBio Revio runs of HG002, each typed from a single file.

| Sample | Run | Instrument | HLA&nbsp;coverage | Genes&nbsp;typed | Alleles&nbsp;called | 1-field | 2-field | 3-field | 4-field |
|---|---|---|---:|---|---|---|---|---|---|
| HG002 | 231005_s1 | Revio | 26.3× | 8/8 | 16/16 | 16/16 | 16/16 | 16/16 | 14/15 |
| HG002 | 230928_s3 | Revio | 23.4× | 8/8 | 16/16 | 16/16 | 16/16 | 16/16 | 15/15 |

**231005_s1**

```
https://ftp-trace.ncbi.nlm.nih.gov/ReferenceSamples/giab/data/AshkenazimTrio/HG002_NA24385_son/PacBio_HiFi-Revio_20231031/HG002_PacBio-Revio_m84039_231005_222902_s1.hifi_reads.bam
```

**230928_s3**

```
https://ftp-trace.ncbi.nlm.nih.gov/ReferenceSamples/giab/data/AshkenazimTrio/HG002_NA24385_son/PacBio_HiFi-Revio_20231031/HG002_PacBio-Revio_m84039_230928_213653_s3.hifi_reads.bam
```


## Human Pangenome Reference Consortium

Empirical coverage across the eight classical HLA genes was 35.2x (median 37.0x, range 12.4–55.9x).

### Call rate

| Metric | Value |
|---|---:|
| Total possible allele calls (39 samples × 16) | 624 |
| Alleles called | 620 |
| Call rate | 99.4% |


### Concordance

| Resolution | Concordance | Concordant alleles | Alleles evaluated |
|---|---:|---:|---:|
| 1-field | 100.0% | 619 | 619 |
| 2-field | 100.0% | 619 | 619 |
| 3-field | 100.0% | 610 | 610 |
| 4-field | 92.8% | 555 | 598 |


### Results by sample

| Sample | Instrument | HLA&nbsp;coverage | Genes&nbsp;typed | Alleles&nbsp;called | 1-field | 2-field | 3-field | 4-field |
|---|---|---:|---|---|---|---|---|---|
| HG002 | Revio | 32.2× | 8/8 | 16/16 | 16/16 | 16/16 | 16/16 | 15/15 |
| HG00438 | Sequel II | 40.9× | 8/8 | 16/16 | 16/16 | 16/16 | 16/16 | 16/16 |
| HG00673 | Sequel II | 40.0× | 8/8 | 16/16 | 16/16 | 16/16 | 16/16 | 14/16 |
| HG00733 | Sequel II | 13.0× | 7/8 | 14/16 | 14/14 | 14/14 | 14/14 | 11/14 |
| HG00735 | Sequel II | 36.4× | 8/8 | 16/16 | 16/16 | 16/16 | 16/16 | 14/16 |
| HG00741 | Sequel II | 35.9× | 8/8 | 16/16 | 16/16 | 16/16 | 16/16 | 14/16 |
| HG01071 | Revio | 12.4× | 7/8 | 14/16 | 14/14 | 14/14 | 14/14 | 13/14 |
| HG01106 | Sequel II | 43.9× | 8/8 | 16/16 | 16/16 | 16/16 | 16/16 | 15/16 |
| HG01109 | Sequel II | 38.1× | 8/8 | 16/16 | 16/16 | 16/16 | 16/16 | 15/16 |
| HG01175 | Sequel II | 36.0× | 8/8 | 16/16 | 16/16 | 16/16 | 16/16 | 15/16 |
| HG01243 | Sequel II | 35.6× | 8/8 | 16/16 | 16/16 | 16/16 | 16/16 | 15/16 |
| HG01258 | Revio | 18.8× | 8/8 | 16/16 | 16/16 | 16/16 | 15/15 | 14/15 |
| HG01358 | Sequel II | 37.0× | 8/8 | 16/16 | 15/15 | 15/15 | 15/15 | 14/14 |
| HG01361 | Revio | 38.7× | 8/8 | 16/16 | 16/16 | 16/16 | 16/16 | 14/15 |
| HG01891 | Revio | 24.5× | 8/8 | 16/16 | 16/16 | 16/16 | 16/16 | 15/16 |
| HG01928 | Revio | 20.8× | 8/8 | 16/16 | 16/16 | 16/16 | 15/15 | 14/15 |
| HG01952 | Sequel II | 41.3× | 8/8 | 16/16 | 16/16 | 16/16 | 16/16 | 13/14 |
| HG01978 | Sequel II | 33.1× | 8/8 | 16/16 | 16/16 | 16/16 | 16/16 | 14/15 |
| HG02055 | Revio | 34.0× | 8/8 | 16/16 | 16/16 | 16/16 | 16/16 | 15/16 |
| HG02109 | Revio | 35.4× | 8/8 | 16/16 | 16/16 | 16/16 | 14/14 | 14/14 |
| HG02145 | Revio | 42.7× | 8/8 | 16/16 | 16/16 | 16/16 | 16/16 | 14/16 |
| HG02148 | Sequel II | 36.3× | 8/8 | 16/16 | 16/16 | 16/16 | 16/16 | 14/16 |
| HG02257 | Sequel II | 37.0× | 8/8 | 16/16 | 16/16 | 16/16 | 16/16 | 14/15 |
| HG02572 | Sequel II | 19.8× | 8/8 | 16/16 | 16/16 | 16/16 | 15/15 | 12/14 |
| HG02622 | Sequel II | 42.3× | 8/8 | 16/16 | 16/16 | 16/16 | 16/16 | 12/16 |
| HG02630 | Sequel II | 49.5× | 8/8 | 16/16 | 16/16 | 16/16 | 16/16 | 15/16 |
| HG02717 | Sequel II | 43.6× | 8/8 | 16/16 | 16/16 | 16/16 | 15/15 | 11/15 |
| HG02723 | Revio | 34.4× | 8/8 | 16/16 | 16/16 | 16/16 | 15/15 | 14/15 |
| HG02818 | Sequel II | 37.1× | 8/8 | 16/16 | 16/16 | 16/16 | 16/16 | 15/16 |
| HG02886 | Sequel II | 48.2× | 8/8 | 16/16 | 16/16 | 16/16 | 16/16 | 14/16 |
| HG03098 | Sequel II | 35.0× | 8/8 | 16/16 | 16/16 | 16/16 | 15/15 | 15/15 |
| HG03486 | Sequel II | 39.5× | 8/8 | 16/16 | 16/16 | 16/16 | 16/16 | 16/16 |
| HG03516 | Revio | 24.4× | 8/8 | 16/16 | 16/16 | 16/16 | 16/16 | 16/16 |
| HG03540 | Sequel II | 55.9× | 8/8 | 16/16 | 16/16 | 16/16 | 15/15 | 14/14 |
| HG03579 | Sequel II | 15.2× | 8/8 | 16/16 | 16/16 | 16/16 | 16/16 | 15/16 |
| NA18906 | Sequel II | 46.8× | 8/8 | 16/16 | 16/16 | 16/16 | 16/16 | 15/15 |
| NA19240 | Sequel II | 41.8× | 8/8 | 16/16 | 16/16 | 16/16 | 16/16 | 15/16 |
| NA20129 | Sequel II | 38.2× | 8/8 | 16/16 | 16/16 | 16/16 | 16/16 | 15/15 |
| NA21309 | Sequel II | 37.2× | 8/8 | 16/16 | 16/16 | 16/16 | 16/16 | 15/15 |
| **Total** | | **35.2× mean** | **310/312** | **620/624** | **619/619** | **619/619** | **610/610** | **555/598** |


### Input files

Each sample was typed from a single unaligned BAM produced by combining the files listed below. All files are in the HPRC repository on the AWS Open Data registry. Build a path by substituting the **Cohort**, **Sample** and **Files** values from the table into:

```
s3://human-pangenomics/working/<COHORT>/<SAMPLE>/raw_data/PacBio_HiFi/<FILE>
```

For example, the first file listed for HG00438 (cohort `HPRC`) is at:

```
s3://human-pangenomics/working/HPRC/HG00438/raw_data/PacBio_HiFi/m64043_200710_174426.ccs.bam
```

Files are public and need no credentials:

```bash
aws s3 cp --no-sign-request \
  s3://human-pangenomics/working/HPRC/HG00438/raw_data/PacBio_HiFi/m64043_200710_174426.ccs.bam .
```

> [!NOTE]
> HG002's file is under an additional `wMods/` subdirectory: `s3://human-pangenomics/working/HPRC_PLUS/HG002/raw_data/PacBio_HiFi/wMods/m84011_220902_175841_s1.hifi_reads.bam`

| Sample | Cohort | Instrument | Files |
|---|---|---|---|
| HG002 | HPRC_PLUS | Revio | `m84011_220902_175841_s1.hifi_reads.bam` |
| HG00438 | HPRC | Sequel II | `m64043_200710_174426.ccs.bam`<br>`m64043_200711_235708.ccs.bam`<br>`m64043_200713_062240.ccs.bam`<br>`m64043_200714_124814.ccs.bam` |
| HG00673 | HPRC | Sequel II | `m64043_200716_182902.ccs.bam`<br>`m64043_200718_004213.ccs.bam`<br>`m64043_200719_070806.ccs.bam`<br>`m64043_200720_133355.ccs.bam` |
| HG00733 | HPRC_PLUS | Sequel II | `m64076_211214_012715.hifi_reads.bam` |
| HG00735 | HPRC | Sequel II | `m64043_200702_173033.ccs.bam`<br>`m64043_200703_234328.ccs.bam`<br>`m64043_200705_060840.ccs.bam`<br>`m64043_200706_123423.ccs.bam` |
| HG00741 | HPRC | Sequel II | `m64136_200625_174949.ccs.bam`<br>`m64136_200627_000247.ccs.bam`<br>`m64136_200628_062837.ccs.bam`<br>`m64136_200629_125431.ccs.bam` |
| HG01071 | HPRC | Revio | `m84081_231222_081401_s1.hifi_reads.bc2026.bam`<br>`m84081_231222_101332_s2.hifi_reads.bc2026.bam`<br>`m84081_231222_121301_s3.hifi_reads.bc2026.bam`<br>`m84081_231222_141231_s4.hifi_reads.bc2026.bam` |
| HG01106 | HPRC | Sequel II | `m64043_200625_174853.ccs.bam`<br>`m64043_200627_000137.ccs.bam`<br>`m64043_200628_062711.ccs.bam`<br>`m64043_200629_125238.ccs.bam` |
| HG01109 | HPRC_PLUS | Sequel II | `m64043_200827_191459.ccs.bam`<br>`m64043_200829_012836.ccs.bam`<br>`m64043_200830_075523.ccs.bam` |
| HG01175 | HPRC | Sequel II | `m64043_200618_201934.ccs.bam`<br>`m64043_200620_173220.ccs.bam`<br>`m64043_200621_234442.ccs.bam`<br>`m64043_200623_060946.ccs.bam` |
| HG01243 | HPRC_PLUS | Sequel II | `m64136_200827_191603.ccs.bam`<br>`m64136_200829_012933.ccs.bam`<br>`m64136_200830_075556.ccs.bam` |
| HG01258 | HPRC | Revio | `m84046_231202_090949_s3.hifi_reads.bc2054.bam` |
| HG01358 | HPRC | Sequel II | `m64076_200201_051547.ccs.bam`<br>`m64076_200203_181219.ccs.bam`<br>`m64076_200206_215943.ccs.bam`<br>`m64076_200208_041234.ccs.bam` |
| HG01361 | HPRC | Revio | `m84046_230724_203319_s4.hifi_reads.bc2056.bam`<br>`m84046_231202_071034_s1.hifi_reads.bc2052.bam` |
| HG01891 | HPRC | Revio | `m84046_230712_224626_s1.hifi_reads.bc2066.bam` |
| HG01928 | HPRC | Revio | `m84081_231222_081401_s1.hifi_reads.bc2022.bam`<br>`m84081_231222_101332_s2.hifi_reads.bc2022.bam`<br>`m84081_231222_121301_s3.hifi_reads.bc2022.bam`<br>`m84081_231222_141231_s4.hifi_reads.bc2022.bam` |
| HG01952 | HPRC | Sequel II | `m64136_200612_201033.ccs.bam`<br>`m64136_200614_192134.ccs.bam`<br>`m64136_200616_013426.ccs.bam`<br>`m64136_200617_080019.ccs.bam` |
| HG01978 | HPRC | Sequel II | `m64136_200530_164818.ccs.bam`<br>`m64136_200602_153012.ccs.bam`<br>`m64136_200603_214308.ccs.bam`<br>`m64136_200605_040848.ccs.bam` |
| HG02055 | HPRC_PLUS | Revio | `m84039_230303_005138_s2.hifi_reads.bc2010.bam`<br>`m84039_230311_024907_s1.hifi_reads.bc2010.bam`<br>`m84039_230314_210032_s1.hifi_reads.bc2010.bam`<br>`m84039_230316_185945_s1.hifi_reads.bc2010.bam` |
| HG02109 | HPRC_PLUS | Revio | `m84036_230317_175945_s2.hifi_reads.bc2013.bam`<br>`m84039_230303_012244_s3.hifi_reads.bc2013.bam`<br>`m84039_230314_213047_s2.hifi_reads.bc2013.bam`<br>`m84039_230316_193003_s2.hifi_reads.bc2013.bam` |
| HG02145 | HPRC_PLUS | Revio | `m84036_230317_175945_s2.hifi_reads.bc2014.bam`<br>`m84039_230303_012244_s3.hifi_reads.bc2014.bam`<br>`m84039_230314_213047_s2.hifi_reads.bc2014.bam`<br>`m84039_230316_193003_s2.hifi_reads.bc2014.bam` |
| HG02148 | HPRC | Sequel II | `m64136_200618_202033.ccs.bam`<br>`m64136_200620_173618.ccs.bam`<br>`m64136_200621_234916.ccs.bam`<br>`m64136_200623_061529.ccs.bam` |
| HG02257 | HPRC | Sequel II | `m64076_200125_231256.ccs.bam`<br>`m64076_200127_180545.ccs.bam`<br>`m64076_200129_001835.ccs.bam`<br>`m64076_200130_064345.ccs.bam` |
| HG02572 | HPRC | Sequel II | `m64076_200313_161705.ccs.bam`<br>`m64076_200317_041201.ccs.bam`<br>`m64076_200318_103811.ccs.bam`<br>`m64076_210502_044702.hifi_reads.bam` |
| HG02622 | HPRC | Sequel II | `m64043_200530_164723.ccs.bam`<br>`m64043_200601_191521.ccs.bam`<br>`m64043_200603_012738.ccs.bam`<br>`m64043_200604_075218.ccs.bam` |
| HG02630 | HPRC | Sequel II | `m64043_200501_162248.ccs.bam`<br>`m64043_200502_223511.ccs.bam`<br>`m64043_200504_050026.ccs.bam`<br>`m64043_200505_112554.ccs.bam` |
| HG02717 | HPRC | Sequel II | `m64043_200403_163826.ccs.bam`<br>`m64043_200405_180950.ccs.bam`<br>`m64043_200407_002219.ccs.bam`<br>`m64043_200408_064651.ccs.bam` |
| HG02723 | HPRC_PLUS | Revio | `m84036_230317_175945_s2.hifi_reads.bc2012.bam`<br>`m84039_230303_012244_s3.hifi_reads.bc2012.bam`<br>`m84039_230314_213047_s2.hifi_reads.bc2012.bam`<br>`m84039_230316_193003_s2.hifi_reads.bc2012.bam` |
| HG02818 | HPRC_PLUS | Sequel II | `m64043_200206_173947.ccs.bam`<br>`m64043_200207_235213.ccs.bam`<br>`m64043_200209_061852.ccs.bam`<br>`m64043_200314_004623.ccs.bam`<br>`m64043_200315_071057.ccs.bam`<br>`m64043_200316_214923.ccs.bam`<br>`m64043_200318_040100.ccs.bam` |
| HG02886 | HPRC | Sequel II | `m64043_200410_214826.ccs.bam`<br>`m64043_200412_040054.ccs.bam`<br>`m64043_200413_102554.ccs.bam`<br>`m64043_200414_165107.ccs.bam` |
| HG03098 | HPRC_PLUS | Sequel II | `m64043_201123_083343.ccs.bam`<br>`m64043_201128_031055.ccs.bam`<br>`m64043_201203_004011.ccs.bam` |
| HG03486 | HPRC_PLUS | Sequel II | `m64043_200424_162541.ccs.bam`<br>`m64043_200425_223840.ccs.bam`<br>`m64043_200428_155222.ccs.bam`<br>`m64043_200429_220517.ccs.bam` |
| HG03516 | HPRC | Revio | `m84046_231202_110908_s4.hifi_reads.bc2056.bam` |
| HG03540 | HPRC | Sequel II | `m64043_200521_171703.ccs.bam`<br>`m64043_200522_232930.ccs.bam`<br>`m64043_200524_055430.ccs.bam`<br>`m64043_200525_121851.ccs.bam` |
| HG03579 | HPRC | Sequel II | `m64043_200516_230634.ccs.bam` |
| NA18906 | HPRC_PLUS | Sequel II | `m64136_200521_171936.ccs.bam`<br>`m64136_200523_195722.ccs.bam`<br>`m64136_200525_021027.ccs.bam`<br>`m64136_200526_083627.ccs.bam` |
| NA19240 | HPRC_PLUS | Sequel II | `m64043_200128_181438.ccs.bam`<br>`m64043_200130_175216.ccs.bam`<br>`m64043_200201_000449.ccs.bam`<br>`m64043_200202_062937.ccs.bam` |
| NA20129 | HPRC_PLUS | Sequel II | `m64043_191227_185626.ccs.bam`<br>`m64043_191229_010753.ccs.bam`<br>`m64043_191230_073311.ccs.bam`<br>`m64043_200111_140530.ccs.bam`<br>`m64043_200114_192155.ccs.bam` |
| NA21309 | HPRC_PLUS | Sequel II | `m64043_191210_201113.ccs.bam`<br>`m64043_191213_191857.ccs.bam`<br>`m64043_191215_014401.ccs.bam`<br>`m64043_191219_192900.ccs.bam` |
