# HLA-Resolve Documentation

Documentation for HLA-Resolve.

## 1. Preparing phased VCF for vcf2fasta
Following phasing with HiPhase or longphase, the SNV and SV VCFs are merged and normalized with bcftools. The following steps are taken to prepare the VCF records for use with vcf2fasta. 
1. Filter the joint phased VCF by gene. The gene coordinates are defined by the GRCh38 GFF3 records.
2. Check for the edge-case of 1 unphased heterozygous genotype in the gene region (see section 2. Notes on the use of vcf2fasta). 
3. Bin genotypes into three different VCFs.  
    a. Pass-filter genotypes with unphased heterozygous genotypes removed  
    b. Fail-filter genotypes  
    c. Pass-filter unphased heterozygous genotypes. These will be reported as they will not be able to be applied during haplotype reconstruction with vcf2fasta  

## 2. Notes on the use of vcf2fasta. 

The vcf2fasta included in hla_resolve/data was originally written by Santiago Sanchez-Ramirez https://github.com/santiagosnchez/vcf2fasta. Unfortunately, it is no longer maintained. Quanyu Chen (yeeus) forked the repository and made some updates https://github.com/yeeus/vcf2fasta. This is the base version that is included under hla_resolve/data/vcf2fasta/, with a few minor edits of my own.

Currently, when vcf2fasta is parsing GFF3 input, it does not concatenate CDS in the correct order for minus (-) strand genes. See this GitHub issue: https://github.com/santiagosnchez/vcf2fasta/issues/23. The easiest solution, which I implemented here, is to manually reverse sort the CDS for minus strand genes in the order they appear in the transcript (i.e., hy coordinate, high to low) before running vcf2fasta. I have manually reverse sorted the minus strand HLA GFF3 genes and included them in hla_resolve/data/hla_gff. In that directory, hla_<gene>.gff3 is the unaltered annotation from GRCh38.p14. hla_<gene>_cds_sorted.gff3 is an altered, CDS-only reverse-sorted that is used as input for vcf2fasta. 

Another bug is that vcf2fasta checks the first genotype of the overall vcf to determine phase status. So if the first variant is unphased, it treats the whole vcf as unphased and won’t output phased sequences in FASTA format. Instead, it will output one diploid sequence with IUPAC codes representing heterozygous bases. For the purpose of hla_resolve, I have made a minor edit to main() in vcf2fasta.py to force it to operate in "phased" mode. The HLA-Resolve pipeline ensures this is appropriate, as it first ensures each gene is completely spanned a haplotype block, and features special treatment for genes that have incomplete phasing.

```diff
- are genotypes phased
- phased = v2f.getPhased(vcf)
- if not phased:
-     print('No phased genotypes found on first variant. Treating as "unphased"')
- else:
-     print('Phased genotypes found on first variant. Treating as "phased"')
+ # Matt's patch: forcibly treat all input as phased
+ phased = True
+ print('Treating as "phased"')
```

When vcf2fasta treats a VCF as "phased," it treats unphased heterozygous genotypes as phased and arbitrarily assigns each allele to a haplotype. This means that unphased heterozygous genotypes must be removed from the phased VCF prior to running vcf2fasta. However, if there is only one heterozygous genotype in the single-gene VCF and it remains unphased, it should not be filtered before vcf2fasta. Because we care only about phased haplotypes within genes, it doesn’t matter which haplotype gets which allele in this scenario. The two reconstructed haplotypes will be correct. However, the standard pipeline removes unphased heterozygous genotypes from the vcf before vcf2fasta, because when operating in “phased’ mode, if there is a mix of phased and unphased heterozygous genotypes, vcf2fasta randomly assigns alleles from the unphased heterozygous genotypes to haplotypes, resulting in inevitable switch errors. I addressed the edge case of a single (unphased) heterozygous genotype in reconstruct_fasta_methods.py.

Matt note to self: If an updated version of vcf2fasta is used in a future version of hla_resolve, I will need to change the files in hla_resolve/data/hla_gff, including make_coords.py.

### Additional vcf2fasta Notes
1. Multiple transcripts exist for all HLA genes in the hg38 gff3. If you don't select one, vcf2fasta will make FASTA files for all annotated transcripts. 

2. As far as I can tell, vcf2fasta applies all variants it sees when generating the FASTA. It does not check the FILTER status of records. Therefore, it is important to filter the VCF prior to running vcf2fasta. 

3. Before running vcf2fasta in phased mode, you need to ensure that the genes you’re making fastas for are fully contained within a single phase set. Vcf2fasta will apply the variants to the hap1 and hap2 fasta, but if there are two or more phase sets contained within your gene, it could cause a switch error.

4. I made minor edits to `vcf2fasta/v2f/functions.py` to address issues with the `getPloidy()` and `getPhased()` functions. Originally, both functions relied on the first variant record to infer ploidy and phasing status. If that first record contained missing genotypes or was of an unexpected format, the program exited with an error.  

    The updated versions now iterate through the VCF until they encounter a valid genotype, allowing assignment of ploidy and phasing information even when the first variant is missing or incomplete.  

    **Functional Changes**

    **1. `getPloidy()`**

    ```diff
    - def getPloidy(vcf):
    -     var = [ y for x,y in next(vcf.fetch()).samples.items() ]
    -     p = [ len(v.get('GT')) for v in var if v.get('GT')[0] is not None ]
    -     return p[0]
    + def getPloidy(vcf):
    +     for rec in vcf.fetch():
    +         for sample in rec.samples.values():
    +             gt = sample.get('GT')
    +             if gt and all(g is not None for g in gt):
    +                 return len(gt)
    +     raise ValueError("No valid genotypes found to infer ploidy.")
    ```

    **2. `getPhased()`**

    ```diff
    - def getPhased(vcf):
    -     var = [ y for x,y in next(vcf.fetch()).samples.items() ]
    -     p = any([ not v.phased for v in var ])
    -     return not p
    + def getPhased(vcf):
    +     for rec in vcf.fetch():
    +         for sample in rec.samples.values():
    +             gt = sample.get('GT')
    +             if gt and all(g is not None for g in gt):
    +                 return sample.phased
    +     raise ValueError("No phased genotypes found to determine phasing.")
    ```
