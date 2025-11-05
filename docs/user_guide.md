# HLA-Resolve Documentation

Documentation for HLA-Resolve.

## Notes on the use of vcf2fasta. 

The vcf2fasta included in hla_resolve/data was originally written by Santiago Sanchez-Ramirez https://github.com/santiagosnchez/vcf2fasta. Unfortunately, it is no longer maintained. Quanyu Chen (yeeus) forked the repository and made some updates https://github.com/yeeus/vcf2fasta. This is the base version that is included under hla_resolve/data/vcf2fasta/, with a few minor edits of my own.

Currently, when vcf2fasta is parsing GFF3 input, it does not concatenate CDS in the correct order for minus (-) strand genes. See this GitHub issue: https://github.com/santiagosnchez/vcf2fasta/issues/23. The easiest solution, which I implemented here, is to manually reverse sort the CDS for minus strand genes in the order they appear in the transcript (i.e., hy coordinate, high to low) before running vcf2fasta. I have manually reverse sorted the minus strand HLA GFF3 genes and included them in hla_resolve/data/hla_gff. In that directory, hla_<gene>.gff3 is the unaltered annotation from GRCh38.p14. hla_<gene>_cds_sorted.gff3 is an altered, CDS-only reverse-sorted that is used as input for vcf2fasta. 

Another bug is that vcf2fasta checks the first genotype of the overall vcf to determine phase status. So if the first variant is unphased, it treats the whole vcf as unphased and won’t output phased sequences in FASTA format. Instead, it will output one diploid sequence with IUPAC codes representing heterozygous bases. To work around this, hla_resolve/reconstruct_fasta_methods.py filters the input VCF to begin with the first phased variant (see filter_vcf()). 

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
