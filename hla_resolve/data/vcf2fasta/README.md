# vcf2fasta
**Major release. Current version should now work with haploid, diploid, phased, and unphased (IUPAC) outputs.**

vcf2fasta.py is Python program that extracts FASTA alignments from VCF files given a GFF file with feature coordinates and a reference FASTA file.

## Note

**I just editted some codes for my own purpose. It now should work with bed or gff file only. The vcf file I used was generated from MC for pangenome purpose and therefore samples are phased and haploid.**

There should be several bugs in this version, but that worked for me haha. So if anyone has other problems or custom purposes, we could discuss them together or you can fork this repo and edit it anyway.
### update in 2/12/2025
Glad to see people working on this script! As what I said in this [issue](https://github.com/santiagosnchez/vcf2fasta/issues/23), I would update this project in the near future.

## Preprocessing

The reference must be indexed using:

```
samtools faidx ref.fa
```

And the VCF file should be tabix indexed and compressed:

```
(bgzip variants.vcf)
tabix variants.vcf(.gz) (or bcftools index variants.vcf(.gz))
```

For most GFF3 formats, no modification is needed for the GFF file if the structure follow [Ensembl](https://m.ensembl.org/info/website/upload/gff3.html). However, it is important to keep the whole structure of the GFF file, including complete gene features. If CDSs are the focus they should be accompanied by it's corresponding gene or parent feature:

* gene
* CDS/exon

Similarly of the focus are introns:

* gene
* intron

Or transcripts:

* gene
* transcript

etc.. Alternatively, all features can be left on the GFF. However, the `--feat | -e` argument must be used at all times.

**If multiple transcript isoforms are on the GFF, all of them will be fetched.**

## Requirements
* `pysam`
* `art`

```bash
pip3 install pysam art
```

## Options
Run with `-h` option for more details

```
usage: vcf2fasta.py [-h] --fasta GENOME --vcf VCF [--gff GFF] [--bed BED] [--feat FEAT] [--blend] [--inframe] [--out OUT] [--addref] [--skip] [--remove-gap]

        Converts regions/intervals in the genome into FASTA alignments
        provided a VCF file, a GFF file, and FASTA reference.

options:
  -h, --help            show this help message and exit
  --fasta GENOME, -f GENOME
                        FASTA file with the reference genome.
  --vcf VCF, -v VCF     VCF file.
  --gff GFF, -g GFF     GFF file.
  --bed BED             BED file. BED6 is preferred while the fifth column is phase value which can be used for option "--inframe". Header beginning with "#" will be ignored.
  --feat FEAT, -e FEAT  feature/annotation in the GFF file. (i.e. gene, CDS, intron). Only for gff input.
  --blend, -b           concatenate GFF entries of FEAT into a single alignment. Useful for CDS. (default: False)
  --inframe, -i         force the first codon of the sequence to be inframe. Useful for incomplete CDS. (default: False)
  --out OUT, -o OUT     provide a name for the output directory (optional)
  --addref, -r          include the reference sequence in the FASTA alignment (default: False)
  --skip, -s            skips features without variants (default: False)
  --remove-gap, -rg     remove gaps in the final alignments (default: False)

        Before running the code make sure
        that your reference FASTA file is indexed:

        samtools faidx genome.fa

        BGZIP compress and index your VCF file:

        (bgzip variants.vcf)
        tabix variants.vcf(.gz) (or bcftools index variants.vcf(.gz))

        The GFF file does not need to be indexed.

        examples:
        python vcf2fasta.py -f genome.fa -v variants.vcf.gz -g intervals.gff -e five_prime_UTR -r -rg
```

If running on CDS use the `--blend | -b` option to concatenate coding sequences. Otherwise it will spit out FASTA alignments for each CDS.
