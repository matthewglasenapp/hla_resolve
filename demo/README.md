Input data: PacBio Revio HiFi targeted sequencing reads from HG002 (Ashkenazi Son), a sample from the GIAB and HPRC benchmarks.

```
hla_resolve \
--input_file hla_resolve/demo/HG002.hifi_reads.fastq.gz \
--sample_name HG002 \
--platform pacbio \
--scheme targeted \
--output_dir test \
--trim_adapters \
--adapter_file hla_resolve/demo/adapters.fasta \
--threads 6
```
