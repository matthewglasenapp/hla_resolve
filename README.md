# HLA Resolve
HLA typing from long reads

**Authors:** [Matthew Glasenapp](https://github.com/matthewglasenapp), [Alex Symons](https://github.com/FlyingFish800)

```
usage: hla_resolve.py [-h] --input_file INPUT_FILE --sample_name SAMPLE_NAME --platform {pacbio,ont} --output_dir OUTPUT_DIR --aligner
                      {minimap2,vg} [--genotyper {bcftools,clair3,deepvariant}] [--threads THREADS] [--read_group_string READ_GROUP_STRING]

Run HLA-Resolve

options:
  -h, --help            show this help message and exit
  --input_file INPUT_FILE
                        Path to the raw sequencing reads file (default: None)
  --sample_name SAMPLE_NAME
                        Override the parsed sample name (default: None)
  --platform {pacbio,ont}
                        Specify sequencing platform (pacbio, ont) (default: None)
  --output_dir OUTPUT_DIR
                        Output Directory (default: None)
  --aligner {minimap2,vg}
                        Tool for reference genome alignment (default: None)
  --genotyper {bcftools,clair3,deepvariant}
                        Tool for genotyping (default: deepvariant)
  --threads THREADS     Number of threads to use (default: 6)
  --read_group_string READ_GROUP_STRING
                        Override the parsed read group string (default: None)

Examples: python3 script.py --input_file reads.bam --sample_name HG002 --platform pacbio --output_dir out --aligner minimap2 --genotyper
deepvariant --threads 10
```
