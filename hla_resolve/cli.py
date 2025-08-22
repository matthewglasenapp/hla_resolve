import time
import textwrap
import argparse
import sys
from sample import Samples
from config import *
from utils import check_required_commands
from ont_workflow import preprocess_ont_sample
from pacbio_workflow import preprocess_pacbio_sample
from resolve_alleles import resolve_alleles

def main():
    parser = argparse.ArgumentParser(
    description="Run HLA-Resolve",
    formatter_class=argparse.ArgumentDefaultsHelpFormatter,
    epilog=textwrap.dedent("""\
        Examples:
          python3 -m hla_resolve.cli --input_file reads.bam --sample_name HG002 --platform pacbio --output_dir out --aligner minimap2 --genotyper deepvariant --threads 10
    """),
)
    parser.add_argument("--input_file", required=True, help="Path to the raw sequencing reads file")
    parser.add_argument("--sample_name", required=True, help="Override the parsed sample name", default=None)
    parser.add_argument("--platform", choices=["pacbio", "ont"], required=True, help="Specify sequencing platform (pacbio, ont)")
    parser.add_argument("--output_dir", required=True, help="Output Directory", default=None)
    parser.add_argument("--aligner", choices=["minimap2", "vg"], required=True, help="Tool for reference genome alignment", default=None)
    parser.add_argument("--genotyper", choices=["bcftools", "clair3", "deepvariant"], required=False, help="Tool for genotyping", default="deepvariant")
    parser.add_argument("--trim_adapters", action="store_true", help="Enable adapter trimming before processing")
    parser.add_argument("--adapter_file", type=str, required=False, default=None, help="Path to a file with custom adapter sequences (FASTA/FASTQ). If not provided, default adapters will be used.")
    parser.add_argument("--threads", type=int, required=False, help="Number of threads to use", default=6)
    parser.add_argument("--read_group_string", required=False, help="Override the parsed read group string", default=None)
    parser.add_argument("--clean-up", action="store_true", help="Remove intermediate files")
    
    # Show help and exit if no arguments were provided
    if len(sys.argv) == 1:
        parser.print_help()
        parser.exit()

    args = parser.parse_args()

    # Check that all required tools are installed
    check_required_commands()
    
    start_time = time.time()
    
    sample = Samples(input_file=args.input_file, sample_name=args.sample_name, platform=args.platform, output_dir=args.output_dir, aligner=args.aligner, genotyper=args.genotyper, trim_adapters=args.trim_adapters, adapter_file=args.adapter_file, threads=args.threads, read_group_string=args.read_group_string, clean_up=args.clean_up)

    if sample.platform == "PACBIO":
        preprocess_pacbio_sample(
            sample=sample,
            reference_fasta=Samples.reference_fasta,
            vg=vg,
            reference_gbz=reference_gbz,
            ref_paths=ref_paths,
            deepvariant_sif=Samples.deepvariant_sif,
            chr6_bed=Samples.chr6_bed,
            clair3_ont_model_path=clair3_ont_model_path,
            clair3_hifi_model_path=clair3_hifi_model_path,
            sawfish=sawfish,
            pbtrgt_repeat_file=Samples.pbtrgt_repeat_file
        )
    elif sample.platform == "ONT":
        preprocess_ont_sample(
            sample=sample,
            reference_fasta=Samples.reference_fasta,
            vg=vg,
            reference_gbz=reference_gbz,
            ref_paths=ref_paths,
            deepvariant_sif=Samples.deepvariant_sif,
            chr6_bed=Samples.chr6_bed,
            clair3_ont_model_path=clair3_ont_model_path,
            clair3_hifi_model_path=clair3_hifi_model_path,
            longphase=longphase,
            tandem_repeat_bed=Samples.tandem_repeat_bed
        )
    
    resolve_alleles(
        sample=sample,
        mosdepth_regions_file=mosdepth_regions_file,
        depth_thresh=depth_thresh,
        prop_20x_thresh=prop_20x_thresh,
        prop_30x_thresh=prop_30x_thresh,
        mhc_start=mhc_start,
        mhc_stop=mhc_stop,
        genes_bed=genes_bed,
        genes_of_interest=genes_of_interest_extended,
        hla_genes_regions_file=hla_genes_regions_file,
        vcf2fasta_script=vcf2fasta_script,
        reference_genome=reference_genome,
        DNA_bases=DNA_bases,
        stop_codons=stop_codons,
        IMGT_XML=IMGT_XML
    )

    end_time = time.time()
    elapsed_time = end_time - start_time
    minutes, seconds = divmod(elapsed_time,60)
    print(f"Processed sampled in {int(minutes)}:{seconds:.2f}!")

if __name__ == "__main__":
    main()
