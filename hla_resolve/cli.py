# This software is Copyright ©2026. The Regents of the University of California
# ("Regents"). All Rights Reserved.
#
# See LICENSE.txt for license details.

import textwrap
import argparse
import sys
from importlib.metadata import version

def main():
    # `hla_resolve setup`: one-time download/build of every external dependency.
    # Run once after install so later array jobs find everything present and never
    # race to fetch or rebuild it. Handled before the main parser, which requires
    # --input_file etc. for an actual typing run.
    if len(sys.argv) >= 2 and sys.argv[1] == "setup":
        from . import config
        config.run_setup()
        return

    parser = argparse.ArgumentParser(
    description="Run HLA-Resolve",
    formatter_class=argparse.ArgumentDefaultsHelpFormatter,
    epilog=textwrap.dedent("""\
        One-time setup (downloads references, binaries, and images):
          hla_resolve setup

        Example run:
          hla_resolve --input_file reads.bam --sample_name HG002 --platform pacbio --scheme hybrid_capture --output_dir out --threads 10

        HLA-Resolve is pre-release software intended for research use only
        and not for use in diagnostic procedures.
    """),
)
    parser.add_argument("--version", action="version", version=f"%(prog)s {version('hla_resolve')}")
    parser.add_argument("--input_file", required=True, help="Path to the raw sequencing reads file")
    parser.add_argument("--sample_name", required=True, help="Override the parsed sample name", default=None)
    parser.add_argument("--platform", choices=["pacbio", "ont"], required=True, help="Specify sequencing platform (pacbio, ont)")
    parser.add_argument("--scheme", choices=["WGS", "WES", "hybrid_capture", "amplicon"], required=True, help="Sequencing scheme")
    parser.add_argument("--output_dir", required=True, help="Output Directory", default=None)
    parser.add_argument("--trim_adapters", action="store_true", help="Enable adapter trimming before processing")
    parser.add_argument("--adapter_file", type=str, required=False, default=None, help="Path to a file with custom adapter sequences (FASTA/FASTQ). If not provided, default adapters will be used.")
    parser.add_argument("--threads", type=int, required=False, help="Number of threads to use", default=6)
    parser.add_argument("--read_group_string", required=False, help="Override the parsed read group string", default=None)
    parser.add_argument("--clean-up", action="store_true", help="Remove intermediate files")
    parser.add_argument("--clair3_model", type=str, required=False, default=None, help="Clair3 model name (bundled in SIF). Defaults to r1041_e82_400bps_sup_v500 for ONT and hifi_revio for PacBio.")
    parser.add_argument("--verbose", action="store_true", help="Print detailed per-variant diagnostic output (overlap suppression, RefCall rescue, unphased het records, CDS sanity check)")

    # Show help and exit if no arguments were provided
    if len(sys.argv) == 1:
        parser.print_help()
        parser.exit()

    args = parser.parse_args()

    if args.platform == "ont":
        parser.error("ONT support is not yet available; only --platform pacbio is supported.")

    # Defer heavy imports until after argument parsing so that
    # `hla_resolve` (no args) prints help instantly.
    import time
    import os
    from . import config
    from .sample_manager import Samples, build_workflow_config
    from .utils import check_required_commands, setup_logging
    from .ont_pipeline import preprocess_ont_sample
    from .pacbio_pipeline import preprocess_pacbio_sample
    from .resolve_alleles_pipeline import resolve_alleles
    from .cleanup import cleanup_intermediate_files

    config.VERBOSE = args.verbose
    
    args.aligner = "rammap"
    if args.platform == "ont":
        args.snp_caller = "clair3"
        args.indel_caller = "clair3"
        args.rescue_refcalls = False
    else:
        args.snp_caller = "bcftools"
        args.indel_caller = "deepvariant"
        args.rescue_refcalls = True

    setup_logging(output_dir=args.output_dir, sample_name=args.sample_name)

    # Check that all required tools are installed
    check_required_commands()
    
    start_time = time.time()
    
    sample = Samples(input_file=args.input_file, sample_name=args.sample_name, platform=args.platform, output_dir=args.output_dir, aligner=args.aligner, snp_caller=args.snp_caller, indel_caller=args.indel_caller, trim_adapters=args.trim_adapters, adapter_file=args.adapter_file, threads=args.threads, read_group_string=args.read_group_string, clean_up=args.clean_up, scheme=args.scheme, clair3_model=args.clair3_model, rescue_refcalls=args.rescue_refcalls)

    # Build workflow configuration from sample object
    workflow_config = build_workflow_config(sample)
    
    if workflow_config['platform'] == "PACBIO":
        preprocess_pacbio_sample(config=workflow_config)
    elif workflow_config['platform'] == "ONT":
        preprocess_ont_sample(config=workflow_config)
    
    # Check if variant calling was successful before proceeding to HLA resolution
    if os.path.exists(workflow_config['snv_vcf']):
        resolve_alleles(config=workflow_config)
    else:
        print(f"Skipping HLA allele resolution for {workflow_config['sample_ID']} due to insufficient reads for variant calling")
    
    # Clean up intermediate files if requested
    cleanup_intermediate_files(config=workflow_config)

    end_time = time.time()
    elapsed_time = end_time - start_time
    minutes, seconds = divmod(elapsed_time,60)
    print(f"Processed sample in {int(minutes)}:{seconds:.2f}!")

if __name__ == "__main__":
    main()
