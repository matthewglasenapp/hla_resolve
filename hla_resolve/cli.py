# This software is Copyright ©2026. The Regents of the University of California
# ("Regents"). All Rights Reserved.
#
# See LICENSE.txt for license details.

import textwrap
import argparse
import shlex
import sys
from importlib.metadata import version

class _HelpFormatter(argparse.ArgumentDefaultsHelpFormatter, argparse.RawDescriptionHelpFormatter):
    # Show argument defaults AND keep the epilog's line breaks (raw), so the
    # setup/example blocks don't get rewrapped into one run-on paragraph.
    pass

def main():
    # `hla_resolve setup`: one-time download/build of every external dependency.
    # Run once after install so later array jobs find everything present and never
    # race to fetch or rebuild it. Handled before the main parser, which requires
    # --input_file etc. for an actual typing run.
    if len(sys.argv) >= 2 and sys.argv[1] == "setup":
        from . import config
        config.run_setup()
        return

    # Needed before the parser is built so --help shows the real default.
    # utils imports nothing heavy at module level.
    from .utils import detect_cpus
    available_cpus, cpus_known = detect_cpus()

    parser = argparse.ArgumentParser(
    description="Run HLA-Resolve",
    formatter_class=_HelpFormatter,
    epilog=textwrap.dedent("""\
        One-time setup (downloads references, binaries, and images):
          hla_resolve setup

        Example run:
          hla_resolve --input_file reads.bam --sample_name HG002 --platform pacbio --scheme hybrid_capture --output_dir out

        HLA-Resolve is pre-release software intended for research use only
        and not for use in diagnostic procedures.
    """),
)
    parser.add_argument("--version", action="version", version=f"%(prog)s {version('hla_resolve')}")
    parser.add_argument("--input_file", required=True, help="Path to the raw sequencing reads file")
    parser.add_argument("--sample_name", required=True, help="Name for this sample. Used for output filenames and the read group")
    parser.add_argument("--platform", choices=["pacbio", "ont"], required=True, help="Sequencing platform. Only pacbio is supported; ont is not yet available")
    parser.add_argument("--scheme", choices=["WGS", "WES", "hybrid_capture", "amplicon"], required=True, help="Sequencing scheme")
    parser.add_argument("--output_dir", required=True, help="Output directory. Results are written to <output_dir>/<sample_name>/")
    parser.add_argument("--trim_adapters", action="store_true", help="Enable adapter trimming before processing")
    parser.add_argument("--adapter_file", type=str, required=False, default=None, help="Path to a file with custom adapter sequences (FASTA/FASTQ). If not provided, fastplong auto-detection will be used.")
    parser.add_argument("--threads", type=int, required=False, help="Number of threads to use, lowered to the CPU count when fewer are available", default=min(6, available_cpus))
    parser.add_argument("--read_group_string", required=False, help="Override the parsed read group string", default=None)
    parser.add_argument("--clean_up", action="store_true", help="Keep only the HLA typing results and the run log. Everything else HLA-Resolve wrote, including the haplotagged BAM, the VCFs, and the FASTA haplotypes, is removed at the end of the run. Files it did not write are reported and left in place")
    parser.add_argument("--keep_all_intermediates", action="store_true", help="Retain every intermediate file, including the read copies and superseded BAMs that are otherwise removed as soon as a stage finishes with them. For debugging. Uses several times the storage of a default run, so it is a poor choice for a large cohort")
    parser.add_argument("--keep_full_bam", action="store_true", help="Keep the reference-genome BAM. It is deleted by default once reads are filtered to the MHC")
    parser.add_argument("--clair3_model", type=str, required=False, default=None, help="Clair3 model name (bundled in SIF). Defaults to r1041_e82_400bps_sup_v500 for ONT and hifi_revio for PacBio.")
    parser.add_argument("--verbose", action="store_true", help="Print intermediate file paths and detailed per-variant diagnostic output (overlap suppression, RefCall rescue, unphased het records, CDS sanity check)")
    parser.add_argument("--quiet", action="store_true", help="Print only stage headers, warnings, the final results tables, and the output file paths. The full log is still written to the log file")

    # Show help and exit if no arguments were provided
    if len(sys.argv) == 1:
        parser.print_help()
        parser.exit()

    args = parser.parse_args()

    if args.platform == "ont":
        parser.error("ONT support is not yet available; only --platform pacbio is supported.")

    if args.threads < 1:
        parser.error("--threads must be at least 1")

    if args.clean_up and args.keep_all_intermediates:
        parser.error("--clean_up and --keep_all_intermediates ask for opposite things. Pass one or neither.")

    # Defer heavy imports until after argument parsing so that
    # `hla_resolve` (no args) prints help instantly.
    import time
    import os
    import subprocess
    from . import config
    from .sample_manager import InsufficientReads, Samples, build_workflow_config
    from .utils import announce, check_required_commands, setup_logging, version_string
    from .ont_pipeline import preprocess_ont_sample
    from .pacbio_pipeline import preprocess_pacbio_sample
    from .resolve_alleles_pipeline import resolve_alleles
    from .cleanup import cleanup_intermediate_files, prune_empty_dirs, report_discarded, set_policy

    config.VERBOSE = args.verbose
    config.QUIET = args.quiet and not args.verbose
    config.ACTIVE_STAGES = config.active_stages(args.scheme)

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

    # Provenance header. Repeated on the Finished line so one grep over a cohort
    # of logs shows the version beside each sample.
    version_text = version_string()
    print(version_text)
    print(f"Command: {shlex.join([os.path.basename(sys.argv[0])] + sys.argv[1:])}")
    print()

    # Oversubscription is the user's call, so warn and continue.
    if cpus_known and args.threads > available_cpus:
        announce(f"WARNING: --threads {args.threads} but only {available_cpus} CPUs are allocated. Continuing anyway.")

    # Check that all required tools are installed
    check_required_commands()

    # Set the retention policy before anything writes. Several stages read the
    # user's input file in place rather than copying it, so it is marked
    # protected and no discard can reach it.
    set_policy(keep_all_intermediates=args.keep_all_intermediates,
               protected_paths=[args.input_file])
    
    start_time = time.time()
    
    # A bad or too-small input is an expected outcome, not a crash. Report it the
    # same way as a run that fails later: a status line and exit 1.
    try:
        sample = Samples(input_file=args.input_file, sample_name=args.sample_name, platform=args.platform, output_dir=args.output_dir, aligner=args.aligner, snp_caller=args.snp_caller, indel_caller=args.indel_caller, trim_adapters=args.trim_adapters, adapter_file=args.adapter_file, threads=args.threads, read_group_string=args.read_group_string, clean_up=args.clean_up, scheme=args.scheme, clair3_model=args.clair3_model, rescue_refcalls=args.rescue_refcalls, keep_full_bam=args.keep_full_bam, keep_all_intermediates=args.keep_all_intermediates)
    except (InsufficientReads, FileNotFoundError, OSError, ValueError) as err:
        status = "insufficient_reads" if isinstance(err, InsufficientReads) else "input_error"
        announce(f"ERROR: {err}")
        announce(f"Finished {args.sample_name} (status: {status}) [{version_text}]")
        sys.exit(1)

    # Build workflow configuration from sample object
    workflow_config = build_workflow_config(sample)
    
    # An external tool that dies takes the whole run with it. Report which command
    # failed rather than a traceback, and never leave a partial result looking ok.
    try:
        reads = None
        if workflow_config['platform'] == "PACBIO":
            reads = preprocess_pacbio_sample(config=workflow_config)
        elif workflow_config['platform'] == "ONT":
            preprocess_ont_sample(config=workflow_config)

        # Check if variant calling was successful before proceeding to HLA resolution
        state = None
        if os.path.exists(workflow_config['snv_vcf']):
            state = resolve_alleles(config=workflow_config)
            status = "ok"
        else:
            print(f"Skipping HLA allele resolution for {workflow_config['sample_ID']} due to insufficient reads for variant calling")
            status = "insufficient_reads"
    except subprocess.CalledProcessError as err:
        announce(f"ERROR: a command exited with code {err.returncode}")
        if err.returncode in (137, -9):
            announce("That code means the process was killed. The usual cause is the out-of-memory killer, so raise the job's memory or lower --threads.")
        announce(f"  {str(err.cmd).replace('set -o pipefail; ', '', 1)}")
        announce(f"Finished {workflow_config['sample_ID']} (status: tool_failed) [{version_text}]")
        sys.exit(1)

    # Directories emptied by the incremental discards above serve no purpose.
    prune_empty_dirs(config=workflow_config)
    report_discarded()

    # --clean_up: keep the typing results and nothing else
    cleanup_intermediate_files(config=workflow_config)

    elapsed_time = time.time() - start_time

    minutes, seconds = divmod(elapsed_time, 60)
    announce(f"Finished {workflow_config['sample_ID']} in {int(minutes)}m {seconds:.0f}s (status: {status}) [{version_text}]")

    if status != "ok":
        sys.exit(1)

if __name__ == "__main__":
    main()
