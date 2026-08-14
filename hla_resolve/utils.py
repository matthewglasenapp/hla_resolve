# This software is Copyright ©2026. The Regents of the University of California
# ("Regents"). All Rights Reserved.
#
# See LICENSE.txt for license details.

import os
import shutil
import subprocess
import sys


def run_quiet(cmd):
    """Run a shell command with check=True. Stderr/stdout are discarded on
    success; on non-zero exit the captured streams are re-emitted (stdout to
    sys.stdout, stderr to sys.stderr — both of which go through TeeStream into
    the log file) and CalledProcessError is re-raised. In --verbose mode,
    streams pass through unfiltered so tool output reaches the terminal."""
    from . import config
    if config.VERBOSE:
        return subprocess.run(cmd, shell=True, check=True)
    try:
        return subprocess.run(cmd, shell=True, check=True, capture_output=True, text=True)
    except subprocess.CalledProcessError as e:
        if e.stdout:
            print(e.stdout, end="")
        if e.stderr:
            print(e.stderr, end="", file=sys.stderr)
        raise


class TeeStream:
    def __init__(self, *streams):
        self.streams = streams

    def write(self, message):
        for stream in self.streams:
            stream.write(message)
        for stream in self.streams:
            stream.flush()

    def flush(self):
        for stream in self.streams:
            stream.flush()

    def isatty(self):
        return any(getattr(stream, "isatty", lambda: False)() for stream in self.streams)

    @property
    def encoding(self):
        for stream in self.streams:
            encoding = getattr(stream, "encoding", None)
            if encoding:
                return encoding
        return "utf-8"


def announce(*args, **kwargs):
    # Reaches the terminal even under --quiet, and always the log. Without
    # --quiet sys.stdout already tees both, so one print is enough.
    from . import config
    print(*args, **kwargs)
    if config.QUIET and hasattr(sys, "_hla_resolve_stdout"):
        print(*args, **kwargs, file=sys._hla_resolve_stdout)


def stage(name):
    # Numbered from config.STAGES so the count matches the documented workflow.
    from . import config
    total = len(config.STAGES)
    number = config.STAGES.index(name) + 1
    announce(f"\n[{number}/{total}] {name}")


def setup_logging(output_dir, sample_name=None):
    os.makedirs(output_dir, exist_ok=True)
    log_basename = f"{sample_name}.hla_resolve.log" if sample_name else "hla_resolve.log"
    log_path = os.path.join(output_dir, log_basename)

    if getattr(sys, "_hla_resolve_log_path", None) == log_path:
        return log_path

    log_file = open(log_path, "w", buffering=1)

    if not hasattr(sys, "_hla_resolve_stdout"):
        sys._hla_resolve_stdout = sys.stdout
        sys._hla_resolve_stderr = sys.stderr

    from . import config
    if config.QUIET:
        # Everything still reaches the log. Only stage headers, warnings and the
        # results table reach the terminal, via announce().
        sys.stdout = TeeStream(log_file)
    else:
        sys.stdout = TeeStream(sys._hla_resolve_stdout, log_file)
    sys.stderr = TeeStream(sys._hla_resolve_stderr, log_file)

    sys._hla_resolve_log_path = log_path
    sys._hla_resolve_log_file = log_file

    print(f"Logging to {log_path}")
    print("\n")
    return log_path

def check_required_commands():    
    """Check that all required bioinformatics tools are installed and executable"""
    print("Checking the installation status of the required bioinformatics tools!")

    required_commands = [
        "bam2fastq",
        "bcftools",
        "bgzip",
        "freebayes",
        "cutadapt",
        "fastplong",
        # FASTQC only for dev
        #"fastqc",
        "java",
        "hiphase",
        "pbmarkdup",
        "pbsv",
        "pigz",
        "samtools",
        "singularity",
        "sniffles",
        "tabix",
        "trgt",
    ]

    missing_commands = []
    for command in required_commands:
        if shutil.which(command) is None:
            missing_commands.append(command)
    
    if len(missing_commands) != 0:
        print(f"Error: Missing the following commands: {', '.join(missing_commands)}")
        print("Run `hla_resolve setup` to download and build the external dependencies.")
        print("If setup has already been run and tools are still missing, check that the")
        print("hla_resolve conda environment is active, then open an issue at")
        print("https://github.com/matthewglasenapp/hla_resolve/issues")
        sys.exit(1)
    else:
        print("All tools required are installed!")
        print("\n")
