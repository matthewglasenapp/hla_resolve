# This software is Copyright ©2026. The Regents of the University of California
# ("Regents"). All Rights Reserved.
#
# See LICENSE.txt for license details.

import os
import shutil
import subprocess
import sys
import time


def run_quiet(cmd):
    """Run a shell command with check=True. Stderr/stdout are discarded on
    success; on non-zero exit the captured streams are re-emitted (stdout to
    sys.stdout, stderr to sys.stderr — both of which go through TeeStream into
    the log file) and CalledProcessError is re-raised. In --verbose mode,
    streams pass through unfiltered so tool output reaches the terminal."""
    from . import config
    # pipefail, run under bash. /bin/sh returns only the last command's status, so
    # a process killed mid-pipe (OOM, usually) hands truncated output to the next
    # stage and the whole command still reports success.
    cmd = f"set -o pipefail; {cmd}"
    if config.VERBOSE:
        return subprocess.run(cmd, shell=True, check=True, executable="/bin/bash")
    try:
        return subprocess.run(cmd, shell=True, check=True, capture_output=True, text=True,
                              executable="/bin/bash")
    except subprocess.CalledProcessError as e:
        if e.stdout:
            print(e.stdout, end="")
        if e.stderr:
            print(e.stderr, end="", file=sys.stderr)
        raise


def version_string():
    # importlib.metadata reports the version recorded at install time, which goes
    # stale after a bare git pull. Append the commit when the package sits in a
    # checkout, so the log records the code that actually ran.
    from importlib.metadata import version
    text = f"hla_resolve {version('hla_resolve')}"
    repo = os.path.dirname(os.path.dirname(os.path.abspath(__file__)))
    try:
        commit = subprocess.run(["git", "-C", repo, "rev-parse", "--short", "HEAD"],
                                capture_output=True, text=True, timeout=5).stdout.strip()
        if commit:
            # Untracked files are excluded, or a stray output file would mark
            # every run dirty.
            modified = subprocess.run(["git", "-C", repo, "status", "--porcelain", "-uno"],
                                      capture_output=True, text=True, timeout=5).stdout.strip()
            text += f" ({commit}-dirty)" if modified else f" ({commit})"
    except Exception:
        pass
    return text


def detect_cpus():
    # Returns (count, confident). Only Slurm's own numbers are trusted enough to
    # warn on. The affinity mask reports the whole node wherever cores are not
    # constrained by cgroups, and cpu_count always does.
    for variable in ("SLURM_CPUS_PER_TASK", "SLURM_CPUS_ON_NODE"):
        value = os.environ.get(variable, "")
        if value.isdigit() and int(value) > 0:
            return int(value), True
    if hasattr(os, "sched_getaffinity"):
        return len(os.sched_getaffinity(0)), False
    return os.cpu_count() or 1, False


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


def detail(*args, **kwargs):
    # Plumbing that belongs in the log but not on a user's terminal: paths of
    # intermediate files, tool input arguments. --verbose puts it back on screen.
    from . import config
    if config.VERBOSE:
        print(*args, **kwargs)
        return
    log_file = getattr(sys, "_hla_resolve_log_file", None)
    if log_file is not None:
        print(*args, **kwargs, file=log_file)
    elif not config.QUIET:
        print(*args, **kwargs)


_stage_start = None


def stage(name):
    # Numbered from config.ACTIVE_STAGES so the count matches the stages this run
    # performs. A scheme that skips one gets a shorter list, never a gap.
    from . import config
    finish_stage()
    stages = config.ACTIVE_STAGES or config.STAGES
    total = len(stages)
    number = stages.index(name) + 1 if name in stages else config.STAGES.index(name) + 1
    announce(f"\n[{number}/{total}] {name}")
    global _stage_start
    _stage_start = time.time()


def finish_stage():
    """Report how long the running stage took. Called when the next stage starts,
    and once by the pipeline when the last stage ends."""
    global _stage_start
    if _stage_start is None:
        return
    elapsed = time.time() - _stage_start
    _stage_start = None
    minutes, seconds = divmod(elapsed, 60)
    if minutes:
        announce(f"Finished in {int(minutes)}m {seconds:.0f}s")
    elif seconds < 0.1:
        # Rounding a fast stage to 0.0s reads as a stopped clock.
        announce("Finished in <0.1s")
    else:
        announce(f"Finished in {seconds:.1f}s")


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
    print()
    return log_path

def check_required_commands():    
    """Check that all required bioinformatics tools are installed and executable"""
    print("Checking for the required bioinformatics tools...")

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
        print(f"ERROR: Missing the following commands: {', '.join(missing_commands)}")
        print("Run `hla_resolve setup` to download and build the external dependencies.")
        print("If setup has already been run and tools are still missing, check that the")
        print("hla_resolve conda environment is active, then open an issue at")
        print("https://github.com/matthewglasenapp/hla_resolve/issues")
        sys.exit(1)
    else:
        print("All required tools are installed")
        print()
