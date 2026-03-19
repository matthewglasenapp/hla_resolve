# This software is Copyright ©2026. The Regents of the University of California
# ("Regents"). All Rights Reserved.
#
# Licensed under the UC Santa Cruz Noncommercial License (the "License").
# You may not use this file except in compliance with the License.
# A copy of the License is included in this repository as LICENSE.txt.

import os
import shutil
import sys


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


def setup_logging(output_dir, sample_name=None):
    os.makedirs(output_dir, exist_ok=True)
    log_basename = f"{sample_name}.hla_resolve.log" if sample_name else "hla_resolve.log"
    log_path = os.path.join(output_dir, log_basename)

    if getattr(sys, "_hla_resolve_log_path", None) == log_path:
        return log_path

    log_file = open(log_path, "a", buffering=1)

    if not hasattr(sys, "_hla_resolve_stdout"):
        sys._hla_resolve_stdout = sys.stdout
        sys._hla_resolve_stderr = sys.stderr

    sys.stdout = TeeStream(sys._hla_resolve_stdout, log_file)
    sys.stderr = TeeStream(sys._hla_resolve_stderr, log_file)

    sys._hla_resolve_log_path = log_path
    sys._hla_resolve_log_file = log_file

    print(f"Logging to {log_path}")
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
        sys.exit(1)
    else:
        print("All tools required are installed!")
        print("\n\n")
