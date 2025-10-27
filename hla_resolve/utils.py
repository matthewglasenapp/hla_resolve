import os
import shutil
import sys

def check_required_commands():    
    """Check that all required bioinformatics tools are installed and executable"""
    print("Checking the installation status of the required bioinformatics tools!")

    required_commands = [
        "bam2fastq",
        "bcftools",
        "bgzip",
        "cutadapt",
        "fastplong",
        # FASTQC only for dev
        #"fastqc",
        "java",
        "hiphase",
        "pbmarkdup",
        # pbmm2 deprecated for minimap2
        #"pbmm2",
        # pbsv depreceated for sawfish
        #"pbsv",
        "pigz",
        "samtools",
        # Only required for DeepVariant in dev mode
        #"singularity",
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
