# Configuration constants and paths for HLA-Resolve

# HLA genes of interest
genes_of_interest = ("HLA-A", "HLA-B", "HLA-C", "HLA-DRB1", "HLA-DQA1", "HLA-DQB1", "HLA-DPA1", "HLA-DPB1")

# Minimum reads per sample
# DeepVariant is stalling and not exiting for samples with very few BAM records (e.g., HG01891: 35 mapped reads to chr6)
# Set mapped chr6 reads threshold at which variant calling should not proceed
min_reads_sample = 100

# This program is for long-read data only. 
# Require that mean read length is at least 300 bp or higher
min_read_length = 300

# Pangenome graph reference info 
reference_gbz = "/hb/scratch/ogarci12/hybridcapture_pangenome/ref/hprc-v1.0-mc-grch38-minaf.0.1.gbz"
#reference_gbz = "/hb/scratch/mglasena/graph/hprc-v2.0-mc-grch38.gbz"
ref_paths = "/hb/scratch/ogarci12/hybridcapture_pangenome/ref/hprc-v1.0-mc-grch38-minaf.0.1.dict"
#ref_paths = "/hb/scratch/mglasena/graph/grch38_primary.dict"
vg = "/hb/scratch/ogarci12/hybridcapture_pangenome/vg"
mosdepth_regions_file = "/hb/scratch/mglasena/test_hla_resolve/hla_resolve/hla_resolve/hla_genes.bed"
#mosdepth_regions_file = "/hb/scratch/mglasena/test_hla_resolve/hla_resolve/hla_resolve/hla_genes.shifted.bed"


# Transposase mosaic end binding sequence
# The TE sequence (and its reverse complement) introduced during tagmentation still needs to be removed
# Adapters and barcodes were removed by PacBio with lima
me = "AGATGTGTATAAGAGACAG"
me_rc = "CTGTCTCTTATACACATCT"

# Tool paths
longphase = "/hb/home/mglasena/software/longphase/longphase_linux-x64"
prowler_trimmer = "/hb/home/mglasena/software/ProwlerTrimmer/TrimmerLarge.py"
sawfish = "/hb/home/mglasena/software/sawfish-v2.0.3-x86_64-unknown-linux-gnu/bin/sawfish"
clair3_ont_model_path = "/hb/home/mglasena/.conda/envs/clair3/bin/models/r941_prom_sup_g5014"
clair3_hifi_model_path = "/hb/home/mglasena/.conda/envs/clair3/bin/models/hifi_revio"

# Coverage Thresholds
# Might have to relax if you didn't get the flanking regions of the gene (i.e., UTR)
depth_thresh = 10
prop_20x_thresh = 0.0
prop_30x_thresh = 0.0

# Additional paths from reconstruct_fasta_methods.py
vcf2fasta_script = "/hb/scratch/mglasena/vcf2fasta/vcf2fasta.py"
# Reference genome paths - set based on aligner
#reference_genome = "/hb/groups/cornejo_lab/matt/hla_capture/input_data/reference/GCA_000001405.15_GRCh38_no_alt_analysis_set.fa"
#reference_genome_minimap2 = "/hb/scratch/mglasena/test_hla_resolve/hla_resolve/hla_resolve/data/reference/augmented_hg38.fa"
reference_genome_minimap2 = "/hb/scratch/mglasena/test_hla_resolve/hla_resolve/hla_resolve/data/reference/augmented_hg38_drb_alt.fa"
#reference_genome_minimap2 = "/hb/scratch/mglasena/test_hla_resolve/hla_resolve/hla_resolve/data/reference/augmented_hg38_with_long_drb1.fa"
reference_genome_vg = "/hb/groups/cornejo_lab/matt/hla_capture/input_data/reference/hprc-v1.0-chr-renamed.fa"

gff_dir = "/hb/groups/cornejo_lab/matt/hla_capture/input_data/hla_gff"
#gff_dir = "/hb/groups/cornejo_lab/matt/hla_capture/input_data/hla_gff/coord_shift/"
hla_genes_regions_file = "/hb/scratch/mglasena/test_hla_resolve/hla_resolve/hla_resolve/hla_genes.bed"
#hla_genes_regions_file = "/hb/scratch/mglasena/test_hla_resolve/hla_resolve/hla_resolve/hla_genes.shifted.bed"


# Additional paths and constants from investigate_haploblocks_methods.py
genes_bed = "/hb/groups/cornejo_lab/matt/hla_capture/input_data/reference/parse_haploblocks_bed.bed"
#genes_bed = "/hb/groups/cornejo_lab/matt/hla_capture/input_data/reference/parse_haploblocks_bed.shifted.bed"
genes_of_interest_extended = ("HLA-A", "HLA-B", "HLA-C", "HLA-DRB1", "HLA-DQA1", "HLA-DQA2", "HLA-DQB1", "HLA-DQB2", "HLA-DPA1", "HLA-DPB1")

# Extended MHC coordinates
mhc_start = 29555628
mhc_stop = 33409896

# DNA bases and stop codons
DNA_bases = {"A", "T", "G", "C"}
stop_codons = ["TAA", "TAG", "TGA"]

# IPD/IMGT HLA XML file
IMGT_XML = "/hb/scratch/mglasena/test_hla_resolve/hla_resolve/hla_resolve/hla.xml"
