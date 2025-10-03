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
reference_genome_minimap2 = "/hb/scratch/mglasena/test_hla_resolve/hla_resolve/hla_resolve/data/reference/augmented_hg38.fa"
#reference_genome_minimap2 = "/hb/scratch/mglasena/test_hla_resolve/hla_resolve/hla_resolve/data/reference/augmented_hg38_drb_alt.fa"
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

# GRCh38 HLA gene antigen recognition sequence coordinates
# BED (0-based coordinates)
ARS_dict = {
	"HLA-A": "chr6:29942756-29943543",
	"HLA-B": "chr6:31356166-31356957",
	"HLA-C": "chr6:31271072-31271868",
	"HLA-DPA1": "chr6:33069640-33069886",
	"HLA-DPB1": "chr6:33080671-33080935",
	"HLA-DQA1": "chr6:32641309-32641558",
	"HLA-DQB1": "chr6:32664797-32665067",
	"HLA-DRB1": "chr6:32584108-32584378"
}

gene_dict = {
	"HLA-A":    (29941260, 29949572),
	"HLA-B":    (31353872, 31367067),
	"HLA-C":    (31268749, 31272130),
	"HLA-DPA1": (33064569, 33080775),
	"HLA-DPB1": (33075990, 33089696),
	"HLA-DQA1": (32628179, 32647062),
	"HLA-DQB1": (32659467, 32668383),
	"HLA-DRB1": (32577902, 32589848)
}

CDS_dict = {
	'HLA-A': [[29942554, 29942626], [29942757, 29943026], [29943268, 29943543], [29944122, 29944397], [29944500, 29944616], [29945059, 29945091], [29945234, 29945281], [29945451, 29945455]],
	'HLA-B': [[31354483, 31354526], [31354633, 31354665], [31355107, 31355223], [31355317, 31355592], [31356167, 31356442], [31356688, 31356957], [31357086, 31357158]],
	'HLA-C': [[31269169, 31269173], [31269338, 31269385], [31269493, 31269525], [31269966, 31270085], [31270210, 31270485], [31271073, 31271348], [31271599, 31271868], [31271999, 31272071]],
	'HLA-DPA1': [[33068650, 33068804], [33069019, 33069300], [33069641, 33069886], [33073471, 33073570]], 
	'HLA-DPB1': [[33076042, 33076141], [33080672, 33080935], [33084950, 33085231], [33085779, 33085889], [33086219, 33086238]], 
	'HLA-DQA1': [[32637459, 32637540], [32641310, 32641558], [32641972, 32642253], [32642610, 32642764]], 
	'HLA-DQB1': [[32660236, 32660249], [32660859, 32660882], [32661347, 32661457], [32661967, 32662248], [32664798, 32665067], [32666499, 32666607]], 
	'HLA-DRB1': [[32579091, 32579104], [32580247, 32580270], [32580746, 32580856], [32581557, 32581838], [32584109, 32584378], [32589643, 32589742]]
}