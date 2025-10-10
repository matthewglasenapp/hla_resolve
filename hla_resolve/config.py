# Configuration constants and paths for HLA-Resolve
import os

# Get the data directory relative to this config file
_data_dir = os.path.join(os.path.dirname(os.path.abspath(__file__)), "data")

# HLA genes of interest for HLA typing
genes_of_interest = ("HLA-A", "HLA-B", "HLA-C", "HLA-DPA1", "HLA-DPB1", "HLA-DQA1", "HLA-DQB1", "HLA-DRB1")

# Small dummy reference used to identify raw sequencing reads originating from DRB3 and DRB4 genes
# This reference is used in the bait_DRB_paralogs() function of preprocess_methods.py
# This reference consists of a contig for each of DRB1, DRB3 and DRB4
# The contigs represent exon2 +/- 2kb, with the 270 bases of exon 2 hardmasked with "N"
dummy_reference = os.path.join(_data_dir, "reference/DRB_1_3_4.fa")

# Paths to several software tools installed from source
longphase = "/hb/home/mglasena/software/longphase/longphase_linux-x64"
prowler_trimmer = "/hb/home/mglasena/software/ProwlerTrimmer/TrimmerLarge.py"
sawfish = "/hb/home/mglasena/software/sawfish-v2.0.3-x86_64-unknown-linux-gnu/bin/sawfish"
vcf2fasta_script = "/hb/scratch/mglasena/vcf2fasta/vcf2fasta.py"
vg = "/hb/scratch/ogarci12/hybridcapture_pangenome/vg"
picard = "/hb/home/mglasena/software/picard/picard.jar"

# Paths to CLAR3 models for ONT and HiFi data
# Consequence of installing CLAR3 into a conda environment
clair3_ont_model_path = "/hb/home/mglasena/.conda/envs/clair3/bin/models/r941_prom_sup_g5014"
clair3_hifi_model_path = "/hb/home/mglasena/.conda/envs/clair3/bin/models/hifi_revio"

# Paths to reference genome for reference genome alignment with minimap2

# Reference fasta with added HLA-OLI/HLA-Y contig for baiting out reads originating from HLA-Y
reference_genome_minimap2 = os.path.join(_data_dir, "reference/augmented_hg38_chr6.fa")
#reference_genome_minimap2 = os.path.join(_data_dir, "reference/augmented_hg38_drb_alt.fa")
#reference_genome_minimap2 = os.path.join(_data_dir, "reference/augmented_hg38_with_long_drb1.fa")

# Paths to reference genome files for genome alignment to HPRC pangenome graph with vg giraffe
# This approach is experimental and is not implemented in the current release
# These files are only used when aligner == "vg"
# Headers renamed using the following command: 
# sed 's/^>GRCh38\.chr/>chr/' /hb/scratch/ogarci12/hybridcapture_pangenome/ref/hprc-v1.0-mc-grch38-minaf.0.1.fa > /hb/groups/cornejo_lab/matt/hla_capture/input_data/reference/hprc-v1.0-chr-renamed.fa
reference_genome_vg = "/hb/groups/cornejo_lab/matt/hla_capture/input_data/reference/hprc-v1.0-chr-renamed.fa"
reference_genome_vg_gbz = "/hb/scratch/ogarci12/hybridcapture_pangenome/ref/hprc-v1.0-mc-grch38-minaf.0.1.gbz"
#reference_genome_vg_gbz = "/hb/scratch/mglasena/graph/hprc-v2.0-mc-grch38.gbz"
reference_genome_vg_paths = "/hb/scratch/ogarci12/hybridcapture_pangenome/ref/hprc-v1.0-mc-grch38-minaf.0.1.dict"
#reference_genome_vg_paths = "/hb/scratch/mglasena/graph/grch38_primary.dict"

# DeepVariant SIF file for variant calling with DeepVariant
deepvariant_sif = os.path.join(_data_dir, "deepvariant_sif/deepvariant.sif")

# GRCh38 tandem repeat mask file for pbsv
# Downloaded from https://github.com/PacificBiosciences/pbsv/blob/master/annotations/human_GRCh38_no_alt_analysis_set.trf.bed
tandem_repeat_bed = os.path.join(_data_dir, "repeats_bed/human_GRCh38_no_alt_analysis_set.trf.bed")

# BED file of coordinates for chr6
# Used to constrain variant calling with Clair3 and Sniffles to chr6
# The file is simply chr6\t28000000\t34000000
chr6_bed = os.path.join(_data_dir, "reference/chr6.bed")

# GRCh38 tandem repeat definition file for pbtrgt
# Downloaded from https://zenodo.org/records/8329210
pbtrgt_repeat_file = os.path.join(_data_dir, "repeats_bed/chr6_polymorphic_repeats.hg38.bed")

# GFF files of the HLA genes of interest for FASTA reconstruction with vcf2fasta
# Used in reconstruct_fasta_methods.py
gff_dir = os.path.join(_data_dir, "hla_gff")
#gff_dir = os.path.join(_data_dir, "hla_gff/coord_shift")

# BED file of coordinates for the 8 HLA genes of interest
# The coordinates were taken from the GRCh38 GFF3 file
# Used for both mosdepth coverage analysis and VCF filtering
# HLA-B and HLA-DQA1 coordinates were slightly modified from the raw GFF3 file to exclude exons that are not part of the MANE Select transcript
# The second, commented-out bed file contains shifted coordinates following an experimental modification to HLA-DRB1 that changed the reference genome length
hla_genes_regions_file = os.path.join(_data_dir, "mosdepth/hla_genes.bed")
#hla_genes_regions_file = os.path.join(_data_dir, "mosdepth/hla_genes.shifted.bed")

# BED file for parsing haploblocks in the extended MHC region
# Used in evaluate_gene_haploblocks() function of investigate_haploblocks_methods.py
genes_bed = os.path.join(_data_dir, "reference/parse_haploblocks_bed.bed")
#genes_bed = os.path.join(_data_dir, "reference/parse_haploblocks_bed.shifted.bed")
genes_of_interest_extended = ("HLA-A", "HLA-B", "HLA-C", "HLA-DRB1", "HLA-DQA1", "HLA-DQA2", "HLA-DQB1", "HLA-DQB2", "HLA-DPA1", "HLA-DPB1")

# Paramters

# Minimum reads per sample
# DeepVariant stalls when a sample has very few BAM records (e.g., HG01891 had only 35 mapped reads to chr6)
# Set threshold at which variant calling should not proceed
# chr6_read_count is returned by filter_reads() function
# Variant calling requires that chr6_read_count >= min_reads_sample:
min_reads_sample = 100

# This program is for long-read data only. Require that mean read length is at least 300 bp or higher
# is returned by run_fastplong(), called within parse_input_file() of sample_manager.py
# The program exits with a ValueError if mean_read_length falls below the specified threshold.
min_read_length = 300

# Coverage depth thresholds to proceed with HLA typing
# Used in parse_mosdepth() function of preprocess_methods.py
# Mean gene-wide depth
depth_thresh = 10
# Proportion of bases with mean depth >= 20x
prop_20x_thresh = 0.0
# Proportion of bases with mean depth >= 30x
prop_30x_thresh = 0.0
# Extended MHC coordinates
mhc_start = 29555628
mhc_stop = 33409896

# DNA bases and stop codons
# Used in parse_fastas() function of reconstruct_fasta_methods.py for vcf2fasta sanity checking
DNA_bases = {"A", "T", "G", "C"}
stop_codons = ["TAA", "TAG", "TGA"]

# IPD/IMGT HLA XML file for HLA allele classification
# Downloaded from https://github.com/ANHIG/IMGTHLA
# Used in hla_typer.py for HLA typing
IMGT_XML = os.path.join(_data_dir, "IPD_IMGT_XML/hla.xml")

# GRCh38 HLA gene antigen recognition sequence coordinates (1-based coordinates, GFF format)
# Used in evaluate_gene_haploblocks() function of investigate_haploblocks_methods.py 
# For HLA-C/B/C, the ARS is CDS 2 and 3
# For HLA-DPA1/DPB1/DQA1/DQB1/DRB1, the ARS is CDS 2
ARS_dict = {
	"HLA-A": "chr6:29942757-29943543",
	"HLA-B": "chr6:31356167-31356957",
	"HLA-C": "chr6:31271073-31271868",
	"HLA-DPA1": "chr6:33069641-33069886",
	"HLA-DPB1": "chr6:33080672-33080935",
	"HLA-DQA1": "chr6:32641310-32641558",
	"HLA-DQB1": "chr6:32664798-32665067",
	"HLA-DRB1": "chr6:32584109-32584378"
}

# Gene coordinates (1-based coordinates, GFF format)
# Used in parse_fastas() function of reconstruct_fasta_methods.py to clamp haploblock coordinates to gene coordinates
# HLA-B and HLA-DQA1 coordinates were slightly modified from the raw GFF3 file to exclude exons that are not part of the MANE Select transcript
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

# CDS coordinates for each HLA gene of interest for HLA typing (1-based coordinates, GFF format)
# Used in parse_fastas() function of reconstruct_fasta_methods.py to determine whether each CDS is fully contained in the haploblock
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