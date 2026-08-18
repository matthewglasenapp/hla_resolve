# This software is Copyright ©2026. The Regents of the University of California
# ("Regents"). All Rights Reserved.
#
# See LICENSE.txt for license details.

import glob
import os
import shutil

def remove_stale_sort_temps(config):
	"""Remove samtools sort spill files left behind by an interrupted run."""
	# A sort that finishes removes its own temps, and this runs before alignment
	# starts, so anything matching here is from an earlier job that was killed.
	stale = []
	for bam in (config['hg38_bam'], config['hg38_bam_drb']):
		stale.extend(glob.glob(bam + ".tmp.*.bam"))

	freed = 0
	removed = 0
	for path in stale:
		try:
			size = os.path.getsize(path)
			os.remove(path)
		except OSError:
			continue
		freed += size
		removed += 1

	if removed:
		print(f"Removed {removed} sort temp file(s) from an interrupted run, freed {freed / 1e9:.1f} GB")
		print()

def discard_full_genome_bam(config):
	"""Remove the whole-genome BAM once reads are filtered to chr6."""
	if config['scheme'] not in ("WGS", "WES") or config.get('keep_full_bam', False):
		return

	bam = config['hg38_bam']
	freed = 0
	for path in (bam, bam + ".bai", bam + ".csi"):
		if os.path.exists(path):
			freed += os.path.getsize(path)
			os.remove(path)

	if freed:
		print(f"Discarded whole-genome BAM, freed {freed / 1e9:.1f} GB: {bam}")
		print("Pass --keep_full_bam to retain it.")
		print()

def cleanup_intermediate_files(config):
	"""
	Clean up intermediate files and directories after HLA typing is complete.
	This is a separate concern from allele resolution.
	"""
	if not config.get('clean_up', False):
		return
	
	print("Cleaning up intermediate files...")
	
	# Define directories to clean up
	directories_to_clean = [
		config['fastq_raw_dir'],
		config['fastq_trimmed_dir'],
		config['mapped_bam_dir'],
		config['genotypes_dir'],
		config['sv_dir'],
		config['phased_vcf_dir'],
		config['mosdepth_dir'],
		config['parsed_haploblock_dir'],
		config['filtered_vcf_dir'],
		config['vcf2fasta_out_dir']
	]
	
	for directory in directories_to_clean:
		if os.path.exists(directory) and directory != config['hla_typing_dir']:
			print(f"Removing: {directory}")
			shutil.rmtree(directory)
	
	print("Cleanup completed")
