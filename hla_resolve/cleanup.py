# This software is Copyright ©2026. The Regents of the University of California
# ("Regents"). All Rights Reserved.
#
# See LICENSE.txt for license details.

import os
import shutil

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
	
	print("Cleanup completed!")
