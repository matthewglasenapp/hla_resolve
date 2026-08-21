# This software is Copyright ©2026. The Regents of the University of California
# ("Regents"). All Rights Reserved.
#
# See LICENSE.txt for license details.

"""File retention policy.

A run keeps what a user needs to interpret or re-analyze the calls. That is the
haplotagged MHC BAM, the chr6 and per-gene VCFs, the coverage and phasing
reports, the reconstructed FASTA haplotypes, and the typing results. Everything
else is a working copy of the reads that some later stage supersedes.

Those working copies are removed as soon as their last consumer finishes rather
than at the end of the run, so peak disk stays near the size of the input
instead of several times it. A cohort of a thousand samples is the case that
matters here. A discarded file makes the sample non-resumable, which is the
intended trade.

--keep_all_intermediates turns every discard below into a no-op.
"""

import glob
import os
import shutil

from .utils import announce, detail

# Index and companion files travel with the file they describe, so a caller
# names the BAM or VCF and these go with it.
SIDECAR_SUFFIXES = (".bai", ".csi", ".crai", ".tbi", ".pbi", ".fai")

_KEEP_ALL = False
_PROTECTED = set()
_FREED_BYTES = 0


def set_policy(keep_all_intermediates=False, protected_paths=()):
	"""Set the retention policy for the run. Called once, before any stage runs.

	protected_paths names files the tool must never remove whatever the policy,
	which in practice means the user's own input file. Several stages read the
	input in place instead of copying it, so an unguarded discard of "the reads
	this stage consumed" could otherwise delete it.
	"""
	global _KEEP_ALL, _PROTECTED, _FREED_BYTES
	_KEEP_ALL = bool(keep_all_intermediates)
	_PROTECTED = {os.path.realpath(path) for path in protected_paths if path}
	_FREED_BYTES = 0


def _human(size):
	if size >= 1e9:
		return f"{size / 1e9:.1f} GB"
	if size >= 1e6:
		return f"{size / 1e6:.0f} MB"
	return f"{size / 1e3:.0f} KB"


def is_protected(path):
	return bool(path) and os.path.realpath(path) in _PROTECTED


def _remove(path):
	"""Remove one file and report its size. Returns 0 if it was not there."""
	try:
		size = os.path.getsize(path)
		os.remove(path)
	except OSError:
		return 0
	return size


def discard(paths, label):
	"""Remove intermediate files whose last consumer has finished.

	Sidecar index files are removed alongside each named path. Missing files are
	not an error: a stage that did not run leaves nothing to clean up.
	"""
	global _FREED_BYTES

	if _KEEP_ALL:
		return 0

	freed = 0
	removed = []
	for path in paths:
		if not path or is_protected(path):
			continue
		for candidate in (path,) + tuple(path + suffix for suffix in SIDECAR_SUFFIXES):
			size = _remove(candidate)
			if size:
				freed += size
				removed.append(candidate)

	if removed:
		detail(f"Discarded {label}, freed {_human(freed)}:")
		for path in removed:
			detail(f"  {path}")

	_FREED_BYTES += freed
	return freed


def discard_temp(*paths):
	"""Remove a tool's own scratch files. Same policy, no log line.

	Used for files the pipeline never names in its output, such as the SA-stripped
	pbsv input copy or the uncompressed VCF that only exists to be bgzipped.
	"""
	global _FREED_BYTES

	if _KEEP_ALL:
		return 0

	freed = 0
	for path in paths:
		if not path or is_protected(path):
			continue
		for candidate in (path,) + tuple(path + suffix for suffix in SIDECAR_SUFFIXES):
			freed += _remove(candidate)

	_FREED_BYTES += freed
	return freed


def prune_empty_dirs(config):
	"""Remove per-sample directories left empty once their contents are discarded."""
	if _KEEP_ALL:
		return

	for key in ('fastq_raw_dir', 'fastq_trimmed_dir', 'mapped_bam_dir', 'sv_dir',
	            'pbtrgt_dir', 'genotypes_dir', 'vcf2fasta_out_dir', 'mosdepth_dir',
	            'parsed_haploblock_dir', 'filtered_vcf_dir', 'phased_vcf_dir'):
		directory = config.get(key)
		if not directory:
			continue
		try:
			os.rmdir(directory)
		except OSError:
			continue


def report_discarded():
	"""One line at the end of the run with the total reclaimed."""
	if _KEEP_ALL or _FREED_BYTES < 1e8:
		return
	announce(f"Discarded intermediate files, freed {_human(_FREED_BYTES)}. "
	         "Pass --keep_all_intermediates to retain them.")


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
		size = _remove(path)
		if not size:
			continue
		freed += size
		removed += 1

	if removed:
		print(f"Removed {removed} sort temp file(s) from an interrupted run, freed {_human(freed)}")
		print()


def discard_full_genome_bam(config):
	"""Remove the reference-genome BAM once reads are filtered to the MHC.

	Every stage after read filtering works on the MHC window, so the full
	alignment is dead weight. It is the largest single file a run produces.
	"""
	if config.get('keep_full_bam', False) or _KEEP_ALL:
		return

	bam = config['hg38_bam']
	freed = discard_temp(bam)

	if freed:
		announce(f"Discarded the reference-genome BAM, freed {_human(freed)}")
		detail(f"  {bam}")
		detail("Pass --keep_full_bam to retain it.")


def cleanup_intermediate_files(config):
	"""--clean_up: keep the HLA typing results and nothing else.

	For cohorts where only the calls are wanted. The run log lives beside the
	sample directory rather than inside it, so it survives this.
	"""
	if not config.get('clean_up', False):
		return

	print("Removing everything except the HLA typing results...")

	keep = os.path.realpath(config['hla_typing_dir'])

	for entry in os.scandir(config['output_dir']):
		if os.path.realpath(entry.path) == keep:
			continue
		if entry.is_dir(follow_symlinks=False):
			shutil.rmtree(entry.path, ignore_errors=True)
		else:
			_remove(entry.path)

	print(f"HLA typing results retained in {config['hla_typing_dir']}")
	print()
