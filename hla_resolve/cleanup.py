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

# Every directory hla_resolve creates under the sample directory, minus the
# typing results, which --clean_up keeps. A recursive delete should only ever
# reach what this tool wrote, so anything else found there is reported and left
# alone. A directory added in a future version has to be added here too; the
# warning makes that visible rather than silently deleting a user's data.
ARTIFACT_DIR_KEYS = (
	'fastq_raw_dir', 'fastq_trimmed_dir', 'mapped_bam_dir', 'parsed_haploblock_dir',
	'genotypes_dir', 'sv_dir', 'filtered_vcf_dir', 'vcf2fasta_out_dir',
	'hla_fasta_dir', 'mosdepth_dir', 'phased_vcf_dir', 'pbtrgt_dir',
)

_KEEP_ALL = False
_PROTECTED = set()       # realpaths
_PROTECTED_IDS = set()   # (st_dev, st_ino), so a hard link cannot slip past
_FREED_BYTES = 0


def set_policy(keep_all_intermediates=False, protected_paths=()):
	"""Set the retention policy for the run. Called once, before any stage runs.

	protected_paths names files the tool must never remove whatever the policy,
	which in practice means the user's own input file. Several stages read the
	input in place instead of copying it, so an unguarded discard of "the reads
	this stage consumed" could otherwise delete it.
	"""
	global _KEEP_ALL, _PROTECTED, _PROTECTED_IDS, _FREED_BYTES
	_KEEP_ALL = bool(keep_all_intermediates)
	_PROTECTED = {os.path.realpath(path) for path in protected_paths if path}

	# Match on the file itself, not on how it was spelled. A hard link to the
	# input has a different path and the same inode, so a realpath comparison
	# alone would not recognize it.
	_PROTECTED_IDS = set()
	for path in _PROTECTED:
		try:
			info = os.stat(path)
		except OSError:
			continue
		_PROTECTED_IDS.add((info.st_dev, info.st_ino))

	_FREED_BYTES = 0


def _human(size):
	if size >= 1e9:
		return f"{size / 1e9:.1f} GB"
	if size >= 1e6:
		return f"{size / 1e6:.0f} MB"
	return f"{size / 1e3:.0f} KB"


def keep_all():
	"""True when the run was asked to retain every intermediate.

	Lets a stage skip producing a debug-only artifact rather than writing it and
	deleting it again.
	"""
	return _KEEP_ALL


def is_protected(path):
	"""True for a file the tool must never remove, however it was reached."""
	if not path:
		return False
	if os.path.realpath(path) in _PROTECTED:
		return True
	try:
		info = os.stat(path)
	except OSError:
		return False
	return (info.st_dev, info.st_ino) in _PROTECTED_IDS


def holds_protected(directory):
	"""True when a protected file sits anywhere inside this directory.

	Guards the recursive delete, which would otherwise reach a protected file by
	removing the tree above it rather than by naming it.
	"""
	directory = os.path.realpath(directory)
	for path in _PROTECTED:
		try:
			if os.path.commonpath([path, directory]) == directory:
				return True
		except ValueError:
			continue
	return False


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
	ours = set()
	for key in ARTIFACT_DIR_KEYS:
		directory = config.get(key)
		if directory:
			ours.add(os.path.realpath(directory))

	unknown = []
	for entry in os.scandir(config['output_dir']):
		real = os.path.realpath(entry.path)
		if real == keep:
			continue

		# Only remove what this tool wrote.
		if real not in ours:
			unknown.append(entry.path)
			continue

		# Never take a tree down over the top of an input file.
		if holds_protected(entry.path):
			announce(f"WARNING: keeping {entry.path}, it holds an input file")
			continue

		if entry.is_dir(follow_symlinks=False):
			shutil.rmtree(entry.path, ignore_errors=True)
		else:
			_remove(entry.path)

	for path in unknown:
		announce(f"WARNING: keeping {path}, hla_resolve did not write it")

	print(f"HLA typing results retained in {config['hla_typing_dir']}")
	print()
