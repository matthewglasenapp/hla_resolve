# This software is Copyright ©2026. The Regents of the University of California
# ("Regents"). All Rights Reserved.
#
# See LICENSE.txt for license details.

"""Recover HLA-A reads that long whole-genome reads place at a neighboring locus."""

import os
import subprocess

import edlib
import pysam

from . import config
from .cleanup import discard_temp
from .utils import detail

FLANK = 500
MAX_ERROR = 0.08
WINDOW_PAD = 200000


def _region(region):
	chrom, span = region.split(":")
	start, stop = span.split("-")
	return chrom, int(start), int(stop)


def _read_fasta(path):
	seqs, name = {}, None
	with open(path) as fh:
		for line in fh:
			line = line.strip()
			if line.startswith(">"):
				name = line[1:].split()[0]
				seqs[name] = []
			elif name:
				seqs[name].append(line.upper())
	return {k: "".join(v) for k, v in seqs.items()}


def _revcomp(seq):
	return seq[::-1].translate(str.maketrans("ACGTN", "TGCAN"))


def _error(probe, seq):
	"""Edit distance of probe placed inside seq, either strand, per probe base."""
	return min(edlib.align(p, seq, mode="HW")["editDistance"] for p in (probe, _revcomp(probe))) / len(probe)


def rescue_hla_a_copies(mhc_bam, reference_fasta, panel_fasta, threads):
	"""Add recovered HLA-A reads to mhc_bam in place. Returns the count."""
	chrom, h_start, h_stop = _region(config.hla_h_region)
	_, _, a_stop = _region(config.hla_a_region)
	alleles = {n: s for n, s in _read_fasta(panel_fasta).items() if n.startswith("A*")}
	with pysam.FastaFile(reference_fasta) as fa:
		hla_h = fa.fetch(chrom, h_start - 1, h_stop).upper()
		win_start = max(0, h_start - 1 - WINDOW_PAD)
		window = fa.fetch(chrom, win_start, a_stop + WINDOW_PAD).upper()

	copies = []
	with pysam.AlignmentFile(mhc_bam) as bam:
		header = bam.header
		done = {r.query_name[:-5] for r in bam.fetch(chrom, h_start - 1, a_stop) if r.query_name.endswith("_hlaA")}
		for read in bam.fetch(chrom, h_start - 1, h_stop):
			if read.is_secondary or read.is_supplementary or read.is_unmapped:
				continue
			if read.query_name.endswith("_hlaA") or read.query_name in done:
				continue
			if read.reference_start > h_start - 1 or read.reference_end < h_stop:
				continue
			qs = qe = None
			for q, p in read.get_aligned_pairs(matches_only=True):
				if p >= h_start - 1 and qs is None:
					qs = q
				if p <= h_stop - 1:
					qe = q
			if qs is None or qe is None:
				continue
			lo, hi = max(0, qs - FLANK), qe + 1 + FLANK
			segment = read.query_sequence[lo:hi]
			h_err = _error(hla_h, segment)
			a_err, allele = min((_error(s, segment), n) for n, s in alleles.items())
			if a_err < h_err and a_err <= MAX_ERROR:
				quals = read.query_qualities
				qual = pysam.qualities_to_qualitystring(quals[lo:hi]) if quals is not None else "?" * len(segment)
				rg = read.get_tag("RG") if read.has_tag("RG") else None
				copies.append((read.query_name + "_hlaA", segment, qual, rg))
				detail(f"{read.query_name}: {allele} ({a_err:.3f} against {h_err:.3f})")
	if not copies:
		return 0

	work = mhc_bam + ".hla_a_rescue"
	window_fa, reads_fq, sam = work + ".window.fa", work + ".fq", work + ".sam"
	added, merged = work + ".added.bam", work + ".merged.bam"
	try:
		with open(window_fa, "w") as fh:
			fh.write(">window\n")
			for i in range(0, len(window), 80):
				fh.write(window[i:i + 80] + "\n")
		with open(reads_fq, "w") as fh:
			for name, seq, qual, _ in copies:
				fh.write(f"@{name}\n{seq}\n+\n{qual}\n")
		subprocess.run(f"set -o pipefail; {config.rammap} -Y -t {threads} -ax map-ont {window_fa} {reads_fq} > {sam}",
		               shell=True, check=True, executable="/bin/bash",
		               stdout=subprocess.DEVNULL, stderr=subprocess.DEVNULL)

		rgs = {name: rg for name, _, _, rg in copies}
		tid = header.get_tid(chrom)
		records = []
		with pysam.AlignmentFile(sam) as src:
			for rec in src:
				if rec.is_secondary or rec.is_supplementary or rec.is_unmapped:
					continue
				out = pysam.AlignedSegment(header)
				out.query_name = rec.query_name
				out.flag = rec.flag & 0x10
				out.reference_id = tid
				out.reference_start = rec.reference_start + win_start
				out.mapping_quality = rec.mapping_quality
				out.cigartuples = rec.cigartuples
				out.query_sequence = rec.query_sequence
				out.query_qualities = rec.query_qualities
				if rgs.get(rec.query_name):
					out.set_tag("RG", rgs[rec.query_name])
				records.append(out)
				detail(f"{rec.query_name} realigned to {chrom}:{out.reference_start + 1} MAPQ {out.mapping_quality}")
		if not records:
			return 0
		records.sort(key=lambda r: r.reference_start)
		with pysam.AlignmentFile(added, "wb", header=header) as dst:
			for rec in records:
				dst.write(rec)
		subprocess.run(f"samtools merge -f -c -p {merged} {mhc_bam} {added}", shell=True, check=True,
		               stdout=subprocess.DEVNULL, stderr=subprocess.DEVNULL)
		os.replace(merged, mhc_bam)
		subprocess.run(f"samtools index {mhc_bam}", shell=True, check=True)
		return len(records)
	finally:
		discard_temp(window_fa, reads_fq, sam, added, merged)
