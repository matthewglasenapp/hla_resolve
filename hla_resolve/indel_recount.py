# This software is Copyright ©2026. The Regents of the University of California
# ("Regents"). All Rights Reserved.
#
# See LICENSE.txt for license details.

"""Restore the second allele of a heterozygous indel that Clair3 calls homozygous.

Where the two haplotypes carry different insertions or deletions at the same
position (a 1/2 site), Clair3 on ONT reads can report only the more common one
as 1/1, with an allele fraction near 0.5. One haplotype then gets the wrong
sequence. Each such call is recounted from the reads, and when a second indel
allele has enough support the record is rewritten as 1/2.
"""

from collections import Counter

import pysam

MAX_AF = 0.75        # a 1/1 call above this is taken as a true homozygote
MIN_FRACTION = 0.25  # each of the two alleles needs this share of the spanning reads
MIN_READS = 5        # and at least this many reads


def _indel_alleles(bam, chrom, pos0, ref_base):
	"""Count the indel allele each read carries right after 0-based pos0.

	An insertion is written as ref_base + inserted bases and a deletion as
	ref_base + deleted bases with the ALT ref_base, which is how both appear in
	a VCF anchored at pos0. Reads with no indel there count as the reference.
	"""
	counts = Counter()
	for col in bam.pileup(chrom, pos0, pos0 + 1, truncate=True, stepper="nofilter",
	                      ignore_overlaps=False, min_base_quality=0):
		for pr in col.pileups:
			if pr.is_refskip:
				continue
			if pr.is_del:
				continue
			if pr.indel > 0:
				q = pr.alignment.query_sequence
				ins = q[pr.query_position + 1: pr.query_position + 1 + pr.indel]
				counts[("ins", ref_base + ins.upper())] += 1
			elif pr.indel < 0:
				counts[("del", -pr.indel)] += 1
			else:
				counts["ref"] += 1
	return counts


def split_collapsed_indels(input_vcf, output_vcf, bam_path):
	"""Rewrite 1/1 indels with a well-supported second allele as 1/2.

	Returns the number of records rewritten.
	"""
	vcf = pysam.VariantFile(input_vcf)
	out = pysam.VariantFile(output_vcf, "wz", header=vcf.header)
	bam = pysam.AlignmentFile(bam_path)
	sample = list(vcf.header.samples)[0]
	rewritten = 0

	for rec in vcf:
		s = rec.samples[sample]
		gt = s.get("GT")
		af = s.get("AF") if "AF" in rec.format else None
		if isinstance(af, tuple):
			af = af[0]
		if (gt != (1, 1) or rec.alts is None or len(rec.alts) != 1 or af is None
				or af >= MAX_AF or len(rec.ref) == len(rec.alts[0])):
			out.write(rec)
			continue

		# Only simple anchored records: an insertion with a 1 bp REF or a
		# deletion with a 1 bp ALT, so the read counts line up with the call.
		alt = rec.alts[0]
		if len(rec.ref) == 1:
			called = ("ins", alt)
		elif len(alt) == 1:
			called = ("del", len(rec.ref) - 1)
		else:
			out.write(rec)
			continue
		counts = _indel_alleles(bam, rec.chrom, rec.pos - 1, rec.ref[0])
		total = sum(counts.values())
		others = [(k, n) for k, n in counts.most_common() if k != "ref" and k != called]
		if not others or total == 0:
			out.write(rec)
			continue
		key, n = others[0]
		c = counts[called]
		# Both alleles need read support. A call the reads do not show at this
		# position is a different spelling of a nearby event, not a collapsed
		# 1/2 site.
		if (n < MIN_READS or n / total < MIN_FRACTION
				or c < MIN_READS or c / total < MIN_FRACTION):
			out.write(rec)
			continue

		# Spell the second allele against this record's REF. A deletion longer
		# than REF needs a longer REF, which the 1/2 record cannot share, so it
		# is left alone.
		if key[0] == "ins":
			second = key[1] + rec.ref[1:]
		else:
			if key[1] >= len(rec.ref):
				out.write(rec)
				continue
			second = rec.ref[0] + rec.ref[1 + key[1]:]
		if second == alt or second == rec.ref:
			out.write(rec)
			continue

		rec.alts = (alt, second)
		rec.samples[sample]["GT"] = (1, 2)
		out.write(rec)
		rewritten += 1
		print(f"  {rec.chrom}:{rec.pos} {rec.ref}>{alt} 1/1 (AF {af:.2f}) -> {alt},{second} 1/2 "
		      f"({counts[called]} and {n} of {total} reads)")

	out.close()
	vcf.close()
	bam.close()
	pysam.tabix_index(output_vcf, preset="vcf", force=True)
	return rewritten
