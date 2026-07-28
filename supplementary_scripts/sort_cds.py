#!/usr/bin/env python3
# This software is Copyright ©2026. The Regents of the University of California
# ("Regents"). All Rights Reserved.
#
# See LICENSE.txt for license details.

import os
import re

gff_dir = os.path.join(os.path.dirname(__file__), "..", "hla_resolve", "data", "hla_gff")

# Some GRCh38 gene models are shorter than the IPD-IMGT/HLA reference sequences
# they are matched against, so the reconstruction runs out before the reference
# allele does and the trailing bases have nothing to align to. These genes get
# their gene interval padded.
#
# Values are (low_pad, high_pad) in reference bases: low_pad is added below the
# start coordinate, high_pad above the end. Padding is asymmetric because the
# shortfall usually is. For a minus-strand gene the 3' end of the transcript is
# the LOW coordinate, so a 3' shortfall is corrected with low_pad.
#
#   HLA-C    symmetric 1 kb; IPD sequences carry substantially more UTR at both
#            ends than the Ensembl annotation.
#   HLA-B    minus strand, systematic ~510 bp gap at the 3' (low) end across 45
#            haplotypes; padded 602 to clear it.
#   HLA-DPA1 minus strand, systematic ~250 bp gap at the 3' (low) end across 28
#            haplotypes; padded 299 to clear it.
#
# HLA-A, HLA-DPB1 and HLA-DRB1 show no systematic trailing gap and are unpadded.
GENE_PADDING = {
	"HLA-C":    (1000, 1000),
	"HLA-B":    (602, 0),
	"HLA-DPA1": (299, 0),
}

def get_raw_gff_files():
	"""Return only raw .gff3 files (exclude *_cds_sorted.gff3 and *_gene.gff3)."""
	gff_files = []
	for fname in os.listdir(gff_dir):
		if fname.endswith(".gff3") and not re.search(r"_(cds_sorted|gene)\.gff3$", fname):
			gff_files.append(os.path.join(gff_dir, fname))
	return gff_files

# Run only once
def sort_cds(gff_file):
	meta_lines = []
	cds_lines = []
	gene_line = []
	strand = None

	with open(gff_file, "r") as f:
		for line in f:
			if line.startswith("#"):
				meta_lines.append(line)
			else:
				fields = line.strip().split("\t")
				if fields[0] == "6":
					fields[0] = "chr6"
				if fields[2] == "CDS":
					start = int(fields[3])
					strand = fields[6]
					cds_lines.append((start, fields))
				elif fields[2] == "gene":
					gene_line.append(fields)


	sorted_cds = sorted(cds_lines, key=lambda line: line[0])

	outfile_cds = gff_file.replace(".gff3", "_cds_sorted.gff3")
	outfile_gene = gff_file.replace(".gff3", "_gene.gff3")

	with open(outfile_cds, "w") as out:
		for line in meta_lines:
			out.write(line)
		for line in sorted_cds:
			fields = line[1]
			out.write("\t".join(fields) + "\n")

	# Apply padding to gene coordinates if configured
	gene_fields = gene_line[0]
	attrs = gene_fields[8] if len(gene_fields) > 8 else ""
	name_match = re.search(r"Name=([^;]+)", attrs)
	gene_name = name_match.group(1) if name_match else None
	if gene_name in GENE_PADDING:
		low_pad, high_pad = GENE_PADDING[gene_name]
		gene_fields[3] = str(max(1, int(gene_fields[3]) - low_pad))
		gene_fields[4] = str(int(gene_fields[4]) + high_pad)

	with open(outfile_gene, "w") as out2:
		for line in meta_lines:
			out2.write(line)
		out2.write("\t".join(gene_fields) + "\n")

	print(f"Wrote: {outfile_gene}")
	print(f"Wrote: {outfile_cds}")

if __name__ == "__main__":
	gff_files = get_raw_gff_files()
	for gff_file in gff_files:
		sort_cds(gff_file)
