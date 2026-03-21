#!/usr/bin/env python3
# This software is Copyright ©2026. The Regents of the University of California
# ("Regents"). All Rights Reserved.
#
# See LICENSE.txt for license details.

import os
import re

gff_dir = os.path.join(os.path.dirname(__file__), "..", "hla_resolve", "data", "hla_gff")

# HLA-C gene model requires 1kb padding on each end because most IPD-IMGT/HLA
# reference sequences include substantially more UTR than the GRCh38 annotation.
GENE_PADDING = {
	"HLA-C": 1000,
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
		pad = GENE_PADDING[gene_name]
		gene_fields[3] = str(int(gene_fields[3]) - pad)
		gene_fields[4] = str(int(gene_fields[4]) + pad)

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
