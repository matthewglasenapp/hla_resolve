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


	if strand == "-":
		sorted_cds = sorted(cds_lines, key=lambda line: line[0], reverse=True)
	else: 
		sorted_cds = sorted(cds_lines, key=lambda line: line[0])

	outfile_cds = gff_file.replace(".gff3", "_cds_sorted.gff3")
	outfile_gene = gff_file.replace(".gff3", "_gene.gff3")

	with open(outfile_cds, "w") as out:
		for line in meta_lines:
			out.write(line)
		for line in sorted_cds:
			fields = line[1]
			out.write("\t".join(fields) + "\n")

	with open(outfile_gene, "w") as out2:
		for line in meta_lines:
			out2.write(line)
		out2.write("\t".join(gene_line[0]) + "\n")

	print(f"Wrote: {outfile_gene}")
	print(f"Wrote: {outfile_cds}")

sort_cds("/hb/groups/cornejo_lab/matt/hla_capture/input_data/hla_gff/hla_dqb1.gff3")
