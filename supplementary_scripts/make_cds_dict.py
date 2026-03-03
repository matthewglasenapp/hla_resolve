#!/usr/bin/env python3
import os

indir = os.path.join(os.path.dirname(__file__), "..", "hla_resolve", "data", "hla_gff")
CDS_dict = {}

def parse_cds_coords(gff_path):
    cds_coords = []
    with open(gff_path) as f:
        for line in f:
            if line.startswith("#"):
                continue
            parts = line.strip().split("\t")
            if len(parts) < 9:
                continue
            ftype = parts[2]
            if ftype == "CDS":
                start = int(parts[3])  # GFF3 start = 1-based inclusive
                end = int(parts[4])    # GFF3 end   = 1-based inclusive
                cds_coords.append([start, end])
    # make sure they're in genomic order
    cds_coords.sort(key=lambda x: x[0])
    return cds_coords

for fname in os.listdir(indir):
    if fname.endswith("_cds_sorted.gff3"):
        gene_name = fname.split("_cds_sorted.gff3")[0].upper().replace("HLA_", "HLA-")
        fpath = os.path.join(indir, fname)
        CDS_dict[gene_name] = parse_cds_coords(fpath)

# Print nicely
print(CDS_dict)
