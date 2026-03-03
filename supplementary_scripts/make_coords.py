#!/usr/bin/env python3
import os

indir = os.path.join(os.path.dirname(__file__), "..", "hla_resolve", "data", "hla_gff")

def parse_gff3(filename, feature_type):
    coords = []
    strand = None
    with open(filename) as f:
        for line in f:
            if line.startswith("#"):
                continue
            parts = line.strip().split("\t")
            if len(parts) < 9:
                continue
            ftype, start, end = parts[2], int(parts[3]), int(parts[4])
            if ftype == feature_type:
                coords.append((start, end))
                strand = parts[6]  # "+" or "-"
    return coords, strand

for fname in os.listdir(indir):
    if not (fname.endswith("_cds_sorted.gff3") or fname.endswith("_gene.gff3")):
        continue

    fpath = os.path.join(indir, fname)
    feature_type = "CDS" if "_cds_sorted.gff3" in fname else "gene"
    coords, strand = parse_gff3(fpath, feature_type)

    # Sort intervals according to transcript orientation
    if strand == "+":
        coords.sort(key=lambda x: x[0])
    elif strand == "-":
        coords.sort(key=lambda x: x[0], reverse=True)

    outname = fname.replace(".gff3", "_coords.txt")
    outpath = os.path.join(indir, outname)

    with open(outpath, "w") as out:
        for start, end in coords:
            if strand == "+":
                rng = range(start, end + 1)
            else:  # reverse strand
                rng = range(end, start - 1, -1)
            for pos in rng:
                out.write(f"{pos}\n")

    print(f"Wrote {feature_type}-expanded coords: {outpath} (strand {strand})")
