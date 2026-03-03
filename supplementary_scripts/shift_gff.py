#!/usr/bin/env python3
import os

SHIFT = 4959
SKIP_GENES = {"hla_a", "hla_b", "hla_c"}  # skip these gene families

gff_dir = os.path.join(os.path.dirname(__file__), "..", "hla_resolve", "data", "hla_gff")
outdir = os.path.join(gff_dir, "coord_shift")

os.makedirs(outdir, exist_ok=True)

for fname in os.listdir(gff_dir):
    if not fname.endswith(".gff3"):
        continue

    # skip hla-a/b/c related files
    if any(fname.lower().startswith(g) for g in SKIP_GENES):
        print(f"Skipping {fname}")
        continue

    infile = os.path.join(gff_dir, fname)
    outfile = os.path.join(outdir, fname)

    with open(infile) as fin, open(outfile, "w") as fout:
        for line in fin:
            if line.startswith("#") or not line.strip():
                fout.write(line)
                continue

            fields = line.rstrip("\n").split("\t")
            if len(fields) < 9:
                fout.write(line)
                continue

            try:
                start = int(fields[3]) + SHIFT
                end   = int(fields[4]) + SHIFT
                fields[3] = str(start)
                fields[4] = str(end)
                fout.write("\t".join(fields) + "\n")
            except ValueError:
                fout.write(line)

    print(f"Shifted {fname} -> {outfile}")

print(f"\nAll done. Shifted files written to {outdir}")
