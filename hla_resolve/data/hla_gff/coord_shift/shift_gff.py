#!/usr/bin/env python3
import os
import shutil

SHIFT = 4959
SKIP_GENES = {"hla_a", "hla_b", "hla_c"}  # skip these gene families

indir = "."
outdir = "tmp_shifted"

os.makedirs(outdir, exist_ok=True)

for fname in os.listdir(indir):
    if not fname.endswith(".gff3"):
        continue

    # skip hla-a/b/c related files
    if any(fname.lower().startswith(g) for g in SKIP_GENES):
        print(f"Skipping {fname}")
        continue

    infile = os.path.join(indir, fname)
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

# overwrite originals
for fname in os.listdir(outdir):
    shutil.move(os.path.join(outdir, fname), os.path.join(indir, fname))
os.rmdir(outdir)
print("All done. Shift applied.")
