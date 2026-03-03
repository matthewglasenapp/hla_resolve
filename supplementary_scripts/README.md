# Supplementary Scripts

One-time data preparation scripts used to generate the HLA GFF3-derived files in `hla_resolve/data/hla_gff/`.

## Pipeline

The scripts should be run in this order:

### 1. `sort_cds.py`

Takes the raw GFF3 files (e.g., `hla_a.gff3`) downloaded from GRCh38.p14 Ensembl and produces two derived files per gene:

- `*_cds_sorted.gff3` — CDS entries sorted by coordinate in transcript order. For minus-strand genes (e.g., HLA-B), CDS are reverse-sorted (high-to-low) to work around a [vcf2fasta bug](https://github.com/santiagosnchez/vcf2fasta/issues/23) where CDS are not concatenated in the correct order for minus-strand genes. Plus-strand genes are sorted ascending (no change needed).
- `*_gene.gff3` — Contains only the gene feature record, used for gene boundary extraction.

### 2. `make_coords.py`

Takes the `*_cds_sorted.gff3` and `*_gene.gff3` files produced by `sort_cds.py` and expands each interval into per-position coordinate files:

- `*_cds_sorted_coords.txt` — One genomic position per line for all CDS exons, ordered by transcript orientation.
- `*_gene_coords.txt` — One genomic position per line spanning the full gene interval.

These coordinate files are used in `reconstruct_fasta_methods.py` for coordinate clamping, ARS extraction, and mapping FASTA positions to genomic coordinates.

### 3. `make_cds_dict.py`

Parses the `*_cds_sorted.gff3` files to produce the `CDS_dict` dictionary defined in `hla_resolve/config.py` (line 335). This dictionary maps each HLA gene to its list of CDS interval coordinates (1-based, inclusive) in genomic order.

## Other scripts

### `shift_gff.py`

Not currently used. Shifts all GFF3 coordinates by +4959 bp for HLA class II genes, skipping class I genes (HLA-A, B, C). This was used to correct a coordinate offset between annotation sources for the class II region. Writes shifted files to `hla_gff/coord_shift/`.

### `shift_bed_after_drb1.py`

Not currently used. Shifts BED coordinates by +4959 bp for entries downstream of HLA-DRB1 (start >= 32,589,742) on chr6. This offset is needed when using a GRCh38 reference augmented with DRB1 ALT scaffolds. Reads `parse_haploblocks_bed.bed` from `glasenapp_2026_manuscript/data/bed_files/parse_haploblocks/` and writes `parse_haploblocks_bed.shifted.bed` to the same directory.
