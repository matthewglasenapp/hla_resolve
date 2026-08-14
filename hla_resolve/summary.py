# This software is Copyright ©2026. The Regents of the University of California
# ("Regents"). All Rights Reserved.
#
# See LICENSE.txt for license details.

import json
import os
from datetime import datetime, timezone

# Pass-3 result tuple, as built in hla_typer.pass_3_classification:
#   (allele, distance, match_length, identity, mismatch_identity, tiebreak_used, equidistant)
ALLELE, DISTANCE, TIEBREAK, EQUIDISTANT = 0, 1, 5, -1


def gene_and_index(record_name):
	# HG002_HLA-A_1 -> ("HLA-A", "1"). Returns (None, None) if it does not parse.
	name = record_name[:-len("_incomplete")] if record_name.endswith("_incomplete") else record_name
	head, sep, index = name.rpartition("_")
	if not sep or index not in ("1", "2"):
		return None, None
	hla = head[head.find("HLA-"):] if "HLA-" in head else None
	return hla, index


def reconstruction_of(gene, cds_rescued_genes):
	# DR/DQ re-consensus is not represented here. Whether its accept gate fired
	# is decided inside reconsensus_drdq and is not returned to the pipeline.
	if gene in cds_rescued_genes:
		return cds_rescued_genes[gene].get("tier", "cds_rescue")
	return "full_gene"


def build(config, classifications, coverage_stats, phased_genes, cds_rescued_genes,
		  reads, runtime_seconds, status, version):
	genes = {}
	for gene in config['genes_of_interest']:
		genes[gene] = {
			"calls": [None, None],
			"distance": [None, None],
			"equidistant": [0, 0],
			"phasing": "fully_phased" if gene in phased_genes else "partial",
			"reconstruction": "none",
			"mean_depth": None,
			"prop_20x": None,
		}
		stats = (coverage_stats or {}).get(gene)
		if stats:
			genes[gene]["mean_depth"] = round(stats["depth"], 1)
			genes[gene]["prop_20x"] = round(stats["prop_20x"], 4)

	for record_name, result in (classifications or {}).items():
		gene, index = gene_and_index(record_name)
		if gene is None or gene not in genes:
			continue
		slot = int(index) - 1
		genes[gene]["calls"][slot] = result[ALLELE]
		genes[gene]["distance"][slot] = result[DISTANCE]
		equidistant = result[EQUIDISTANT]
		genes[gene]["equidistant"][slot] = len(equidistant) if equidistant else 1
		genes[gene]["reconstruction"] = reconstruction_of(gene, cds_rescued_genes)

	alleles_called = sum(1 for g in genes.values() for c in g["calls"] if c)
	genes_typed = sum(1 for g in genes.values() if any(g["calls"]))

	return {
		"sample": config['sample_ID'],
		"hla_resolve_version": version,
		"completed": datetime.now(timezone.utc).strftime("%Y-%m-%dT%H:%M:%SZ"),
		"status": status,
		"platform": config['platform'].lower(),
		"scheme": config.get('scheme'),
		"input_file": config.get('input_file'),
		"reads": reads,
		"runtime_seconds": runtime_seconds,
		"genes_typed": genes_typed,
		"genes_total": len(genes),
		"alleles_called": alleles_called,
		"genes": genes,
	}


def write(path, payload):
	with open(path, "w") as f:
		json.dump(payload, f, indent=2)
		f.write("\n")
	return path
