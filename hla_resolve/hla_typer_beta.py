"""
Beta version of hla_typer.py for metric analysis.
Outputs detailed CIGAR-derived metrics for pass 3 candidates to help
determine the optimal classification metric.
"""

from Bio import SeqIO
import xml.etree.ElementTree as ET
import pandas as pd
import edlib
import time
import math
import re
import csv
import argparse
import os

NUM_SAMPLES = None

# Function decorator to print time taken to run a function
def print_time_taken(fun):
    def wrapper(self, *args, **kwargs):
        start_time = time.time()
        result = fun(self, *args, **kwargs)
        end_time = time.time()
        time_taken = end_time-start_time
        seconds = time_taken % 60
        minutes = math.floor(time_taken / 60)
        print(f"INFO: {fun.__name__} took {minutes}m {seconds:2.0f}s")
        return result
    return wrapper

# Get Gene from sample name
def get_gene(sample_name, astrisk=True):
    hla_gene_index = sample_name[sample_name.find("HLA"):]
    hla_gene = hla_gene_index.split("_")[0]
    hla_gene = hla_gene[hla_gene.index("-")+1:]
    if astrisk: hla_gene = hla_gene + "*"
    return hla_gene

# Get sample ID from sample name
def get_sampleid(sample_name):
    return sample_name[:sample_name.find("HLA")-1]

def get_haplotype_index(sample_name):
    """Extract haplotype index (1 or 2) from sample name like 'IHW09021_HLA-A_1'"""
    return sample_name.split("_")[-1]

# =============================================================================
# ENHANCED DISTANCE FUNCTION - Returns detailed metrics
# =============================================================================
def get_detailed_metrics(sequence, sequence_data):
    """
    Get detailed alignment metrics from CIGAR string.
    Returns: dict with matches, mismatches, gap_events, total_gap_bp,
             edit_distance, gap_compressed_distance
    """
    # Ensure first argument is shorter sequence for semi-global alignment
    if len(sequence) <= len(sequence_data):
        results = edlib.align(sequence, sequence_data, task="path", mode="HW")
    else:
        results = edlib.align(sequence_data, sequence, task="path", mode="HW")

    cigar = results["cigar"]
    raw_edit_distance = results["editDistance"]

    # Parse CIGAR string for detailed metrics
    # = : match
    # X : mismatch
    # I : insertion in query
    # D : deletion in query (insertion in reference)

    matches = 0
    mismatches = 0
    insertions = []  # list of insertion lengths
    deletions = []   # list of deletion lengths

    # Parse all CIGAR operations
    for m in re.finditer(r"(\d+)([=XIDMNS])", cigar):
        length = int(m.group(1))
        op = m.group(2)

        if op == "=":
            matches += length
        elif op == "X":
            mismatches += length
        elif op == "I":
            insertions.append(length)
        elif op == "D":
            deletions.append(length)
        elif op == "M":
            # M can be match or mismatch in basic CIGAR
            # edlib with task="path" should use extended CIGAR (=, X)
            # but handle M just in case
            matches += length

    gap_events = len(insertions) + len(deletions)
    total_insertion_bp = sum(insertions)
    total_deletion_bp = sum(deletions)
    total_gap_bp = total_insertion_bp + total_deletion_bp

    # Gap-compressed distance: mismatches + gap_events (not gap bp)
    gap_compressed_distance = mismatches + gap_events

    # Alignment length (for reference)
    alignment_length = matches + mismatches + total_gap_bp

    return {
        "matches": matches,
        "mismatches": mismatches,
        "gap_events": gap_events,
        "total_gap_bp": total_gap_bp,
        "insertions": len(insertions),
        "insertion_bp": total_insertion_bp,
        "deletions": len(deletions),
        "deletion_bp": total_deletion_bp,
        "raw_edit_distance": raw_edit_distance,
        "gap_compressed_distance": gap_compressed_distance,
        "alignment_length": alignment_length,
        "cigar": cigar
    }

def compute_identity_metrics(metrics):
    """
    Compute various identity formulas from detailed metrics.
    Returns dict with different identity calculations.
    """
    m = metrics["matches"]
    mm = metrics["mismatches"]
    ge = metrics["gap_events"]
    gb = metrics["total_gap_bp"]
    gc_dist = metrics["gap_compressed_distance"]
    raw_dist = metrics["raw_edit_distance"]

    results = {
        "match_length": m,
        "edit_distance": raw_dist,
        "gap_compressed_distance": gc_dist,
    }

    # Original identity formula: 1 - (gap_compressed_dist / matches)
    if m > 0:
        results["identity_original"] = 1 - (gc_dist / m)
    else:
        results["identity_original"] = 0

    # New identity formula: matches / (matches + gap_compressed_dist)
    denom = m + gc_dist
    if denom > 0:
        results["identity_v2"] = m / denom
    else:
        results["identity_v2"] = 0

    # Mismatch-only identity: matches / (matches + mismatches)
    # Ignores gaps entirely - only scores where both have data
    denom = m + mm
    if denom > 0:
        results["identity_mismatch_only"] = m / denom
    else:
        results["identity_mismatch_only"] = 0

    # Weighted identity: mismatches count 10x more than gap events
    weighted_penalty = mm * 10 + ge
    denom = m + weighted_penalty
    if denom > 0:
        results["identity_weighted_10x"] = m / denom
    else:
        results["identity_weighted_10x"] = 0

    # Weighted identity: mismatches count 5x more than gap events
    weighted_penalty = mm * 5 + ge
    denom = m + weighted_penalty
    if denom > 0:
        results["identity_weighted_5x"] = m / denom
    else:
        results["identity_weighted_5x"] = 0

    return results

# =============================================================================
# STANDARD FUNCTIONS (mostly unchanged from original)
# =============================================================================

def get_distance(sequence, sequence_data, get_length=False, gap_compressed=True):
    """Original distance function for passes 1 and 2"""
    task = "distance"
    if gap_compressed or get_length:
        task = "path"

    if len(sequence) <= len(sequence_data):
        results = edlib.align(sequence, sequence_data, task=task, mode="HW")
    else:
        results = edlib.align(sequence_data, sequence, task=task, mode="HW")

    edit_distance = results["editDistance"]

    if gap_compressed:
        pattern = r"(\d+)(I|D)"
        for match in re.finditer(pattern, results["cigar"]):
            length = int(match.group(1))
            if length > 1:
                edit_distance -= length - 1

    if get_length:
        match_length = 0
        pattern = r"(\d+)(=)"
        for match in re.finditer(pattern, results["cigar"]):
            length = int(match.group(1))
            match_length += length
        return edit_distance, match_length

    return edit_distance

@print_time_taken
def build_g_group_dict(xml_file, ignore_unconfirmed=False, ignore_incomplete=False):
    """Build G-group dictionary from XML (unchanged)"""
    print("INFO: Parsing metadata XML file")
    tree = ET.parse(xml_file)
    root = tree.getroot()

    version = root.attrib['version']
    print("INFO: XML version:", version)
    group_accessor = "name" if version > "3.61.0" else "status"

    g_groups = dict()
    p_groups = dict()
    sequence_data = dict()

    tag = lambda x:'{http://hla.alleles.org/xml}' + x

    alleles = root.find(tag('alleles'))
    for allele in alleles.findall(tag('allele')):
        if ignore_unconfirmed and allele.find(tag('releaseversions')).attrib['confirmed'] == "Unconfirmed":
            continue

        class_type = allele.find(tag('locus')).attrib['class']

        sequence = allele.find(tag('sequence'))
        if sequence == None:
            continue

        dna_sequence = sequence.find(tag('nucsequence'))
        features = sequence.findall(tag('feature'))
        if dna_sequence == None or features == None:
            continue

        feature_orders = sorted([int(feature.attrib['order']) for feature in features if 'order' in feature.attrib.keys()])
        complete_order_list = list(range(1, len(feature_orders)+1))
        if ignore_incomplete and feature_orders != complete_order_list:
            continue

        name = allele.attrib["name"]

        sequence_data[name] = dict()
        for i, feature in enumerate(features):
            feature_type = feature.attrib["featuretype"]

            if feature_type not in ["UTR", "Exon", "Intron"]:
                continue

            if feature_type not in sequence_data[name].keys():
                sequence_data[name][feature_type] = []

            locations = feature.find(tag('SequenceCoordinates'))
            start, end = int(locations.attrib["start"]), int(locations.attrib["end"])
            cut_sequence = dna_sequence.text[start-1:end]

            sequence_data[name][feature_type].append((i, cut_sequence))

            if feature.attrib["name"] == "Exon 2":
                sequence_data[name]["peptide_binding_domain"] = [(2, cut_sequence)]

            if feature.attrib["name"] == "Exon 3" and class_type == "I":
                if sequence_data[name].get("peptide_binding_domain") == None:
                    sequence_data[name]["peptide_binding_domain"] = [(3, cut_sequence)]
                    continue
                sequence_data[name]["peptide_binding_domain"].append((3, cut_sequence))

        g_group_xml = allele.find(tag('hla_g_group'))
        if g_group_xml != None:
            g_group = g_group_xml.attrib[group_accessor]
            g_groups[name] = g_group

        p_group_xml = allele.find(tag('hla_p_group'))
        if p_group_xml != None:
            p_group = p_group_xml.attrib[group_accessor]
            p_groups[name] = p_group

    return g_groups, p_groups, sequence_data

def generate_allele_dict(allele_to_g_groups):
    g_group_to_allele = dict()
    for k, v in allele_to_g_groups.items():
        if v not in g_group_to_allele.keys():
            g_group_to_allele[v] = []
        g_group_to_allele[v].append(k)
    return g_group_to_allele

@print_time_taken
def get_g_group_exons(allele_to_g_groups, sequence_data):
    g_group_to_allele = generate_allele_dict(allele_to_g_groups)

    print("INFO: Finding common g group sequences")
    common_sequence = dict()
    for g_group, alleles in g_group_to_allele.items():
        if g_group == "None":
            continue

        peptide_domain_sequence = sequence_data[alleles[0]]["peptide_binding_domain"]

        for allele in alleles[1:]:
            if len(sequence_data[allele]["peptide_binding_domain"]) != len(peptide_domain_sequence):
                print("WARN: Number of peptide binding domain related exons not consistent across g group")

            for entry in sequence_data[allele]["peptide_binding_domain"]:
                corresponding_seq = [s[1] for s in peptide_domain_sequence if s[0] == entry[0]][0]
                if get_distance(entry[1], corresponding_seq) != 0:
                    print(f"ERR: Allele {allele} contains sequence not found in {g_group}")
                    return None

            common_sequence[g_group] = peptide_domain_sequence

    return common_sequence

def load_test_data(sample_file):
    samples = {}
    for record in SeqIO.parse(sample_file, "fasta"):
        samples[record.id] = record.seq
    return samples

def assign_classification_to_sample(common_sequenes, sequence, sample_name, logfile=None):
    best = (None, None, 0)
    same_dist = []

    hla_gene = get_gene(sample_name)

    for class_name, exons in common_sequenes.items():
        if hla_gene not in class_name:
            continue

        distance = 0
        for _, exon_sequence in exons:
            distance += get_distance(sequence, exon_sequence)

        if distance == best[1]:
            same_dist += [class_name]

        if best[1] == None or distance < best[1]:
            best = (class_name, distance)
            same_dist = [class_name]

    best = (*best, len(same_dist))

    if best[1] == 0:
        if logfile != None:
            logfile.write(f"For {sample_name}, assigned {best[0]} dist {best[1]}\n")
        return best

    if logfile != None:
        logfile.write(f"For {sample_name}, found no perfect match. (Best: {best[0]} with distance {best[1]})\n")

    return (None, -1, 0)

def assign_classification_to_sample_full_seq(full_sequence, sequence, full_sample_name, logfile=None, norm_distance=False, eval_metric="edit_distance"):
    best = (None, None, None, None)
    same_dist = []

    hla_gene = get_gene(full_sample_name)

    for class_name, sequence_data in full_sequence.items():
        if hla_gene not in class_name:
            continue

        distance, match_len = get_distance(sequence, sequence_data, get_length=True)

        if norm_distance:
            distance = distance/match_len

        seq_identity = match_len / (match_len + distance)

        metrics = {"edit_distance": distance, "match_length": match_len, "identity": seq_identity}
        metric_idxs = {"edit_distance": 1, "match_length": 2, "identity": 3}

        metric_idx = metric_idxs[eval_metric]
        best_metric = best[metric_idx]
        metric = metrics[eval_metric]

        maximize = lambda : metric > best_metric
        minimize = lambda : metric < best_metric
        cmp_fn = minimize if eval_metric == "edit_distance" else maximize

        if best_metric == None or cmp_fn():
            best = (class_name, distance, match_len, seq_identity)
            same_dist = [class_name]
        elif metric == best_metric:
            same_dist += [class_name]

    if logfile != None:
        logfile.write(f"For {full_sample_name}, assigned {best[0]} dist {best[1]} len {best[2]} id {best[3]} using {eval_metric}\n")

    return (*best, same_dist)

@print_time_taken
def pass_1_classification(common_sequenes, samples, g_groups_dict, truth_data=None):
    print("INFO: Beginning classification pass 1...")

    results = {}
    perfect = 0
    current = 0

    # Track results by gene
    gene_totals = {}
    gene_perfect = {}

    pass_1_logfile = open("g_group_assignment.log", "w")

    for sample_name, sample_sequence in samples.items():
        current += 1
        if current % 20 == 0:
            print(f"INFO: Processing {current}/{len(samples)} ({(current/len(samples))*100:.2f}%)")

        result = assign_classification_to_sample(common_sequenes, sample_sequence, sample_name, logfile=pass_1_logfile)
        results[sample_name] = result

        # Track by gene
        gene = get_gene(sample_name, astrisk=False)
        if gene not in gene_totals:
            gene_totals[gene] = 0
            gene_perfect[gene] = 0
        gene_totals[gene] += 1

        if result[1] == 0:
            perfect += 1
            gene_perfect[gene] += 1

    if pass_1_logfile != None:
        pass_1_logfile.close()

    print(f"\nINFO: 0 distance g group assignments: {perfect}/{len(samples)}")
    print("\n" + "="*60)
    print("PASS 1 G-GROUP ASSIGNMENT BY GENE:")
    print("="*60)
    for gene in sorted(gene_totals.keys()):
        total = gene_totals[gene]
        perf = gene_perfect[gene]
        pct = (perf / total * 100) if total > 0 else 0
        print(f"  {gene:6s}: {perf:3d}/{total:3d} = {pct:5.1f}%")
    print("="*60 + "\n")

    return results

def produce_allele_seq_db(sequence_data, selected_alleles=None, exon_only=True):
    allele_sequence_db = {}

    for allele, sequences in sequence_data.items():
        if selected_alleles != None and allele not in selected_alleles:
            continue

        segments = []
        for feature, seq in sequences.items():
            if feature == "peptide_binding_domain":
                continue

            if exon_only and feature != "Exon":
                continue

            segments += seq

        segments.sort()
        full_sequence = "".join(list(map(lambda x:x[1], segments)))
        allele_sequence_db[allele] = full_sequence

    return allele_sequence_db

@print_time_taken
def pass_2_classification(sequence_data, allele_to_g_groups, results_dict, samples, truth_data=None, exon_only=True, metric="edit_distance"):
    print("INFO: Beginning classification pass 2...")

    g_group_to_allele = generate_allele_dict(allele_to_g_groups)
    all_allele_sequence_db = produce_allele_seq_db(sequence_data, exon_only=exon_only)

    pass_2_logfile = open("3_field_allele_assignment.log", "w")
    perfect = 0
    current = 0

    results = {}
    for sample_name, classification in results_dict.items():
        current += 1
        if current % 20 == 0:
            print(f"INFO: Processing {current}/{len(results_dict)} ({(current/len(results_dict))*100:.2f}%)")

        g_group = classification[0]

        if g_group != None and ';' in g_group and \
            classification[1] == 0 and classification[2] == 1:

            print("INFO: Differing allele found in first pass")
            results[sample_name] = g_group.split(";")[1]
            perfect += 1
            continue

        search_all = g_group not in g_group_to_allele.keys()

        if not search_all:
            alleles = g_group_to_allele[g_group]
            allele_sequence_db = produce_allele_seq_db(sequence_data, selected_alleles=alleles, exon_only=exon_only)
        else:
            allele_sequence_db = all_allele_sequence_db

        result = assign_classification_to_sample_full_seq(allele_sequence_db, samples[sample_name], sample_name, logfile=pass_2_logfile, eval_metric=metric)
        equidistant = result[-1]
        if len(equidistant) > 1:
            equidistant = sorted(equidistant)
            result = (equidistant[0], *result[1:])

        if result[0] == None:
            print(f"WARN: Ecountered allele that couldn't be classified in pass 2:", sample_name)

        if result[1] == 0:
            perfect += 1
        results[sample_name] = result

    if pass_2_logfile != None:
        pass_2_logfile.close()

    print(f"INFO: 0 distance allele assignments: {perfect}/{len(samples)}. See logfile for details")

    return results

# =============================================================================
# PASS 3 BETA - Outputs detailed metrics for analysis
# =============================================================================

def load_truth_data_beta(truth_file):
    """
    Load truth data from CSV file.
    Returns: {sample_id: {gene_hap: [list of acceptable alleles]}}

    Handles:
    - Multiple acceptable alleles separated by '/'
    - Incomplete alleles (e.g., 'A*02' matches any A*02:xx:xx:xx)
    - NA entries
    """
    truth = {}

    with open(truth_file, 'r', encoding='utf-8-sig') as f:
        reader = csv.DictReader(f)
        for row in reader:
            sample_id = row['sample']
            truth[sample_id] = {}

            # Process each gene column
            for col in row.keys():
                if col == 'sample':
                    continue
                if row[col] in ['NA', '', None]:
                    continue

                # Parse gene and haplotype from column name (e.g., "HLA-DPA1_1")
                # Column format: HLA-GENE_N or just GENE_N
                gene_hap = col  # e.g., "HLA-A_1"

                # Handle multiple acceptable alleles (separated by /)
                allele_options = row[col].split('/')

                # Normalize allele names
                normalized = []
                for allele in allele_options:
                    allele = allele.strip()
                    # Add HLA- prefix if missing
                    if not allele.startswith('HLA-'):
                        allele = 'HLA-' + allele
                    normalized.append(allele)

                truth[sample_id][gene_hap] = normalized

    return truth

def allele_matches_truth(candidate, truth_alleles):
    """
    Check if candidate allele matches any of the truth alleles.
    Handles partial matches (e.g., truth 'HLA-A*02' matches 'HLA-A*02:01:01:01')
    """
    for truth in truth_alleles:
        # Exact match
        if candidate == truth:
            return True

        # Truth is partial (fewer fields) - check if candidate starts with it
        # e.g., truth = "HLA-A*02:01:01" should match "HLA-A*02:01:01:01"
        truth_fields = truth.split(':')
        cand_fields = candidate.split(':')

        if len(truth_fields) <= len(cand_fields):
            # Compare up to the number of fields in truth
            truth_prefix = ':'.join(truth_fields)
            cand_prefix = ':'.join(cand_fields[:len(truth_fields)])
            if truth_prefix == cand_prefix:
                return True

    return False

@print_time_taken
def pass_3_beta(sequence_data, results_dict, samples, truth_data, output_csv="pass3_metrics.csv"):
    """
    Beta pass 3 that outputs detailed metrics for each candidate.
    Compares against truth data for analysis.
    """
    print("INFO: Beginning classification pass 3 (BETA - metric analysis)...")

    all_allele_sequence_db = produce_allele_seq_db(sequence_data, exon_only=False)

    # CSV output for analysis
    csv_fields = [
        "sample_name", "sample_id", "gene", "haplotype",
        "candidate_allele", "is_truth", "is_best_match_length", "is_best_identity_v2",
        "is_best_mismatch_only", "is_best_weighted_5x", "is_best_weighted_10x",
        "matches", "mismatches", "gap_events", "total_gap_bp",
        "insertions", "insertion_bp", "deletions", "deletion_bp",
        "raw_edit_distance", "gap_compressed_distance",
        "match_length", "identity_original", "identity_v2",
        "identity_mismatch_only", "identity_weighted_5x", "identity_weighted_10x",
        "truth_alleles"
    ]

    csv_file = open(output_csv, "w", newline='')
    csv_writer = csv.DictWriter(csv_file, fieldnames=csv_fields)
    csv_writer.writeheader()

    current = 0
    for sample_name, classification in results_dict.items():
        current += 1
        if current % 20 == 0:
            print(f"INFO: Processing {current}/{len(results_dict)} ({(current/len(results_dict))*100:.2f}%)")

        classified_allele = classification[0]

        if classified_allele == None:
            print(f"ERR: Encountered allele that couldn't be classified:", sample_name)
            continue

        fields = classified_allele.split(":")

        if len(fields) > 3:
            fields = fields[:-1]
        else:
            continue  # Skip samples without fourth field potential

        exon_match = ":".join(fields)

        # Build list of all alleles that match the exon substring
        all_exon_matches = []
        for allele in sequence_data.keys():
            dist = get_distance(exon_match, allele)
            if dist == 0:
                all_exon_matches.append(allele)

        if not all_exon_matches:
            print(f"WARN: No matching alleles found for {sample_name}")
            continue

        # Get sample info
        sample_id = get_sampleid(sample_name)
        gene = get_gene(sample_name, astrisk=False)
        haplotype = get_haplotype_index(sample_name)

        # Get truth alleles for this sample/gene/haplotype
        gene_hap_key = f"HLA-{gene}_{haplotype}"
        truth_alleles = []
        if sample_id in truth_data and gene_hap_key in truth_data[sample_id]:
            truth_alleles = truth_data[sample_id][gene_hap_key]

        # Collect metrics for all candidates
        candidate_metrics = []

        allele_sequence_db = produce_allele_seq_db(sequence_data, selected_alleles=all_exon_matches, exon_only=False)

        for candidate_allele, ref_sequence in allele_sequence_db.items():
            sample_sequence = samples[sample_name]

            # Get detailed metrics
            detailed = get_detailed_metrics(str(sample_sequence), ref_sequence)
            identities = compute_identity_metrics(detailed)

            # Check if this is a truth allele
            is_truth = allele_matches_truth(candidate_allele, truth_alleles)

            candidate_metrics.append({
                "allele": candidate_allele,
                "is_truth": is_truth,
                "detailed": detailed,
                "identities": identities
            })

        # Determine best by each metric
        if candidate_metrics:
            best_match_length = max(candidate_metrics, key=lambda x: x["identities"]["match_length"])["allele"]
            best_identity_v2 = max(candidate_metrics, key=lambda x: x["identities"]["identity_v2"])["allele"]
            best_mismatch_only = max(candidate_metrics, key=lambda x: x["identities"]["identity_mismatch_only"])["allele"]
            best_weighted_5x = max(candidate_metrics, key=lambda x: x["identities"]["identity_weighted_5x"])["allele"]
            best_weighted_10x = max(candidate_metrics, key=lambda x: x["identities"]["identity_weighted_10x"])["allele"]

        # Write all candidates to CSV
        for cm in candidate_metrics:
            row = {
                "sample_name": sample_name,
                "sample_id": sample_id,
                "gene": gene,
                "haplotype": haplotype,
                "candidate_allele": cm["allele"],
                "is_truth": cm["is_truth"],
                "is_best_match_length": cm["allele"] == best_match_length,
                "is_best_identity_v2": cm["allele"] == best_identity_v2,
                "is_best_mismatch_only": cm["allele"] == best_mismatch_only,
                "is_best_weighted_5x": cm["allele"] == best_weighted_5x,
                "is_best_weighted_10x": cm["allele"] == best_weighted_10x,
                "matches": cm["detailed"]["matches"],
                "mismatches": cm["detailed"]["mismatches"],
                "gap_events": cm["detailed"]["gap_events"],
                "total_gap_bp": cm["detailed"]["total_gap_bp"],
                "insertions": cm["detailed"]["insertions"],
                "insertion_bp": cm["detailed"]["insertion_bp"],
                "deletions": cm["detailed"]["deletions"],
                "deletion_bp": cm["detailed"]["deletion_bp"],
                "raw_edit_distance": cm["detailed"]["raw_edit_distance"],
                "gap_compressed_distance": cm["detailed"]["gap_compressed_distance"],
                "match_length": cm["identities"]["match_length"],
                "identity_original": cm["identities"]["identity_original"],
                "identity_v2": cm["identities"]["identity_v2"],
                "identity_mismatch_only": cm["identities"]["identity_mismatch_only"],
                "identity_weighted_5x": cm["identities"]["identity_weighted_5x"],
                "identity_weighted_10x": cm["identities"]["identity_weighted_10x"],
                "truth_alleles": "|".join(truth_alleles) if truth_alleles else "UNKNOWN"
            }
            csv_writer.writerow(row)

    csv_file.close()
    print(f"INFO: Detailed metrics written to {output_csv}")

# =============================================================================
# MAIN FUNCTIONS
# =============================================================================

def run_classification_beta(reference_xml_file, samples_file, full_sample_file, truth_file,
                            pass2_metric="match_length", output_csv="pass3_metrics.csv"):
    """Run classification with beta pass 3 for metric analysis"""

    truth_data = load_truth_data_beta(truth_file)
    print(f"INFO: Loaded truth data for {len(truth_data)} samples")

    g_group_dict, p_group_dict, sequence_data = build_g_group_dict(reference_xml_file)

    g_group_common_sequences = get_g_group_exons(g_group_dict, sequence_data)
    if g_group_common_sequences == None:
        exit(1)

    samples = load_test_data(samples_file)
    full_samples = load_test_data(full_sample_file)

    print("INFO: Classifying samples to G group")
    g_group_classifications = pass_1_classification(g_group_common_sequences, samples, g_group_dict)

    print("INFO: Classifying the samples to an allele (pass 2)")
    allele_classifications = pass_2_classification(sequence_data, g_group_dict, g_group_classifications, samples, metric=pass2_metric)

    print("INFO: Running beta pass 3 with detailed metrics output")
    pass_3_beta(sequence_data, allele_classifications, full_samples, truth_data, output_csv=output_csv)

    print("INFO: Beta classification complete. Analyze", output_csv, "to determine optimal metric.")

def run_pass1_only(reference_xml_file, samples_file):
    """Run only pass 1 classification and output G-group assignment by gene"""

    g_group_dict, p_group_dict, sequence_data = build_g_group_dict(reference_xml_file)

    g_group_common_sequences = get_g_group_exons(g_group_dict, sequence_data)
    if g_group_common_sequences == None:
        exit(1)

    samples = load_test_data(samples_file)

    print("INFO: Classifying samples to G group (PASS 1 ONLY)")
    g_group_classifications = pass_1_classification(g_group_common_sequences, samples, g_group_dict)

    print("INFO: Pass 1 complete. Exiting (--pass1-only mode).")


# CLI interface
if __name__ == "__main__":
    parser = argparse.ArgumentParser(description='Beta HLA Typer - Metric Analysis')
    parser.add_argument('--hla-xml', required=True, help='Input hla.xml reference file')
    parser.add_argument('--samples', required=True, help='Input FASTA file with CDS sequences')
    parser.add_argument('--full-sequence', required=False, help='Input FASTA file with full gene sequences')
    parser.add_argument('--truth', required=False, help='Input CSV file containing truth data')
    parser.add_argument("--pass2-metric", default="match_length", help="Metric for pass 2")
    parser.add_argument("--output", default="pass3_metrics.csv", help="Output CSV for metric analysis")
    parser.add_argument("--pass1-only", action='store_true', help="Only run pass 1 (G-group assignment) and exit")

    args = parser.parse_args()

    if args.pass1_only:
        run_pass1_only(
            reference_xml_file=args.hla_xml,
            samples_file=args.samples
        )
    else:
        if not args.full_sequence or not args.truth:
            print("ERROR: --full-sequence and --truth are required unless using --pass1-only")
            exit(1)
        run_classification_beta(
            reference_xml_file=args.hla_xml,
            samples_file=args.samples,
            full_sample_file=args.full_sequence,
            truth_file=args.truth,
            pass2_metric=args.pass2_metric,
            output_csv=args.output
        )

# Pipeline entrypoint
def main(reference_xml_file, hla_fasta_dir, sample_ID, truth_file,
         pass2_metric="match_length", output_csv="pass3_metrics.csv"):
    """Entrypoint for pipeline integration"""
    samples_file = os.path.join(hla_fasta_dir, str(sample_ID) + "_HLA_haplotypes_CDS.fasta")
    full_sample_file = os.path.join(hla_fasta_dir, str(sample_ID) + "_HLA_haplotypes_gene.fasta")

    run_classification_beta(
        reference_xml_file=reference_xml_file,
        samples_file=samples_file,
        full_sample_file=full_sample_file,
        truth_file=truth_file,
        pass2_metric=pass2_metric,
        output_csv=output_csv
    )
