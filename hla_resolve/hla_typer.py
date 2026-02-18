from Bio import SeqIO
import xml.etree.ElementTree as ET
import pandas as pd
import edlib
import time
import math
import re

import argparse
import os

NUM_SAMPLES = None

# Known problematic entries whose exons 2 (or 3 if applicable) do not
# match with other alleles in their G group
non_matching_exons = {"HLA-C*07:02:01:17N": "C*07:02:01G",
                      "HLA-C*15:02:01:08N": "C*15:02:01G",
                      "HLA-DRB4*01:03:01:02N": "DRB4*01:01:01G",
                      "HLA-DRB4*01:03:01:13N": "DRB4*01:01:01G",
                      "HLA-DRB4*01:03:01:26N": "DRB4*01:01:01G",
                      "HLA-DRB4*01:03:01:28N": "DRB4*01:01:01G",
                      "HLA-DRB4*01:03:41N": "DRB4*01:01:01G",
                      "HLA-DMB*01:05": "DMB*01:01:01G"}

# Get Gene from sample name
def get_gene(sample_name, astrisk = True):
    # Get HLA-<gene>_<idx> portion, then get portion between '_' and '-'
    hla_gene_index = sample_name[sample_name.find("HLA"):]
    hla_gene = hla_gene_index.split("_")[0]
    hla_gene = hla_gene[hla_gene.index("-")+1:]
    if astrisk: hla_gene = hla_gene + "*"
    return hla_gene

# Get sample ID from sample name
def get_sampleid(sample_name):
    # Get everything before 'HLA'
    return sample_name[:sample_name.find("HLA")-1]

# Function decorator to print time taken to run a function
def print_time_taken(fun):
    def wrapper(self, *args, **kwargs):
        # Get start time
        start_time = time.time()

        # Run decorated function
        result = fun(self, *args, **kwargs)

        # Get end time, and print minutes and seconds taken
        end_time = time.time()
        time_taken = end_time-start_time
        seconds = time_taken % 60
        minutes = math.floor(time_taken / 60)
        if NUM_SAMPLES == None:
            print(f"INFO: {fun.__name__} took {minutes}m {seconds:2.0f}s")
        else:
            time_per_sample = time_taken / NUM_SAMPLES
            seconds_per_sample = time_per_sample % 60
            minutes_per_sample = math.floor(time_per_sample / 60)
            print(f"INFO: {fun.__name__} took {minutes}m {seconds:2.0f}s")
            print(f"INFO: {fun.__name__} took {minutes_per_sample}m {seconds_per_sample:2.0f}s per sample")
        
        # Return result of decorated function
        return result
    
    return wrapper

# Pull hla xml database from web. Put in ./tmp directory
def build_g_group_dict_from_web():
    # These imports are only needed here
    from zipfile import ZipFile
    import os
    import requests
    import datetime

    db_url = "https://ftp.ebi.ac.uk/pub/databases/ipd/imgt/hla/xml/hla.xml.zip"

    file_exists = os.path.isfile("tmp/hla.xml")

    # If file exists, check if up to date
    if file_exists:
        # Get file last modified timestamp as datetime object
        file_modified = os.path.getmtime("tmp/hla.xml")
        file_modified = datetime.datetime.fromtimestamp(file_modified)

        # Get database last modified timestamp as datetime object
        response = requests.head(db_url)
        response.raise_for_status()
        db_modified = response.headers.get('Last-Modified')
        db_modified = datetime.datetime.strptime(db_modified, '%a, %d %b %Y %H:%M:%S %Z')

        # If file is more recent than latest database update, no need to get new file
        if file_modified > db_modified:
            print("INFO: Cached xml pulled from web is up-to-date")
            return
        else:
            print("INFO: Cached xml pulled from web is not up-to-date. Redownloading...")
            os.system("rm tmp/hla.xml.zip")
    else:
        print("INFO: No cached xml database. Pulling from web...")

    # Get the zip file
    command = "wget -P ./tmp " + db_url 
    os.system("mkdir -p tmp")
    os.system(command)
    
    # Extract using zipfile library, unzip command not available
    zip = ZipFile('tmp/hla.xml.zip')
    zip.extractall('tmp/')
    zip.close()

# Builds dictionary databases of allele's g group associations and code
# Input: xml_file: path to XML file with referebnce database info
#        ignore_unconfirmed: (bool) ignore XML database entries marked as 'unconfirmed'
#        ignore_incomplete: (bool) ignore XML entries missing ANY features
# Outputs: g_groups: {allele: g_group}
#          sequence_data: {allele: {feature: sequence}}
@print_time_taken
def build_g_group_dict(xml_file, ignore_unconfirmed=False, ignore_incomplete=False):
    # If None, pull from web
    # if xml_file == None:
    #     build_g_group_dict_from_web()
    #     xml_file = "./tmp/hla.xml"

    # Parse metadata XML file
    print("INFO: Parsing metadata XML file")
    tree = ET.parse(xml_file)
    root = tree.getroot()

    # Handle version specific changes
    version = root.attrib['version']
    print("INFO: XML version:", version)
    group_accessor = "name" if version > "3.61.0" else "status"

    # Dictionary databases
    g_groups = dict()
    p_groups = dict()
    sequence_data = dict()

    # All tags prefixed with this string
    tag = lambda x:'{http://hla.alleles.org/xml}' + x

    # From alleles, get every allele's g group
    alleles = root.find(tag('alleles'))
    for allele in alleles.findall(tag('allele')):
        # Skip unconfirmed entries
        if ignore_unconfirmed and allele.find(tag('releaseversions')).attrib['confirmed'] == "Unconfirmed":
            continue

        # Class (I or II)
        class_type = allele.find(tag('locus')).attrib['class']
        
        # Sequence section of file structure
        sequence = allele.find(tag('sequence'))
        if sequence == None:
            continue

        # DNA sequence, and feature annotation (exon, intron, etc.)
        dna_sequence = sequence.find(tag('nucsequence'))
        features = sequence.findall(tag('feature'))
        if dna_sequence == None or features == None:
            continue

        # Check all features are present using their orders
        feature_orders = sorted([int(feature.attrib['order']) for feature in features if 'order' in feature.attrib.keys()])
        complete_order_list = list(range(1, len(feature_orders)+1))
        if ignore_incomplete and feature_orders != complete_order_list:
            continue

        # Get full allele name and g group
        name = allele.attrib["name"]

        # Collect sequences by type
        sequence_data[name] = dict()
        for i, feature in enumerate(features):
            feature_type = feature.attrib["featuretype"]

            # Skip features that don't denote region
            if feature_type not in ["UTR", "Exon", "Intron"]:
                continue

            if feature_type not in sequence_data[name].keys():
                sequence_data[name][feature_type] = []

            # Get snippet for this feature
            locations = feature.find(tag('SequenceCoordinates'))
            start, end = int(locations.attrib["start"]), int(locations.attrib["end"])
            cut_sequence = dna_sequence.text[start-1:end]

            sequence_data[name][feature_type].append((i, cut_sequence))

            # Add sequences relevant to determining G group here too
            if feature.attrib["name"] == "Exon 2":
                sequence_data[name]["peptide_binding_domain"] = [(2, cut_sequence)]

            if feature.attrib["name"] == "Exon 3" and class_type == "I":
                # Case where class 1 entry doesn't have Exon 2? Just HLA-E*01:04
                if sequence_data[name].get("peptide_binding_domain") == None:
                    sequence_data[name]["peptide_binding_domain"] = [(3, cut_sequence)]
                    continue

                sequence_data[name]["peptide_binding_domain"].append((3, cut_sequence))

        # Add G group if possible
        g_group_xml = allele.find(tag('hla_g_group'))
        if g_group_xml != None:
            g_group = g_group_xml.attrib[group_accessor]
            g_groups[name] = g_group

        # Add G group if possible
        p_group_xml = allele.find(tag('hla_p_group'))
        if p_group_xml != None:
            p_group = p_group_xml.attrib[group_accessor]
            p_groups[name] = p_group

    return g_groups, p_groups, sequence_data

# Create dictionary containing all alleles within a g group
# Input: allele_to_g_groups: {allele: g_group}
# Output: g_group_to_allele: {g_group: alleles}
def generate_allele_dict(allele_to_g_groups):
    g_group_to_allele = dict()
    for k, v in allele_to_g_groups.items():
        if v not in g_group_to_allele.keys():
            g_group_to_allele[v] = []

        g_group_to_allele[v].append(k)

    return g_group_to_allele

# Get the edit distance between two samples
# Inputs:   sequence: sample sequence (str)
#           sequence_data: reference sequence (str)
#           get_length: return length (bool) (optional)
# Output:   edit distance (int), match_length (int) (optional)
def get_distance(sequence, sequence_data, get_length=False, gap_compressed=True, get_detailed=False):
    # Don't spend extra time calculating locations unless necessary
    task = "distance"
    
    # If gap compressed enabled, or match length, we need to calculate full path
    if gap_compressed or get_length:
        task = "path"

    # Ensure first argument is shorter sequence
    if len(sequence) <= len(sequence_data):
        results = edlib.align(sequence, sequence_data, task=task, mode="HW")
    else:
        results = edlib.align(sequence_data, sequence, task=task, mode="HW")

    edit_distance = results["editDistance"]
    
    # If using gap compressed distance, adjust to remove extra distance
    if gap_compressed or get_distance:
        # Find all instances of insertion or deletion in cigar
        pattern = r"(\d+)(I|D)"
        for match in re.finditer(pattern, results["cigar"]):

            # For any insertion/deletion longer than 1 bp,
            # remove it from the edit distance leaving only
            # one behind. All insertions/deletions count as 1
            length = int(match.group(1))
            if length > 1:
                edit_distance -= length - 1

    # If we want to get match length, calculate from cigar and return it too
    if get_length or get_detailed:
        match_length = 0
        mis_match_length = 0

        # Sum all matches
        pattern = r"(\d+)(=)"
        for match in re.finditer(pattern, results["cigar"]):
            length = int(match.group(1))
            match_length += length

        # Sum all matches
        pattern = r"(\d+)(X)"
        for match in re.finditer(pattern, results["cigar"]):
            length = int(match.group(1))
            mis_match_length += length

        if get_detailed:
            prop_mismatch = match_length / (match_length + mis_match_length) if (match_length + mis_match_length) > 0 else 0
            return (results["cigar"], results["locations"][0][0], results["locations"][0][1], results["editDistance"], edit_distance, prop_mismatch, match_length)

        return edit_distance, match_length, mis_match_length
    
    # Otherwise just return edit distance
    return edit_distance

# Get common peptide binding exons for each g group
# Inputs: allele_to_g_groups: {allele: g_group}
#         sequence_data: {allele: {feature: sequence}}
# Ouput: common_sequence: {g_group: sequence}
#        orphan_alleles: {allele: sequence}
@print_time_taken
def get_g_group_exons(allele_to_g_groups, sequence_data):
    # Create dictionary containing all alleles within a g group
    g_group_to_allele = generate_allele_dict(allele_to_g_groups)

    print("INFO: Finding common g group sequences")
    common_sequence = dict()
    # Check that peptide binding domain is same for every allele in g group
    for g_group, alleles in g_group_to_allele.items():
        if g_group == "None":
            continue

        # Grab first one as comparison reference
        peptide_domain_sequence = sequence_data[alleles[0]]["peptide_binding_domain"]

        # Iterate over other alleles to compare
        for allele in alleles[1:]:
            # Skip known issue alleles
            # TODO: This is nto required. Talk to Matt about removal
            # if allele in non_matching_exons.keys():
            #     print("INFO: Skipping known problematic allele", allele)
                
            #     # Add second entry of this allele for matching algorithm
            #     differing_ars_sequence = sequence_data[allele]["peptide_binding_domain"]

            #     # We need a unique entry in the dictionary for this differing sequence.
            #     # Use <g_group>;<allele>
            #     g_group_entry_name = f"{g_group};{allele}"

            #     common_sequence[g_group_entry_name] = differing_ars_sequence

            #     continue

            if len(sequence_data[allele]["peptide_binding_domain"]) != len(peptide_domain_sequence):
                print("WARN: Number of peptide binding domain related exons not consistent across g group")

            # Check each entry in the peptide binding domain sequence
            # Will contiain Exon 2 and 3 for type I, Exon 2 for type II
            for entry in sequence_data[allele]["peptide_binding_domain"]:

                # If entries not in common, return
                corresponding_seq = [s[1] for s in peptide_domain_sequence if s[0] == entry[0]][0]
                if get_distance(entry[1], corresponding_seq) != 0:
                    print(f"ERR: Allele {allele} contains sequence not found in {g_group}")
                    print("INFO: Test sequence:", entry)
                    print("INFO: Corresponding G Group sequence:", corresponding_seq)
                    print("INFO: Calculated distance:", get_distance(entry[1], corresponding_seq))
                    return None
                
            common_sequence[g_group] = peptide_domain_sequence

    return common_sequence #, orphan_alleles

# Load samples from a fasta file
# Input: Path to file with samples
# Output: Dictionary containing {sample_name: sequence}
def load_test_data(sample_file):
    samples = {}
    for record in SeqIO.parse(sample_file, "fasta"):
        samples[record.id] = record.seq

    return samples

# Do first pass of classification (to g group level)
# This function operates on concatinated exon sequence
# Input: common_sequences: {class: (exon, sequence)} Common G group ARS
#        sequence: (str) Sample exon sequence
#        sample_name: (str) Name of sample
#        logfile: (file | None) (optional)
# Output: (class, distance, num_equidistant)

def assign_classification_to_sample(common_sequenes, sequence, sample_name, logfile=None):
    best = (None, None, 0) # (class_name, distance, num_matches)
    same_dist = []

    hla_gene = get_gene(sample_name)

    # Go through each reference sequence and find closest match
    for class_name, exons in common_sequenes.items():
        if hla_gene not in class_name:
            continue

        # Check all exons for reference sample
        distance = 0
        for _, exon_sequence in exons: # one or two exons
            distance += get_distance(sequence, exon_sequence)

        if distance == best[1]:
            same_dist += [class_name]

        # If first, or closest sample, set as best
        if best[1] == None or distance < best[1]:
            best = (class_name, distance)
            same_dist = [class_name]

    best = (*best, len(same_dist))

    # Return found g group if edit distance zero
    if best[1] == 0:
        if logfile != None:
            logfile.write(f"For {sample_name}, assigned {best[0]} dist {best[1]}\n")
            if len(same_dist) > 1:
                logfile.write(f"Equidistant for {sample_name}: {', '.join(same_dist)}\n")
        return best
    
    # No exact g group found, do unrestricted search during pass 2
    if logfile != None:
        logfile.write(f"For {sample_name}, found no perfect match. (Best: {best[0]} with distance {best[1]})\n")

    return (None, -1, 0)

# Do full sequence based classification. Used in pass 2 and 3
# Input: full_sequence: {class: sequence} Full sequences of samples to check
#        sequence: (str) Sample sequence to classify
#        full_sample_name: (str) Name of sample
#        logfile: (file | None) (optional)
#        norm_distance: (bool)
#        eval_metric: (str: "edit_distance" | "match_length" | "identity")
# Output: (class, distance, match_length, num_equidistant)
def assign_classification_to_sample_full_seq(full_sequence, sequence, full_sample_name, logfile=None, norm_distance=False, eval_metric="edit_distance"):
    best = (None, None, None, None, None) # (class_name, distance, match_length, identity, mismatch_identity)
    same_dist = []

    hla_gene = get_gene(full_sample_name)

    for class_name, sequence_data in full_sequence.items():
        # Skip if not the same gene
        if hla_gene not in class_name:
            continue

        # Get distance metrics between test and sample sequences
        distance, match_len, mismatch_len = get_distance(sequence, sequence_data, get_length=True)

        # Normalize edit distance by match length if this is enabled
        if norm_distance:
            distance = distance/match_len

        # Sequence identity, 1 - (edit distance/match length)
        seq_identity = 1 - (distance / match_len)

        mismatch_identity = match_len / (match_len + mismatch_len) if (match_len+mismatch_len) > 0 else 0

        # LUT for evaluation metric based on provided choice
        metrics = {"edit_distance": distance, "match_length": match_len, "identity": seq_identity, "mismatch_identity": mismatch_identity}
        metric_idxs = {"edit_distance": 1, "match_length": 2, "identity": 3, "mismatch_identity": 4} # Index into saved best match

        # Get index into best match based on chosen evaluation metric
        metric_idx = metric_idxs[eval_metric]
        
        # Grab desired metric from saved best match
        best_metric = best[metric_idx]

        # Get desired metric for the current test sample
        metric = metrics[eval_metric]

        # Either choose to maximize or minimize metric. Min(edit dist), or Max(match len OR identity)
        maximize = lambda : metric > best_metric
        minimize = lambda : metric < best_metric
        cmp_fn = minimize if eval_metric == "edit_distance" else maximize

        # If first sample, or closest sample, set as best
        if best_metric == None or cmp_fn():
            best = (class_name, distance, match_len, seq_identity, mismatch_identity)
            same_dist = [class_name]
        elif metric == best_metric:
            if eval_metric == "mismatch_identity":
                if match_len > best[2]:
                    best = (class_name, distance, match_len, seq_identity, mismatch_identity)
                    same_dist = [class_name]
                elif match_len == best[2]:
                    same_dist += [class_name]
            else: same_dist += [class_name]
    
    if logfile != None:
        logfile.write(f"For {full_sample_name}, assigned {best[0]} dist {best[1]} len {best[2]} id {best[3]} mismatch {best[4]} using {eval_metric}\n")
        if len(same_dist) > 1:
            logfile.write(f"Equidistant for {full_sample_name}: {', '.join(same_dist)}\n")

    return (*best, same_dist)

# Do first pass of classifiaction (to g group level)
# Input: common_sequence: {g_group: sequence} Common ARS of each G group
#        samples: {sample_name: ARS sequence} ARS of all samples to be classified
#        g_groups_dict: {allele: g_group}
#        truth_data: {sample_name: {<gene>_<num>: allele}} | None
# Output: {sample_name: (g_group, distance)}

@print_time_taken
def pass_1_classification(common_sequenes, samples, g_groups_dict, truth_data = None):
    print("INFO: Beginning classification pass 1...")

    results = {}
    perfect = 0
    current = 0

    # pass_1_logfile = None
    pass_1_logfile = open("g_group_assignment.log", "w")

    for sample_name, sample_sequence in samples.items():
        current += 1
        if current % 20 == 0:
            print(f"INFO: Processing {current}/{len(samples)} ({(current/len(samples))*100:.2f}%)")

        result = assign_classification_to_sample(common_sequenes, sample_sequence, sample_name, logfile=pass_1_logfile)
        results[sample_name] = result
        if result[1] == 0:
            perfect += 1

        # Get just unique string corresponding to sample
        sample_id = sample_name.split("_")[0]

        # If truth data was given and has an entry of this sample
        if truth_data != None and sample_id in truth_data.keys():
            dist_to_truth_g_group(sample_name, truth_data, common_sequenes, sample_sequence, pass_1_logfile, g_groups_dict)

    if pass_1_logfile != None:
        pass_1_logfile.close()

    print(f"INFO: 0 distance g group assignments: {perfect}/{len(samples)}. See logfile for details")
    
    return results

# Builds a database of relevant allele sequences
# Inputs:   sequence_data: {allele: {feature: (position, sequence)}}
#           selected_alleles: [allele] | None
# Outputs:  {allele: sequence}
def produce_allele_seq_db(sequence_data, selected_alleles = None, exon_only = True):
    allele_sequence_db = {}

    # Iterate over the sequence data for all alleles
    for allele, sequences in sequence_data.items():
        # Skip alleles we don't want to seach
        if selected_alleles != None and allele not in selected_alleles:
            continue

        # Start building up the full sequence for this alleles
        segments = []
        for feature, seq in sequences.items():
            # Skip duplicate entries for Exon 2/3 used in first pass
            if feature == "peptide_binding_domain":
                continue

            # Skip non-exons if we are classifying concatinated exon data
            if exon_only and feature != "Exon":
                continue

            # Add entry of form (index, sequence)
            segments += seq

        # Sort based on index (first element in each tuple)
        segments.sort()

        # Assemble full sequence (remove index and concatinate)
        full_sequence = "".join(list(map(lambda x:x[1], segments)))
        allele_sequence_db[allele] = full_sequence

    return allele_sequence_db

# General function for matching partially incomplete name to database of proper names
# Matches based on keys of db. Used in truth distance checking functions to resolve
# incomplete allele names so they can be used to access dictionary
# Inputs:   x:  (str)
#           db: {allele: (value never used/accessed)}
# Outputs:  (str)
def match_partial(x, db):
    match = x    
    
    # Some alleles are not full. Ex: HLA-C*03:04:01. This is one of the HLA-C*03:04:01:* genes
    # Use our existing matching functions for allele names.

    # If exact match for allele name not found in our seq database
    if not x in db.keys():
        # Find best match for the allele name
        best = (None, None)
        for allele in db.keys():
            dist = get_distance(x, allele, gap_compressed=False)
            if best[1] == None or dist < best[1]:
                best = (allele, dist)

        match = best[0]
        
    return match

# Function to grab the allele from the truth data corresponding to the sample
# This function handles reformatting, ensuring leading zeros, etc.
# Inputs:   sample_name: (str)
#           truth_data: {sample_id: {<gene>_<index>: allele}}
#           index: (int) | None
# Outputs:  allele (str)
def get_allele(sample_name, truth_data, index=None):
    # Get just unique string corresponding to sample
    sample_id = sample_name.split("_")[0]

    # This will run twice for each sample/gene pair, so grab the index
    # Get gene/index with form <gene>_<index>. Ex: A_1
    hla_gene = get_gene(sample_name)

    # If index is None, grab from sample name. Otherwise use provided
    if index == None:
        # TODO: Is this reliable way to grab index?
        gene_with_index = f"{hla_gene[0]}_{sample_name[-1]}"
    else:
        gene_with_index = f"{hla_gene[0]}_{index}"

    # Get correct alleles for this sample
    correct_allele = truth_data[sample_id][gene_with_index]

    fields = correct_allele.split("*")[1].split(":")

    # first field is one digit, zero padding required
    if len(fields[0]) < 2:
        # Prepend zero for padding
        fields[0] = "0"+fields[0]

        # Reconstruct with zero padded first field
        correct_allele = f"{hla_gene}*{':'.join(fields)}"

    # Our database contains HLA- prefix
    correct_allele = "HLA-"+correct_allele

    return correct_allele, sample_id, gene_with_index

# Get the distance of a sample from the true allele
# Input:    full_sample_name: (str)
#           truth_data: {sample_name: {<gene>_<num>: allele}}
#           sequence_db: {g_group: ars}
#           sequence: (str)
#           logfile: (file)
#           allele_to_g_group: {allele: g_group}
#           index: (int) | None
# Output:   None
def dist_to_truth_g_group(sample_name, truth_data, sequence_db, sequence, logfile, allele_to_g_group, index=None):
    # Repeat for both alleles. Pick best to report (of id index set, do that)
    for allele_index in (1,2):
        # Reformat to standard format
        correct_allele, sample_id, gene_with_index = get_allele(sample_name, truth_data, index=allele_index)

        # Acceptable alleles are keys in the allele to g group dictionary
        correct_allele = match_partial(correct_allele, allele_to_g_group)

        # Convert to G group. Sequence db should b G group ARS.
        correct_g_group = allele_to_g_group[correct_allele]

        # Handle case where g group is none
        if correct_g_group == "None":
            logfile.writelines(f"Sample {sample_name} with truth table entry {truth_data[sample_id][gene_with_index]}")
            logfile.writelines(f" (determined using the {gene_with_index}: {correct_allele} allele) belongs to no G group.\n")
            return

        # Get sequence, concatinate all exons. Db has form [(exon_num: seq), ...]
        correct_sequence = sequence_db[correct_g_group]
        correct_sequence = ''.join(map(lambda x:x[1], correct_sequence))

        # Get distance and log it
        correct_g_group_dist = get_distance(sequence, correct_sequence)

        # Best 0: edit distance 1: g group 2: sample id (just unique sample string) 3: <gene>_<index> (ex: A_1)
        logfile.writelines(f"Sample {sample_name} with truth table entry {truth_data[sample_id][gene_with_index]}")
        logfile.writelines(f" (calculated distance to {gene_with_index}: {correct_g_group}) had distance {correct_g_group_dist}\n")

# Get the distance of a sample from the true allele
# Input:    full_sample_name: (str)
#           truth_data: {sample_name: {<gene>_<num>: allele}}
#           sequence_db: {allele: sequence}
#           sequence: (str)
#           logfile: (file)
# Output:   None
def dist_to_truth_allele(sample_name, truth_data, sequence_db, sequence, logfile):
    # Repeat for both alleles. Pick best to report (of id index set, do that)
    for allele_index in (1,2):
        # Get just unique string corresponding to sample
        correct_allele, sample_id, gene_with_index = get_allele(sample_name, truth_data, index=allele_index)

        # Some alleles are not full. Ex: HLA-C*03:04:01. This is one of the HLA-C*03:04:01:* genes
        # Check if allele found. If so, continue as normal.
        # If not, find partial match. Missing fields are wildcards. Then proceede with first match
        correct_allele = match_partial(correct_allele, sequence_db)

        # Get sequence and calculate distance
        correct_sequence = sequence_db[correct_allele]
        correct_allele_dist, leng = get_distance(sequence, correct_sequence, get_length=True)

        logfile.writelines(f"Sample {sample_name} with truth table entry {truth_data[sample_id][gene_with_index]}")
        logfile.writelines(f" (calculated distance to {gene_with_index}: {correct_allele}) had distance {correct_allele_dist} and match length {leng}\n")

# Do second pass classification to assign allele from g group to sample
# This pass is the exon matching stage
# Input: sequence_data: {allele: {feature: sequence}}
#        allele_to_g_group: {allele: g_group}
#        results_dict: {sample_name: (g_group, distance)}
#        samples: {sample_name: sequence}
#        truth_data: {sample_name: {<gene>_<num>: allele}} | None
#        exon_only: bool
# Output: {sample_name: (g_group, distance)}
@print_time_taken
def pass_2_classification(sequence_data, allele_to_g_groups, results_dict, samples, truth_data=None, exon_only=True, metric="edit_distance"):
    print("INFO: Beginning classification pass 2...")

    # {g_group: allele}
    g_group_to_allele = generate_allele_dict(allele_to_g_groups)

    # Precompute database of all sequences
    all_allele_sequence_db = produce_allele_seq_db(sequence_data, exon_only=exon_only)

    # pass_2_logfile = None
    pass_2_logfile = open("3_field_allele_assignment.log", "w")
    perfect = 0
    current = 0

    results = {}
    for sample_name, classification in results_dict.items():
        current += 1
        if current % 20 == 0:
            print(f"INFO: Processing {current}/{len(results_dict)} ({(current/len(results_dict))*100:.2f}%)")

        # TODO: Handle here g groups of form <g_group>;<allele>
        g_group = classification[0]

        # Perfect differing allele found, and it was only match. No need to search everywhere
        # If it is not perfect, search_all will be set
        # Full search required because we don't know what next closest allele would be
        if g_group != None and ';' in g_group and \
            classification[1] == 0 and classification[2] == 1:
            
            print("INFO: Differing allele found in first pass")
            results[sample_name] = g_group.split(";")[1]
            perfect += 1
            continue

        # Handle: G group known (g_group is g group (string))
        #         Orphaned allele (g_group is "None" (string). Still a key in g_group_to_allele)
        #         G group cannot be determined (g_group is None (nonetype))

        search_all = g_group not in g_group_to_allele.keys()

        # Get all allele sequences in g group

        # Classify within known G group
        if not search_all:
            # Alleles to build database are those found within the g group
            alleles = g_group_to_allele[g_group]
            allele_sequence_db = produce_allele_seq_db(sequence_data, selected_alleles=alleles, exon_only=exon_only)
        else:
            # If we should search all, use precomputed db
            allele_sequence_db = all_allele_sequence_db

        # Re-use classification function
        result = assign_classification_to_sample_full_seq(allele_sequence_db, samples[sample_name], sample_name, logfile=pass_2_logfile, eval_metric=metric)
        equidistant = result[-1]
        if len(equidistant) > 1:
            equidistant = sorted(equidistant)
            old = result[0]
            result = (equidistant[0], *result[1:])
            
            pass_2_logfile.writelines(f"{sample_name} picked lowest {old} -> {result}\n")

        # If class is None, print warning. Will error in pass 3
        if result[0] == None:
            print(f"WARN: Ecountered allele that couldn't be classified in pass 2:", sample_name)

        # If perfect match, add to score
        if result[1] == 0:
            perfect += 1
        results[sample_name] = result

        # Get just unique string corresponding to sample
        sample_id = sample_name.split("_")[0]

        # If truth data was given and has an entry of this sample
        if truth_data != None and sample_id in truth_data.keys():
            dist_to_truth_allele(sample_name, truth_data, all_allele_sequence_db, samples[sample_name], pass_2_logfile)

    if pass_2_logfile != None:
        pass_2_logfile.close()

    print(f"INFO: 0 distance allele assignments: {perfect}/{len(samples)}. See logfile for details")

    return results

# Do a third pass, to deteremine 4th field based on introns/UTRs
# Strip fourth field (if found) and classify based on subset where 
# fourth field is only difference.
# Input: sequence_data: {allele: {feature: sequence}}
#        results_dict: {sample_name: (allele, distance)}
#        samples: {sample_name: sequence}
#        truth_data: {sample_name: {<gene>_<num>: allele}} | None
#        exon_only: bool
# Output: {sample_name: (g_group, distance)}
@print_time_taken
def pass_3_classification(sequence_data, results_dict, samples, truth_data=None, metric="identity", generate_query_ref_comp=False):
    print("INFO: Beginning classification pass 3...")

    headers = ["sample", "ref_allele_name", "CIGAR", "alignment_path_start", "alignment_path_stop", "raw_edit", "gc_edit", "prop_mismatch", "match_length"]
    entries = {header:[] for header in headers}

    # Precompute database of all sequences
    all_allele_sequence_db = produce_allele_seq_db(sequence_data, exon_only=False)

    # pass_3_logfile = None
    pass_3_logfile = open("allele_assignment.log", "w")
    perfect = 0
    current = 0

    results = {}
    for sample_name, classification in results_dict.items():
        current += 1
        if current % 20 == 0:
            print(f"INFO: Processing {current}/{len(results_dict)} ({(current/len(results_dict))*100:.2f}%)")

        # Strip last field from the allele, to get substring that matches the exon reigon
        classified_allele = classification[0]

        # Log if classification was None, and skip
        if classified_allele == None:
            print(f"ERR: Ecountered allele that couldn't be classified:", sample_name)
            if pass_3_logfile != None:
                pass_3_logfile.writelines(f"{sample_name} was assigned None in pass 2")
            continue

        fields = classified_allele.split(":")

        # For four field alleles, trim of last field so we can replace it with wildcard
        if len(fields) > 3:
            fields = fields[:-1]
        else:
            # If there is no fourth field, no need for wildcard search
            pass_3_logfile.writelines(f"{sample_name} does not have fourth field. Short circuit to {classified_allele}\n")
            results[sample_name] = (classified_allele,0,1)
            continue

        # NOTE: Idea: Maybe when we have equidistant matches, pick the longest match?

        # Construct sequence specifying all coding reigons we want to match
        exon_match = ":".join(fields)
        
        # Build list of all alleles that match the exon substring
        all_exon_matches = []
        for allele in sequence_data.keys():
            dist = get_distance(exon_match, allele)
            if dist == 0: # Allele names will match perfectly!
                all_exon_matches.append(allele)

        pass_3_logfile.writelines(f"{classified_allele} -> {exon_match} ({len(all_exon_matches)})     {','.join(all_exon_matches)}\n")

        # Database of allele data from alleles found to have matching exons
        allele_sequence_db = produce_allele_seq_db(sequence_data, selected_alleles=all_exon_matches, exon_only=False)

        # Same as pass 2, but this time, classify full sequence input against full sequence db
        result = assign_classification_to_sample_full_seq(allele_sequence_db, samples[sample_name], sample_name, logfile=pass_3_logfile, eval_metric=metric)
        if result[1] == 0:
            perfect += 1
        results[sample_name] = result


        pass_3_logfile.writelines(f"Sample: {samples[sample_name]}\n")
        pass_3_logfile.writelines(f"Result: {allele_sequence_db[result[0]]}\n")

        # Get just unique string corresponding to sample
        sample_id = get_sampleid(sample_name)

        # If truth data was given and has an entry of this sample
        if truth_data != None and sample_id in truth_data.keys():
            dist_to_truth_allele(sample_name, truth_data, all_allele_sequence_db, samples[sample_name], pass_3_logfile)

        if generate_query_ref_comp:
            detailed_metrics = get_distance(allele_sequence_db[result[0]], samples[sample_name], get_detailed=True)
            entries["sample"].append(sample_name)
            entries["ref_allele_name"].append(result[0])
            entries["CIGAR"].append(detailed_metrics[0])
            entries["alignment_path_start"].append(detailed_metrics[1])
            entries["alignment_path_stop"].append(detailed_metrics[2])
            entries["raw_edit"].append(detailed_metrics[3])
            entries["gc_edit"].append(detailed_metrics[4])
            entries["prop_mismatch"].append(detailed_metrics[5])
            entries["match_length"].append(detailed_metrics[6])

    if pass_3_logfile != None:
        pass_3_logfile.close()

    print(f"INFO: 0 distance refined allele assignments: {perfect}/{len(samples)}. See logfile for details")

    if generate_query_ref_comp:
        pd.DataFrame(entries, columns=headers).to_csv("sample_ref_comp.csv", index=False)

    return results

# Write results to csv file
# Input: results: {sample_name: (g_group, distance)}
#        file_path: path to output file
#        equidistant_file_path: Path to file which should contain all equidistant options
# Output: None
@print_time_taken
def output_results(results, file_path, equidistant_file_path=None):
    # Helper functions to get allele
    def get_allele (x):
        if "incomplete" in x:
            return "HLA-"+get_gene(x, astrisk=False)+"_"+x.split("_")[-2]+"_incomplete"
        else:
            return "HLA-"+get_gene(x, astrisk=False)+"_"+x.split("_")[-1]

    # Get list of all samples and alleles
    entry_names = results.keys()
    unique_samples = sorted(list(set(map(get_sampleid, entry_names))))
    unique_alleles = sorted(list(set(map(get_allele, entry_names))))

    # Array containing sample id, and empty slots for allele classification
    data = [[sample] + [None] * len(unique_alleles) for sample in unique_samples]

    if equidistant_file_path:
        data_equidist = [[sample] + [None] * len(unique_alleles) for sample in unique_samples]
    
    # Go through each entry in the results, and add it to the data array
    for entry, classification in results.items():
        sample_name = get_sampleid(entry)
        sample_allele = get_allele(entry)

        # Insert classification into correct slot in array
        name_index = unique_samples.index(sample_name)
        allele_index = unique_alleles.index(sample_allele)
        data[name_index][allele_index+1] = classification[0]

        if equidistant_file_path:
            data_equidist[name_index][allele_index+1] = ";".join(classification[-1])

    # Swap alleles 1 and 2 when 2 > 1
    for entry in range(len(data)):
        for i in range(1,len(unique_alleles),2):
            # Get alleles HLA-* 1 and 2
            allele_1 = data[entry][i]
            allele_2 = data[entry][i+1]

            # If 2 alphabetically greater than 1, swap
            if allele_1 != None and allele_2 != None:
                if allele_2 < allele_1:
                    data[entry][i] = allele_2
                    data[entry][i+1] = allele_1
            
    # Turn into dataframe and output to csv file
    df = pd.DataFrame(data, columns=["sample"]+unique_alleles)
    df.to_csv(file_path, index=False)

    if equidistant_file_path != None:
        df_equidist = pd.DataFrame(data_equidist, columns=["sample"]+unique_alleles)
        df_equidist.to_csv(equidistant_file_path, index=False)

# Load a dictionary with truth data from the file if provided
# Inputs:   path: (str | None)
#           source: (str)
# Outputs:  {sample_name: {"A_1": a_1_allele, "A_2": a_2_allele, ...}} | None
@print_time_taken
def load_truth_data(path, source = "IHW"):
    # Skip if no truth data
    if path == None:
        return None

    truth_df = pd.read_csv(path)

    # If file contains info from multiple typers, select only one
    if "source" in truth_df.columns:
        truth_df = truth_df[truth_df["source"] == source]
        
    # Copy from col 1 into col 2 for empty entires
    for col in ["A", "B", "C"]:
        null_entries = truth_df[col+"_2"].isna()
        truth_df.loc[null_entries, col+"_2"] = truth_df[col+"_1"][null_entries]

    # Ensure all alleles in truth set begin with (A|B|C)*
    for col in ["A_1","A_2","B_1","B_2","C_1","C_2"]:
        letter = col[0]
        missing_prefix = ~truth_df[col].str.contains(rf"{letter}\*", regex=True)
        truth_df[col][missing_prefix] = f"{letter}*"+truth_df[col][missing_prefix] 

    # Build dictionary of allele labels for samples
    sample_labels = {}
    for sample in truth_df["sample"]:
        sample_labels[sample] = {}

        # Grab sample data, and insert into dictionary categorized by gene
        sample_data = truth_df[truth_df["sample"] == sample]
        for gene in ["A_1", "A_2", "B_1", "B_2", "C_1", "C_2"]:
            sample_labels[sample][gene] = sample_data[gene].item()

    return sample_labels

# Write out json file containing sequence data as JSON and G group info
# Input:    reference_xml_file: (str) XML file with allele reference sequences
#           samples_file: (str) fasta file full of input sample concatinated exons
# Output:   (None) Writes to database_data.json
@print_time_taken
def write_json(sequence_data, g_group_dict):
    import json

    filename = "database_data.json"

    g_group_to_allele = generate_allele_dict(g_group_dict)

    # Reformat data. Currently in form [allele][feature_group] -> [(feature index, seq), ...]
    # Want in form [allele][feature_group][index] -> seq
    reformatted_sequence_data = {}
    for allele in sequence_data.keys():
        # [allele] entry
        reformatted_sequence_data[allele] = {}

        all_feature_data = sequence_data[allele]
        for feature_name in all_feature_data.keys():
            # [feature_group] entry
            reformatted_sequence_data[allele][feature_name] = {}

            # [(feature index, seq), ...] list
            features = all_feature_data[feature_name]

            # Create a mapping from the feature index to 0 based sequential index
            # feature_indicies = {x:int(i) for i, x in enumerate(list(map(lambda x:x[0], features)))}

            reformatted_sequence_data[allele][feature_name] = list(map(lambda x:x[1], features))

            # For each of those features, create a [allele][feature_group][index] entry and assign sequence
            # for feature in features:
            #     feature_index = feature[0]
            #     feature_sequence = feature[1]
            #     reformatted_sequence_data[allele][feature_name][feature_indicies[feature_index]] = feature_sequence

    with open(filename, "w") as f:
        json.dump({"sequence_data": reformatted_sequence_data, "g_group_alleles": g_group_to_allele}, f, indent=4)

# Runs all classification passes and writes results
# Called by main() in the pipeline
# Input:    reference_xml_file: (str) XML file with allele reference sequences
#           samples_file: (str) fasta file full of input sample concatinated exons
#           full_sample_file: (str | None) fasta file with full sample sequences for fourth field
#           truth_file: (str | None) csv file with true allele calls, for tracing
#           pass2_metric: (str) assignment metric to use during classification pass 2
#           pass3_metric: (str) assignment metric to use during classification pass 3
#           ignore_unconfirmed: (bool) ignore XML database entries marked as 'unconfirmed'
#           ignore_incomplete: (bool) ignore XML entries missing ANY features
# Output:   (None) Writes to assignement.log and output.csv files for each stage
def run_classification(reference_xml_file, samples_file, full_sample_file=None, truth_file=None, 
                       pass2_metric="edit_distance", pass3_metric="mismatch_identity", ignore_unconfirmed=False,
                       ignore_incomplete=False, write_full=False, generate_query_ref_comp=False):

    truth_data = load_truth_data(truth_file)

    g_group_dict, p_group_dict, sequence_data = build_g_group_dict(reference_xml_file, ignore_unconfirmed, ignore_incomplete)

    g_group_common_sequences = get_g_group_exons(g_group_dict, sequence_data)
    if g_group_common_sequences == None:
        exit(1)

    write_json(sequence_data, g_group_dict)
        
    samples = load_test_data(samples_file)
    NUM_SAMPLES = len(samples)

    if full_sample_file != None:
        full_samples = load_test_data(full_sample_file)
        if NUM_SAMPLES != len(full_samples):
            print(f"ERR: Different samples in exon only and full sequence sample files.")
            print(f"{NUM_SAMPLES} (in exon only) vs {len(full_sample_file)} (in full sequence)")
            exit(1)

    print("INFO: Classifying samples to G group")
    g_group_classifications = pass_1_classification(g_group_common_sequences, samples, g_group_dict, truth_data=truth_data)

    print("INFO: Writing g group results")
    output_results(g_group_classifications, "g_group_output.csv", "g_group_output_full.csv" if write_full else None)

    print("INFO: Classifying the samples to an allele")
    allele_classifications = pass_2_classification(sequence_data, g_group_dict, g_group_classifications, samples, truth_data=truth_data, metric=pass2_metric)

    print("INFO: Writing allele results")
    output_results(allele_classifications, "3_field_allele_output.csv", "3_field_allele_output_full.csv" if write_full else None)
    
    if full_sample_file == None:
        exit(0)

    print("INFO: Refining allele classifications based on non-coding regions", pass3_metric)
    refined_classifications = pass_3_classification(sequence_data, allele_classifications, full_samples, truth_data=truth_data, metric=pass3_metric, generate_query_ref_comp=generate_query_ref_comp)

    print("INFO: Writing refined allele results")
    output_results(refined_classifications, "allele_output.csv", "allele_output_full.csv" if write_full else None)

# CLI interface used for testing
if __name__ == "__main__":
    parser = argparse.ArgumentParser(description='Two-Step HLA Allele Identification System')
    parser.add_argument('--hla-xml', required=False, default=None, help='Input hla.xml reference file')
    parser.add_argument('--samples', required=False, default="../data/HLA_Class_I_haplotypes.fa", help='Input FASTA file with full sequences')
    parser.add_argument('--truth', required=False, default=None, help='Input csv file containg truth data, for testing purposes')
    parser.add_argument("--full-sequence", required=False, default=None, help="To enable the third intron/UTR classification stage, supply full sequence data here")
    parser.add_argument("--pass2-metric", required=False, default="edit_distance", help="Metric used to assign fourth field, 'edit_distance' (default), 'match_length', 'identity' or 'mismatch_identity'")
    parser.add_argument("--pass3-metric", required=False, default="mismatch_identity", help="Metric used to assign fourth field, 'edit_distance', 'match_length', 'identity' or 'mismatch_identity' (default)")
    parser.add_argument("--ignore-unconfirmed", action='store_true', help="Do not consider 'uncomfirmed' database entries")
    parser.add_argument("--ignore-incomplete", action='store_true', help="Do not consider database entries that are missing any features")
    parser.add_argument("--write-full", action='store_true', help="Write all equidistant options to a file that ends with ..._full.csv")
    parser.add_argument("--generate-query-ref-comp", action='store_true', help="Generate CSV containing query-reference comparisons. File is named 'sample_ref_comp.csv'")
    
    args = parser.parse_args()

    xml_file = args.hla_xml
    samples_file = args.samples
    truth_file = args.truth
    full_sample_file = args.full_sequence
    pass3_metric = args.pass3_metric
    pass2_metric = args.pass2_metric
    ignore_unconfirmed = args.ignore_unconfirmed
    ignore_incomplete = args.ignore_incomplete
    write_full = args.write_full
    generate_query_ref_comp = args.generate_query_ref_comp

    run_classification(xml_file, samples_file, full_sample_file, truth_file, pass2_metric, pass3_metric, ignore_unconfirmed, ignore_incomplete, write_full, generate_query_ref_comp)

# Main entrypoint from pipeline
# Input:    reference_xml_file: (str) Path to XML file with allele reference sequences
#           hla_fasta_dir: (str) Path to directory containing pipeline output fasta files
#           sample_ID: (str | int) ID number of the sample to be processed
#           pass2_metric: (str) assignment metric to use during classification pass 2
#           pass3_metric: (str) assignment metric to use during classification pass 3
#           ignore_unconfirmed: (bool) ignore XML database entries marked as 'unconfirmed'
#           ignore_incomplete: (bool) ignore XML entries missing ANY features
# Output:   (None) Writes to assignement.log and output.csv files for each stage
def main(reference_xml_file, hla_fasta_dir, sample_ID, pass2_metric = "edit_distance",
         pass3_metric = "mismatch_identity", ignore_unconfirmed = False, ignore_incomplete = False,
         generate_query_ref_comp = False):
    samples_file = os.path.join(hla_fasta_dir, str(sample_ID) + "_HLA_haplotypes_CDS.fasta")
    full_sample_file = os.path.join(hla_fasta_dir, str(sample_ID) + "_HLA_haplotypes_gene.fasta")

    run_classification(reference_xml_file, samples_file, full_sample_file, pass2_metric=pass2_metric,
                       pass3_metric=pass3_metric, ignore_unconfirmed=ignore_unconfirmed, ignore_incomplete=ignore_incomplete,
                       generate_query_ref_comp=generate_query_ref_comp)