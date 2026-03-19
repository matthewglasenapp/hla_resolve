# This software is Copyright ©2026. The Regents of the University of California
# ("Regents"). All Rights Reserved.
#
# Licensed under the UC Santa Cruz Noncommercial License (the "License").
# You may not use this file except in compliance with the License.
# A copy of the License is included in this repository as LICENSE.txt.

import pandas as pd
import re
import argparse
from hla_typer import build_g_group_dict

problematic_samples = ["IHW09117"] #["IHW09117", "IHW09125"]
problematic_genes = [] #["IHW09049 DQB1_1", "IHW09049 DQB1_2"]

# Patter to match uncalled alleles
NO_CALL_PAT = re.compile(r"(?:^|\b)(?:no\s*call|na|none)$", re.IGNORECASE)

# Helper function to identify if allele was called
def is_called(x: str) -> bool:
    return isinstance(x, str) and "*" in x and not NO_CALL_PAT.search(x)

# Split multiple options into list
def split_truth_options(a: str):
    if not isinstance(a, str): return []

    parts = re.split(r"[|/]", a)

    return [p.strip() for p in parts]

# Truncate allele to only n fields
def truncate_to_fields(allele: str, n: int):
    # If not called, return none
    if not is_called(allele):
        return None

    head, tail = allele.split("*", 1)
    fields = tail.split(":")

    # Return none if insufficient fields
    if len(fields) < n:
        return None

    # Zero-pad numeric parts to two digits
    padded = []
    for f in fields[:n]:
        match = re.match(r"(\d+)(.*)", f)
        if match:
            padded.append(f"{int(match.group(1)):02d}{match.group(2)}")
        else:
            padded.append(f)

    return f"{head}*{':'.join(padded)}"

def log_f (file, message):
    if file == None:
        return

    file.writelines(message + "\n")

def get_percentages(truth_df, class_df, genes, genes_numbered, tool, platform, g_group_dict=None, quiet = False):
    n_fields = 5 if g_group_dict == None else 4

    # Filter truth set by tool and platform if those columns exist
    if "tool" in class_df.columns:
        if not quiet: print("INFO: Filtering truth set by tool:", tool)
        class_df = class_df[class_df["tool"] == tool]

    if "platform" in class_df.columns:
        if not quiet: print("INFO: Filtering truth set by platform:", platform)
        class_df = class_df[class_df["platform"] == platform]

    # No allele is no call
    for col in genes:
        null_entries = truth_df[col+"_1"].isna()
        truth_df.loc[null_entries, col+"_1"] = "no call"

    # Apply to classifier output too
    for col in genes:
        null_entries = class_df[col+"_1"].isna()
        class_df.loc[null_entries, col+"_1"] = f"no call"

    # Remove HLA- prefix if present in classifier output
    for col in genes_numbered:
        letter = col.split("_")[0]
        has_prefix = ~class_df[col].str.contains(f"HLA-{letter}\*", regex=True, na=True)
        class_df[col][has_prefix] = class_df[col][has_prefix].str.removeprefix('HLA-')

    # Get samples found in input data and reference set
    common_samples = set(truth_df["sample"]) & set(class_df["sample"])

    # Problematic samples
    for problematic_sample in set(problematic_samples) & set(truth_df["sample"]):
        if not quiet: print(f"INFO: Skipping problematic sample {problematic_sample}")
        common_samples.discard(problematic_sample)

    # Only consider entries in both truth set and classified data
    common_truth = truth_df[truth_df["sample"].isin(common_samples)]
    common_class = class_df[class_df["sample"].isin(common_samples)]

    # Create dictionaries to track concordance for every gene at every field level
    concordance = {f"{n}-field":{g:0 for g in genes} for n in range(1,n_fields)}
    fields = {f"{n}-field":{g:0 for g in genes} for n in range(1,n_fields)}

    # Select IHW as truth source
    if "source" in common_truth.columns:
        src_common_truth = common_truth[common_truth["source"] == args.source]
    else:
        src_common_truth = common_truth

    logfile = "truth_check.log"
    if quiet:
        f = open(logfile, "w")
    else:
        f = None

    # This is so we print a single warning message for this case. Otherwise, ouput is spammed
    found_no_g_groups_for_a_sample = False

    # Calculate feild condorances
    for sample in common_samples:
        # Get data for just this sample
        truth = src_common_truth[src_common_truth["sample"] == sample].iloc[0]
        classification = common_class[common_class["sample"] == sample].iloc[0]

        # Check each gene
        for gene in genes:
            gene_name = gene.split("HLA-")[1]
            truth_alleles = [truth.get(f"{gene}_{i}","NA") for i in [1,2] if f"{sample} {gene_name}_{i}" not in problematic_genes]
            test_alleles  = [classification.get(f"{gene}_{i}","NA") for i in [1,2] if f"{sample} {gene_name}_{i}" not in problematic_genes]
            truth_alleles = [a for a in truth_alleles if is_called(a)]
            test_alleles  = [a for a in test_alleles  if is_called(a)]
            if not truth_alleles: continue

            # When checking g groups, convert alleles to g groups
            if g_group_dict != None:
                truth_alleles = list(map(lambda x: g_group_dict.get('HLA-'+x, None), truth_alleles))
                truth_alleles = list(filter(lambda x: x != None, truth_alleles))

                skip_because_not_g_group = False
                for allele in test_alleles:
                    if not allele.endswith('G'):
                        # print("Must supply G groups when in G group checking mode. Problematic allele:", allele)
                        skip_because_not_g_group = True

                if skip_because_not_g_group:
                    if not found_no_g_groups_for_a_sample:
                        print("WARN: Skipping some samples because non-g-group found in g-group comparison.")
                    found_no_g_groups_for_a_sample = True
                    continue

                if truth_alleles == [] or test_alleles == []:
                    continue

            # Check each field
            for n in range(1,n_fields):
                used_tests = set() # Set to avoid reusing same match

                for truth_allele in truth_alleles:
                    log_f(f, f"Entry: {sample} {gene}\nTruth: ")
                    log_f(f, ",".join(truth_alleles))
                    log_f(f, "\nClassification: ")
                    log_f(f, ",".join(test_alleles))
                    log_f(f, "\n")

                    if not is_called(truth_allele):
                        log_f(f, "Truth not called\n")
                        continue

                    if len(truth_allele.split("*",1)[1].split(":")) < n:
                        log_f(f, f"Truth allele not resolved to {n} fields\n")
                        continue

                    # Get every option, truncated to n fields
                    truth_opts = split_truth_options(truth_allele) or [truth_allele]
                    truth_truncs = [truncate_to_fields(a, n) for a in truth_opts if truncate_to_fields(a, n)]
                    if not truth_truncs:
                        log_f(f, f"No truncated truth options available\n")
                        continue

                    fields[f"{n}-field"][gene] += 1

                    if not test_alleles: continue

                    test_truncs = [truncate_to_fields(a, n) for a in test_alleles if truncate_to_fields(a, n)]

                    # Try matching both alleles
                    matched = False
                    for i, test_allele in enumerate(test_alleles):
                        # Truncate test allele to n fields. Skip if non-existant
                        test_trunc = truncate_to_fields(test_allele, n)
                        if not test_trunc:
                            continue

                        # Allow both truth alleles to use same test allele if both sets are homozygous
                        truth_options = list(map(lambda x:truncate_to_fields(x, n), truth_alleles))
                        if truth_options.count(test_trunc) == 2 and test_truncs.count(test_trunc) == 2: # Check by seeing if this allele shows up twice
                            log_f(f, f"Multiple {test_trunc} options found\n{test_alleles} {test_allele}\n")
                            if test_trunc in truth_truncs:
                                log_f(f, f"Match!\n")
                                matched = True
                                break

                        elif i not in used_tests and test_trunc in truth_truncs:
                            log_f(f, f"NOT in used tests!\n")
                            # Otherwise prevent re-use
                            matched = True
                            used_tests.add(i)
                            break

                    if matched:
                        log_f(f, f"Matched at {n} fields\n")
                        concordance[f"{n}-field"][gene] += 1
                    else:
                        log_f(f, f"NOT Matched at {n} fields\n")

    return concordance, fields

if __name__ == "__main__":
    parser = argparse.ArgumentParser(description='Two-Step HLA Allele Identification System')
    parser.add_argument('--hla-xml', required=False, default="../data/hla.xml", help='Input hla.xml reference file')
    parser.add_argument('--reference', required=False, default="../data/ihw_truth.csv", help='Input csv file with base classifications')
    parser.add_argument('--samples', required=False, default="allele_output.csv", help='Input csv file with full classifications')
    parser.add_argument('--source', required=False, default="IHW", help='For reference files with multiple sources, select the source here. Ignored if no source column')
    parser.add_argument('--tool', required=False, default="hla_resolve", help='HLA typing tool name to select from reference file. Ignored if no tool column')
    parser.add_argument('--platform', required=False, default="pacbio", help='Sequencing platform to select from reference file. Ignored if no platform column')
    parser.add_argument("--generate-all", required=False, default=None, help="Prefix to output all percentage concordances as csv. Disabled unless specified. Output files will be named <prefix>_1-field.csv, <prefix>_2-field.csv, etc.")
    parser.add_argument("--g-groups", action='store_true', default=None, help="Print g group concordance")

    args = parser.parse_args()
    xml_file = args.hla_xml #"../data/hla.xml"
    gen_all = args.generate_all

    truth_df = pd.read_csv(args.reference) #"../data/ihw_truth.csv")
    class_df = pd.read_csv(args.samples) #"allele_output.csv")

    # Set default values for g_group
    if args.g_groups == None:
        # Default behaviour for tool = 'hla_la' is compare g groups
        if args.tool == 'hla_la': g_group_mode = True
        else: g_group_mode = False
    # Allow g_group level comparison to be explicitly set
    else: g_group_mode = args.g_groups

    n_fields = 4 if g_group_mode else 5

    if g_group_mode: g_group_dict, p_group_dict, sequence_data = build_g_group_dict(xml_file)

    # Get genes from input data file
    genes_numbered = sorted(list(set(truth_df.columns) & set(class_df.columns)))
    genes_numbered.remove("sample")
    if "source" in genes_numbered:
        genes_numbered.remove("source")
    genes = sorted(list(set(map(lambda x:x.split("_")[0], genes_numbered))))

    concordance, fields = get_percentages(truth_df, class_df, genes, genes_numbered, args.tool, args.platform, g_group_dict=g_group_dict if g_group_mode else None)

    # Print concordance by field, and gene
    for n in range(1,n_fields):
        field = f"{n}-field"
        print(f"{field} concordance")
        for a in genes:
            if fields[field][a] == 0:
                print(f"{a}: NaN ({concordance[field][a]}/{fields[field][a]})")
                continue

            print(f"{a}: {(concordance[field][a]/fields[field][a])*100}% ({concordance[field][a]}/{fields[field][a]})")

    # If set, put data into dataframes
    if gen_all != None:
        print("INFO: ===GENERATING ALL PERCENTAGES===")
        if (type(class_df.get("tool")) == type(None)) or (type(class_df.get("platform")) == type(None)):
            print("ERR: Did not generate information on all tools/platforms, one of those columns was not found.")
            exit(1)

        tools = set(class_df.get("tool"))
        platforms = set(class_df.get("platform"))

        per_field_data_pct = {i:[] for i in range(n_fields)}
        per_field_data_num = {i:[] for i in range(n_fields)}

        # Calculate concordances for each tool/platform pair
        for tool in tools:
            for plat in platforms:
                # Get concordance data silently
                all_concordance, all_fields = get_percentages(truth_df, class_df, genes, genes_numbered, tool, plat,
                                                              g_group_dict=g_group_dict if g_group_mode else None, quiet=True)

                # Calculate concordance percentage and insert into apropriate dataframe
                for field in range (1,n_fields):
                    concordances = all_concordance[f"{field}-field"]

                    # Calculate concordance percentages per gene
                    concordance_pcts = []
                    concordance_nums = []
                    for gene in genes:
                        if all_fields[f"{field}-field"][gene] == 0:
                            print(f"WARNING: No entries for {tool} {plat} {field}-field {gene}, setting concordance to NaN")
                            concordance_pcts.append(float('nan'))
                            continue

                        pct = (concordances[gene]/all_fields[f"{field}-field"][gene])*100
                        pct = round(pct)
                        concordance_pcts.append(pct)
                        concordance_nums.append(concordances[gene])

                    # Create new row contianing calculated concordance percentages
                    per_field_data_pct[field].append([tool, plat] + concordance_pcts)
                    per_field_data_num[field].append([tool, plat] + concordance_nums)

        # Output each dataframe to csv
        for field in range(1,n_fields):
            out_file_pct = f"{gen_all}_percent_{field}-field.csv"
            out_file_num = f"{gen_all}_count_{field}-field.csv"
            if g_group_mode:
                out_file_pct = f"{gen_all}_percent_{field}-field_g-group.csv"
                out_file_num = f"{gen_all}_count_{field}-field_g-group.csv"

            print(f"INFO: Writing all percentages for {field}-field to {out_file_pct}")
            field_pct_df = pd.DataFrame(per_field_data_pct[field], columns=["Tool", "Platform"] + genes)
            field_pct_df = field_pct_df.sort_values(by=["Tool", "Platform"])
            field_pct_df.to_csv(out_file_pct, index=False)

            print(f"INFO: Writing all concordance counts for {field}-field to {out_file_num}")
            field_num_df = pd.DataFrame(per_field_data_num[field], columns=["Tool", "Platform"] + genes)
            field_num_df = field_num_df.sort_values(by=["Tool", "Platform"])
            field_num_df.to_csv(out_file_num, index=False)
