# This software is Copyright ©2026. The Regents of the University of California
# ("Regents"). All Rights Reserved.
#
# See LICENSE.txt for license details.

# 4th-field read re-consensus refinement for HLA-DRB1/DQA1/DQB1 only.
# Gated by config.reconsensus_drdq; wired in for PacBio via resolve_alleles_pipeline.

import os
import re
import gzip
import shutil
import tempfile
import subprocess

import edlib

DRDQ_GENES = ("HLA-DRB1", "HLA-DQA1", "HLA-DQB1")

_WILDCARD_EQUALITIES = [("N", b) for b in "ACGT"] + [(b, "N") for b in "ACGT"]
_BIG = 10 ** 9


def _run(cmd):
    subprocess.run(cmd, shell=True, check=True, capture_output=True, text=True)


def _trunc_to_3_fields(allele):
    fields = allele.split(":")
    return ":".join(fields[:3]) if len(fields) >= 4 else allele


def _num_fields(allele):
    if "*" not in allele:
        return 0
    return len(allele.split("*", 1)[1].split(":"))


def _fourth_field_num(allele):
    if "*" not in allele:
        return _BIG
    fields = allele.split("*", 1)[1].split(":")
    if len(fields) < 4:
        return _BIG
    match = re.match(r"(\d+)", fields[3])
    return int(match.group(1)) if match else _BIG


def _core_sequence(sequence_data, allele):
    segments = []
    feats = sequence_data.get(allele)
    if not feats:
        return ""
    for feature_type, entries in feats.items():
        if feature_type in ("Exon", "Intron"):
            segments += entries
    segments.sort()
    return "".join(seq for _, seq in segments)


def _full_sequence(sequence_data, allele):
    segments = []
    feats = sequence_data.get(allele)
    if not feats:
        return ""
    for feature_type, entries in feats.items():
        if feature_type in ("UTR", "Exon", "Intron"):
            segments += entries
    segments.sort()
    return "".join(seq for _, seq in segments)


def _edit_distance(query, ref):
    if not query or not ref:
        return _BIG
    a, b = (query, ref) if len(query) <= len(ref) else (ref, query)
    return edlib.align(a, b, mode="HW", task="distance",
                       additionalEqualities=_WILDCARD_EQUALITIES)["editDistance"]


def _read_fasta_seq(path):
    if not os.path.isfile(path):
        return ""
    seq = []
    with open(path) as fh:
        for line in fh:
            if not line.startswith(">"):
                seq.append(line.strip())
    return "".join(seq).upper()


def _gene_region(gene, phased_vcf, gene_dict, pad=3000):
    coords = gene_dict.get(gene) if gene_dict else None
    chrom = "chr6"
    positions = []
    if phased_vcf and os.path.isfile(phased_vcf):
        try:
            opener = gzip.open if phased_vcf.endswith(".gz") else open
            with opener(phased_vcf, "rt") as fh:
                for line in fh:
                    if line.startswith("#"):
                        continue
                    parts = line.split("\t")
                    if len(parts) < 2:
                        continue
                    chrom = parts[0]
                    try:
                        positions.append(int(parts[1]))
                    except ValueError:
                        continue
        except OSError:
            positions = []
    if positions:
        start, stop = min(positions) - pad, max(positions) + pad
    elif coords:
        start, stop = coords[0] - pad, coords[1] + pad
    else:
        return None
    return f"{chrom}:{max(1, start)}-{stop}"


def _has_het_genotype(phased_vcf):
    if not phased_vcf or not os.path.isfile(phased_vcf):
        return False
    opener = gzip.open if phased_vcf.endswith(".gz") else open
    try:
        with opener(phased_vcf, "rt") as fh:
            for line in fh:
                if line.startswith("#"):
                    continue
                parts = line.rstrip("\n").split("\t")
                if len(parts) < 10:
                    continue
                gt = parts[9].split(":")[0]
                present = {a for a in re.split(r"[|/]", gt) if a not in (".", "")}
                if len(present) > 1:
                    return True
    except OSError:
        return False
    return False


def _extract_reads_fastq(bam, region, out_fq, hp=None):
    hp_flag = f"-d HP:{hp} " if hp is not None else ""
    _run(f"samtools view -b -F 0x900 {hp_flag}{bam} {region} | samtools fastq - > {out_fq}")
    return os.path.isfile(out_fq) and os.path.getsize(out_fq) > 0


def _consensus_on_scaffold(scaffold_fa, reads_fq, workdir, tag):
    aln = os.path.join(workdir, f"{tag}.aln.bam")
    cons = os.path.join(workdir, f"{tag}.cons.fa")
    _run(f"minimap2 -a -x map-hifi {scaffold_fa} {reads_fq} | "
         f"samtools sort -o {aln} - && samtools index {aln}")
    _run(f"samtools consensus -f fasta {aln} > {cons}")
    return _read_fasta_seq(cons)


def _consensus_self_sort(scaffold_seq_1, scaffold_seq_2, reads_fq, workdir):
    ref = os.path.join(workdir, "two_scaffold.fa")
    with open(ref, "w") as fh:
        fh.write(f">sc1\n{scaffold_seq_1}\n>sc2\n{scaffold_seq_2}\n")
    aln = os.path.join(workdir, "sort.aln.bam")
    prim = os.path.join(workdir, "sort.prim.bam")
    _run(f"minimap2 -a -x map-hifi {ref} {reads_fq} | "
         f"samtools sort -o {aln} - && samtools index {aln}")
    _run(f"samtools view -b -F 0x900 {aln} > {prim} && samtools index {prim}")
    cons1 = os.path.join(workdir, "sort.cons1.fa")
    cons2 = os.path.join(workdir, "sort.cons2.fa")
    _run(f"samtools consensus -f fasta -r sc1 {prim} > {cons1}")
    _run(f"samtools consensus -f fasta -r sc2 {prim} > {cons2}")
    return _read_fasta_seq(cons1), _read_fasta_seq(cons2)


def _align_and_count(reads_fq, ref_fa, workdir, tag):
    aln = os.path.join(workdir, f"{tag}.aln.bam")
    prim = os.path.join(workdir, f"{tag}.prim.bam")
    _run(f"minimap2 -a -x map-hifi {ref_fa} {reads_fq} | "
         f"samtools sort -o {aln} - && samtools index {aln}")
    _run(f"samtools view -b -F 0x900 {aln} > {prim} && samtools index {prim}")
    result = subprocess.run(f"samtools idxstats {prim}", shell=True, check=True,
                            capture_output=True, text=True)
    counts = {}
    for line in result.stdout.splitlines():
        parts = line.split("\t")
        if len(parts) >= 3 and parts[0] != "*":
            counts[parts[0]] = int(parts[2])
    return prim, counts


def _consensus_region(prim_bam, contig, workdir, tag):
    cons = os.path.join(workdir, f"{tag}.cons.fa")
    _run(f"samtools consensus -f fasta -r {contig} {prim_bam} > {cons}")
    return _read_fasta_seq(cons)


def _best_in_lineage(consensus, lineage, gene, sequence_data, core_cache):
    if not consensus:
        return None
    cores = []
    for name in sequence_data.keys():
        if not name.startswith(gene + "*") or _num_fields(name) < 4:
            continue
        if _trunc_to_3_fields(name) != lineage:
            continue
        if name not in core_cache:
            core_cache[name] = _core_sequence(sequence_data, name)
        cores.append((name, core_cache[name]))
    if not cores:
        return None
    scored = [(_edit_distance(consensus, core), name) for name, core in cores]
    min_dist = min(d for d, _ in scored)
    tied = sorted(name for d, name in scored if d == min_dist)
    best = min(tied, key=_fourth_field_num)
    return (best, min_dist, tied)


def _match_metrics(consensus, core):
    # (raw_distance, match_len, mismatch_len, seq_identity, mismatch_identity)
    if not consensus or not core:
        return (_BIG, 0, 0, 0.0, 0.0)
    a, b = (consensus, core) if len(consensus) <= len(core) else (core, consensus)
    result = edlib.align(a, b, mode="HW", task="path",
                         additionalEqualities=_WILDCARD_EQUALITIES)
    cigar = result["cigar"] or ""
    match = sum(int(x) for x in re.findall(r"(\d+)=", cigar))
    mism = sum(int(x) for x in re.findall(r"(\d+)X", cigar))
    dist = result["editDistance"]
    seq_id = 1 - dist / match if match else 0.0
    mm_id = match / (match + mism) if (match + mism) else 0.0
    return (dist, match, mism, seq_id, mm_id)


def _log_hap(logfile, name, tup, changed):
    if logfile is None or not tup:
        return
    allele, dist, mlen, seq_id, mm_id = tup[0], tup[1], tup[2], tup[3], tup[4]
    tie = tup[6] if len(tup) > 6 else [allele]
    verb = "refined to" if changed else "kept"
    logfile.writelines(f"For {name}, re-consensus {verb} {allele} dist {dist} len {mlen} "
                       f"id {seq_id} mismatch {mm_id} using edit_distance (Exon+Intron core)\n")
    if isinstance(tie, (list, tuple)) and len(tie) > 1:
        logfile.writelines(f"Equidistant (re-consensus) for {name}: {', '.join(tie)}\n")


def _refine_one_gene(sample_ID, gene, name1, name2, results, query_seqs,
                     sequence_data, ctx, core_cache, logfile=None):
    changed = set()
    try:
        best_guess_1 = results[name1][0]
        best_guess_2 = results[name2][0]
        if not best_guess_1 or not best_guess_2:
            return

        lineage_1 = _trunc_to_3_fields(best_guess_1)
        lineage_2 = _trunc_to_3_fields(best_guess_2)

        phased_vcf = ctx.get("gene_vcfs", {}).get(gene)
        region = _gene_region(gene, phased_vcf, ctx.get("gene_dict"))
        if region is None:
            return

        same_lineage = lineage_1 == lineage_2
        pooled = same_lineage or not _has_het_genotype(phased_vcf)

        bam = ctx["bam"]
        mode = ctx.get("mode", "hp_tag")

        workdir = tempfile.mkdtemp(prefix=f"reconsensus_{sample_ID}_{gene.replace('HLA-', '')}_")
        try:
            refined = {}
            cons_of = {}
            if pooled:
                reads_fq = os.path.join(workdir, "pooled.fq")
                if not _extract_reads_fastq(bam, region, reads_fq):
                    return
                scaffold_allele = best_guess_1 if _num_fields(best_guess_1) >= 4 else best_guess_2
                scaffold_seq = _full_sequence(sequence_data, scaffold_allele)
                if not scaffold_seq:
                    return
                scaffold_fa = os.path.join(workdir, "scaffold.fa")
                with open(scaffold_fa, "w") as fh:
                    fh.write(f">scaffold\n{scaffold_seq}\n")
                cons = _consensus_on_scaffold(scaffold_fa, reads_fq, workdir, "pooled")
                for slot, guess, lineage in ((1, best_guess_1, lineage_1), (2, best_guess_2, lineage_2)):
                    if _num_fields(guess) >= 4:
                        cons_of[slot] = cons
                        refined[slot] = _best_in_lineage(cons, lineage, gene, sequence_data, core_cache)
            elif mode == "self_sort":
                reads_fq = os.path.join(workdir, "all.fq")
                if not _extract_reads_fastq(bam, region, reads_fq):
                    return
                scaffold_seq_1 = _full_sequence(sequence_data, best_guess_1)
                scaffold_seq_2 = _full_sequence(sequence_data, best_guess_2)
                if not scaffold_seq_1 or not scaffold_seq_2:
                    return
                cons1, cons2 = _consensus_self_sort(scaffold_seq_1, scaffold_seq_2, reads_fq, workdir)
                for slot, guess, lineage, cons in ((1, best_guess_1, lineage_1, cons1), (2, best_guess_2, lineage_2, cons2)):
                    if _num_fields(guess) >= 4:
                        cons_of[slot] = cons
                        refined[slot] = _best_in_lineage(cons, lineage, gene, sequence_data, core_cache)
            else:
                # hp_tag: assign each HP tag to a scaffold by primary-alignment fraction,
                # then consensus each tag on its assigned scaffold for the 4th-field call.
                scaffold_seq_1 = _full_sequence(sequence_data, best_guess_1)
                scaffold_seq_2 = _full_sequence(sequence_data, best_guess_2)
                if not scaffold_seq_1 or not scaffold_seq_2:
                    return
                ref_fa = os.path.join(workdir, "two_scaffold.fa")
                with open(ref_fa, "w") as fh:
                    fh.write(f">sc1\n{scaffold_seq_1}\n>sc2\n{scaffold_seq_2}\n")
                prim = {}
                frac = {}
                for hp in (1, 2):
                    fq = os.path.join(workdir, f"hp{hp}.fq")
                    if not _extract_reads_fastq(bam, region, fq, hp=hp):
                        prim[hp] = None
                        frac[(hp, "sc1")] = frac[(hp, "sc2")] = 0.0
                        continue
                    prim_bam, counts = _align_and_count(fq, ref_fa, workdir, f"hp{hp}")
                    prim[hp] = prim_bam
                    total = counts.get("sc1", 0) + counts.get("sc2", 0)
                    frac[(hp, "sc1")] = counts.get("sc1", 0) / total if total else 0.0
                    frac[(hp, "sc2")] = counts.get("sc2", 0) / total if total else 0.0
                straight = frac[(1, "sc1")] + frac[(2, "sc2")]
                swap = frac[(1, "sc2")] + frac[(2, "sc1")]
                # tag assigned to lineage 1 (sc1) and lineage 2 (sc2)
                tag_by_slot = {1: 1, 2: 2} if straight >= swap else {1: 2, 2: 1}
                slot_contig = {1: "sc1", 2: "sc2"}
                slot_lineage = {1: lineage_1, 2: lineage_2}
                slot_guess = {1: best_guess_1, 2: best_guess_2}
                for slot in (1, 2):
                    if _num_fields(slot_guess[slot]) < 4:
                        continue
                    hp = tag_by_slot[slot]
                    if prim[hp] is None:
                        continue
                    cons = _consensus_region(prim[hp], slot_contig[slot], workdir, f"slot{slot}")
                    cons_of[slot] = cons
                    refined[slot] = _best_in_lineage(cons, slot_lineage[slot], gene, sequence_data, core_cache)

            present = [i for i in (1, 2) if refined.get(i)]
            if not present:
                return

            name_by_idx = {1: name1, 2: name2}
            guess_by_idx = {1: best_guess_1, 2: best_guess_2}

            # assign[slot] = (allele, consensus, tie_set) for the slots to overwrite
            assign = {}
            if same_lineage:
                winner = min(present, key=lambda i: refined[i][1])
                allele, _, tie = refined[winner]
                assign[1] = assign[2] = (allele, cons_of[winner], tie)
            else:
                for i in present:
                    assign[i] = (refined[i][0], cons_of[i], refined[i][2])

            refined_q = sum(refined[i][1] for i in present)
            original_q = 0
            for i in present:
                query = query_seqs.get(name_by_idx[i])
                query = str(query) if query is not None else ""
                original_q += _edit_distance(query, _core_sequence(sequence_data, guess_by_idx[i]))

            if refined_q <= original_q:
                for i, (allele, cons, tie) in assign.items():
                    dist, mlen, _, seq_id, mm_id = _match_metrics(cons, _core_sequence(sequence_data, allele))
                    used_tb = len(tie) > 1
                    results[name_by_idx[i]] = (allele, dist, mlen, seq_id, mm_id, used_tb, tie)
                    changed.add(name_by_idx[i])
        finally:
            shutil.rmtree(workdir, ignore_errors=True)
    finally:
        _log_hap(logfile, name1, results.get(name1), name1 in changed)
        _log_hap(logfile, name2, results.get(name2), name2 in changed)


def refine_drdq(results, query_seqs, sequence_data, ctx, logfile=None):
    # Refine the DRDQ result tuples (allele + metrics + tie set) in place.
    if not ctx or not ctx.get("bam") or not os.path.isfile(ctx["bam"]):
        return

    core_cache = {}
    groups = {}
    for name in results.keys():
        for gene in DRDQ_GENES:
            token = "_" + gene + "_"
            if token in name:
                sample_ID = name.split(token)[0]
                groups.setdefault((sample_ID, gene), []).append(name)
                break

    for (sample_ID, gene), names in groups.items():
        if len(names) != 2:
            continue
        names.sort(key=lambda n: n.rsplit("_", 1)[-1])
        try:
            _refine_one_gene(sample_ID, gene, names[0], names[1],
                            results, query_seqs, sequence_data, ctx, core_cache, logfile=logfile)
        except Exception:
            continue
