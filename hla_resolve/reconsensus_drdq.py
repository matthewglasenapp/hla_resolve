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

from . import config

DRDQ_GENES = ("HLA-DRB1", "HLA-DQA1", "HLA-DQB1")

_WILDCARD_EQUALITIES = [("N", b) for b in "ACGT"] + [(b, "N") for b in "ACGT"]
_BIG = 10 ** 9

# A het DR/DQ read whose competitive primary lands on the OTHER haplotype's
# scaffold with MAPQ at or above this is treated as mis-phased (decisively
# wrong-haplotype, e.g. reads HiPhase tagged to the wrong HP) and dropped before
# that slot's consensus. Reads on their own scaffold, or with low MAPQ over
# regions identical between the two alleles, are retained so Force HP keeps its
# no-starvation behavior.
_PHASING_DROP_MAPQ = 30


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


def _utr_parts(sequence_data, allele):
    # (5' UTR, 3' UTR) of an allele. Returns ("", "") unless the record carries
    # both, so a partially-annotated entry contributes no UTR to the shared span
    # rather than an ambiguously-oriented one.
    feats = sequence_data.get(allele)
    if not feats:
        return "", ""
    utrs = sorted(feats.get("UTR", []))
    if len(utrs) != 2:
        return "", ""
    return utrs[0][1], utrs[-1][1]


def _full_with_exon_intervals(sequence_data, allele):
    # Full (UTR+Exon+Intron) sequence plus the [start, stop) offsets of each
    # exon within it. Sorted identically to _full_sequence so offsets line up.
    feats = sequence_data.get(allele)
    if not feats:
        return "", []
    tagged = []
    for feature_type, entries in feats.items():
        if feature_type in ("UTR", "Exon", "Intron"):
            for pos, seq in entries:
                tagged.append((pos, seq, feature_type))
    tagged.sort(key=lambda x: (x[0], x[1]))
    parts = []
    intervals = []
    offset = 0
    for _, seq, feature_type in tagged:
        if feature_type == "Exon":
            intervals.append((offset, offset + len(seq)))
        parts.append(seq)
        offset += len(seq)
    return "".join(parts), intervals


def _project_cds(consensus, allele, sequence_data):
    # Extract the CDS (exon concatenation) of the re-consensus sequence by
    # aligning the consensus to the chosen allele's full sequence and pulling
    # the consensus bases that fall within each exon interval.
    full, intervals = _full_with_exon_intervals(sequence_data, allele)
    if not consensus or not full or not intervals:
        return ""
    result = edlib.align(consensus, full, mode="NW", task="path",
                         additionalEqualities=_WILDCARD_EQUALITIES)
    cigar = result.get("cigar") or ""
    if not cigar:
        return ""
    # ref_q[r] = consensus index aligned at full-sequence position r
    ref_q = [0] * (len(full) + 1)
    r = q = 0
    for count, op in re.findall(r"(\d+)([=XID])", cigar):
        count = int(count)
        if op in ("=", "X"):
            for _ in range(count):
                ref_q[r] = q
                r += 1
                q += 1
        elif op == "D":  # base in full only (gap in consensus)
            for _ in range(count):
                ref_q[r] = q
                r += 1
        else:            # "I": base in consensus only
            q += count
    ref_q[len(full)] = q
    return "".join(consensus[ref_q[s]:ref_q[e]] for s, e in intervals)


def _iter_fasta_records(path):
    header, seq = None, []
    with open(path) as fh:
        for line in fh:
            if line.startswith(">"):
                if header is not None:
                    yield header, "".join(seq)
                header, seq = line[1:].strip(), []
            else:
                seq.append(line.strip())
    if header is not None:
        yield header, "".join(seq)


def _rewrite_fasta(path, replacements):
    # Overwrite the sequence of any record whose id is in replacements; leave
    # all other records untouched. Wrap at 60 cols to match SeqIO output.
    if not path or not os.path.isfile(path) or not replacements:
        return
    records = []
    for header, seq in _iter_fasta_records(path):
        rid = header.split()[0]
        if rid in replacements:
            seq = replacements[rid]
        records.append((rid, seq))
    with open(path, "w") as fh:
        for rid, seq in records:
            fh.write(f">{rid}\n")
            for i in range(0, len(seq), 60):
                fh.write(seq[i:i + 60] + "\n")


def _edit_distance(query, ref):
    if not query or not ref:
        return _BIG
    a, b = (query, ref) if len(query) <= len(ref) else (ref, query)
    return edlib.align(a, b, mode="HW", task="distance",
                       additionalEqualities=_WILDCARD_EQUALITIES)["editDistance"]


def _penalized_distance(query, ref):
    # Edit distance with N treated as a real mismatch (NO wildcard) so an
    # ambiguous/blurry consensus (Ns at unresolved het positions) scores as a
    # WORSE match than a clean one. Used only by the un-pooling accept gate: the
    # pooled consensus of a true het carries Ns exactly at the positions that
    # distinguish the two 4th fields, and must not win them for free.
    if not query or not ref:
        return _BIG
    a, b = (query, ref) if len(query) <= len(ref) else (ref, query)
    return edlib.align(a, b, mode="HW", task="distance")["editDistance"]


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


def _has_pass_phased_het_genotype(phased_vcf):
    # True only if the gene VCF carries at least one PASS-filter, phased,
    # heterozygous genotype. This is the trigger for un-pooling: a same-lineage
    # call is split into two haplotypes only when there is a real phased het to
    # split on. Unphased hets (HiPhase could not phase them) and non-PASS records
    # (DR-region reconstruction artifacts) do not qualify, so those loci stay
    # pooled and are called homozygous.
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
                if parts[6] not in ("PASS", "."):
                    continue
                gt = parts[9].split(":")[0]
                if "|" not in gt:
                    continue
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
    _run(f"{config.rammap} -a -x map-hifi {scaffold_fa} {reads_fq} | "
         f"samtools sort -o {aln} - && samtools index {aln}")
    _run(f"samtools consensus -f fasta {aln} > {cons}")
    return _read_fasta_seq(cons)


def _consensus_self_sort(scaffold_seq_1, scaffold_seq_2, reads_fq, workdir):
    ref = os.path.join(workdir, "two_scaffold.fa")
    with open(ref, "w") as fh:
        fh.write(f">sc1\n{scaffold_seq_1}\n>sc2\n{scaffold_seq_2}\n")
    aln = os.path.join(workdir, "sort.aln.bam")
    prim = os.path.join(workdir, "sort.prim.bam")
    _run(f"{config.rammap} -a -x map-hifi {ref} {reads_fq} | "
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
    _run(f"{config.rammap} -a -x map-hifi {ref_fa} {reads_fq} | "
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


def _retained_reads_fastq(prim_bam, own_contig, out_fq, min_mapq=_PHASING_DROP_MAPQ):
    # From the competitive (two-scaffold) primary BAM for one HP tag, keep every
    # read EXCEPT those whose primary aligns to the OTHER haplotype's scaffold
    # with MAPQ >= min_mapq (decisively mis-phased). Reads on their own contig,
    # unmapped reads, or low-MAPQ reads over identical regions are retained.
    # Returns True if the retained fastq is non-empty.
    expr = f'(rname == "{own_contig}") || (mapq < {min_mapq})'
    _run(f"samtools view -b -e '{expr}' {prim_bam} | samtools fastq - > {out_fq}")
    return os.path.isfile(out_fq) and os.path.getsize(out_fq) > 0


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

    if len(tied) > 1:
        # The exon+intron core could not separate these candidates. In most
        # such sets it never can: the tied alleles are byte-identical across
        # every exon and intron, and the only sequence that distinguishes them
        # is UTR.
        #
        # Comparing full sequences would be unfair, because the tied records
        # carry different amounts of UTR and an allele could win simply by
        # having more of it annotated. So each candidate is trimmed to the UTR
        # span that ALL of them share -- the innermost min(5' UTR) and
        # min(3' UTR) bases, UTRs extending outward from the coding sequence --
        # and every candidate is scored over an identical window.
        #
        # This runs only on sets the core left tied, so it cannot disturb a
        # call the core already decided. Where it still cannot separate them,
        # the lowest fourth field decides as before.
        utrs = {name: _utr_parts(sequence_data, name) for name in tied}
        keep5 = min(len(utrs[name][0]) for name in tied)
        keep3 = min(len(utrs[name][1]) for name in tied)
        if keep5 or keep3:
            def _shared_span(name):
                utr5, utr3 = utrs[name]
                return ((utr5[-keep5:] if keep5 else "")
                        + core_cache[name]
                        + (utr3[:keep3] if keep3 else ""))

            rescored = [(_edit_distance(consensus, _shared_span(name)), name)
                        for name in tied]
            best_shared = min(d for d, _ in rescored)
            narrowed = sorted(name for d, name in rescored if d == best_shared)
            if len(narrowed) < len(tied):
                tied = narrowed

    best = min(tied, key=_fourth_field_num)
    return (best, min_dist, tied)


def _pooled_alternative(bam, region, scaffold_seq, lineage, gene,
                        sequence_data, core_cache, workdir):
    # Build the single pooled (homozygous) consensus on the shared scaffold from
    # ALL reads in the region — the alternative to splitting by HP tag — and
    # return (consensus, best-in-lineage tuple). Used only by the same-lineage
    # un-pooling accept gate.
    fq = os.path.join(workdir, "pooled_alt.fq")
    if not _extract_reads_fastq(bam, region, fq):
        return "", None
    fa = os.path.join(workdir, "pooled_alt.scaffold.fa")
    with open(fa, "w") as fh:
        fh.write(f">scaffold\n{scaffold_seq}\n")
    cons = _consensus_on_scaffold(fa, fq, workdir, "pooled_alt")
    ref = _best_in_lineage(cons, lineage, gene, sequence_data, core_cache)
    return cons, ref


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
    if changed:
        outcome = f"refined to {allele}"
    elif _num_fields(allele) < 4:
        outcome = f"not refined (guess {allele} has no 4th field to resolve); kept prior call"
    else:
        outcome = f"not refined (re-consensus not accepted); kept prior call {allele}"
    if dist is None:
        logfile.writelines(f"For {name}, re-consensus: {outcome} (no full-sequence match measured)\n")
    else:
        logfile.writelines(f"For {name}, re-consensus: {outcome} dist {dist} len {mlen} "
                           f"id {seq_id} mismatch {mm_id} using edit_distance "
                           f"(Exon+Intron core; ties broken on shared UTR)\n")
    if isinstance(tie, (list, tuple)) and len(tie) > 1:
        logfile.writelines(f"Equidistant (re-consensus) for {name}: {', '.join(tie)}\n")


def _refine_one_gene(sample_ID, gene, name1, name2, results, query_seqs,
                     sequence_data, ctx, core_cache, logfile=None, overrides=None):
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
        # Pool (call homozygous) only when there is no real phased het to split
        # on. A same-lineage locus that carries a PASS-phased het is NOT pooled
        # up front; it goes through the HP split and the un-pooling accept gate
        # below, which decides whether the two haplotypes truly differ at the
        # 4th field or the locus is homozygous with artifact hets.
        pooled = not _has_pass_phased_het_genotype(phased_vcf)

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
                if same_lineage:
                    # Un-pooled same-lineage het: both haplotypes share ONE
                    # 4th-field scaffold so neither per-HP consensus is biased
                    # toward a per-slot seed (the seed-echo that would otherwise
                    # decide the 4th field). The competitive two-scaffold step
                    # below then drops nothing (identical contigs), so the split
                    # is purely by HP tag.
                    shared = scaffold_seq_1 if _num_fields(best_guess_1) >= 4 else scaffold_seq_2
                    scaffold_seq_1 = scaffold_seq_2 = shared
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
                slot_scaffold_seq = {1: scaffold_seq_1, 2: scaffold_seq_2}
                for slot in (1, 2):
                    if _num_fields(slot_guess[slot]) < 4:
                        continue
                    hp = tag_by_slot[slot]
                    if prim[hp] is None:
                        continue
                    # Force HP: build this slot's consensus from the assigned HP tag's
                    # reads mapped to the slot's scaffold ALONE. The two-scaffold step
                    # above only decides which tag pairs with which allele; it must not
                    # also re-drop reads whose primary lands on the other scaffold over
                    # regions identical between the two alleles, which starved a slot
                    # under the competitive consensus and produced runs of N.
                    #
                    # Margin gate: before building the consensus, drop only reads whose
                    # competitive primary decisively prefers the OTHER slot's scaffold
                    # (MAPQ >= _PHASING_DROP_MAPQ) — mis-phased reads HiPhase tagged to
                    # the wrong HP. Reads on their own scaffold or ambiguous over
                    # identical regions (low MAPQ) are retained, so no slot is starved.
                    slot_scaffold_fa = os.path.join(workdir, f"slot{slot}.scaffold.fa")
                    with open(slot_scaffold_fa, "w") as fh:
                        fh.write(f">{slot_contig[slot]}\n{slot_scaffold_seq[slot]}\n")
                    gated_fq = os.path.join(workdir, f"slot{slot}.gated.fq")
                    if _retained_reads_fastq(prim[hp], slot_contig[slot], gated_fq):
                        reads_fq = gated_fq
                    else:
                        reads_fq = os.path.join(workdir, f"hp{hp}.fq")
                    cons = _consensus_on_scaffold(slot_scaffold_fa, reads_fq, workdir, f"slot{slot}")
                    cons_of[slot] = cons
                    refined[slot] = _best_in_lineage(cons, slot_lineage[slot], gene, sequence_data, core_cache)

            present = [i for i in (1, 2) if refined.get(i)]
            if not present:
                return

            name_by_idx = {1: name1, 2: name2}
            guess_by_idx = {1: best_guess_1, 2: best_guess_2}

            # assign[slot] = (allele, consensus, tie_set) for the slots to overwrite
            assign = {}
            refined_q = sum(refined[i][1] for i in present)
            if pooled and same_lineage:
                # Genuinely homozygous locus (no phased het): one pooled
                # consensus assigned to both slots.
                winner = min(present, key=lambda i: refined[i][1])
                allele, _, tie = refined[winner]
                assign[1] = assign[2] = (allele, cons_of[winner], tie)
            elif same_lineage and len(present) == 2:
                # Un-pooled same-lineage het: accept the HP split into two
                # distinct 4th fields only if the two per-haplotype consensuses
                # explain the reads better (N-penalized) than a single pooled
                # homozygous consensus. The N penalty stops a blurry pooled blend
                # from matching an allele for free at the very positions that
                # distinguish the two 4th fields. On a true homozygote the pooled
                # consensus (full depth) is cleaner than either half-depth split,
                # so the split loses and the call stays homozygous.
                shared_seq = _full_sequence(
                    sequence_data,
                    best_guess_1 if _num_fields(best_guess_1) >= 4 else best_guess_2)
                pooled_cons, pooled_ref = _pooled_alternative(
                    bam, region, shared_seq, lineage_1, gene, sequence_data,
                    core_cache, workdir)
                split_pen = sum(
                    _penalized_distance(cons_of[i], _core_sequence(sequence_data, refined[i][0]))
                    for i in present)
                pooled_pen = (2 * _penalized_distance(pooled_cons, _core_sequence(sequence_data, pooled_ref[0]))
                              if pooled_ref and pooled_cons else _BIG)
                if split_pen < pooled_pen:
                    for i in present:
                        assign[i] = (refined[i][0], cons_of[i], refined[i][2])
                    if logfile is not None:
                        logfile.writelines(
                            f"For {sample_ID} {gene}, un-pooled het accepted: split "
                            f"residual {split_pen} < pooled {pooled_pen}\n")
                else:
                    allele, _, tie = pooled_ref
                    assign[1] = assign[2] = (allele, pooled_cons, tie)
                    refined_q = 2 * pooled_ref[1]
                    if logfile is not None:
                        logfile.writelines(
                            f"For {sample_ID} {gene}, un-pooled het rejected (kept "
                            f"homozygous): split residual {split_pen} >= pooled {pooled_pen}\n")
            else:
                for i in present:
                    assign[i] = (refined[i][0], cons_of[i], refined[i][2])

            original_q = 0
            for i in present:
                query = query_seqs.get(name_by_idx[i])
                query = str(query) if query is not None else ""
                original_q += _edit_distance(query, _core_sequence(sequence_data, guess_by_idx[i]))

            if refined_q <= original_q:
                for i, (allele, cons, tie) in assign.items():
                    dist, mlen, _, seq_id, mm_id = _match_metrics(cons, _core_sequence(sequence_data, allele))
                    used_tb = len(tie) > 1
                    name = name_by_idx[i]
                    results[name] = (allele, dist, mlen, seq_id, mm_id, used_tb, tie)
                    changed.add(name)
                    # Override the deposited haplotype with the re-consensus
                    # sequence: full gene = consensus, CDS = exon projection.
                    if overrides is not None and cons:
                        overrides[name] = {
                            "gene": cons,
                            "cds": _project_cds(cons, allele, sequence_data),
                        }
                        if query_seqs is not None:
                            query_seqs[name] = cons
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
    overrides = {}
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
                            results, query_seqs, sequence_data, ctx, core_cache,
                            logfile=logfile, overrides=overrides)
        except Exception:
            continue

    # Deposit the re-consensus sequence into the final haplotype FASTAs for the
    # refined DR/DQ haplotypes only; declined haplotypes keep vcf2fasta output.
    if overrides:
        gene_reps = {name: ov["gene"] for name, ov in overrides.items() if ov.get("gene")}
        cds_reps = {name: ov["cds"] for name, ov in overrides.items() if ov.get("cds")}
        _rewrite_fasta(ctx.get("gene_fasta"), gene_reps)
        _rewrite_fasta(ctx.get("cds_fasta"), cds_reps)
