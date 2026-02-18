#!/usr/bin/env python3
"""
Diagnostic script to instrument v2f.getSequences and trace
exactly what happens with pbsv variant records during execution.

Run from the directory containing the VCF, GFF, and reference files:
    python3 debug_getSequences.py \
        --vcf IHW09364_HLA-DRB1_PASS_phased.vcf.gz \
        --fasta /path/to/reference.fa \
        --gff /path/to/hla_drb1.gff3
"""

import sys
import os
import argparse
import pysam
import collections

# Add the vcf2fasta directory to the path so we can import v2f
# Adjust this path as needed
SCRIPT_DIR = os.path.dirname(os.path.abspath(__file__))
V2F_DIR = os.path.join(SCRIPT_DIR, "data", "vcf2fasta")
sys.path.insert(0, V2F_DIR)

import v2f.functions as v2f

def main():
    parser = argparse.ArgumentParser(description="Debug getSequences for pbsv variants")
    parser.add_argument("--vcf", required=True, help="Input VCF file")
    parser.add_argument("--fasta", required=True, help="Reference FASTA")
    parser.add_argument("--gff", required=True, help="GFF3 file")
    parser.add_argument("--feat", default="gene", help="Feature type (default: gene)")
    args = parser.parse_args()

    # ===== STEP 0: Check which code is actually running =====
    print("=" * 70)
    print("STEP 0: Verify which code is being executed")
    print("=" * 70)
    func_file = v2f.__file__
    print(f"v2f module loaded from: {func_file}")

    # Check if it's a .pyc file (cached bytecode)
    if func_file.endswith('.pyc'):
        print("WARNING: Running from compiled .pyc bytecache!")
        print("The .py source may differ from what's actually executing.")
        print("Consider deleting __pycache__ and re-running.")

    # Check for __pycache__ alongside the source
    v2f_src_dir = os.path.dirname(func_file.replace('.pyc', '.py'))
    pycache_dir = os.path.join(v2f_src_dir, '__pycache__')
    if os.path.exists(pycache_dir):
        print(f"__pycache__ exists at: {pycache_dir}")
        import glob
        pyc_files = glob.glob(os.path.join(pycache_dir, 'functions*.pyc'))
        for pf in pyc_files:
            mtime = os.path.getmtime(pf)
            print(f"  {pf} (mtime: {mtime})")
        src_file = os.path.join(v2f_src_dir, 'functions.py')
        if os.path.exists(src_file):
            src_mtime = os.path.getmtime(src_file)
            print(f"  functions.py mtime: {src_mtime}")
            for pf in pyc_files:
                pyc_mtime = os.path.getmtime(pf)
                if pyc_mtime > src_mtime:
                    print(f"  WARNING: {os.path.basename(pf)} is NEWER than functions.py!")
                    print(f"           Python may be using STALE cached bytecode!")

    # Inspect the actual getSequences source to check which version is running
    import inspect
    src = inspect.getsource(v2f.getSequences)
    if 'reverse=True' in src:
        print("getSequences contains 'reverse=True' -> NEW (reverse-sort) code")
    elif 'posadd' in src:
        print("getSequences contains 'posadd' -> OLD (forward+posadd) code")
    else:
        print("getSequences: UNKNOWN version - neither reverse nor posadd detected")

    # Print the critical section of getSequences
    print("\n--- Critical section of getSequences source ---")
    lines = src.split('\n')
    for i, line in enumerate(lines):
        if any(kw in line for kw in ['vrec_list', 'posadd', 'reverse', 'vrec.pos', 'UpdateSeq', 'getAlleles', 'vcf.fetch']):
            print(f"  {i}: {line}")
    print("--- End critical section ---\n")

    # ===== STEP 1: Load data =====
    print("=" * 70)
    print("STEP 1: Loading data")
    print("=" * 70)
    vcf = pysam.VariantFile(args.vcf)
    ref = pysam.FastaFile(args.fasta)
    samples = list(vcf.header.samples)
    print(f"Samples: {samples}")

    ploidy = v2f.getPloidy(vcf)
    print(f"Ploidy: {ploidy}")

    # Read GFF
    gff = v2f.ReadGFF(args.gff, parser)
    intervals = v2f.filterFeatureInGFF(gff, args.feat)
    genes = list(intervals.keys())
    print(f"Genes: {genes}")
    for gene in genes:
        print(f"  {gene}: {len(intervals[gene])} intervals")
        for rec in intervals[gene]:
            print(f"    {rec[0]}:{rec[3]}-{rec[4]} strand={rec[6]}")

    # ===== STEP 2: Check what vcf.fetch returns =====
    print("\n" + "=" * 70)
    print("STEP 2: Check VCF records from fetch")
    print("=" * 70)
    gene = genes[0]
    rec0 = intervals[gene][0]
    chrom = rec0[0]
    start = int(rec0[3]) - 1
    end = int(rec0[4])
    print(f"Fetching region: {chrom}:{start}-{end}")

    all_records = list(vcf.fetch(chrom, start, end))
    print(f"Total records fetched: {len(all_records)}")

    # Find pbsv records
    pbsv_records = []
    for i, vrec in enumerate(all_records):
        rec_id = vrec.id or ""
        if "pbsv" in rec_id.lower():
            pbsv_records.append((i, vrec))
            gt_str = "/".join(str(a) for a in list(vrec.samples.values())[0].get("GT"))
            phased = list(vrec.samples.values())[0].phased
            print(f"  PBSV record [{i}]: {rec_id} pos={vrec.pos} ref={vrec.ref} alt={vrec.alts[0][:30]}... "
                  f"GT={gt_str} phased={phased} ref_len={len(vrec.ref)} alt_len={len(vrec.alts[0])}")

    # Check for other variants at the same positions as pbsv records
    for idx, pbsv_rec in pbsv_records:
        pbsv_pos = pbsv_rec.pos
        same_pos = [(i, r) for i, r in enumerate(all_records) if r.pos == pbsv_pos and i != idx]
        if same_pos:
            print(f"\n  WARNING: Other variants at same position as {pbsv_rec.id} (pos={pbsv_pos}):")
            for j, r in same_pos:
                print(f"    [{j}] id={r.id} ref={r.ref} alt={r.alts[0] if r.alts else 'None'}")

    # ===== STEP 3: Monkeypatch UpdateSeq to trace pbsv processing =====
    print("\n" + "=" * 70)
    print("STEP 3: Run getSequences with instrumented UpdateSeq")
    print("=" * 70)

    # Store pbsv positions for detection
    pbsv_positions = set()
    for _, pr in pbsv_records:
        pbsv_positions.add(pr.pos)

    # Monkeypatch UpdateSeq
    original_UpdateSeq = v2f.UpdateSeq
    update_seq_call_count = [0]
    pbsv_update_count = [0]

    def traced_UpdateSeq(alleles, samp, pos, ref_len, seq):
        update_seq_call_count[0] += 1
        result = original_UpdateSeq(alleles, samp, pos, ref_len, seq)

        # Check if this might be a pbsv variant (large allele)
        allele_len = len(alleles[samp])
        if allele_len > 10 and samp.endswith("_0"):
            pbsv_update_count[0] += 1
            print(f"\n  >>> LARGE ALLELE UpdateSeq call #{update_seq_call_count[0]}:")
            print(f"      sample={samp}, pos={pos}, ref_len={ref_len}, allele_len={allele_len}")
            print(f"      seq_len_before={len(seq)}, seq_len_after={len(result)}")
            print(f"      allele[0:20]={alleles[samp][:20]}...")
            print(f"      seq[pos-3:pos+3] before = {seq[max(0,pos-3):pos+3]}")
            print(f"      result[pos-3:pos+allele_len+3] = {result[max(0,pos-3):pos+allele_len+3][:50]}...")
            # Verify the allele was actually inserted
            if alleles[samp] in result:
                print(f"      VERIFY: allele IS present in result ✓")
            else:
                print(f"      VERIFY: allele NOT found in result ✗")
                # Check if first 20 chars are present
                if alleles[samp][:20] in result:
                    print(f"      VERIFY: first 20 chars found ✓ (might be partially overwritten later)")
                else:
                    print(f"      VERIFY: first 20 chars NOT found ✗")

        return result

    v2f.UpdateSeq = traced_UpdateSeq

    # Also monkeypatch getAlleles to trace pbsv records
    original_getAlleles = v2f.getAlleles

    def traced_getAlleles(rec, ploidy, phased, addref):
        result, max_len = original_getAlleles(rec, ploidy, phased, addref)
        rec_id = rec.id or ""
        if "pbsv" in rec_id.lower():
            print(f"\n  >>> getAlleles for {rec_id} (pos={rec.pos}, ref={rec.ref}, alt={rec.alts[0][:30]}...):")
            print(f"      max_len={max_len}")
            for samp, allele in result.items():
                print(f"      {samp}: allele_len={len(allele)}, allele[:30]={allele[:30]}...")
        return result, max_len

    v2f.getAlleles = traced_getAlleles

    # Create mock args
    class MockArgs:
        def __init__(self):
            self.blend = True
            self.feat = args.feat
            self.bed = None
            self.addref = False
            self.inframe = False
            self.skip = False
            self.remove_gap = False

    mock_args = MockArgs()
    phased = True

    # Detect auto-blend
    single = [len(intervals[i]) == 1 for i in genes]
    if all(single):
        mock_args.blend = True
        print("Auto-blend: True (all genes have single records)")

    # Run getSequences
    print(f"\nCalling v2f.getSequences for gene={gene}...")
    sequences, varsites = v2f.getSequences(
        intervals, gene, ref, vcf, ploidy, phased, samples, mock_args
    )
    print(f"\nTotal UpdateSeq calls: {update_seq_call_count[0]}")
    print(f"Large allele UpdateSeq calls: {pbsv_update_count[0]}")
    print(f"Variant sites: {varsites}")

    # Restore originals
    v2f.UpdateSeq = original_UpdateSeq
    v2f.getAlleles = original_getAlleles

    # ===== STEP 4: Check output sequences for insertion content =====
    print("\n" + "=" * 70)
    print("STEP 4: Check output sequences")
    print("=" * 70)

    for featname in sequences.keys():
        print(f"\nFeature: {featname}")
        for samp in sequences[featname].keys():
            seq = sequences[featname][samp]
            gap_count = seq.count('-')
            print(f"  {samp}: length={len(seq)}, gaps={gap_count}")

        # Check for pbsv insertion sequences
        for _, pbsv_rec in pbsv_records:
            alt = pbsv_rec.alts[0]
            rec_id = pbsv_rec.id
            print(f"\n  Checking for {rec_id} (ALT len={len(alt)}):")

            for samp in sequences[featname].keys():
                seq = sequences[featname][samp]
                # Full ALT
                if alt in seq:
                    pos_found = seq.index(alt)
                    print(f"    {samp}: FULL ALT found at position {pos_found} ✓")
                else:
                    print(f"    {samp}: FULL ALT NOT found ✗")
                    # Try substrings
                    for substr_len in [40, 20, 10]:
                        if len(alt) > substr_len:
                            # Try middle
                            mid = len(alt) // 2
                            substr = alt[mid:mid+substr_len]
                            if substr in seq:
                                pos_found = seq.index(substr)
                                print(f"    {samp}: middle {substr_len}bp found at {pos_found} ✓")
                                break
                    else:
                        print(f"    {samp}: No substrings found either ✗")

    # ===== STEP 5: Compare with manual reimplementation =====
    print("\n" + "=" * 70)
    print("STEP 5: Manual reimplementation for comparison")
    print("=" * 70)

    rec0 = intervals[gene][0]
    chrom = rec0[0]
    start = int(rec0[3]) - 1
    end = int(rec0[4])

    refseq = ref.fetch(chrom, start, end).upper()
    print(f"Reference length: {len(refseq)}")

    # Initialize sequences
    manual_seqs = {}
    for sample in samples:
        for i in range(ploidy):
            manual_seqs[f"{sample}_{i}"] = refseq

    # Fetch and sort variants
    vcf2 = pysam.VariantFile(args.vcf)
    vrec_list = list(vcf2.fetch(chrom, start, end))
    vrec_list.sort(key=lambda r: r.pos, reverse=True)
    print(f"Variant records: {len(vrec_list)}")

    for vrec in vrec_list:
        pos = vrec.pos - start - 1  # Same formula as v2f code
        alleles, max_len = original_getAlleles(vrec, ploidy, True, False)
        ref_len = len(vrec.ref)

        rec_id = vrec.id or ""
        if "pbsv" in rec_id.lower():
            print(f"\n  Applying {rec_id}: pos={pos}, ref_len={ref_len}, max_len={max_len}")
            for samp in alleles:
                before_len = len(manual_seqs.get(samp, ""))
                print(f"    {samp}: allele={alleles[samp][:30]}... len={len(alleles[samp])}")

        for samp in alleles:
            if samp in manual_seqs:
                manual_seqs[samp] = original_UpdateSeq(alleles, samp, pos, ref_len, manual_seqs[samp])

        if "pbsv" in rec_id.lower():
            for samp in alleles:
                if samp in manual_seqs:
                    print(f"    {samp} after: len={len(manual_seqs[samp])}")

    # Check manual results
    print("\nManual reimplementation results:")
    for samp in sorted(manual_seqs.keys()):
        seq = manual_seqs[samp]
        gap_count = seq.count('-')
        print(f"  {samp}: length={len(seq)}, gaps={gap_count}")

    for _, pbsv_rec in pbsv_records:
        alt = pbsv_rec.alts[0]
        rec_id = pbsv_rec.id
        print(f"\n  Checking for {rec_id} (ALT len={len(alt)}):")
        for samp in sorted(manual_seqs.keys()):
            seq = manual_seqs[samp]
            if alt in seq:
                print(f"    {samp}: FULL ALT found ✓")
            else:
                print(f"    {samp}: FULL ALT NOT found ✗")

    # ===== STEP 6: Diff the outputs =====
    print("\n" + "=" * 70)
    print("STEP 6: Compare v2f vs manual output")
    print("=" * 70)
    for samp in sorted(manual_seqs.keys()):
        v2f_seq = None
        for featname in sequences:
            if samp in sequences[featname]:
                v2f_seq = sequences[featname][samp]
                break

        manual_seq = manual_seqs[samp]
        if v2f_seq is None:
            print(f"  {samp}: NOT in v2f output")
            continue

        if v2f_seq == manual_seq:
            print(f"  {samp}: IDENTICAL (len={len(v2f_seq)})")
        else:
            print(f"  {samp}: DIFFERENT!")
            print(f"    v2f:    len={len(v2f_seq)}, gaps={v2f_seq.count('-')}")
            print(f"    manual: len={len(manual_seq)}, gaps={manual_seq.count('-')}")
            # Find first difference
            min_len = min(len(v2f_seq), len(manual_seq))
            for i in range(min_len):
                if v2f_seq[i] != manual_seq[i]:
                    print(f"    First diff at position {i}:")
                    print(f"      v2f[{i-5}:{i+10}]    = {v2f_seq[max(0,i-5):i+10]}")
                    print(f"      manual[{i-5}:{i+10}] = {manual_seq[max(0,i-5):i+10]}")
                    break


if __name__ == "__main__":
    main()
