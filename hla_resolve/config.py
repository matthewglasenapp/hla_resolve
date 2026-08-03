# This software is Copyright ©2026. The Regents of the University of California
# ("Regents"). All Rights Reserved.
#
# See LICENSE.txt for license details.

# Configuration constants and paths for HLA-Resolve
import os
import re
import time
import fcntl
import subprocess
from pathlib import Path
from zipfile import ZipFile
from contextlib import contextmanager

# Runtime verbosity flag. Set to True by cli.py when --verbose is passed.
# Consumers should `from . import config` and check `config.VERBOSE` so the
# runtime mutation is visible (a `from .config import VERBOSE` import would
# bind the False value at import time).
VERBOSE = False

# Get the data directory relative to this config file
_data_dir = os.path.join(os.path.dirname(os.path.abspath(__file__)), "data")

# Pinned versions for downloaded tools. The version is baked into the on-disk
# filename, so the exists() check in each ensure_* below is version-aware:
# bumping a version here changes the target filename, which makes the check fail
# for the old artifact and re-download the new one. Without this, exists() only
# saw a generic name (e.g. picard.jar) and a version bump would silently keep
# the stale artifact on existing installs.
PICARD_VERSION = "2.27.4"
LONGPHASE_VERSION = "v2.0"
RAMMAP_VERSION = "v1.0.0"
DEEPVARIANT_VERSION = "1.6.1"
# WARNING: the official google/deepvariant:1.6.1 image self-reports "1.6.0" via
# `run_deepvariant --version` — upstream left `ARG VERSION=1.6.0` in the r1.6.1
# Dockerfile (google/deepvariant issue #830). The image IS 1.6.1; do NOT "correct"
# the pin to 1.6.0. Verify identity by manifest digest (below), not by --version.
DEEPVARIANT_DIGEST = "sha256:ccab95548e6c3ec28c75232987f31209ff1392027d67732435ce1ba3d0b55c68"
CLAIR3_VERSION = "latest"  # NOTE: a moving tag — filename can't detect upstream :latest updates
# IPD-IMGT/HLA database release. Keep in sync with the active (uncommented) db_url
# in ensure_hla_xml() below; verify_download_versions.sh checks the installed
# hla.xml's <release version="..."> header against this.
IMGT_RELEASE = "3.64.0"

@contextmanager
def _setup_lock(resource_dir):
    """Serialize first-run setup across processes with an flock on .setup.lock.

    Array jobs launched right after install would otherwise download and build the
    same shared artifacts concurrently. Unix-only (fcntl).
    """
    resource_dir = Path(resource_dir)
    resource_dir.mkdir(parents=True, exist_ok=True)
    lock_path = resource_dir / ".setup.lock"
    with open(lock_path, "w") as lock:
        fcntl.flock(lock, fcntl.LOCK_EX)
        try:
            yield
        finally:
            fcntl.flock(lock, fcntl.LOCK_UN)

def _wget_atomic(url, dest, executable=False):
    """Download `url` to a temp file in dest's directory, then atomically rename.

    os.replace() is atomic on a single filesystem, so other processes never
    observe a partially-downloaded `dest`: they see either no file or the
    complete one. Pair with _setup_lock so a crashed download leaves only a
    .tmp stub (ignored by the exists() check), never a half-written `dest`.
    """
    dest = Path(dest)
    tmp = dest.with_name(dest.name + ".tmp")
    subprocess.run(["wget", url, "-O", str(tmp)], check=True)
    if executable:
        subprocess.run(["chmod", "+x", str(tmp)], check=True)
    os.replace(tmp, dest)

# 50-base spot-check windows inside each hard-mask region.
_DRB5_MASK_CHECK = "chr6:32517400-32517450"
_DRB6_MASK_CHECK = "chr6:32556800-32556850"

def _region_all_N(augmented_file, region, retries=3):
    """Return True if every base of `region` is N, False if any base is real.

    Raises RuntimeError if the region cannot be read (missing index, samtools
    failure, transient network-FS blip) after `retries` attempts. The point is to
    NEVER silently mistake "could not read it" for "not masked": a transient read
    failure under a cold N-way array launch on a shared filesystem was what let a
    good reference get deleted and rebuilt, cascading the whole array into a race.
    """
    last_err = None
    for _ in range(retries):
        try:
            result = subprocess.run(
                ["samtools", "faidx", str(augmented_file), region],
                capture_output=True, text=True, check=True,
            )
            seq = "".join(line for line in result.stdout.splitlines() if not line.startswith(">"))
            if not seq:
                raise RuntimeError(f"empty read for {region}")
            return all(c.upper() == "N" for c in seq)
        except (subprocess.CalledProcessError, FileNotFoundError, RuntimeError) as e:
            last_err = e
            time.sleep(2)
    raise RuntimeError(f"could not read {region} from {augmented_file}: {last_err}")

def _drb5_is_masked(augmented_file):
    """Fast-path check: is HLA-DRB5 hard-masked? Swallows read errors (returns
    False) so a flake just sends us into the lock, where the rebuild decision is
    made strictly via _region_all_N. Never deletes anything on its own.
    """
    if not augmented_file.exists():
        return False
    try:
        return _region_all_N(augmented_file, _DRB5_MASK_CHECK, retries=1)
    except RuntimeError:
        return False

def _drb6_is_masked(augmented_file):
    """Fast-path check: is HLA-DRB6 hard-masked? Mirrors _drb5_is_masked (swallows
    read errors so a flake sends us into the lock, never deletes on its own).
    """
    if not augmented_file.exists():
        return False
    try:
        return _region_all_N(augmented_file, _DRB6_MASK_CHECK, retries=1)
    except RuntimeError:
        return False

def ensure_reference_genome():
    """Build augmented_hg38.fa: GRCh38 plus the HLA-Y/OLI scaffold, with HLA-DRB5
    (chr6:32517353-32530287) and HLA-DRB6 (chr6:32552713-32560022) hard-masked so
    divergent DRB1 alleles do not lose MAPQ to those paralogs. Real DRB5 and DRB6
    reads are removed downstream by the DRB bait. DRB9 is left intact.

    An existing build with either region unmasked is removed and rebuilt.
    """
    ref_dir = Path(_data_dir) / "reference"
    grch38_file = ref_dir / "GCA_000001405.15_GRCh38_no_alt_analysis_set.fna"
    augmented_file = ref_dir / "augmented_hg38.fa"
    hla_y_file = ref_dir / "hla_y_scaffold.fasta"

    # Fast path (no lock): a properly-masked augmented reference already exists.
    if _drb5_is_masked(augmented_file) and _drb6_is_masked(augmented_file):
        return

    with _setup_lock(ref_dir):
        # Re-check inside the lock: another process may have built it while we
        # were waiting, in which case there is nothing left to do.
        if _drb5_is_masked(augmented_file) and _drb6_is_masked(augmented_file):
            return

        # Otherwise rebuild — either no file at all, or a v0.1.0 unmasked file,
        # or a DRB5-only-masked reference from before the DRB6 mask was added.
        # But a file that exists is only deleted if we can CONFIRM it is genuinely
        # missing a mask. The fast-path check above swallows read errors, so an
        # existing-but-not-confirmed file might just have flaked under FS load;
        # verify strictly here and refuse to destroy a reference we cannot read.
        if augmented_file.exists():
            try:
                confirmed_masked = (
                    _region_all_N(augmented_file, _DRB5_MASK_CHECK)
                    and _region_all_N(augmented_file, _DRB6_MASK_CHECK)
                )
            except RuntimeError as e:
                raise RuntimeError(
                    f"Could not verify HLA-DRB5/DRB6 masks on existing {augmented_file} "
                    f"({e}). Refusing to delete or rebuild it, so a transient filesystem "
                    f"error cannot destroy a good reference. Re-run 'hla_resolve setup' "
                    f"when the filesystem is quiet, or remove the file manually if it is "
                    f"truly stale."
                )
            if confirmed_masked:
                # The reference is fine; the fast-path check merely flaked earlier.
                return
            print("Existing augmented_hg38.fa is missing a required mask (DRB5/DRB6). Rebuilding...")
            augmented_file.unlink()
            fai = augmented_file.with_suffix(".fa.fai")
            if fai.exists():
                fai.unlink()
        # Clean up the v0.1.0 sidecar that's no longer used by v0.2.0+
        for stale in (ref_dir / "augmented_hg38_drb_masked.fa",
                      ref_dir / "augmented_hg38_drb_masked.fa.fai"):
            if stale.exists():
                stale.unlink()

        original_cwd = os.getcwd()
        os.chdir(ref_dir)

        try:
            # Download base GRCh38 if missing. Download+decompress to temp names
            # and atomically rename, so an interrupted fetch never leaves a
            # partial .fna that a later run would mistake for complete.
            if not grch38_file.exists():
                print("Downloading GRCh38 reference genome...")
                tmp_gz = ref_dir / "grch38.tmp.fna.gz"
                tmp_fna = ref_dir / "grch38.tmp.fna"
                subprocess.run([
                    "wget",
                    "https://ftp.ncbi.nlm.nih.gov/genomes/all/GCA/000/001/405/GCA_000001405.15_GRCh38/seqs_for_alignment_pipelines.ucsc_ids/GCA_000001405.15_GRCh38_no_alt_analysis_set.fna.gz",
                    "-O", str(tmp_gz)
                ], check=True)
                subprocess.run(["gunzip", str(tmp_gz)], check=True)
                os.replace(tmp_fna, grch38_file)
                subprocess.run(["samtools", "faidx", str(grch38_file)], check=True)

            if not hla_y_file.exists():
                raise FileNotFoundError(f"HLA-Y/OLI scaffold not found: {hla_y_file}")

            # Concatenate GRCh38 + HLA-Y/OLI scaffold into a temporary unmasked file
            unmasked = ref_dir / "augmented_hg38.unmasked.tmp.fa"
            print("Augmenting GRCh38 with HLA-Y/OLI scaffold...")
            subprocess.run(["bash", "-c", f"cat {grch38_file} {hla_y_file} > {unmasked}"], check=True)

            # Hard-mask HLA-DRB5 and HLA-DRB6.
            # GFF3 is 1-based inclusive; BED is 0-based half-open -> start-1, end unchanged.
            # DRB6 exon 2 competes with DRB1 for divergent DR4 reads (e.g. DRB1*04:12),
            # costing them MAPQ and coverage. Real DRB6 reads are still removed by the
            # DRB panel bait (DRB6_GRCh38), as with DRB5.
            print("Hard-masking HLA-DRB5 (chr6:32517353-32530287) and HLA-DRB6 (chr6:32552713-32560022)...")
            mask_bed = ref_dir / "drb5_mask.bed"
            with open(mask_bed, "w") as f:
                f.write("chr6\t32517352\t32530287\n")
                f.write("chr6\t32552712\t32560022\n")

            # Build into temp names, then atomically publish the reference and its
            # index together, so the augmented_file only appears fully built and
            # indexed — the fast-path check above never sees a half-masked file.
            masked_tmp = ref_dir / "augmented_hg38.masked.tmp.fa"
            subprocess.run([
                "bedtools", "maskfasta",
                "-fi", str(unmasked),
                "-bed", str(mask_bed),
                "-fo", str(masked_tmp)
            ], check=True)
            subprocess.run(["samtools", "faidx", str(masked_tmp)], check=True)

            os.replace(masked_tmp, augmented_file)
            os.replace(str(masked_tmp) + ".fai", str(augmented_file) + ".fai")

            unmasked.unlink()
            print(f"Reference ready: {augmented_file} (HLA-Y/OLI added, HLA-DRB5 and HLA-DRB6 hard-masked)")
        finally:
            os.chdir(original_cwd)

def ensure_longphase():
    """Download and extract longphase if not present"""
    longphase_dir = Path(_data_dir) / "longphase"
    longphase_bin = longphase_dir / f"longphase_linux-x64_{LONGPHASE_VERSION}"
    tar_file = longphase_dir / f"longphase_linux-x64_{LONGPHASE_VERSION}.tar.xz"

    if longphase_bin.exists():
        return str(longphase_bin)

    with _setup_lock(longphase_dir):
        if not longphase_bin.exists():
            print(f"Longphase {LONGPHASE_VERSION} not found! Downloading longphase...")

            # Download
            subprocess.run([
                "wget",
                f"https://github.com/twolinin/longphase/releases/download/{LONGPHASE_VERSION}/longphase_linux-x64.tar.xz",
                "-O", str(tar_file)
            ], check=True)

            # Extract into a temp dir, then atomically move the binary into place,
            # so a crash mid-extract can't leave a partial binary that a later run
            # mistakes for complete.
            print("Extracting longphase...")
            extract_tmp = longphase_dir / ".extract.tmp"
            subprocess.run(["rm", "-rf", str(extract_tmp)], check=True)
            extract_tmp.mkdir(parents=True)
            subprocess.run([
                "tar", "-xJf", str(tar_file), "-C", str(extract_tmp)
            ], check=True)

            extracted_bin = extract_tmp / "longphase_linux-x64"
            subprocess.run(["chmod", "+x", str(extracted_bin)], check=True)
            os.replace(extracted_bin, longphase_bin)

            # Clean up the tarball and temp extraction dir
            tar_file.unlink()
            subprocess.run(["rm", "-rf", str(extract_tmp)], check=True)
            print("Longphase download complete!")

    return str(longphase_bin)

def ensure_rammap():
    """Download rammap binary if not present"""
    rammap_dir = Path(_data_dir) / "rammap"
    rammap_bin = rammap_dir / f"rammap_{RAMMAP_VERSION}"

    if rammap_bin.exists():
        return str(rammap_bin)

    with _setup_lock(rammap_dir):
        if not rammap_bin.exists():
            print(f"Rammap {RAMMAP_VERSION} not found! Downloading rammap...")
            _wget_atomic(
                f"https://github.com/jwanglab/rammap/releases/download/{RAMMAP_VERSION}/rammap_x86_64-unknown-linux-gnu_{RAMMAP_VERSION}",
                rammap_bin,
                executable=True,
            )
            print("Rammap download complete!")

    return str(rammap_bin)

def ensure_picard():
    """Download Picard if not present"""
    picard_dir = Path(_data_dir) / "picard"
    picard_jar = picard_dir / f"picard_{PICARD_VERSION}.jar"

    if picard_jar.exists():
        return str(picard_jar)

    with _setup_lock(picard_dir):
        if not picard_jar.exists():
            print(f"Picard {PICARD_VERSION} not found! Downloading Picard")
            _wget_atomic(
                f"https://github.com/broadinstitute/picard/releases/download/{PICARD_VERSION}/picard.jar",
                picard_jar,
            )
            print("Picard download complete!")

    return str(picard_jar)

def ensure_deepvariant_sif():
    """Pull DeepVariant Singularity image if not present"""
    sif_dir = Path(_data_dir) / "deepvariant_sif"
    sif_file = sif_dir / f"deepvariant_{DEEPVARIANT_VERSION}.sif"

    if sif_file.exists():
        return str(sif_file)

    with _setup_lock(sif_dir):
        if not sif_file.exists():
            print(f"DeepVariant SIF {DEEPVARIANT_VERSION} not found! Pulling from Docker Hub...")
            tmp_sif = sif_file.with_name(sif_file.name + ".tmp")
            # Pin by DIGEST, not tag, so we get exactly the verified 1.6.1 image
            # (immune to tag mutation), and --disable-cache so the pull never
            # reuses a previously-cached image. _setup_lock already bounds this to
            # a single pull per install (other array tasks wait and reuse the
            # resulting .sif), so disabling the cache costs at most one ~2.7 GB
            # download, not one per task.
            subprocess.run([
                "singularity", "pull", "--force", "--disable-cache",
                str(tmp_sif),
                f"docker://google/deepvariant@{DEEPVARIANT_DIGEST}"
            ], check=True)
            os.replace(tmp_sif, sif_file)
            print("DeepVariant SIF download complete!")

    return str(sif_file)

def ensure_clair3_sif():
    """Pull Clair3 Singularity image if not present"""
    sif_dir = Path(_data_dir) / "clair3_sif"
    sif_file = sif_dir / f"clair3_{CLAIR3_VERSION}.sif"

    if sif_file.exists():
        return str(sif_file)

    with _setup_lock(sif_dir):
        if not sif_file.exists():
            print(f"Clair3 SIF {CLAIR3_VERSION} not found! Pulling from Docker Hub...")
            tmp_sif = sif_file.with_name(sif_file.name + ".tmp")
            subprocess.run([
                "singularity", "pull", "--force",
                str(tmp_sif),
                f"docker://hkubal/clair3:{CLAIR3_VERSION}"
            ], check=True)
            os.replace(tmp_sif, sif_file)
            print("Clair3 SIF download complete!")

    return str(sif_file)

def ensure_hla_xml():
    """Download HLA XML database if not present"""
    xml_dir = Path(_data_dir) / "IPD_IMGT_XML"
    xml_file = xml_dir / f"hla_{IMGT_RELEASE}.xml"
    zip_file = xml_dir / "hla.xml.zip"
    # IPD-IMGT/HLA Release 3.60.0
    #db_url = "https://raw.githubusercontent.com/ANHIG/IMGTHLA/652dbe954426f117a9f3619826fc4e3687713d90/xml/hla.xml.zip"
    # IPD-IMGT/HLA Release 3.61.0
    #db_url = "https://raw.githubusercontent.com/ANHIG/IMGTHLA/93c70bcfe271a737bc75b7ca7f5f9844bf65136d/xml/hla.xml.zip"
    # IPD-IMGT/HLA Release 3.64.0 (v3.64.0-alpha) — pinned; developed and validated against this release
    db_url = "https://raw.githubusercontent.com/ANHIG/IMGTHLA/bdc6233a79b1d0323140c25ef6e2e6d4e7187c26/xml/hla.xml.zip"
    # Latest (uncomment to always pull the newest IPD-IMGT/HLA release; NOT reproducible)
    #db_url = "https://raw.githubusercontent.com/ANHIG/IMGTHLA/Latest/xml/hla.xml.zip"

    # Fast path (no lock): the database is already present.
    if xml_file.exists():
        return

    with _setup_lock(xml_dir):
        if xml_file.exists():
            return
        print(f"INFO: Downloading HLA XML database (IPD-IMGT/HLA {IMGT_RELEASE})")

        # Download the zip file
        subprocess.run([
            "wget",
            db_url,
            "-O", str(zip_file)
        ], check=True)

        # Extract into a temp dir, then atomically move the versioned hla.xml into
        # place, so a crash mid-extract can't leave a partial file that a later run
        # treats as complete.
        extract_tmp = xml_dir / ".extract.tmp"
        subprocess.run(["rm", "-rf", str(extract_tmp)], check=True)
        extract_tmp.mkdir(parents=True)
        with ZipFile(zip_file) as zip_ref:
            zip_ref.extractall(extract_tmp)
        extracted = extract_tmp / "hla.xml"

        # Self-check: the downloaded file's <release version="..."> header must
        # match the pinned IMGT_RELEASE, guarding against the db_url commit and the
        # IMGT_RELEASE constant drifting out of sync. (hla.xml is ISO-8859-1.)
        with open(extracted, encoding="latin-1") as fh:
            header = fh.read(4096)
        match = re.search(r'<release version="([0-9.]+)"', header)
        found = match.group(1) if match else None
        if found != IMGT_RELEASE:
            raise ValueError(
                f"Downloaded hla.xml is IPD-IMGT/HLA release {found}, but IMGT_RELEASE "
                f"is {IMGT_RELEASE} -- the db_url commit and IMGT_RELEASE are out of sync."
            )

        os.replace(extracted, xml_file)

        # Clean up the zip and temp extraction dir
        zip_file.unlink()
        subprocess.run(["rm", "-rf", str(extract_tmp)], check=True)
        print(f"HLA XML database download complete! (IPD-IMGT/HLA {IMGT_RELEASE})")

def run_setup():
    """Download and build every external dependency once, up front.

    Intended to be run a single time right after install (`hla_resolve setup`),
    so later runs — including array jobs launched simultaneously — find every
    artifact already present and never race to download or rebuild it. Idempotent:
    safe to re-run, each step no-ops when its artifact is already there.
    """
    ensure_reference_genome()
    ensure_picard()
    ensure_longphase()
    ensure_rammap()
    ensure_hla_xml()
    ensure_deepvariant_sif()
    # ensure_clair3_sif()  # re-enable when ONT support lands
    print("hla_resolve setup complete: all dependencies present.")

# Download reference genome, Picard, longphase, rammap, HLA XML database, DeepVariant SIF, and Clair3 SIF on first import
ensure_reference_genome()
picard = ensure_picard()
longphase = ensure_longphase()
rammap = ensure_rammap()
ensure_hla_xml()
deepvariant_sif = ensure_deepvariant_sif()
# clair3_sif = ensure_clair3_sif()  # TODO: re-enable when running ONT
clair3_sif = None

# HLA genes of interest for HLA typing
genes_of_interest = ("HLA-A", "HLA-B", "HLA-C", "HLA-DPA1", "HLA-DPB1", "HLA-DQA1", "HLA-DQB1", "HLA-DRB1")

# Class I HLA genes (ARS = CDS exon 2 + exon 3, with intron 2 between them)
CLASS_I_GENES = {"HLA-A", "HLA-B", "HLA-C"}

# Multi-allele DRB reference for competitive read classification.
# 33 entries: 13 DRB1 + 3 DRB3 + 1 DRB4 (IPD-IMGT/HLA), GRCh38-extracted
# DRB5/DRB6/DRB9, DRB7/DRB8 from the SSTO haplotype, one CM089004.1 element, and
# 10 paralog/intergenic blocks lifted from HPRC assemblies (most with their own
# DRB1 copy masked, so they compete only for paralog reads). Reads whose
# competitive primary is anything other than DRB1*XX are flagged for removal by
# filter_reads(). Used in classify_DRB_reads() function of preprocess_methods.py.
drb_multiallele_reference = os.path.join(_data_dir, "reference/DRB_reference.fa")

# GRCh38 extent of the DR sub-region (HLA-DRA .. HLA-DRB1), spanning the DRB
# paralog cluster (DRB5/DRB6/DRB1). On every PacBio scheme, competitive DRB
# classification is restricted to primary reads already placed in this window,
# rather than the whole read set: we are reclassifying reads that mapped into
# the DR region, not rehoming unmapped reads. Removes genome-wide homology noise
# from the kill-list. Reads are selected by overlap, so a read starting inside
# the window and extending past HLA-DRB1's 3' end is still classified.
drb_region = "chr6:32439878-32589848"


# Clair3 model names — bundled inside the Clair3 SIF at /opt/models/
# Users can override the ONT model via --clair3_model at the command line
# Available models: https://github.com/HKU-BAL/Clair3#pre-trained-models
clair3_ont_model = "r1041_e82_400bps_sup_v500"   # R10.4.1 SUP (default)
clair3_hifi_model = "hifi_revio"

# Reference fasta: GRCh38 + HLA-OLI/HLA-Y scaffold, with HLA-DRB5 and HLA-DRB6
# hard-masked. Masking forces divergent DRB1 reads to anchor at DRB1 instead of
# misrouting to a paralog locus: DRB5 draws *07/*09, DRB6 draws divergent DR4
# alleles such as DRB1*04:12. Built by ensure_reference_genome().
reference_genome_minimap2 = os.path.join(_data_dir, "reference/augmented_hg38.fa")

# DeepVariant SIF file path — populated by ensure_deepvariant_sif() above

# GRCh38 tandem repeat mask file for pbsv
# Downloaded from https://github.com/PacificBiosciences/pbsv/blob/master/annotations/human_GRCh38_no_alt_analysis_set.trf.bed
tandem_repeat_bed = os.path.join(_data_dir, "repeats_bed/human_GRCh38_no_alt_analysis_set.trf.bed")

# BED file of coordinates for chr6
# Used to constrain variant calling with Clair3 and Sniffles to chr6
# The file is simply chr6\t28000000\t34000000
chr6_bed = os.path.join(_data_dir, "reference/chr6.bed")

# GRCh38 tandem repeat definition file for pbtrgt
# Downloaded from https://zenodo.org/records/8329210
pbtrgt_repeat_file = os.path.join(_data_dir, "repeats_bed/chr6_polymorphic_repeats.hg38.bed")

# GFF files of the HLA genes of interest for FASTA reconstruction with vcf2fasta
# Used in reconstruct_fasta_methods.py
gff_dir = os.path.join(_data_dir, "hla_gff")

# BED file of coordinates for the 8 HLA genes of interest
# The coordinates were taken from the GRCh38 GFF3 file
# Used for both mosdepth coverage analysis and VCF filtering
# HLA-B and HLA-DQA1 coordinates were slightly modified from the raw GFF3 file to exclude exons that are not part of the MANE Select transcript
hla_genes_regions_file = os.path.join(_data_dir, "mosdepth/hla_genes.bed")

# BED file for parsing haploblocks in the extended MHC region
# Used in evaluate_gene_haploblocks() function of investigate_haploblocks_methods.py
genes_bed = os.path.join(_data_dir, "reference/parse_haploblocks_bed.bed")
genes_of_interest_extended = ("HLA-A", "HLA-B", "HLA-C", "HLA-DRB1", "HLA-DQA1", "HLA-DQA2", "HLA-DQB1", "HLA-DQB2", "HLA-DPA1", "HLA-DPB1")

# Parameters

# Minimum reads per sample
# DeepVariant stalls when a sample has very few BAM records (e.g., HG01891 had only 35 mapped reads to chr6)
# Set threshold at which variant calling should not proceed
# chr6_read_count is returned by filter_reads() function
# Variant calling requires that chr6_read_count >= min_reads_sample:
min_reads_sample = 1000

# This program is for long-read data only. Require that mean read length is at least 300 bp or higher
# is returned by run_fastplong(), called within parse_input_file() of sample_manager.py
# The program exits with a ValueError if mean_read_length falls below the specified threshold.
min_read_length = 300

# Coverage depth thresholds to proceed with HLA typing
# Used in parse_mosdepth() function of preprocess_methods.py
# CDS thresholds (pooled across all coding exons of a gene)
# Mean depth over concatenated CDS
cds_depth_thresh = 0
# Proportion of CDS bases with depth >= 20x
cds_prop_20x_thresh = 0.0
# Proportion of CDS bases with depth >= 30x
cds_prop_30x_thresh = 0.0
# ARS (antigen recognition site) thresholds
# Mean ARS depth
ars_depth_thresh = 8
# Proportion of ARS bases with depth >= 20x
ars_prop_20x_thresh = 0
# Proportion of ARS bases with depth >= 30x
ars_prop_30x_thresh = 0
# Extended MHC coordinates
mhc_start = 29555628
mhc_stop = 33409896

# DeepVariant calling region. Kept tighter than the extended MHC to avoid a homopolymer
# RefCall artifact near chr6 31354430 that the wider MHC bounds trigger for HG002 HLA-B.
# Still spans all eight typed genes with margin.
dv_region_start = 29900000
dv_region_stop = 33150000

# DNA bases and stop codons
# Used in parse_fastas() function of reconstruct_fasta_methods.py for vcf2fasta sanity checking
DNA_bases = {"A", "T", "G", "C"}
stop_codons = ["TAA", "TAG", "TGA"]

# Read re-consensus refinement of the 4th field for HLA-DRB1/DQA1/DQB1 only.
# Rewrites the 4th field of those genes' calls in place; never changes the
# 3-field lineage and never touches HLA-A/B/C/DPA1/DPB1. Set False to disable.
# Used in hla_typer.py (threaded from resolve_alleles_pipeline.py, PacBio only).
reconsensus_drdq = True

# Read-assignment mode for the DR/DQ re-consensus: "hp_tag" (default) splits reads
# by their HP:i:1 / HP:i:2 tag; "self_sort" sorts reads by primary alignment to a
# 2-contig scaffold reference instead.
reconsensus_read_assignment = "hp_tag"

# IPD/IMGT HLA XML file for HLA allele classification
# Downloaded from https://github.com/ANHIG/IMGTHLA
# Used in hla_typer.py for HLA typing
IMGT_XML = os.path.join(_data_dir, f"IPD_IMGT_XML/hla_{IMGT_RELEASE}.xml")

# GRCh38 HLA gene antigen recognition sequence coordinates (1-based coordinates, GFF format)
# Used in evaluate_gene_haploblocks() function of investigate_haploblocks_methods.py 
# For HLA-C/B/C, the ARS is CDS 2 and 3
# For HLA-DPA1/DPB1/DQA1/DQB1/DRB1, the ARS is CDS 2
ARS_dict = {
	"HLA-A": "chr6:29942757-29943543",
	"HLA-B": "chr6:31356167-31356957",
	"HLA-C": "chr6:31271073-31271868",
	"HLA-DPA1": "chr6:33069641-33069886",
	"HLA-DPB1": "chr6:33080672-33080935",
	"HLA-DQA1": "chr6:32641310-32641558",
	"HLA-DQB1": "chr6:32664798-32665067",
	"HLA-DRB1": "chr6:32584109-32584378"
}

# Gene coordinates (1-based coordinates, GFF format)
# Used in parse_fastas() function of reconstruct_fasta_methods.py to clamp haploblock coordinates to gene coordinates
# These match the gene feature records in data/hla_gff/<gene>_gene.gff3, which are
# the intervals vcf2fasta actually reconstructs. HLA-C is padded 1 kb on each side
# relative to the Ensembl gene (see GENE_PADDING in supplementary_scripts/sort_cds.py).
gene_dict = {
	"HLA-A":    (29941260, 29949572),
	"HLA-B":    (31353270, 31367067),
	"HLA-C":    (31267749, 31273130),
	"HLA-DPA1": (33064270, 33080775),
	"HLA-DPB1": (33075990, 33089696),
	"HLA-DQA1": (32628179, 32647062),
	"HLA-DQB1": (32659467, 32668383),
	"HLA-DRB1": (32577902, 32589848)
}

# CDS coordinates for each HLA gene of interest for HLA typing (1-based coordinates, GFF format)
# Used in parse_fastas() function of reconstruct_fasta_methods.py to determine whether each CDS is fully contained in the haploblock
CDS_dict = {
	'HLA-A': [[29942554, 29942626], [29942757, 29943026], [29943268, 29943543], [29944122, 29944397], [29944500, 29944616], [29945059, 29945091], [29945234, 29945281], [29945451, 29945455]],
	'HLA-B': [[31354483, 31354526], [31354633, 31354665], [31355107, 31355223], [31355317, 31355592], [31356167, 31356442], [31356688, 31356957], [31357086, 31357158]],
	'HLA-C': [[31269169, 31269173], [31269338, 31269385], [31269493, 31269525], [31269966, 31270085], [31270210, 31270485], [31271073, 31271348], [31271599, 31271868], [31271999, 31272071]],
	'HLA-DPA1': [[33068650, 33068804], [33069019, 33069300], [33069641, 33069886], [33073471, 33073570]], 
	'HLA-DPB1': [[33076042, 33076141], [33080672, 33080935], [33084950, 33085231], [33085779, 33085889], [33086219, 33086238]], 
	'HLA-DQA1': [[32637459, 32637540], [32641310, 32641558], [32641972, 32642253], [32642610, 32642764]], 
	'HLA-DQB1': [[32660236, 32660249], [32660859, 32660882], [32661347, 32661457], [32661967, 32662248], [32664798, 32665067], [32666499, 32666607]], 
	'HLA-DRB1': [[32579091, 32579104], [32580247, 32580270], [32580746, 32580856], [32581557, 32581838], [32584109, 32584378], [32589643, 32589742]]
}