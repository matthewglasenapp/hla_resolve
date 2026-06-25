#!/usr/bin/env bash
# verify_download_versions.sh
#
# Confirm each downloaded dependency actually IS the version it claims — don't
# trust the filename. Expected versions/digests are read from config.py (the
# source of truth), and each download is checked against what it reports about
# itself:
#   - picard      : version embedded in the jar manifest
#   - longphase   : --version
#   - rammap      : --version
#   - deepvariant : MANIFEST DIGEST of the pinned docker tag, NOT --version.
#                   The google/deepvariant:1.6.1 image self-reports "1.6.0"
#                   (upstream ARG VERSION bug, issue #830), so --version is
#                   useless here; the digest is the only reliable identity.
#
# Usage:  bash hla_resolve/verify_download_versions.sh
set -uo pipefail

# Locate the installed package's data dir + config.py. Importing the bare
# package is safe; importing hla_resolve.config would trigger the downloads.
PKG=$(python -c "import hla_resolve, os; print(os.path.dirname(hla_resolve.__file__))")
DATA="$PKG/data"; CONFIG="$PKG/config.py"

# Pull EXPECTED values from config.py constants (not from filenames).
ver() { grep -oE "^$1 = \"[^\"]+\"" "$CONFIG" | sed -E 's/.*"([^"]+)".*/\1/'; }
PICARD_V=$(ver PICARD_VERSION);     LONGPHASE_V=$(ver LONGPHASE_VERSION)
RAMMAP_V=$(ver RAMMAP_VERSION);     DV_V=$(ver DEEPVARIANT_VERSION)
DV_DIGEST=$(ver DEEPVARIANT_DIGEST)

pass=0; fail=0
ok()  { echo "  PASS  $1"; pass=$((pass+1)); }
bad() { echo "  FAIL  $1"; fail=$((fail+1)); }
check_substr() { # name expected output
  if grep -qF "$2" <<<"$3"; then ok "$1: reports '$2'"
  else bad "$1: expected '$2'; tool said: $(head -1 <<<"$3")"; fi
}

echo "== Verifying downloads against config.py (querying each download's own metadata) =="

# Picard — version embedded in the jar manifest (independent of filename)
JAR="$DATA/picard/picard_${PICARD_V}.jar"
if [[ -f "$JAR" ]]; then
  check_substr picard "$PICARD_V" "$(unzip -p "$JAR" META-INF/MANIFEST.MF 2>/dev/null | tr -d '\r')"
else bad "picard: MISSING $JAR"; fi

# longphase / rammap — best-effort --version (strip leading 'v'; binaries print '2.0')
LP="$DATA/longphase/longphase_linux-x64_${LONGPHASE_V}"
if [[ -x "$LP" ]]; then check_substr longphase "${LONGPHASE_V#v}" "$("$LP" --version 2>&1 || true)"
else bad "longphase: MISSING $LP"; fi
RM="$DATA/rammap/rammap_${RAMMAP_V}"
if [[ -x "$RM" ]]; then check_substr rammap "${RAMMAP_V#v}" "$("$RM" --version 2>&1 || "$RM" -h 2>&1 || true)"
else bad "rammap: MISSING $RM"; fi

# DeepVariant — resolve the pinned tag's manifest digest from Docker Hub and
# compare to the expected digest in config.py. (Do NOT use --version; see header.)
SIF="$DATA/deepvariant_sif/deepvariant_${DV_V}.sif"
if [[ ! -f "$SIF" ]]; then
  bad "deepvariant: MISSING $SIF"
else
  repo=google/deepvariant
  tok=$(curl -s "https://auth.docker.io/token?service=registry.docker.io&scope=repository:${repo}:pull" \
        | python -c "import sys,json;print(json.load(sys.stdin)['token'])" 2>/dev/null)
  got=$(curl -sI -H "Authorization: Bearer $tok" \
             -H "Accept: application/vnd.oci.image.index.v1+json" \
             -H "Accept: application/vnd.docker.distribution.manifest.list.v2+json" \
             "https://registry-1.docker.io/v2/${repo}/manifests/${DV_V}" \
        | awk 'tolower($1)=="docker-content-digest:"{print $2}' | tr -d '\r')
  if [[ -z "$got" ]]; then
    bad "deepvariant: could not resolve digest for :${DV_V} (no network to Docker Hub?)"
  elif [[ "$got" == "$DV_DIGEST" ]]; then
    ok "deepvariant: tag :${DV_V} resolves to expected digest ${DV_DIGEST}"
  else
    bad "deepvariant: tag :${DV_V} digest is $got, expected $DV_DIGEST"
  fi
fi

echo "== $pass passed, $fail failed/missing =="
[[ $fail -eq 0 ]]
