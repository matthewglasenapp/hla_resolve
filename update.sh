#!/usr/bin/env bash
# update.sh — update an existing hla_resolve clone to the latest version.
#
# Safe to run from any working directory: it cd's to its own location (the repo
# root) first, so `bash update.sh`, `bash ./update.sh`, and
# `bash /abs/path/to/hla_resolve/update.sh` all behave identically.
set -euo pipefail

cd "$(dirname "$0")"

echo "Updating hla_resolve in: $(pwd)"
git fetch origin --tags
git checkout main
git pull --ff-only

# Re-run the editable install so the recorded version, entry points, and deps
# are refreshed. Without this, `hla_resolve --version` keeps reporting the
# version captured at the previous install.
pip install -e .

echo
echo "Now on:        $(hla_resolve --version)"
echo "Latest tagged: $(git tag --sort=-v:refname | head -1)"
