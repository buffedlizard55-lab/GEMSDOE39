#!/usr/bin/env bash
# Placeholder. Official DrivenData downloads require a competition account and
# cannot be fetched without authentication (DrivenData ToS prohibits robots).
#
# To get competition data locally:
#   Option A (authenticated browser): log into
#     https://www.drivendata.org/competitions/306/competition-doe-gems/data/
#     and download training_features.tif, labels.tif, sample_submission.tif into data/.
#   Option B (owner-mirror via gh, SHA-256 pinned):
#     python scripts/restore_data.py --group core
#     python scripts/restore_data.py --group external
#
# Option B is used in CI; it requires gh to be authenticated against
# github.com/buffedlizard55-lab (already configured in this sandbox), and
# every byte is SHA-256 verified against registry/data_manifest.json.
set -euo pipefail
cd "$(dirname "$0")/.."
echo "Use: python scripts/restore_data.py --group core && python scripts/restore_data.py --group external"
echo "     (authenticated gh; SHA-256 pinned owner mirrors). See README."
