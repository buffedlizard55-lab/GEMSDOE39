#!/usr/bin/env bash
# Restore the data needed by the local pipeline from SHA-256-pinned owner mirrors.
# The official DrivenData data page is login-gated; this script does not bypass
# that login and does not claim the GitHub mirrors are organizer-authenticated.
set -euo pipefail
ROOT="$(cd "$(dirname "$0")/.." && pwd)"
cd "$ROOT"
if [[ -n "${PYTHON:-}" ]]; then
  PYTHON_BIN="$PYTHON"
elif [[ -x "$ROOT/.venv/bin/python" ]]; then
  PYTHON_BIN="$ROOT/.venv/bin/python"
else
  PYTHON_BIN="python3"
fi
command -v gh >/dev/null 2>&1 || {
  echo "ERROR: GitHub CLI (gh) is required for the pinned owner mirrors." >&2
  exit 1
}
if ! gh auth status >/dev/null 2>&1; then
  echo "ERROR: gh is not authenticated. Use the official competition download after login, or connect GitHub in Arena." >&2
  exit 1
fi
"$PYTHON_BIN" scripts/restore_data.py --group core
"$PYTHON_BIN" scripts/restore_data.py --group external
"$PYTHON_BIN" scripts/prepare_data.py
