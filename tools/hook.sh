#!/usr/bin/env bash
# Runs one tool of the toolkit environment for a pre-commit hook, without activating the environment.
# Environment: BULKRNASEQ_ENV (prefix) if set; otherwise the folder of the Python running pre-commit
# (/proc/$PPID/exe, Linux), which is the environment pre-commit was installed from; otherwise PATH.
# Usage: tools/hook.sh <tool> [args...]
set -euo pipefail

if [[ -n "${BULKRNASEQ_ENV:-}" ]]; then
  PATH="${BULKRNASEQ_ENV}/bin:${PATH}"
elif [[ -e "/proc/${PPID}/exe" ]]; then
  parent_bin="$(dirname "$(readlink -f "/proc/${PPID}/exe")")"
  if [[ -x "${parent_bin}/Rscript" || -x "${parent_bin}/ruff" ]]; then
    PATH="${parent_bin}:${PATH}"
  fi
fi

tool="$1"
shift
if ! command -v "${tool}" > /dev/null; then
  echo "hook: '${tool}' not found; set BULKRNASEQ_ENV to the toolkit environment prefix" >&2
  exit 1
fi
exec "${tool}" "$@"
