"""Compares two results folders file by file (md5).

Files that differ between two runs by design are skipped: run_timestamps.json and the benchmarks/ folder
(wall-clock times), and .snakemake/. With --code-changed (results of different toolkit code, e.g. a refactor
checked against a baseline) toolkit.json is skipped too and run_manifest.json is compared as JSON without its
"toolkit" block, the only parts that legitimately depend on the code state.

Prints one line per file that is new, missing or changed, and exits with 1 if there is any.

Usage: python tools/compare_results.py [--code-changed] <baseline> <results>
"""

from __future__ import annotations

import hashlib
import json
import sys
from pathlib import Path

TIMING_FILES = {"run_timestamps.json"}
TIMING_DIRS = {"benchmarks", ".snakemake"}


def result_files(root: Path, code_changed: bool) -> dict[str, Path]:
    out = {}
    for path in sorted(root.rglob("*")):
        rel = path.relative_to(root)
        if not path.is_file() or TIMING_DIRS.intersection(rel.parts[:-1]) or path.name in TIMING_FILES:
            continue
        if code_changed and path.name == "toolkit.json":
            continue
        out[rel.as_posix()] = path
    return out


def md5(path: Path) -> str:
    return hashlib.md5(path.read_bytes()).hexdigest()


def manifest_without_toolkit(path: Path) -> dict:
    data = json.loads(path.read_text(encoding="utf-8"))
    data.pop("toolkit", None)
    return data


def compare(baseline: Path, results: Path, code_changed: bool = False) -> list[str]:
    a, b = result_files(baseline, code_changed), result_files(results, code_changed)
    problems = []
    for rel in sorted(set(a) | set(b)):
        if rel not in a:
            problems.append(f"new      {rel}")
        elif rel not in b:
            problems.append(f"missing  {rel}")
        elif code_changed and Path(rel).name == "run_manifest.json":
            ma, mb = manifest_without_toolkit(a[rel]), manifest_without_toolkit(b[rel])
            if ma != mb:
                keys = sorted(k for k in set(ma) | set(mb) if ma.get(k) != mb.get(k))
                problems.append(f"changed  {rel} (outside 'toolkit': {', '.join(keys)})")
        elif md5(a[rel]) != md5(b[rel]):
            problems.append(f"changed  {rel}")
    return problems


def main(argv: list[str]) -> int:
    code_changed = "--code-changed" in argv
    paths = [a for a in argv if not a.startswith("--")]
    if len(paths) != 2:
        print(__doc__.strip().splitlines()[-1], file=sys.stderr)
        return 2
    baseline, results = (Path(p) for p in paths)
    problems = compare(baseline, results, code_changed)
    for p in problems:
        print(p)
    if not problems:
        print(f"identical: {len(result_files(baseline, code_changed))} files")
    return 1 if problems else 0


if __name__ == "__main__":
    sys.exit(main(sys.argv[1:]))
