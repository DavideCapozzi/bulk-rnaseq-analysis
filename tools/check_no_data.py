"""Repository data guard.

Fails when a file that is about to be committed (or, with --all, any tracked file) looks like real data:
an absolute path to a user or mount folder, a data file outside the allowed folders, a file larger than the
size limit, or a project-specific pattern. Project patterns (e.g. sample names) are read from a file outside
the repository, so that the patterns themselves are never published; its path comes from the
BULKRNASEQ_FORBIDDEN_PATTERNS environment variable or, failing that, from the local git setting
`git config bulkrnaseq.forbiddenPatterns <path>` (stored in .git/config, never pushed). Without either,
only the generic checks run.

--all checks every tracked file plus every untracked file that is not ignored (what `git add -A` would stage).

Usage: python tools/check_no_data.py [--all] [FILE ...]
"""

from __future__ import annotations

import os
import re
import subprocess
import sys
import zipfile
from pathlib import Path

MAX_BYTES = 500 * 1024
DATA_EXTENSIONS = {
    ".xlsx",
    ".xls",
    ".csv",
    ".tsv",
    ".tab",
    ".gz",
    ".zip",
    ".bg",
    ".bedgraph",
    ".bam",
    ".sam",
    ".cram",
    ".bai",
    ".fastq",
    ".fq",
    ".gtf",
    ".gff",
    ".gff3",
    ".rds",
    ".rda",
    ".rdata",
    ".h5",
    ".hdf5",
    ".mtx",
}
# Folders where data-like files are expected: synthetic fixtures, generic gene lists, example sample sheets.
ALLOWED_DATA_DIRS = ("tests/fixtures/", "inst/extdata/gene_lists/", "config/examples/")
ABSOLUTE_PATH = re.compile(
    r"(?<![\w.])(/home/[A-Za-z0-9_.-]+|/Users/[A-Za-z0-9_.-]+|/mnt/[a-z]/|/media/[A-Za-z0-9_.-]+|"
    r"[A-Za-z]:\\\\(Users|Documents))"
)
# This file documents the patterns above; it is the only file allowed to contain them literally.
SELF = "tools/check_no_data.py"


def patterns_path(root: Path) -> str:
    path = os.environ.get("BULKRNASEQ_FORBIDDEN_PATTERNS", "")
    if path:
        return path
    out = subprocess.run(
        ["git", "-C", str(root), "config", "--get", "bulkrnaseq.forbiddenPatterns"],
        capture_output=True,
        text=True,
    )
    return out.stdout.strip()


def load_project_patterns(root: Path) -> list[re.Pattern[str]]:
    path = patterns_path(root)
    if not path or not Path(path).is_file():
        print("data guard: no project patterns configured; generic checks only", file=sys.stderr)
        return []
    lines = Path(path).read_text(encoding="utf-8").splitlines()
    return [re.compile(x.strip()) for x in lines if x.strip() and not x.lstrip().startswith("#")]


def text_of(path: Path) -> str:
    if zipfile.is_zipfile(path):
        parts = []
        with zipfile.ZipFile(path) as zf:
            for name in sorted(zf.namelist()):
                parts.append(zf.read(name).decode("utf-8", errors="replace"))
        return "\n".join(parts)
    data = path.read_bytes()
    if path.suffix == ".gz":
        import gzip

        data = gzip.decompress(data)
    return data.decode("utf-8", errors="replace")


def check_file(rel: str, root: Path, patterns: list[re.Pattern[str]]) -> list[str]:
    path = root / rel
    if not path.is_file() or path.is_symlink():
        return []
    problems = []
    size = path.stat().st_size
    if size > MAX_BYTES:
        problems.append(f"{rel}: {size} bytes > {MAX_BYTES} (size limit)")
    suffix = path.suffix.lower()
    if suffix in DATA_EXTENSIONS and not rel.startswith(ALLOWED_DATA_DIRS):
        problems.append(f"{rel}: data file extension '{suffix}' outside {', '.join(ALLOWED_DATA_DIRS)}")
    try:
        text = text_of(path)
    except (OSError, zipfile.BadZipFile, EOFError):
        return problems
    if rel != SELF:
        for m in ABSOLUTE_PATH.finditer(text):
            problems.append(f"{rel}: absolute path '{m.group(0)}'")
    for k, pat in enumerate(patterns, start=1):
        if pat.search(text) or pat.search(rel):
            # Report the pattern number only, so CI logs never echo project-specific strings.
            problems.append(f"{rel}: matches project-specific forbidden pattern #{k}")
    return problems


def candidate_files(root: Path) -> list[str]:
    cmd = ["git", "-C", str(root), "ls-files", "--cached", "--others", "--exclude-standard"]
    out = subprocess.run(cmd, check=True, capture_output=True, text=True)
    return sorted({x for x in out.stdout.splitlines() if x})


def main(argv: list[str]) -> int:
    root = Path(
        subprocess.run(
            ["git", "rev-parse", "--show-toplevel"], check=True, capture_output=True, text=True
        ).stdout.strip()
    )
    files = candidate_files(root) if "--all" in argv else [a for a in argv if not a.startswith("--")]
    patterns = load_project_patterns(root)
    problems = [p for f in files for p in check_file(f, root, patterns)]
    for p in problems:
        print(f"data guard: {p}", file=sys.stderr)
    return 1 if problems else 0


if __name__ == "__main__":
    sys.exit(main(sys.argv[1:]))
