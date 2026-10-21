"""Toolkit identity and version guard, shared by the Snakefile and the bulkrnaseq launcher.

A project config pins the toolkit release it was written for (toolkit.version). A run is allowed only when
that version equals the Version field of DESCRIPTION and the checkout has no uncommitted changes, so every
results folder maps to one commit. Development mode (BULKRNASEQ_DEV=1, set by `bulkrnaseq run --dev`) turns
both checks into warnings; the run manifest records it.

Standard library only: the launcher may run under any Python 3.10+.
"""

from __future__ import annotations

import hashlib
import os
import subprocess
from pathlib import Path

DEV_ENV = "BULKRNASEQ_DEV"
# Files whose content determines the results: a change here re-runs the whole workflow.
CODE_PATHS = ("DESCRIPTION", "NAMESPACE", "R", "inst", "workflow")


class ToolkitError(RuntimeError):
    """The checked-out toolkit cannot produce results for this config."""


def package_version(repo: str | Path) -> str:
    for line in Path(repo, "DESCRIPTION").read_text().splitlines():
        if line.startswith("Version:"):
            return line.split(":", 1)[1].strip()
    raise ToolkitError(f"no Version field in {Path(repo, 'DESCRIPTION')}")


def _git(repo: str | Path, *args: str) -> str | None:
    try:
        res = subprocess.run(["git", "-C", str(repo), *args], capture_output=True, text=True, check=False)
    except FileNotFoundError:
        return None
    return res.stdout.strip() if res.returncode == 0 else None


def git_state(repo: str | Path) -> dict:
    """Commit, `git describe` and dirty flag; all None outside a git checkout."""
    commit = _git(repo, "rev-parse", "HEAD")
    if commit is None:
        return {"git_commit": None, "git_describe": None, "git_dirty": None}
    status = _git(repo, "status", "--porcelain")
    return {
        "git_commit": commit,
        "git_describe": _git(repo, "describe", "--tags", "--always", "--dirty"),
        "git_dirty": bool(status),
    }


def code_md5(repo: str | Path) -> str:
    """md5 over the relative paths and contents of CODE_PATHS (caches excluded)."""
    repo = Path(repo)
    files = []
    for entry in CODE_PATHS:
        p = repo / entry
        files.extend([p] if p.is_file() else sorted(f for f in p.rglob("*") if f.is_file()))
    h = hashlib.md5()
    for f in files:
        if "__pycache__" in f.parts:
            continue
        h.update(f.relative_to(repo).as_posix().encode() + b"\0")
        h.update(f.read_bytes() + b"\0")
    return h.hexdigest()


def dev_mode(environ: dict | None = None) -> bool:
    value = (environ if environ is not None else os.environ).get(DEV_ENV, "")
    return value.strip().lower() in {"1", "true", "yes"}


def check_toolkit(pinned: str, repo: str | Path, dev: bool) -> tuple[dict, list[str]]:
    """Return (state, problems); raise ToolkitError when there are problems outside development mode."""
    state = {
        "version": package_version(repo),
        "pinned_version": pinned,
        **git_state(repo),
        "code_md5": code_md5(repo),
        "dev_mode": dev,
    }
    problems = []
    if pinned != state["version"]:
        problems.append(
            f"the config pins toolkit.version {pinned} but the checkout is {state['version']}: "
            f"check out that release (git checkout v{pinned}) or update toolkit.version in the config"
        )
    if state["git_dirty"] is None:
        problems.append(f"{repo} is not a git checkout: results could not be tied to a commit")
    elif state["git_dirty"]:
        problems.append("the toolkit checkout has uncommitted changes (git status --porcelain is not empty)")
    if problems and not dev:
        raise ToolkitError(
            "; ".join(problems) + ". Commit the changes, or use development mode (bulkrnaseq run --dev, "
            "BULKRNASEQ_DEV=1) for runs that are not final results"
        )
    return state, problems
