"""bulkrnaseq launcher: create a project, show what would run, fetch resources, run the workflow.

The toolkit is one checkout; each project is a folder with its own data, config and results:

    <project>/input/                      provider files (read-only)
    <project>/analysis/<analysis>/config/ config.yaml, samples.tsv
    <project>/analysis/<analysis>/resources/, results/

Snakemake is taken from BULKRNASEQ_ENV (prefix of the toolkit environment), then from the folder of the
launcher or of the running Python, then from PATH. The folder that holds snakemake is put first on PATH, so
the workflow's Rscript comes from the same environment.

Standard library only: the launcher may run under any Python 3.10+.
"""

from __future__ import annotations

import argparse
import os
import re
import shutil
import subprocess
import sys
from pathlib import Path

from bulkrnaseq_toolkit import DEV_ENV, git_state, package_version

REPO = Path(__file__).resolve().parents[2]
SNAKEFILE = REPO / "workflow" / "Snakefile"
TEMPLATES = REPO / "workflow" / "templates"
# Scenario -> sample sheet columns (None: the scenario needs no sample sheet).
SCENARIOS = {
    "star": ["sample_id", "input_column", "condition", "label"],
    "excel": ["sample_id", "input_column", "condition"],
    "de_only": None,
}


def find_snakemake(environ: dict, argv0: str, executable: str) -> Path:
    prefix = environ.get("BULKRNASEQ_ENV")
    if prefix:
        exe = Path(prefix).expanduser() / "bin" / "snakemake"
        if not os.access(exe, os.X_OK):
            raise SystemExit(f"bulkrnaseq: BULKRNASEQ_ENV is set but {exe} is not an executable")
        return exe
    for folder in (Path(argv0).parent, Path(executable).parent):
        exe = folder / "snakemake"
        if os.access(exe, os.X_OK):
            return exe
    found = shutil.which("snakemake", path=environ.get("PATH"))
    if found:
        return Path(found)
    raise SystemExit("bulkrnaseq: snakemake not found; set BULKRNASEQ_ENV to the toolkit environment prefix")


def workdir_for(config: Path) -> Path:
    """Snakemake working folder (holds .snakemake/): the analysis folder in the standard layout."""
    return config.parent.parent if config.parent.name == "config" else config.parent


def snakemake_command(
    snakemake: Path, config: Path, targets: list[str], cores: int, dry_run: bool, extra: list[str]
) -> list[str]:
    cmd = [str(snakemake), "-s", str(SNAKEFILE), "--configfile", str(config)]
    cmd += ["--directory", str(workdir_for(config)), "--cores", str(cores)]
    return cmd + (["--dry-run"] if dry_run else []) + targets + extra


def run_workflow(args: argparse.Namespace, targets: list[str]) -> int:
    config = Path(args.config).expanduser().resolve()
    if not config.is_file():
        raise SystemExit(f"bulkrnaseq: config not found: {config}")
    snakemake = find_snakemake(dict(os.environ), sys.argv[0], sys.executable)
    env = dict(os.environ)
    env["PATH"] = os.pathsep.join([str(snakemake.parent), env.get("PATH", "")])
    env[DEV_ENV] = "1" if args.dev else "0"
    extra = [a for a in getattr(args, "snakemake_args", []) if a != "--"]
    cmd = snakemake_command(snakemake, config, targets, args.cores, getattr(args, "dry_run", False), extra)
    return subprocess.run(cmd, env=env, check=False).returncode


def project_name(folder: Path) -> str:
    name = re.sub(r"[^a-z0-9_]+", "_", folder.name.lower()).strip("_")
    return name or "project"


def init_project(project: Path, scenario: str, name: str | None, analysis: str, git: bool) -> list[Path]:
    """Create the standard layout from the scenario templates; never overwrite an existing file."""
    analysis_dir = project / "analysis" / analysis
    config = analysis_dir / "config" / "config.yaml"
    if config.exists():
        raise SystemExit(f"bulkrnaseq: {config} already exists; nothing was changed")
    text = (TEMPLATES / f"config_{scenario}.yaml").read_text()
    text = text.replace("{{name}}", name or project_name(project))
    text = text.replace("{{toolkit_version}}", package_version(REPO))
    for sub in ("config", "resources", "results"):
        (analysis_dir / sub).mkdir(parents=True, exist_ok=True)
    (project / "input").mkdir(exist_ok=True)
    config.write_text(text)
    created = [config]
    if SCENARIOS[scenario]:
        sheet = analysis_dir / "config" / "samples.tsv"
        if not sheet.exists():
            sheet.write_text("\t".join(SCENARIOS[scenario]) + "\n")
            created.append(sheet)
    if git:
        ignore = project / ".gitignore"
        if not ignore.exists():
            shutil.copyfile(TEMPLATES / "project.gitignore", ignore)
            created.append(ignore)
        if not (project / ".git").exists():
            subprocess.run(["git", "init", "-q", str(project)], check=True)
    return created


def cmd_init(args: argparse.Namespace) -> int:
    project = Path(args.project_dir).expanduser().resolve()
    created = init_project(project, args.scenario, args.name, args.analysis, args.git)
    for path in created:
        print(f"created {path}")
    config = created[0]
    print(f"\nnext: put the provider files in {project / 'input'}, edit {config.parent}, then")
    print(f"  bulkrnaseq fetch {config}\n  bulkrnaseq plan  {config}\n  bulkrnaseq run   {config}")
    return 0


def cmd_version(args: argparse.Namespace) -> int:
    state = git_state(REPO)
    dirty = {True: " (uncommitted changes)", False: "", None: " (not a git checkout)"}[state["git_dirty"]]
    print(f"bulkrnaseq {package_version(REPO)}  {state['git_describe'] or ''}{dirty}  [{REPO}]")
    return 0


def build_parser() -> argparse.ArgumentParser:
    parser = argparse.ArgumentParser(prog="bulkrnaseq", description=__doc__.split("\n\n")[0])
    sub = parser.add_subparsers(dest="command", required=True)

    p = sub.add_parser("init", help="create a project folder from a scenario template")
    p.add_argument("project_dir")
    p.add_argument(
        "--scenario",
        required=True,
        choices=sorted(SCENARIOS),
        help="star: STAR/StringTie files; excel: count matrices; de_only: provider DE tables",
    )
    p.add_argument("--name", help="project.name (default: folder name in snake_case)")
    p.add_argument("--analysis", default="pipeline_v1", help="analysis folder name (default: pipeline_v1)")
    p.add_argument("--git", action="store_true", help="also create a local git repository tracking config/")
    p.set_defaults(func=cmd_init)

    for name, target, text in [
        ("plan", ["plan"], "validate the config; show which modules would run and why others are skipped"),
        ("fetch", ["fetch_resources"], "download or copy the resources into the project cache"),
        ("run", [], "run the workflow"),
    ]:
        p = sub.add_parser(name, help=text)
        p.add_argument("config")
        p.add_argument(
            "--dev",
            action="store_true",
            help="development mode: allow uncommitted changes and a different toolkit.version",
        )
        p.add_argument("--cores", type=int, default=1)
        if name == "run":
            p.add_argument("-n", "--dry-run", action="store_true")
            p.add_argument(
                "snakemake_args",
                nargs=argparse.REMAINDER,
                help="after --: extra arguments passed to snakemake",
            )
        p.set_defaults(func=lambda a, t=target: run_workflow(a, t))

    p = sub.add_parser("version", help="toolkit version and git state")
    p.set_defaults(func=cmd_version)
    return parser


def main(argv: list[str] | None = None) -> int:
    args = build_parser().parse_args(argv)
    return args.func(args)
