import subprocess
from pathlib import Path

import pytest
import yaml
from bulkrnaseq_cli import REPO, SCENARIOS, find_snakemake, init_project, project_name, workdir_for
from bulkrnaseq_toolkit import package_version

snakemake_utils = pytest.importorskip("snakemake.utils")
SCHEMA = str(REPO / "config" / "schema.yaml")


@pytest.mark.parametrize("scenario", sorted(SCENARIOS))
def test_init_creates_a_valid_config_pinned_to_this_version(tmp_path, scenario):
    created = init_project(tmp_path / "My-Project", scenario, None, "pipeline_v1", git=False)
    config = tmp_path / "My-Project" / "analysis" / "pipeline_v1" / "config" / "config.yaml"
    assert created[0] == config
    cfg = yaml.safe_load(config.read_text())
    snakemake_utils.validate(cfg, SCHEMA)
    assert cfg["toolkit"]["version"] == package_version(REPO)
    assert cfg["project"]["name"] == "my_project"
    for sub in ("resources", "results"):
        assert (config.parent.parent / sub).is_dir()
    assert (tmp_path / "My-Project" / "input").is_dir()
    sheet = config.parent / "samples.tsv"
    if SCENARIOS[scenario] is None:
        assert not sheet.exists()
    else:
        assert sheet.read_text().rstrip("\n").split("\t")[:2] == ["sample_id", "input_column"]


def test_init_never_overwrites_a_config(tmp_path):
    init_project(tmp_path, "de_only", "p", "pipeline_v1", git=False)
    with pytest.raises(SystemExit, match="already exists"):
        init_project(tmp_path, "star", "p", "pipeline_v1", git=False)


def test_init_git_tracks_config_only(tmp_path):
    init_project(tmp_path, "star", "p", "pipeline_v1", git=True)
    analysis = Path("analysis", "pipeline_v1")
    for f in ("input/x.xlsx", analysis / "results" / "t.tsv", analysis / "resources" / "a.gtf.gz"):
        (tmp_path / f).parent.mkdir(parents=True, exist_ok=True)
        (tmp_path / f).write_text("x")

    def ignored(path):
        return subprocess.run(["git", "-C", str(tmp_path), "check-ignore", "-q", str(path)]).returncode == 0

    assert not ignored(analysis / "config" / "config.yaml")
    assert not ignored(".gitignore")
    assert ignored("input/x.xlsx")
    assert ignored(analysis / "results" / "t.tsv")
    assert ignored(analysis / "resources" / "a.gtf.gz")


def test_examples_pin_the_current_version():
    for f in sorted((REPO / "config" / "examples").glob("*.yaml")):
        cfg = yaml.safe_load(f.read_text())
        snakemake_utils.validate(cfg, SCHEMA)
        assert cfg["toolkit"]["version"] == package_version(REPO), f.name


def test_project_name_is_snake_case():
    assert project_name(Path("Tumour RNA-seq 2")) == "tumour_rna_seq_2"


def test_workdir_is_the_analysis_folder_in_the_standard_layout():
    assert workdir_for(Path("/p/analysis/v1/config/config.yaml")) == Path("/p/analysis/v1")
    assert workdir_for(Path("/p/other/config.yaml")) == Path("/p/other")


def test_find_snakemake_prefers_bulkrnaseq_env(tmp_path):
    exe = tmp_path / "bin" / "snakemake"
    exe.parent.mkdir()
    exe.write_text("#!/bin/sh\n")
    exe.chmod(0o755)
    assert find_snakemake({"BULKRNASEQ_ENV": str(tmp_path)}, "/nowhere/x", "/nowhere/python") == exe
    with pytest.raises(SystemExit, match="BULKRNASEQ_ENV is set"):
        find_snakemake({"BULKRNASEQ_ENV": str(tmp_path / "missing")}, "/nowhere/x", "/nowhere/python")
    assert find_snakemake({"PATH": ""}, str(tmp_path / "bin" / "bulkrnaseq"), "/nowhere/python") == exe
