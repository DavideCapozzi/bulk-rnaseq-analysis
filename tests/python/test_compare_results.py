import json
import sys
from pathlib import Path

REPO = Path(__file__).resolve().parents[2]
sys.path.insert(0, str(REPO / "tools"))
from compare_results import compare, main  # noqa: E402


def make_run(root: Path, value: str, commit: str, seconds: str) -> Path:
    prov = root / "00_provenance"
    (prov / "benchmarks").mkdir(parents=True)
    (root / "02_qc").mkdir()
    (root / "02_qc" / "qc.tsv").write_text(value)
    (prov / "run_timestamps.json").write_text(seconds)
    (prov / "benchmarks" / "02_qc.tsv").write_text(seconds)
    (prov / "toolkit.json").write_text(commit)
    (prov / "run_manifest.json").write_text(json.dumps({"config_md5": "x", "toolkit": {"commit": commit}}))
    return root


def test_timing_files_are_ignored(tmp_path):
    a = make_run(tmp_path / "a", "1", "c1", "10")
    b = make_run(tmp_path / "b", "1", "c1", "99")
    assert compare(a, b) == []


def test_code_changed_ignores_toolkit_identity_only(tmp_path):
    a = make_run(tmp_path / "a", "1", "c1", "10")
    b = make_run(tmp_path / "b", "1", "c2", "10")
    assert compare(a, b) == [
        "changed  00_provenance/run_manifest.json",
        "changed  00_provenance/toolkit.json",
    ]
    assert compare(a, b, code_changed=True) == []
    (b / "02_qc" / "qc.tsv").write_text("2")
    (b / "02_qc" / "new.tsv").write_text("")
    assert compare(a, b, code_changed=True) == ["new      02_qc/new.tsv", "changed  02_qc/qc.tsv"]


def test_manifest_change_outside_toolkit_is_reported(tmp_path):
    a = make_run(tmp_path / "a", "1", "c1", "10")
    b = make_run(tmp_path / "b", "1", "c2", "10")
    (b / "00_provenance" / "run_manifest.json").write_text(json.dumps({"config_md5": "y", "toolkit": {}}))
    assert compare(a, b, code_changed=True) == [
        "changed  00_provenance/run_manifest.json (outside 'toolkit': config_md5)"
    ]


def test_main_exit_codes(tmp_path, capsys):
    a = make_run(tmp_path / "a", "1", "c1", "10")
    b = make_run(tmp_path / "b", "1", "c1", "10")
    assert main([str(a), str(b)]) == 0
    assert "identical: 3 files" in capsys.readouterr().out
    (b / "02_qc" / "qc.tsv").unlink()
    assert main([str(a), str(b)]) == 1
    assert main([str(a)]) == 2
