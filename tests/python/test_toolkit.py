import subprocess

import pytest
from bulkrnaseq_toolkit import ToolkitError, check_toolkit, code_md5, dev_mode, package_version


def git(repo, *args):
    subprocess.run(
        ["git", "-C", str(repo), "-c", "user.name=t", "-c", "user.email=t@example.org", *args],
        check=True,
        capture_output=True,
    )


@pytest.fixture
def repo(tmp_path):
    (tmp_path / "DESCRIPTION").write_text("Package: x\nVersion: 1.2.3\n")
    (tmp_path / "R").mkdir()
    (tmp_path / "R" / "a.R").write_text("f <- function() 1\n")
    git(tmp_path, "init", "-q")
    git(tmp_path, "add", "-A")
    git(tmp_path, "commit", "-q", "-m", "init")
    return tmp_path


def test_package_version_reads_description(repo):
    assert package_version(repo) == "1.2.3"


def test_clean_checkout_with_matching_version_passes(repo):
    state, problems = check_toolkit("1.2.3", repo, dev=False)
    assert problems == []
    assert state["git_dirty"] is False and len(state["git_commit"]) == 40 and state["dev_mode"] is False


def test_dirty_checkout_stops_unless_dev(repo):
    (repo / "R" / "a.R").write_text("f <- function() 2\n")
    with pytest.raises(ToolkitError, match="uncommitted changes"):
        check_toolkit("1.2.3", repo, dev=False)
    state, problems = check_toolkit("1.2.3", repo, dev=True)
    assert state["git_dirty"] is True and len(problems) == 1


def test_untracked_file_counts_as_dirty(repo):
    (repo / "R" / "b.R").write_text("g <- 1\n")
    with pytest.raises(ToolkitError):
        check_toolkit("1.2.3", repo, dev=False)


def test_version_mismatch_names_both_versions(repo):
    with pytest.raises(ToolkitError, match=r"pins toolkit.version 1.0.0 but the checkout is 1.2.3"):
        check_toolkit("1.0.0", repo, dev=False)


def test_not_a_git_checkout_stops(tmp_path):
    (tmp_path / "DESCRIPTION").write_text("Version: 1.2.3\n")
    with pytest.raises(ToolkitError, match="not a git checkout"):
        check_toolkit("1.2.3", tmp_path, dev=False)


def test_code_md5_follows_code_and_ignores_caches(repo):
    before = code_md5(repo)
    (repo / "R" / "__pycache__").mkdir()
    (repo / "R" / "__pycache__" / "x.pyc").write_bytes(b"\0")
    assert code_md5(repo) == before
    (repo / "R" / "a.R").write_text("f <- function() 3\n")
    assert code_md5(repo) != before


@pytest.mark.parametrize("value,expected", [("1", True), ("true", True), ("0", False), ("", False)])
def test_dev_mode_parsing(value, expected):
    assert dev_mode({"BULKRNASEQ_DEV": value}) is expected
