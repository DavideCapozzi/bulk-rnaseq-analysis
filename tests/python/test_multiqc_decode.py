import importlib.util
import json
import random
from pathlib import Path

import pytest

REPO = Path(__file__).resolve().parents[2]
FIXTURE = REPO / "tests" / "fixtures" / "scenario_a" / "multiqc" / "multiqc_report.html"


def load(path, name):
    spec = importlib.util.spec_from_file_location(name, path)
    mod = importlib.util.module_from_spec(spec)
    spec.loader.exec_module(mod)
    return mod


decode = load(REPO / "inst" / "python" / "multiqc_decode.py", "multiqc_decode")
fixture = load(REPO / "tests" / "fixtures" / "make_multiqc_fixture.py", "make_multiqc_fixture")


@pytest.mark.parametrize(
    "text",
    ["a", "hello, hello, hello", '{"x": [1.5, 2, "é", "€"]}', "ab" * 5000],
)
def test_lzstring_round_trip(text):
    assert decode.lz_decompress_from_base64(fixture.lz_compress_to_base64(text)) == text


def test_lzstring_round_trip_random():
    rng = random.Random(1)
    for _ in range(50):
        text = "".join(chr(rng.choice([rng.randint(32, 126), rng.randint(160, 0x2FFF)])) for _ in range(300))
        assert decode.lz_decompress_from_base64(fixture.lz_compress_to_base64(text)) == text


def test_fixture_is_current():
    assert FIXTURE.read_text() == fixture.report_html(fixture.plots())


def test_decode_writes_plot_data_general_stats_and_index(tmp_path):
    decode.main(tmp_path, [("run", str(FIXTURE))])
    plots = json.loads((tmp_path / "run_plotdata.json").read_text())
    assert plots == fixture.plots()
    gs = (tmp_path / "run_general_stats.tsv").read_text().splitlines()
    assert gs[0] == "Sample Name\tReads"
    assert gs[-1] == "lib_C1_R1\t40000"
    index = (tmp_path / "plot_index.tsv").read_text().splitlines()
    assert index[0] == "report\tmultiqc_version\tplot_id\tplot_type\tsamples"
    assert index[1].startswith("run\t1.13\tstar_alignment_plot\tbar_graph\t")
    assert len(index) == 1 + len(plots)


def test_report_without_compressed_data_is_refused(tmp_path):
    report = tmp_path / "new.html"
    report.write_text("<html>MultiQC v1.25</html>")
    with pytest.raises(SystemExit, match="MultiQC version 1.25"):
        decode.main(tmp_path / "out", [("x", str(report))])
