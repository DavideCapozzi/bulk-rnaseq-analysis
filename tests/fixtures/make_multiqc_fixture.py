"""Writes tests/fixtures/scenario_a/multiqc/multiqc_report.html: a small synthetic MultiQC 1.13-style report.

It holds the sections read_multiqc() uses (STAR, RSeQC, FastQC on raw R1) for the scenario (a) libraries plus
a QC-only library (lib_C1, no alignment), LZString-compressed as MultiQC does, and a General Statistics table.
Gene body coverage names the samples "A1", ... (mapped with inputs.qc.sample_aliases).
Run from the repository root before tests/fixtures/make_fixtures.R (which lists it in md5sum.txt):
  python3 -I tests/fixtures/make_multiqc_fixture.py
"""

import json
from pathlib import Path

KEY_B64 = "ABCDEFGHIJKLMNOPQRSTUVWXYZabcdefghijklmnopqrstuvwxyz0123456789+/="


def lz_compress_to_base64(text):
    """Port of LZString.compressToBase64 (pieroxy/lz-string); test-only counterpart of the decoder."""
    dictionary, to_create = {}, set()
    w, enlarge_in, dict_size, num_bits = "", 2, 3, 2
    out, val, pos = [], 0, 0

    def write_bits(n, value):
        nonlocal val, pos
        for _ in range(n):
            val = (val << 1) | (value & 1)
            value >>= 1
            if pos == 5:
                out.append(KEY_B64[val])
                val, pos = 0, 0
            else:
                pos += 1

    def emit(w):
        nonlocal enlarge_in, num_bits
        if w in to_create:
            code = ord(w[0])
            if code < 256:
                write_bits(num_bits, 0)
                write_bits(8, code)
            else:
                write_bits(num_bits, 1)
                write_bits(16, code)
            enlarge_in -= 1
            if enlarge_in == 0:
                enlarge_in, num_bits = 1 << num_bits, num_bits + 1
            to_create.discard(w)
        else:
            write_bits(num_bits, dictionary[w])
        enlarge_in -= 1
        if enlarge_in == 0:
            enlarge_in, num_bits = 1 << num_bits, num_bits + 1

    for c in text:
        if c not in dictionary:
            dictionary[c] = dict_size
            dict_size += 1
            to_create.add(c)
        if w + c in dictionary:
            w += c
        else:
            emit(w)
            dictionary[w + c] = dict_size
            dict_size += 1
            w = c
    if w:
        emit(w)
    write_bits(num_bits, 2)
    while True:
        val <<= 1
        if pos == 5:
            out.append(KEY_B64[val])
            break
        pos += 1
    res = "".join(out)
    return res + "=" * (-len(res) % 4)


ALIGNED = ["lib_A1", "lib_A2", "lib_B1", "lib_B2"]
RAW = ALIGNED + ["lib_C1"]


def bar(samples, cats, datasets=1, labels=None):
    plot = {
        "plot_type": "bar_graph",
        "samples": [samples] * datasets,
        "datasets": [[{"name": k, "data": v} for k, v in cats(d).items()] for d in range(datasets)],
        "config": {},
    }
    if labels:
        plot["config"]["data_labels"] = labels
    return plot


def line(series, config=None):
    return {"plot_type": "xy_line", "datasets": [series], "config": config or {}}


def plots():
    n = len(ALIGNED)
    star = {
        "Uniquely mapped": [6.0e5 + 2e4 * i for i in range(n)],
        "Mapped to multiple loci": [9.0e4 + 1e3 * i for i in range(n)],
        "Mapped to too many loci": [5.0e3 + 1e2 * i for i in range(n)],
        "Unmapped: too short": [2.5e5 - 1e4 * i for i in range(n)],
        "Unmapped: other": [4.0e3 + 50 * i for i in range(n)],
    }

    def gene_counts(d):
        unstranded, same, reverse = (
            [4.0e5, 4.2e5, 4.4e5, 4.6e5],
            [9e3, 9.5e3, 1e4, 1.05e4],
            [3.9e5, 4.1e5, 4.3e5, 4.5e5],
        )
        assigned = [unstranded, same, reverse][d]
        return {
            "Overlapping Genes": assigned,
            "No Feature": [5e4] * n,
            "Ambiguous Features": [3e3] * n,
            "Multimapping": [9e4] * n,
            "Unmapped": [2.5e5] * n,
        }

    read_dist = {
        "CDS_Exons": [3.0e5 + 1e4 * i for i in range(n)],
        "5'UTR_Exons": [2.0e4] * n,
        "3'UTR_Exons": [6.0e4] * n,
        "Introns": [8.0e4 - 5e3 * i for i in range(n)],
        "TSS_up_1kb": [2.0e3] * n,
        "TES_down_1kb": [3.0e3] * n,
        "Other_intergenic": [1.0e4] * n,
    }
    gene_body = [
        {
            "name": s.replace("lib_", ""),
            "data": [[x, round(0.4 + 0.006 * x * (1 + 0.1 * i), 6)] for x in range(1, 101)],
        }
        for i, s in enumerate(ALIGNED)
    ]
    adapter = [
        {
            "name": f"{s}_R1 - illumina_universal_adapter",
            "data": [[x, round(0.02 * x * (1 + i), 6)] for x in range(1, 41)],
        }
        for i, s in enumerate(RAW)
    ]
    seq_q = [
        {"name": f"{s}_R1", "data": [[q, 100 * (q - 20 + i) if q >= 20 else 0] for q in range(2, 42)]}
        for i, s in enumerate(RAW)
    ]
    base_q = [
        {
            "name": f"{s}_R1",
            "data": [[x, round(38 - 0.2 * x * (1 + 1.5 * (s == "lib_C1")), 6)] for x in range(1, 51)],
        }
        for s in RAW
    ]
    dup_cats = [float(k) for k in range(1, 10)] + [">10", ">50", ">100", ">500", ">1k", ">5k", ">10k+"]
    dup = [{"name": f"{s}_R1", "data": [60 - 5 * i] + [2.0] * 10 + [1.0 + i] * 5} for i, s in enumerate(RAW)]
    gc = [
        {"name": f"{s}_R1", "data": [[x, max(0.0, 3.0 - abs(x - 48 - i) * 0.1)] for x in range(0, 101)]}
        for i, s in enumerate(RAW)
    ]
    return {
        "star_alignment_plot": bar(ALIGNED, lambda d: star),
        "star_gene_counts": bar(ALIGNED, gene_counts, 3, ["Unstranded", "Same Stranded", "Reverse Stranded"]),
        "rseqc_infer_experiment_plot": bar(
            [f"{s}.strandness_statistic" for s in ALIGNED],
            lambda d: {
                "Sense": [5.0, 6.0, 7.0, 8.0],
                "Antisense": [93.0, 92.0, 91.0, 90.0],
                "Undetermined": [2.0] * n,
            },
        ),
        "rseqc_read_distribution_plot": bar(ALIGNED, lambda d: read_dist),
        "rseqc_gene_body_coverage_plot": line(gene_body),
        "fastqc_adapter_content_plot": line(adapter),
        "fastqc_per_sequence_quality_scores_plot": line(seq_q),
        "fastqc_per_base_sequence_quality_plot": line(base_q),
        "fastqc_sequence_duplication_levels_plot": line(dup, {"categories": dup_cats}),
        "fastqc_per_sequence_gc_content_plot": line(gc),
    }


def report_html(data):
    reads = [1_000_000 + 250_000 * i for i in range(len(ALIGNED))] + [40_000]  # lib_C1: failed library
    rows = "".join(f"<tr><th>{s}_R1</th><td>{r}</td></tr>" for s, r in zip(RAW, reads, strict=True))
    return (
        "<!DOCTYPE html>\n<html><head><title>MultiQC Report</title></head><body>\n"
        "<p>Report generated by MultiQC v1.13 (synthetic test fixture)</p>\n"
        '<table id="general_stats_table"><thead><tr><th>Sample Name</th><th>Reads</th></tr></thead>'
        f"<tbody>{rows}</tbody></table>\n"
        '<script type="text/javascript" id="mqc_compressed_plotdata">'
        f"{lz_compress_to_base64(json.dumps(data))}</script>\n</body></html>\n"
    )


if __name__ == "__main__":
    out = Path("tests/fixtures/scenario_a/multiqc/multiqc_report.html")
    out.parent.mkdir(parents=True, exist_ok=True)
    out.write_text(report_html(plots()))
