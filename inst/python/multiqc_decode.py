"""Decode MultiQC HTML reports: plot data, General Statistics table and an index of the plots.

MultiQC 1.x reports (tested with 1.12 and 1.13) embed the data of their plots as LZString-compressed JSON
(base64) in <script id="mqc_compressed_plotdata">; the General Statistics table is plain HTML. No R package
decodes LZString, so this script does it with the Python standard library only.

For every report it writes <id>_plotdata.json and <id>_general_stats.tsv, and plot_index.tsv for all reports.

Usage: python3 -I multiqc_decode.py <out_dir> <id>=<report.html> [<id>=<report.html> ...]
"""

import html.parser
import json
import re
import sys
from pathlib import Path

KEY_B64 = "ABCDEFGHIJKLMNOPQRSTUVWXYZabcdefghijklmnopqrstuvwxyz0123456789+/="
B64_VAL = {c: i for i, c in enumerate(KEY_B64)}


def lz_decompress_from_base64(s):
    """Port of LZString.decompressFromBase64 (pieroxy/lz-string)."""
    s = s.strip()
    if not s:
        return ""
    length, reset = len(s), 32
    dictionary = {0: 0, 1: 1, 2: 2}
    enlarge_in, dict_size, num_bits = 4, 4, 3
    result = []
    data_val, data_pos, data_idx = B64_VAL[s[0]], reset, 1

    def read_bits(n):
        nonlocal data_val, data_pos, data_idx
        bits, power, maxpower = 0, 1, 1 << n
        while power != maxpower:
            resb = data_val & data_pos
            data_pos >>= 1
            if data_pos == 0:
                data_pos = reset
                data_val = B64_VAL[s[data_idx]] if data_idx < length else 0
                data_idx += 1
            bits |= (1 if resb > 0 else 0) * power
            power <<= 1
        return bits

    nxt = read_bits(2)
    if nxt == 0:
        c = chr(read_bits(8))
    elif nxt == 1:
        c = chr(read_bits(16))
    else:
        return ""
    dictionary[3] = c
    w = c
    result.append(c)
    while True:
        if data_idx > length:
            return ""
        c = read_bits(num_bits)
        if c == 0:
            dictionary[dict_size] = chr(read_bits(8))
            dict_size += 1
            c = dict_size - 1
            enlarge_in -= 1
        elif c == 1:
            dictionary[dict_size] = chr(read_bits(16))
            dict_size += 1
            c = dict_size - 1
            enlarge_in -= 1
        elif c == 2:
            return "".join(result)
        if enlarge_in == 0:
            enlarge_in = 1 << num_bits
            num_bits += 1
        if c in dictionary:
            entry = dictionary[c]
        elif c == dict_size:
            entry = w + w[0]
        else:
            raise ValueError("LZString: corrupted stream")
        result.append(entry)
        dictionary[dict_size] = w + entry[0]
        dict_size += 1
        enlarge_in -= 1
        w = entry
        if enlarge_in == 0:
            enlarge_in = 1 << num_bits
            num_bits += 1


class GeneralStatsParser(html.parser.HTMLParser):
    """Header and rows of <table id="general_stats_table">, cell text with whitespace collapsed."""

    def __init__(self):
        super().__init__()
        self.in_table = False
        self.in_thead = False
        self.cell = None
        self.row = None
        self.rows = []
        self.header = []

    def handle_starttag(self, tag, attrs):
        if tag == "table" and dict(attrs).get("id") == "general_stats_table":
            self.in_table = True
        if not self.in_table:
            return
        if tag == "thead":
            self.in_thead = True
        if tag == "tr":
            self.row = []
        if tag in ("td", "th"):
            self.cell = []

    def handle_endtag(self, tag):
        if not self.in_table:
            return
        if tag in ("td", "th") and self.cell is not None:
            self.row.append(" ".join("".join(self.cell).split()))
            self.cell = None
        if tag == "tr" and self.row is not None:
            if self.in_thead:
                self.header = self.row
            else:
                self.rows.append(self.row)
            self.row = None
        if tag == "thead":
            self.in_thead = False
        if tag == "table":
            self.in_table = False

    def handle_data(self, data):
        if self.in_table and self.cell is not None:
            self.cell.append(data)


def decode_report(path):
    """Return (MultiQC version, plot data dict, General Statistics header, rows) of one HTML report."""
    txt = Path(path).read_text(encoding="utf-8", errors="replace")
    ver = re.search(r"MultiQC v([0-9][0-9a-z.]*)", txt)
    version = ver.group(1) if ver else "NA"
    m = re.search(r'id="mqc_compressed_plotdata"[^>]*>(.*?)</', txt, re.S)
    if m is None:
        raise SystemExit(
            f"{path}: no compressed plot data (MultiQC version {version}); supported: MultiQC 1.x HTML "
            'reports with <script id="mqc_compressed_plotdata"> (tested with 1.12 and 1.13)'
        )
    plots = json.loads(lz_decompress_from_base64(m.group(1)))
    gs = GeneralStatsParser()
    gs.feed(txt)
    return version, plots, gs.header, gs.rows


def plot_samples(plot):
    """Sample names of a plot, per dataset (bar plots list them; line plots name each series)."""
    return plot.get("samples") or [
        [d.get("name") for d in ds] if isinstance(ds, list) else None for ds in plot.get("datasets", [])
    ]


def main(out, reports):
    out.mkdir(parents=True, exist_ok=True)
    index = []
    for rid, path in reports:
        version, plots, header, rows = decode_report(path)
        (out / f"{rid}_plotdata.json").write_text(json.dumps(plots))
        for pid, p in plots.items():
            index.append((rid, version, pid, p.get("plot_type", "NA"), json.dumps(plot_samples(p))[:300]))
        with open(out / f"{rid}_general_stats.tsv", "w") as fh:
            fh.write("\t".join(header) + "\n")
            for r in rows:
                fh.write("\t".join(r) + "\n")
    with open(out / "plot_index.tsv", "w") as fh:
        fh.write("report\tmultiqc_version\tplot_id\tplot_type\tsamples\n")
        for r in index:
            fh.write("\t".join(r) + "\n")


if __name__ == "__main__":
    if len(sys.argv) < 3:
        raise SystemExit(__doc__)
    main(Path(sys.argv[1]), [a.split("=", 1) for a in sys.argv[2:]])
