"""Writes synthetic reference resources for scenario (a) into tests/fixtures/resources/.

- hpa_immune_cell.tsv.zip: rna_immune_cell_monaco.tsv layout (Gene, Gene name, Immune cell, TPM, pTPM, nTPM);
- hpa_single_cell.tsv.zip: rna_single_cell_type.tsv layout (Gene, Gene name, Cell type, nCPM);
- contaminant_sets.tsv.gz: MSigDB export layout with three synthetic cell-type sets.

Genes come from tests/fixtures/resources/annotation.gtf.gz; values are derived from the md5 of the gene ID, so
the files are deterministic and contain no real measurements. One gene in ten is absent from the references.
Run from the repository root after tests/fixtures/make_fixtures.R:
  python3 -I tests/fixtures/make_reference_fixture.py
"""

import gzip
import hashlib
import io
import re
import zipfile
from pathlib import Path

RES = Path("tests/fixtures/resources")
IMMUNE_CELLS = ["T-reg", "naive CD4 T-cell", "classical monocyte"]
SINGLE_CELLS = ["t-cells", "macrophages", "fibroblasts"]


def genes():
    out = []
    with gzip.open(RES / "annotation.gtf.gz", "rt") as fh:
        for line in fh:
            f = line.rstrip("\n").split("\t")
            if len(f) < 9 or f[2] != "gene":
                continue
            gid = re.search(r'gene_id "([^"]+)"', f[8]).group(1)
            name = re.search(r'gene_name "([^"]+)"', f[8])
            biotype = re.search(r'gene_biotype "([^"]+)"', f[8]).group(1)
            if biotype == "protein_coding" and f[0] not in ("MT", "X", "Y"):
                out.append((gid, name.group(1) if name else ""))
    return out


def value(gid, cell):
    """Deterministic expression in [0, 200), skewed towards low values."""
    h = int(hashlib.md5(f"{gid}|{cell}".encode()).hexdigest()[:8], 16) / 0xFFFFFFFF
    return round(200 * h**3, 1)


def absent(gid):
    return int(hashlib.md5(gid.encode()).hexdigest()[:2], 16) % 10 == 0


def write_zip(path, inner, text):
    info = zipfile.ZipInfo(inner, date_time=(2020, 1, 1, 0, 0, 0))
    info.compress_type = zipfile.ZIP_DEFLATED
    info.external_attr = 0o644 << 16
    with zipfile.ZipFile(path, "w") as zf:
        zf.writestr(info, text)


def write_gz(path, text):
    buf = io.BytesIO()
    with gzip.GzipFile(filename="", mode="wb", fileobj=buf, mtime=0) as gz:
        gz.write(text.encode())
    path.write_bytes(buf.getvalue())


def main():
    g = [x for x in genes() if not absent(x[0])]
    rows = ["Gene\tGene name\tImmune cell\tTPM\tpTPM\tnTPM"]
    for gid, name in g:
        for cell in IMMUNE_CELLS:
            v = value(gid, cell)
            rows.append(f"{gid}\t{name}\t{cell}\t{v}\t{v}\t{v}")
    write_zip(RES / "hpa_immune_cell.tsv.zip", "rna_immune_cell_monaco.tsv", "\n".join(rows) + "\n")
    rows = ["Gene\tGene name\tCell type\tnCPM"]
    for gid, name in g:
        for cell in SINGLE_CELLS:
            rows.append(f"{gid}\t{name}\t{cell}\t{value(gid, cell)}")
    write_zip(RES / "hpa_single_cell.tsv.zip", "rna_single_cell_type.tsv", "\n".join(rows) + "\n")
    # Overlapping sets: genes 0-29 myeloid, 20-49 stromal, 25-34 T/NK (some myeloid/stromal genes are kept).
    allg = genes()
    sets = {
        "SYNTHETIC_MYELOID_CELL": allg[0:30],
        "SYNTHETIC_STROMAL_CELL": allg[20:50],
        "SYNTHETIC_T_NK_CELL": allg[25:35],
    }
    rows = ["gs_collection\tgs_subcollection\tgs_name\tensembl_gene\tgene_symbol\tdb_version"]
    for name, members in sets.items():
        for gid, sym in members:
            rows.append(f"C8\t\t{name}\t{gid}\t{sym}\tsynthetic")
    write_gz(RES / "contaminant_sets.tsv.gz", "\n".join(rows) + "\n")


if __name__ == "__main__":
    main()
