#!/usr/bin/env python3
"""
01_extract_tables.py  --  Smad4 dcKO JDR re-analysis pipeline (step 1/7)

Purpose
-------
Robustly extract the DESeq2 DEG table and the per-sample TPM matrix from the two
project .xlsx workbooks into clean CSVs for downstream analysis.

Why this is not a one-line pandas.read_excel:
  STAR_DEGs_ResOut_01_dcKOvsCON.xlsx declares a malformed worksheet <dimension>
  (ref="A1"), which makes openpyxl's read-only path return an empty frame. We
  therefore parse the sheet XML directly (shared strings + cells), which is exact.

Inputs  (relative to project root /mnt/Disk1/ZXY_2026):
  00.RawData/Bulk-RNAseq/02.TableOutput/STAR_DEGs_ResOut_01_dcKOvsCON.xlsx
  00.RawData/Bulk-RNAseq/02.TableOutput/00.DEGs_Out_TPM.xlsx
Outputs:
  output/DEG_table.csv      (one row per gene; DESeq2 stats + human homolog)
  output/TPM_matrix.csv      (per-sample TPM/log2TPM; header row 1 = the real header)

Run:  python3 scripts/JDR_reanalysis/01_extract_tables.py
"""
import os, re, csv, zipfile
from xml.etree import ElementTree as ET

ROOT = os.path.dirname(os.path.dirname(os.path.dirname(os.path.abspath(__file__))))
OUT  = os.path.join(ROOT, "scripts", "JDR_reanalysis", "output")
os.makedirs(OUT, exist_ok=True)
NS = "{http://schemas.openxmlformats.org/spreadsheetml/2006/main}"


def _shared_strings(z):
    out = []
    try:
        root = ET.fromstring(z.read("xl/sharedStrings.xml"))
    except KeyError:
        return out
    for si in root.iter(f"{NS}si"):
        out.append("".join(t.text or "" for t in si.iter(f"{NS}t")))
    return out


def _first_sheet_target(z):
    wb = ET.fromstring(z.read("xl/workbook.xml"))
    sheet = next(wb.iter(f"{NS}sheet"))
    rid = sheet.get("{http://schemas.openxmlformats.org/officeDocument/2006/relationships}id")
    rels = ET.fromstring(z.read("xl/_rels/workbook.xml.rels"))
    target = {r.get("Id"): r.get("Target") for r in rels}[rid]
    return target if target.startswith("xl/") else "xl/" + target, sheet.get("name")


def _col_idx(cell_ref):
    letters = re.match(r"([A-Z]+)\d+", cell_ref).group(1)
    c = 0
    for ch in letters:
        c = c * 26 + (ord(ch) - 64)
    return c - 1


def extract(xlsx_path, out_csv):
    z = zipfile.ZipFile(xlsx_path)
    sst = _shared_strings(z)
    target, sheet_name = _first_sheet_target(z)
    root = ET.fromstring(z.read(target))
    sd = root.find(f"{NS}sheetData")
    rows = []
    for row in sd.findall(f"{NS}row"):
        cells, maxc = {}, 0
        for c in row.findall(f"{NS}c"):
            ci = _col_idx(c.get("r"))
            v = c.find(f"{NS}v")
            if v is None:
                val = ""
            elif c.get("t") == "s":
                val = sst[int(v.text)]
            else:
                val = v.text
            cells[ci] = val
            maxc = max(maxc, ci)
        rows.append([cells.get(i, "") for i in range(maxc + 1)])
    with open(out_csv, "w", newline="") as fh:
        w = csv.writer(fh)
        for r in rows:
            w.writerow(r)
    print(f"  {os.path.basename(xlsx_path)}  sheet='{sheet_name}'  rows={len(rows)} -> {os.path.basename(out_csv)}")
    return len(rows)


if __name__ == "__main__":
    td = os.path.join(ROOT, "00.RawData", "Bulk-RNAseq", "02.TableOutput")
    print("Extracting workbooks:")
    extract(os.path.join(td, "STAR_DEGs_ResOut_01_dcKOvsCON.xlsx"),
            os.path.join(OUT, "DEG_table.csv"))
    extract(os.path.join(td, "00.DEGs_Out_TPM.xlsx"),
            os.path.join(OUT, "TPM_matrix.csv"))
    print("Done. NOTE: both CSVs have a grouping row 0 ('Basic Information'); the real header is row 1.")
