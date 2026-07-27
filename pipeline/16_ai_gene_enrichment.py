#!/usr/bin/env python3
"""
16_ai_gene_enrichment.py — Mendelian AI-gene enrichment, with an explicit audit of
the ortholog-resolution defect that affected the previously reported result.

WHAT THIS SCRIPT DOES
---------------------
1. Resolves the curated AI panel (`ai_genes.PANEL`) against the DEG table by MGI
   MOUSE symbol, failing loudly if any panel gene cannot be located.
2. Recomputes the one-sided Fisher over-representation of panel genes among
   down- versus up-regulated significant genes, at three thresholds.
3. Reproduces the OLD, defective resolution (BioMart human-ortholog column) side
   by side, so the correction is auditable rather than asserted.

THE DEFECT (found 2026-07-27, pre-submission audit)
---------------------------------------------------
`02_deg_analysis.py` marked risk genes with

    s_risk = s.hs.isin(MENDELIAN)          # hs = hsapiens_homolog_associated_gene_name

`hs` is NaN for 6,209 of the 15,190 tested genes (40.9 %), including seven of the
twenty panel genes: AMBN, AMTN, MMP20, SLC4A4, LAMB3, LAMA3, WNT10A. Because
`isin` evaluates NaN as False, these were not excluded from the test — they were
counted into the NON-risk cells (b and d) of the 2x2 while being genuine AI genes,
all of them significantly DOWN-regulated. The error therefore biased the odds
ratio TOWARD the null: it understated the enrichment it was reporting.

All seven are present in the data under their mouse symbols with strong
statistics (e.g. Ambn padj = 0, Mmp20 padj = 0, Amtn log2FC = -1.09), so this is
an annotation gap in the BioMart ortholog column, not missing data.

Outputs
-------
    output/AI_enrichment_sensitivity.csv   2x2 counts + OR + P, old vs corrected
    output/AI_gene_panel_status.csv        per-gene resolution and statistics
    output/AI_enrichment_report.txt        human-readable audit trail

Run from the project root:
    python3 scripts/JDR_reanalysis/16_ai_gene_enrichment.py
"""

from __future__ import annotations

import os
import sys

import numpy as np
import pandas as pd
from scipy import stats

sys.path.insert(0, os.path.dirname(os.path.abspath(__file__)))
import ai_genes  # noqa: E402  (path shim must precede the import)

ROOT = os.path.dirname(os.path.dirname(os.path.dirname(os.path.abspath(__file__))))
OUT = os.path.join(ROOT, "scripts", "JDR_reanalysis", "output")

PADJ = 0.05
HS_COL = "hsapiens_homolog_associated_gene_name"
MOUSE_COL = "external_gene_name"

# Thresholds reported for sensitivity. The manuscript's headline uses 0.0 (no
# fold-change cutoff) because the core enamel genes carrying the signal have
# |log2FC| ~ 0.4; see the Methods and .wiki/system/lessons-learned.md.
THRESHOLDS = [
    (0.0, "padj<0.05, no LFC cutoff"),
    (0.25, "padj<0.05 & |log2FC|>0.25"),
    (1.0, "padj<0.05 & |log2FC|>1 (DEG rule)"),
]

_lines: list[str] = []


def log(msg: str = "") -> None:
    print(msg)
    _lines.append(msg)


def fisher_2x2(sig: pd.DataFrame, risk: pd.Series) -> tuple[int, int, int, int, float, float]:
    """One-sided Fisher: are panel genes over-represented among DOWN genes?"""
    r = risk.loc[sig.index]
    down = sig.log2FoldChange < 0
    a = int((down & r).sum())  # risk, down
    b = int((down & ~r).sum())  # non-risk, down
    c = int((~down & r).sum())  # risk, up
    d = int((~down & ~r).sum())  # non-risk, up
    odds, p = stats.fisher_exact([[a, b], [c, d]], alternative="greater")
    return a, b, c, d, odds, p


def main() -> None:
    df = pd.read_csv(os.path.join(OUT, "DEG_table.csv"), skiprows=1, low_memory=False)
    tested = df[df.padj.notna()].copy()

    log("=" * 78)
    log("Mendelian AI-gene enrichment — corrected resolution + defect audit")
    log("=" * 78)
    log(f"panel size            : {len(ai_genes.PANEL)} genes")
    log(f"genes tested by DESeq2: {len(tested)}")
    log(f"BioMart human-ortholog column empty for: {tested[HS_COL].isna().sum()} "
        f"({100 * tested[HS_COL].isna().mean():.1f} %)")
    log()

    log(f"panel source          : {ai_genes.PANEL_SOURCE}")
    manual = [g for g in ai_genes.PANEL if g.ortholog_note]
    if manual:
        log(f"manually resolved orthologs ({len(manual)}):")
        for g in manual:
            log(f"    {g.human} -> {g.mouse}")
    log()

    # ---- per-gene resolution status -------------------------------------
    rows = []
    for g in ai_genes.PANEL:
        sub = df[df[MOUSE_COL] == g.mouse]
        hit = sub.iloc[0] if len(sub) else None
        rows.append({
            "human": g.human,
            "mouse": g.mouse,
            "confidence": g.confidence,
            "moi": g.moi,
            "phenotypes": g.phenotypes,
            "in_table": hit is not None,
            "biomart_hs_annotated": bool(hit is not None and pd.notna(hit[HS_COL])),
            "log2FC": None if hit is None else float(hit.log2FoldChange),
            "padj": None if hit is None else (None if pd.isna(hit.padj) else float(hit.padj)),
        })
    status = pd.DataFrame(rows)
    status.to_csv(os.path.join(OUT, "AI_gene_panel_status.csv"), index=False)

    unannotated = status[~status.biomart_hs_annotated & status.in_table]
    log(f"panel genes present in data but MISSING BioMart human-ortholog annotation: "
        f"{len(unannotated)}")
    for _, r in unannotated.iterrows():
        log(f"    {r.human:8} ({r.mouse:8}) log2FC={r.log2FC:+.3f} padj={r.padj:.2e}")
    log("  ^ these were counted as NON-risk by the previous implementation.")
    log()

    # ---- enrichment: corrected vs as-published ---------------------------
    corrected = ai_genes.resolve_in_table(tested, mouse_col=MOUSE_COL, strict=True)
    as_published = tested[HS_COL].isin(ai_genes.human_symbols())

    out = []
    log(f"{'threshold':36} {'resolution':13} {'r_dn':>5}{'nr_dn':>7}{'r_up':>6}"
        f"{'nr_up':>7} {'OR':>8} {'P':>10}")
    log("-" * 100)
    for th, label in THRESHOLDS:
        sig = tested[(tested.padj < PADJ) & (tested.log2FoldChange.abs() > th)]
        for mask, tag in ((as_published, "as-published"), (corrected, "CORRECTED")):
            a, b, c, d, odds, p = fisher_2x2(sig, mask)
            log(f"{label:36} {tag:13} {a:5}{b:7}{c:6}{d:7} {odds:8.2f} {p:10.5f}")
            out.append({"threshold": label, "lfc_cutoff": th, "resolution": tag,
                        "risk_down": a, "nonrisk_down": b, "risk_up": c, "nonrisk_up": d,
                        "odds_ratio": odds, "p_one_sided": p})
        log("-" * 100)

    pd.DataFrame(out).to_csv(os.path.join(OUT, "AI_enrichment_sensitivity.csv"), index=False)

    log()
    # Do not assert a direction here — it depends on the panel in use. With the old
    # 20-gene hand panel the correction raised the OR; with the PanelApp panel it can
    # lower it, because the panel also contains up-regulated genes. Compute and report.
    log("Effect of the ortholog correction (computed, not assumed):")
    df_out = pd.DataFrame(out)
    for th, label in THRESHOLDS:
        sub = df_out[df_out.lfc_cutoff == th].set_index("resolution")
        try:
            old_or, new_or = sub.loc["as-published", "odds_ratio"], sub.loc["CORRECTED", "odds_ratio"]
            old_p, new_p = sub.loc["as-published", "p_one_sided"], sub.loc["CORRECTED", "p_one_sided"]
        except KeyError:
            continue
        direction = ("raises" if new_or > old_or else "lowers" if new_or < old_or else "leaves")
        log(f"  {label:36} OR {old_or:.2f} -> {new_or:.2f} ({direction}); "
            f"P {old_p:.5f} -> {new_p:.5f}")
    ups = [g for g in ai_genes.PANEL
           if len(tested[(tested[MOUSE_COL] == g.mouse) & (tested.padj < PADJ)
                         & (tested.log2FoldChange > 0)])]
    log(f"  panel genes significantly UP-regulated ({len(ups)}): "
        f"{', '.join(sorted(g.human for g in ups)) or 'none'}")

    with open(os.path.join(OUT, "AI_enrichment_report.txt"), "w") as fh:
        fh.write("\n".join(_lines) + "\n")
    log()
    log(f"wrote {OUT}/AI_enrichment_sensitivity.csv, AI_gene_panel_status.csv, "
        f"AI_enrichment_report.txt")


if __name__ == "__main__":
    main()
