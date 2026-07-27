"""
ai_genes.py — SINGLE SOURCE OF TRUTH for the Mendelian amelogenesis-imperfecta (AI)
gene panel used in the enrichment analysis.

The panel is NOT defined in this file. It is loaded from the versioned CSV that
`17_fetch_ai_panel.py` writes after retrieving the Genomics England PanelApp
"Amelogenesis imperfecta" panel (panel 269, GREEN genes only) and mapping it to mouse
orthologs via the MGI homology report. This module is the import point; that script
is the provenance.

WHY IT WORKS THIS WAY
---------------------
The previous panel was a hand-typed table duplicated across an R script and a Python
`set`. A pre-submission audit on 2026-07-27 found its provenance column to be largely
fabricated (18/20 PMIDs cite unrelated papers; 10/20 OMIM IDs wrong, two of them
pointing at entirely different genes) and found no inclusion criterion a reader could
reproduce. It is quarantined under `scripts/_quarantine/`.

Two rules this module enforces, both from .wiki/system/lessons-learned.md:

  §4 One definition, imported everywhere. Nothing re-types this list.
  §1 Fail loudly. `resolve_in_table()` RAISES when a panel gene cannot be located,
     because the previous implementation matched on the DEG table's BioMart
     human-ortholog column, which is NaN for ~41 % of tested genes. Pandas'
     `.isin()` treats NaN as a non-match, so seven genuine AI genes — all strongly
     down-regulated — were silently counted into the NON-risk cells of the 2x2. That
     understated the odds ratio from 16.4 to 9.2. A missing gene must stop the
     analysis, never quietly change the answer.
"""

from __future__ import annotations

import csv
import glob
import os
from dataclasses import dataclass

_OUT = os.path.join(os.path.dirname(os.path.abspath(__file__)), "output")


@dataclass(frozen=True)
class AIGene:
    """One PanelApp GREEN amelogenesis-imperfecta gene, mapped to its mouse ortholog."""

    human: str          # HGNC symbol as listed by PanelApp
    mouse: str          # MGI symbol, as it appears in our Ensembl-116 annotation
    confidence: int     # PanelApp confidence_level (3 = GREEN throughout)
    moi: str            # PanelApp mode of inheritance
    phenotypes: str     # PanelApp phenotype strings
    ortholog_note: str  # non-empty where the ortholog needed manual resolution


class PanelNotBuilt(RuntimeError):
    """The versioned panel CSV is absent — run 17_fetch_ai_panel.py first."""


class OrthologResolutionError(RuntimeError):
    """A panel gene could not be located in the expression table."""


def _load() -> tuple[tuple[AIGene, ...], str]:
    hits = sorted(glob.glob(os.path.join(_OUT, "AI_panel_PanelApp_v*.csv")))
    if not hits:
        raise PanelNotBuilt(
            "No AI_panel_PanelApp_v*.csv in scripts/JDR_reanalysis/output/. "
            "Run `python3 scripts/JDR_reanalysis/17_fetch_ai_panel.py` to retrieve the "
            "pinned PanelApp panel. This module will not fall back to a stored list."
        )
    if len(hits) > 1:
        raise PanelNotBuilt(
            f"Multiple panel versions present: {[os.path.basename(h) for h in hits]}. "
            "Exactly one must be in play so the reported gene set is unambiguous. "
            "Remove the stale file(s)."
        )
    path = hits[0]
    with open(path, newline="") as fh:
        rows = list(csv.DictReader(fh))
    panel = tuple(
        AIGene(
            human=r["human"],
            mouse=r["mouse"],
            confidence=int(r["confidence"]),
            moi=r.get("moi", ""),
            phenotypes=r.get("phenotypes", ""),
            ortholog_note=r.get("ortholog_note", ""),
        )
        for r in rows
    )
    version = os.path.basename(path).replace("AI_panel_PanelApp_v", "").replace(".csv", "")
    return panel, version


PANEL, PANEL_VERSION = _load()
PANEL_SOURCE = (
    f"Genomics England PanelApp, panel 269 'Amelogenesis imperfecta' v{PANEL_VERSION}; "
    f"GREEN (high-evidence) genes only; mouse orthologs from MGI HOM_MouseHumanSequence"
)


def human_symbols() -> set[str]:
    return {g.human for g in PANEL}


def mouse_symbols() -> set[str]:
    return {g.mouse for g in PANEL}


def human_to_mouse() -> dict[str, str]:
    return {g.human: g.mouse for g in PANEL}


def resolve_in_table(df, mouse_col: str = "external_gene_name", strict: bool = True):
    """Boolean mask over `df` marking panel genes, matched on MGI mouse symbol.

    Matching deliberately does NOT use the BioMart human-ortholog column: it is NaN
    for ~41 % of tested genes, and that is the defect this module exists to prevent.

    With `strict=True` a panel gene missing from `df` raises. A gene that is merely
    not expressed still appears with padj = NaN, so a genuine absence means the panel
    and the annotation build have diverged — which invalidates the contingency table
    rather than merely shrinking it.
    """
    present = set(df[mouse_col].dropna().astype(str))
    missing = sorted(m for m in mouse_symbols() if m not in present)
    if missing and strict:
        raise OrthologResolutionError(
            f"{len(missing)} of {len(PANEL)} panel genes are absent from the expression "
            f"table under their MGI symbols: {missing}. Refusing to compute an "
            "enrichment on a silently truncated panel. If these genes are genuinely "
            "not in the annotation, record that fact explicitly and pass strict=False."
        )
    return df[mouse_col].isin(mouse_symbols())


__all__ = [
    "AIGene", "PANEL", "PANEL_VERSION", "PANEL_SOURCE",
    "human_symbols", "mouse_symbols", "human_to_mouse", "resolve_in_table",
    "OrthologResolutionError", "PanelNotBuilt",
]
