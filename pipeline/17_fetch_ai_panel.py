#!/usr/bin/env python3
"""
17_fetch_ai_panel.py — Build the Mendelian amelogenesis-imperfecta (AI) gene panel
from a versioned, citable external source, and map it to mouse orthologs.

WHY THIS REPLACES THE OLD PANEL
-------------------------------
The previous panel lived in `scripts/00_Mendelian_AI_Genes.R` as a hand-typed table
described in the Methods as "curated from OMIM/ClinVar". A pre-submission audit on
2026-07-27 established that this description was not supportable:

  * 18 of 20 PMIDs in its provenance column cite unrelated papers (Atlantic salmon
    population genomics, knee arthroplasty, circadian clocks, congenital heart
    surgery nomenclature); one PMID does not exist in PubMed at all. Verified twice
    independently via NCBI E-utilities.
  * 10 of 20 OMIM IDs are wrong: five are GENE MIMs where PHENOTYPE MIMs are
    required; two point at a different gene entirely (612559 = chronic lymphocytic
    leukaemia susceptibility, not STIM1; 615887 = SLC24A4, not ORAI1); two label the
    right gene with the wrong disease; one is missing.
  * The STIM1/ORAI1 entries pair the enamel phenotype with the dominant
    gain-of-function syndromes (Stormorken, York platelet). The enamel phenotype
    belongs to the RECESSIVE loss-of-function CRAC channelopathy.
  * There was no inclusion criterion a reader could reproduce, and 13 established AI
    genes were absent.

That file is quarantined (see scripts/_quarantine/). Nothing derived from it should
be reused. See .wiki/system/lessons-learned.md for the incident record.

THE REPLACEMENT CRITERION (state this verbatim in the Methods)
--------------------------------------------------------------
Genes are taken from the Genomics England PanelApp "Amelogenesis imperfecta" panel
(panel 269), which is expert-curated, versioned, dated, and publicly retrievable, so
a reader can reproduce the exact gene set. Only GREEN ("high evidence",
confidence_level 3) genes are used — PanelApp's diagnostic-grade rating. AMBER and
RED genes are recorded but excluded, and the exclusion is reported.

The panel VERSION IS PINNED below. Re-running this script reproduces the identical
set; bumping the version is a deliberate, visible edit.

Human -> mouse orthologs come from the MGI vertebrate homology report
(HOM_MouseHumanSequence.rpt), not from string re-casing and not from the BioMart
column in the DEG table (which is empty for ~41 % of tested genes — the defect that
corrupted the previous enrichment; see 16_ai_gene_enrichment.py).

FAILS LOUDLY: live retrieval errors and unresolved GREEN-gene orthologs raise.
Offline reproduction explicitly reads original source records declared in the
data manifest. It never falls back implicitly or substitutes a hand-typed panel.
That rule exists because a hard-coded fallback is what produced the retired
fabricated Figure 7 (.wiki/system/lessons-learned.md §1).

Outputs
-------
    output/AI_panel_PanelApp_v<version>.csv   the pinned, mapped panel
    output/AI_panel_provenance.txt            retrieval record for the Methods

Run via reproduce.py --mode panel, or request --download explicitly.
"""

from __future__ import annotations

import csv
import argparse
import hashlib
import io
import json
import os
import sys
import urllib.request
from pathlib import Path
from collections import defaultdict

ROOT = os.path.dirname(os.path.dirname(os.path.dirname(os.path.abspath(__file__))))
OUT = os.path.join(ROOT, "scripts", "JDR_reanalysis", "output")

PANEL_ID = 269
PANEL_NAME = "Amelogenesis imperfecta"
# Pinned. Changing this changes the reported gene set — do so deliberately.
PANEL_VERSION = "4.37"

PANELAPP = f"https://panelapp.genomicsengland.co.uk/api/v1/panels/{PANEL_ID}/?version={PANEL_VERSION}"
MGI_HOM = "https://www.informatics.jax.org/downloads/reports/HOM_MouseHumanSequence.rpt"

GREEN = 3  # PanelApp confidence_level: 3 = green/high, 2 = amber, 1/0 = red

# ---------------------------------------------------------------------------
# Explicit, justified resolutions for genes MGI's 1:1 filter cannot settle.
# Every entry is a NOMENCLATURE-VERSION mismatch between three sources that each
# moved at a different time: PanelApp (older HGNC symbols), MGI (current HGNC),
# and our DEG table (Ensembl release 116). None is a biological ambiguity.
# Each was checked by hand against MGI and HGNC on 2026-07-27.
# ---------------------------------------------------------------------------
MANUAL_ORTHOLOGS: dict[str, tuple[str, str]] = {
    # human symbol : (mouse symbol, justification)
    "AMELX": (
        "Amelx",
        "MGI homology class 51822737 contains human AMELX and AMELY against a single "
        "mouse gene, so the 1:1 filter rejects it. AMELY is the Y-linked paralogue; "
        "the mouse ortholog of AMELX is unambiguously Amelx.",
    ),
    "C4orf26": (
        "Odaph",
        "PanelApp still lists the retired symbol C4orf26. HGNC renamed it ODAPH "
        "(odontogenesis-associated phosphoprotein); MGI class 51825324 maps ODAPH to "
        "mouse Odaph, which is the symbol present in our annotation.",
    ),
    "FAM83H": (
        "Fam83h",
        "HGNC renamed FAM83H to SACK1H (HGNC:24797) on 2026-02-09, and MGI now lists "
        "the mouse gene as Sack1h. Our DEG table was annotated with Ensembl release "
        "116, which predates the rename and carries Fam83h. We therefore match on "
        "Fam83h and report the gene as SACK1H (formerly FAM83H) in the manuscript.",
    ),
}


class SourceUnavailable(RuntimeError):
    """A required external source could not be retrieved. Never fall back."""


def fetch(url: str, what: str) -> bytes:
    # Both PanelApp and MGI reject the default urllib agent with HTTP 403.
    req = urllib.request.Request(
        url, headers={"User-Agent": "Smad4-dcKO-reanalysis/1.0 (academic; contact via manuscript)",
                      "Accept": "application/json, text/plain, */*"}
    )
    try:
        with urllib.request.urlopen(req, timeout=120) as fh:
            return fh.read()
    except Exception as exc:  # noqa: BLE001 - we intentionally convert everything
        raise SourceUnavailable(
            f"Could not retrieve {what} from {url}: {exc}. "
            "Refusing to substitute a stored copy — rerun when the source is reachable."
        ) from exc


def load_panel(source: bytes | None = None) -> dict:
    data = json.loads((source if source is not None else fetch(PANELAPP, "PanelApp panel")).decode())
    got = str(data.get("version"))
    if got != PANEL_VERSION:
        raise SourceUnavailable(
            f"PanelApp {PANEL_NAME} is now v{got}, but this script pins v{PANEL_VERSION}. "
            "The gene set has changed upstream. Review the diff, then update "
            "PANEL_VERSION and the Methods together — do not silently accept the new set."
        )
    return data


def load_mgi_orthologs(source: bytes | None = None) -> dict[str, str]:
    """Return {human_symbol: mouse_symbol} from the MGI homology classes."""
    txt = (source if source is not None else fetch(MGI_HOM, "MGI homology report")).decode("utf-8", "replace")
    rows = list(csv.DictReader(io.StringIO(txt), delimiter="\t"))
    by_class: dict[str, dict[str, list[str]]] = defaultdict(lambda: defaultdict(list))
    for r in rows:
        key = r.get("DB Class Key", "").strip()
        org = r.get("Common Organism Name", "")
        sym = r.get("Symbol", "").strip()
        if not key or not sym:
            continue
        if org.startswith("human"):
            by_class[key]["human"].append(sym)
        elif org.startswith("mouse"):
            by_class[key]["mouse"].append(sym)
    out: dict[str, str] = {}
    for members in by_class.values():
        # Only accept unambiguous 1:1 classes; anything else is reported, not guessed.
        if len(members["human"]) == 1 and len(members["mouse"]) == 1:
            out[members["human"][0]] = members["mouse"][0]
    return out


def main() -> None:
    ap = argparse.ArgumentParser(description="Build the pinned panel from original source records; no implicit fallback.")
    sources = ap.add_mutually_exclusive_group(required=True)
    sources.add_argument('--sources-dir', type=Path, help='Explicit archived PanelApp JSON and MGI report')
    sources.add_argument('--download', action='store_true', help='Fetch the pinned API version and current MGI report')
    ap.add_argument('--save-sources', type=Path, help='Save original fetched records when using --download')
    ap.add_argument('--output-dir', type=Path, default=Path(__file__).resolve().parent/'output')
    args = ap.parse_args()
    if args.sources_dir:
        panel_raw = (args.sources_dir/'PanelApp_269_v4.37.json').read_bytes()
        mgi_raw = (args.sources_dir/'HOM_MouseHumanSequence.rpt').read_bytes()
    else:
        panel_raw = fetch(PANELAPP, "pinned PanelApp panel")
        mgi_raw = fetch(MGI_HOM, "MGI homology report")
    if args.save_sources:
        if not args.download: ap.error('--save-sources requires --download')
        args.save_sources.mkdir(parents=True, exist_ok=True)
        (args.save_sources/'PanelApp_269_v4.37.json').write_bytes(panel_raw)
        (args.save_sources/'HOM_MouseHumanSequence.rpt').write_bytes(mgi_raw)
    args.output_dir.mkdir(parents=True, exist_ok=True)
    global OUT
    OUT = str(args.output_dir)
    panel = load_panel(panel_raw)
    version = str(panel["version"])
    created = str(panel.get("version_created", ""))[:10]

    genes = panel["genes"]
    green, other = [], []
    for g in genes:
        sym = g["gene_data"]["gene_symbol"]
        rec = {
            "human": sym,
            "confidence": g.get("confidence_level"),
            "moi": (g.get("mode_of_inheritance") or "").strip(),
            "phenotypes": "; ".join(g.get("phenotypes") or [])[:300],
        }
        (green if int(g.get("confidence_level", 0)) == GREEN else other).append(rec)

    print(f"PanelApp panel {PANEL_ID} '{panel['name']}' v{version} ({created})")
    print(f"  total genes {len(genes)}   GREEN {len(green)}   excluded (amber/red) {len(other)}")

    ortho = load_mgi_orthologs(mgi_raw)
    print(f"  MGI 1:1 human->mouse ortholog classes: {len(ortho)}")

    unmapped = [r["human"] for r in green
                if r["human"] not in ortho and r["human"] not in MANUAL_ORTHOLOGS]
    if unmapped:
        raise SourceUnavailable(
            f"{len(unmapped)} GREEN panel genes have no unambiguous 1:1 mouse ortholog "
            f"in MGI: {unmapped}. Resolve each by hand against MGI and record the "
            "decision in MANUAL_ORTHOLOGS before proceeding — do not drop them silently."
        )
    for r in green:
        if r["human"] in MANUAL_ORTHOLOGS:
            r["mouse"], r["ortholog_note"] = MANUAL_ORTHOLOGS[r["human"]]
        else:
            r["mouse"], r["ortholog_note"] = ortho[r["human"]], ""
    if MANUAL_ORTHOLOGS:
        print(f"  manually resolved orthologs: "
              f"{', '.join(g for g in MANUAL_ORTHOLOGS if g in {r['human'] for r in green})}")

    path = os.path.join(OUT, f"AI_panel_PanelApp_v{version}.csv")
    with open(path, "w", newline="") as fh:
        w = csv.DictWriter(fh, fieldnames=["human", "mouse", "confidence", "moi", "phenotypes", "ortholog_note"])
        w.writeheader()
        for r in sorted(green, key=lambda x: x["human"]):
            w.writerow(r)

    prov = os.path.join(OUT, "AI_panel_provenance.txt")
    with open(prov, "w") as fh:
        fh.write(
            f"AI gene panel provenance\n"
            f"========================\n"
            f"Source        : Genomics England PanelApp, panel {PANEL_ID} '{panel['name']}'\n"
            f"Version       : v{version}  (panel version created {created})\n"
            f"URL           : {PANELAPP}\n"
            f"Inclusion     : confidence_level == 3 (GREEN / high evidence) only\n"
            f"GREEN genes   : {len(green)}\n"
            f"Excluded      : {len(other)} amber/red genes -> "
            f"{', '.join(sorted(r['human'] for r in other))}\n"
            f"Ortholog map  : MGI HOM_MouseHumanSequence.rpt, unambiguous 1:1 classes only\n"
            f"                {MGI_HOM}\n"
        )
    print(f"  wrote {path}")
    print(f"  wrote {prov}")
    (args.output_dir/'AI_panel_source_hashes.json').write_text(json.dumps({
        'panel_url':PANELAPP, 'panel_version':PANEL_VERSION, 'mgi_url':MGI_HOM,
        'source_mode':'explicit archived records' if args.sources_dir else 'live download',
        'panel_sha256':hashlib.sha256(panel_raw).hexdigest(),
        'mgi_sha256':hashlib.sha256(mgi_raw).hexdigest()},indent=2))
    print()
    print("  excluded (amber/red): " + ", ".join(sorted(r["human"] for r in other)))


if __name__ == "__main__":
    try:
        main()
    except SourceUnavailable as exc:
        print(f"\nFAILED LOUDLY: {exc}\n", file=sys.stderr)
        sys.exit(1)
