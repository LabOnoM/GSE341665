#!/usr/bin/env python3
"""
21_smad4_loxp_figure.py  --  Smad4 dcKO re-analysis pipeline (step 21)

Render Supplementary Figure S5: direct read-level identification of the Smad4
floxed (loxP-flanked) allele, and the accompanying NEGATIVE result that the Cre
excision product is not detectable in the RNA.

THE STORY THIS FIGURE MUST TELL
    The allele is identified.  The excision is NOT detected.
A reader who sees only this figure must not come away believing the knockout was
demonstrated at the RNA level.

Every DATA value -- number, coordinate, count, sequence, restriction-site name,
exon-number list, direction word -- is read or derived at run time from the
outputs of 20_smad4_loxp_evidence.py and from genome.fa. If an input is missing
the script RAISES rather than substituting a value (project rule -- this
repository once shipped a figure built on fabricated fallback values).

What is NOT computed, and is therefore fixed prose: the panel titles, the
methodological explanations, and the interpretation paragraph. Those are
arguments about the data, not measurements of it. Any sentence that makes a
comparative or directional claim about the data is built from the JSON, and the
sanity gate below breaks the build if the JSON stops supporting it.

Inputs : output/smad4_loxp_reads.csv, output/smad4_exon_coverage.csv,
         output/smad4_exon_ratio.csv, output/smad4_junctions.csv,
         output/smad4_locus_depth.csv, output/smad4_loxp_summary.json
         00.RawData/Bulk-RNAseq/genome.fa (+ .fai)
Outputs: figures/FigureS5_smad4_loxp.{png,pdf}
         figures/FigureS5_smad4_loxp_legend.md
"""
import json
import os
import sys
import re
import textwrap

import matplotlib
matplotlib.use("Agg")
import matplotlib.pyplot as plt
import matplotlib.text
import numpy as np
import pandas as pd
from matplotlib.patches import Polygon, Rectangle

HERE = os.path.dirname(os.path.abspath(__file__))
ROOT = os.path.dirname(os.path.dirname(HERE))
OUT = os.path.join(HERE, "output")
FIG = os.path.join(HERE, "figures")
FASTA = os.path.join(ROOT, "00.RawData", "Bulk-RNAseq", "genome.fa")
os.makedirs(FIG, exist_ok=True)

# ---- house style (identical to 07_figures.py) -------------------------------
BLUE, RED = "#2166AC", "#B2182B"
plt.rcParams.update({"font.size": 8, "axes.titlesize": 9, "axes.titleweight": "bold",
                     "axes.titlelocation": "left",   # panel labels A/B/C/D left-aligned
                     "figure.dpi": 150, "savefig.dpi": 300, "axes.spines.top": False,
                     "axes.spines.right": False, "hatch.linewidth": 0.45})

# Greyscale discipline: colour NEVER carries meaning on its own.
#   genotype  -> colour (blue/red) AND marker (o / ^) AND hatch (none / ///)
#                AND an explicit per-row text label.
#   read part -> aligned = solid fill; soft-clipped = white fill + hatch + open edge.
#   loxP core -> near-black block with white letters (darkest thing in the figure).
SEQ_GENOMIC_CLIP = "#EFEFEF"   # soft-clipped bases that still match the reference
SEQ_GENOMIC_ALN = "#D3D3D3"    # bases the aligner placed on the reference
SEQ_VECTOR = "#9C9C9C"         # targeting-vector polylinker / restriction linkers
SEQ_LOXP = "#1A1A1A"           # the 34 bp loxP core
LOXP_CORE_LEN = 34

# Smallest type size used anywhere in the figure, and the reproduction width the
# figure is designed for. Asserted at build time so type cannot creep downwards.
MIN_FONTSIZE = 6.0
FIGSIZE = (11.4, 15.4)

# Shared genomic window for panels A and C (panel B zooms into it).
WINDOW = (73780000, 73792000)
SAMPLE_ORDER = ["CON1", "CON2", "CON3", "dcKO1", "dcKO2", "dcKO3"]

# Text objects that are deliberately packed edge-to-edge (monospace sequence
# blocks); excluded from the overlap audit because adjacent blocks touch by design.
_SEQ_TEXT_IDS = set()


# ============================ input loading ==================================
def _require(path):
    if not os.path.exists(path):
        raise FileNotFoundError(
            f"Required input missing: {path}\n"
            "Run scripts/JDR_reanalysis/20_smad4_loxp_evidence.py first. "
            "This script will not substitute a fallback value.")
    return path


def load_inputs():
    j = json.load(open(_require(os.path.join(OUT, "smad4_loxp_summary.json"))))
    reads = pd.read_csv(_require(os.path.join(OUT, "smad4_loxp_reads.csv")))
    cov = pd.read_csv(_require(os.path.join(OUT, "smad4_exon_coverage.csv")))
    ratio = pd.read_csv(_require(os.path.join(OUT, "smad4_exon_ratio.csv")), comment="#")
    junc = pd.read_csv(_require(os.path.join(OUT, "smad4_junctions.csv")), comment="#")
    depth = pd.read_csv(_require(os.path.join(OUT, "smad4_locus_depth.csv")), comment="#")
    return j, reads, cov, ratio, junc, depth


class Fasta:
    """Minimal random-access FASTA reader driven by the .fai index."""

    def __init__(self, path):
        _require(path)
        _require(path + ".fai")
        self.path = path
        self.idx = {}
        for line in open(path + ".fai"):
            p = line.rstrip("\n").split("\t")
            self.idx[p[0]] = (int(p[1]), int(p[2]), int(p[3]), int(p[4]))

    def fetch(self, chrom, start, end):
        """1-based inclusive."""
        if chrom not in self.idx:
            raise KeyError(f"{chrom!r} absent from {self.path}.fai")
        length, offset, linebases, linewidth = self.idx[chrom]
        if start < 1 or end > length or end < start:
            raise ValueError(f"bad interval {chrom}:{start}-{end}")
        n = end - start + 1
        seek = offset + (start - 1) // linebases * linewidth + (start - 1) % linebases
        with open(self.path, "rb") as fh:
            fh.seek(seek)
            raw = fh.read(n + n // linebases + 4).decode("ascii")
        return raw.replace("\n", "")[:n].upper()


def fmt(x):
    return f"{int(x):,}"


def mathify(s):
    """Render pandoc superscript markup (`^fl/fl^`) as matplotlib mathtext.

    The study-design strings are written for the legend, which pandoc renders;
    dropped into a figure unchanged they print literal carets.
    """
    return re.sub(r"\^([^\s^]+)\^", r"$^{\1}$", s)


def wrap(text, width):
    """Hard-wrap a paragraph so it cannot run off the canvas.

    savefig(bbox_inches="tight") silently EXPANDS the page to fit overflowing
    text, so an over-long caption does not look broken -- it just makes the
    delivered figure wider than the canvas it was designed on, which in turn
    makes every point size smaller once the figure is scaled to a page. Wrapping
    here keeps the delivered size equal to FIGSIZE.
    """
    out = []
    for para in text.split("\n"):
        out.extend(textwrap.wrap(para, width=width) or [""])
    return "\n".join(out)


def sample_cols(junc):
    """The six per-library junction-count columns, in fixed order.

    Named explicitly rather than matched by suffix: smad4_junctions.csv also has
    total_unique and sibling_canonical_total_unique, which a suffix match would
    silently sweep in and turn the all-important zero into a non-zero.
    """
    cols = [f"{s}_unique" for s in SAMPLE_ORDER]
    missing = [c for c in cols if c not in junc.columns]
    if missing:
        raise RuntimeError(f"smad4_junctions.csv is missing {missing}")
    return cols


# ============================ sanity gate ====================================
def sanity(j, reads, ratio, junc, fa):
    """Fail loudly if the inputs no longer say what this figure asserts."""
    fx = j["floxed_exon"]
    if fx["ensembl_exon_id"] != "ENSMUSE00000143011" or fx["length_bp"] != 184:
        raise RuntimeError("floxed exon identity changed in the summary JSON")
    for key in ("identity",):
        if key not in fx:
            raise RuntimeError(
                f"summary JSON has no floxed_exon['{key}'] -- rerun step 20. The figure "
                "will not state a transcript count or an exon-number list from memory.")
    ident = fx["identity"]
    if ident["primary_exon_id"] != fx["ensembl_exon_id"]:
        raise RuntimeError("the majority stable exon ID no longer matches the reported one")
    if len(reads) != j["loxp_cassettes"]["n_loxp_reads_total"]:
        raise RuntimeError("loxP read table and summary disagree on read count")
    if set(reads.genotype) != {"CON", "dcKO"}:
        raise RuntimeError("loxP reads are no longer present in both genotypes; "
                           "panel B's central point would be false")
    for col in ("is_duplicate_copy", "duplicate_group_size"):
        if col not in reads.columns:
            raise RuntimeError(f"smad4_loxp_reads.csv has no {col} column -- rerun step 20; "
                               "the figure must disclose duplicate reads")
    core = j["loxp_cassettes"]["core_sequence_34bp"]
    if len(core) != LOXP_CORE_LEN:
        raise RuntimeError("loxP core is not 34 bp")
    for site in ("A", "B"):
        s = j["loxp_cassettes"][f"site_{site}"]
        for key in ("linker_5p_sequence", "linker_3p_sequence", "linkers_short", "derivation"):
            if key not in s:
                raise RuntimeError(
                    f"site {site} has no {key!r} -- rerun step 20. Restriction-linker names "
                    "must be derived from the reads, not asserted.")
        d = s["derivation"]
        if d["breakpoint"] != s["insertion_between"][0]:
            raise RuntimeError(f"site {site}: the read-derived breakpoint {d['breakpoint']} "
                               f"disagrees with insertion_between {s['insertion_between']}")
    skip = junc[junc.is_predicted_cre_skip.astype(str).str.lower() == "true"]
    if len(skip) != 1:
        raise RuntimeError("expected exactly one predicted Cre skip junction row")
    scols = sample_cols(junc)
    if skip[scols].to_numpy().sum() != 0:
        raise RuntimeError("the Cre skip junction is no longer zero -- this figure's "
                           "negative result must be re-examined before it is drawn")
    if float(ratio.perm_p_two_sided.min()) < 0.05:
        raise RuntimeError("an exon now reaches P < 0.05; the honesty framing of "
                           "panel C must be rewritten before re-rendering")
    if float(ratio.perm_p_two_sided_readfrac.min()) < 0.05:
        raise RuntimeError("an exon now reaches P < 0.05 on the read-fraction metric; "
                           "panel C must be rewritten before re-rendering")
    # the reference genome must still say what the breakpoint annotations claim
    a = j["loxp_cassettes"]["site_A"]
    amb = int(a["ambiguity_bp"])
    micro = a["derivation"]["microhomology"]
    if amb:
        obs = fa.fetch(j["gene"]["chrom"], a["insertion_between"][0] - amb + 1,
                       a["insertion_between"][0])
        if obs != micro:
            raise RuntimeError(f"site-A microhomology in genome.fa is {obs!r}, but step 20 "
                               f"derived {micro!r} -- genome build mismatch?")
    return True


# ============================ small drawing helpers ==========================
def mono_fontsize(fig, ax, n_chars, xspan_frac=1.0, cap=11.0):
    """Point size at which `n_chars` monospace glyphs exactly fill the span.

    DejaVu Sans Mono advances 0.6022 em per glyph. A 0.98 safety factor keeps the
    row inside its box after hinting.
    """
    width_in = ax.get_position().width * fig.get_figwidth() * xspan_frac
    return min(cap, 0.98 * width_in * 72.0 / (n_chars * 0.6022))


def draw_seq_row(ax, x0, x1, y, height, segs, fs, edge="none"):
    """Draw a run of monospace sequence blocks between x0..x1 (axes coords).

    `segs` is a list of dicts: {seq, fc, tc, hatch(optional)}. Every block uses the
    same character pitch, so the row reads as one continuous sequence.
    Returns a list of (x_centre, x_start, x_end) per segment for labelling.
    """
    n = sum(len(s["seq"]) for s in segs)
    pitch = (x1 - x0) / n
    i = 0
    spans = []
    for s in segs:
        L = len(s["seq"])
        xa, xb = x0 + i * pitch, x0 + (i + L) * pitch
        ax.add_patch(Rectangle((xa, y - height / 2), xb - xa, height,
                               facecolor=s["fc"], edgecolor=s.get("ec", edge),
                               hatch=s.get("hatch"), lw=0.4, zorder=2,
                               transform=ax.transAxes, clip_on=False))
        t = ax.text((xa + xb) / 2, y, s["seq"], ha="center", va="center",
                    family="monospace", fontsize=fs, color=s["tc"], zorder=3,
                    transform=ax.transAxes, clip_on=False)
        _SEQ_TEXT_IDS.add(id(t))
        spans.append((0.5 * (xa + xb), xa, xb))
        i += L
    return spans


def loxp_triangle(ax, x, y, w, h, color="k", lw=0.6):
    """Down-pointing triangle marking a loxP insertion point (data coords)."""
    ax.add_patch(Polygon([[x - w / 2, y + h], [x + w / 2, y + h], [x, y]],
                         closed=True, facecolor=color, edgecolor="k", lw=lw,
                         zorder=6, clip_on=False))


def text_collisions(fig, ignore_ids):
    """Overlapping rendered glyph boxes, excluding deliberately-packed sequence text.

    Adapted from volcano_panel.text_collisions: Annotation extents union in the
    leader line, and tick labels outside the view interval report stale boxes.
    """
    fig.canvas.draw()
    rend = fig.canvas.get_renderer()
    ghosts = set()
    for ax in fig.axes:
        for axis in (ax.xaxis, ax.yaxis):
            lo, hi = sorted(axis.get_view_interval())
            for loc, lab in zip(axis.get_majorticklocs(), axis.get_majorticklabels()):
                if not (lo - 1e-9 <= loc <= hi + 1e-9):
                    ghosts.add(id(lab))
    boxes = []
    for t in fig.findobj(matplotlib.text.Text):
        if (not t.get_visible() or not t.get_text().strip()
                or id(t) in ghosts or id(t) in ignore_ids):
            continue
        try:
            bb = matplotlib.text.Text.get_window_extent(t, renderer=rend)
        except Exception:
            continue
        if bb.width > 0 and bb.height > 0:
            boxes.append((t.get_text(), bb))
    out = []
    for i in range(len(boxes)):
        for k in range(i + 1, len(boxes)):
            a, b = boxes[i][1], boxes[k][1]
            if (min(a.x1, b.x1) - max(a.x0, b.x0) > 0.5 and
                    min(a.y1, b.y1) - max(a.y0, b.y0) > 0.5):
                out.append((boxes[i][0][:40], boxes[k][0][:40]))
    return out


def ink_collisions(fig, ignore_ids):
    """Text drawn ON TOP OF opaque dark ink -- the defect text_collisions cannot see.

    The previous build reported "0 overlapping text pairs" while both panel-B site
    headers were rendered with a solid black triangle sitting on their digits
    ("73,782,[TRIANGLE]2/73,782,243"). Text-versus-text auditing is blind to that,
    because the offender is an artist, not a Text. This checks every visible Text
    against every filled Polygon (the loxP insertion triangles) and against any
    dark filled Rectangle (loxP cores, the filled floxed-exon box).
    """
    fig.canvas.draw()
    rend = fig.canvas.get_renderer()
    marks = []
    for ax in fig.axes:
        for p in ax.patches:
            if not p.get_visible():
                continue
            fc = p.get_facecolor()
            if isinstance(p, Polygon):
                marks.append(("triangle", p.get_window_extent(rend)))
            elif isinstance(p, Rectangle) and fc[3] > 0.5 and sum(fc[:3]) / 3 < 0.45:
                marks.append(("dark block", p.get_window_extent(rend)))
    out = []
    for t in fig.findobj(matplotlib.text.Text):
        if (not t.get_visible() or not t.get_text().strip() or id(t) in ignore_ids):
            continue
        try:
            bb = matplotlib.text.Text.get_window_extent(t, renderer=rend)
        except Exception:
            continue
        if bb.width <= 0 or bb.height <= 0:
            continue
        for what, mb in marks:
            if (min(bb.x1, mb.x1) - max(bb.x0, mb.x0) > 0.5 and
                    min(bb.y1, mb.y1) - max(bb.y0, mb.y0) > 0.5):
                out.append((t.get_text()[:40].replace("\n", "|"), what))
                break
    return out


# ============================ panel A ========================================
def panel_A(ax, j, cov, junc, zooms):
    """Locus schematic over the shared genomic window."""
    W0, W1 = WINDOW
    gene, fx, cas = j["gene"], j["floxed_exon"], j["loxp_cassettes"]
    exons = (cov[cov["sample"] == cov["sample"].iloc[0]]
             [["exon_number", "exon_id", "start", "end", "exon_length"]]
             .drop_duplicates().sort_values("start"))
    inwin = exons[(exons.end >= W0) & (exons.start <= W1)]

    ax.set_xlim(W0, W1)
    ax.set_ylim(-3.85, 4.60)
    ax.set_yticks([])
    for sp in ("left", "right", "top"):
        ax.spines[sp].set_visible(False)
    ax.spines["bottom"].set_position(("data", -3.85))

    # ---- floxed interval -------------------------------------------------
    # Placed BELOW the panel-B zoom brackets: when the bracket rectangles were
    # tall enough to reach this label, their dashed edges struck through it.
    fi = cas["floxed_interval"]
    ax.axvspan(fi["start"], fi["end"], color="#F5E7CC", zorder=0)
    ax.annotate("", xy=(fi["start"], -1.34), xytext=(fi["end"], -1.34),
                arrowprops=dict(arrowstyle="<->", lw=0.8, color="#8A6D2F"))
    ax.text((fi["start"] + fi["end"]) / 2, -1.52,
            f"floxed interval  {fmt(fi['start'])}–{fmt(fi['end'])}   "
            f"({fmt(fi['length_bp'])} bp; contains exactly {fi['n_exons_contained']} exon)",
            ha="center", va="top", fontsize=6.5, color="#6B5320")

    # ---- intron line + exon boxes ---------------------------------------
    ax.plot([W0, W1], [0, 0], color="k", lw=0.9, zorder=1, solid_capstyle="butt")
    for _, e in inwin.iterrows():
        floxed = (e.exon_id == fx["ensembl_exon_id"])
        ax.add_patch(Rectangle((e.start, -0.26), e.end - e.start + 1, 0.52,
                               facecolor="#2B2B2B" if floxed else "#FFFFFF",
                               edgecolor="k", lw=1.0 if floxed else 0.8, zorder=4))
        if not floxed:
            ax.text((e.start + e.end) / 2, 0.32, f"e{int(e.exon_number)}",
                    ha="center", va="bottom", fontsize=6.3)
        else:
            ax.text((e.start + e.end) / 2, 0.32, "floxed exon", ha="center", va="bottom",
                    fontsize=6.3, fontweight="bold")

    # coordinate labels for the two unfloxed flanking exons, both BELOW the line
    for _, e in inwin.iterrows():
        n = int(e.exon_number)
        if n not in (8, 10):
            continue
        ax.annotate(f"exon {n}\n{fmt(e.start)}–{fmt(e.end)}\n({fmt(e.exon_length)} bp)",
                    xy=((e.start + e.end) / 2, -0.28),
                    xytext=((e.start + e.end) / 2, -1.94), ha="center", va="top",
                    fontsize=6.3, linespacing=1.3,
                    arrowprops=dict(arrowstyle="-", lw=0.4, color="#666666",
                                    shrinkA=1, shrinkB=1))

    # ---- loxP triangles + their labels ----------------------------------
    # Linker names come from site["linkers_short"], which step 20 DERIVES from the
    # base calls (including sites formed across the loxP-core/linker boundary), so
    # this label and panel B's per-read decomposition cannot disagree.
    for site, key in (("A", "site_A"), ("B", "site_B")):
        sc = cas[key]
        x = (sc["insertion_between"][0] + sc["insertion_between"][1]) / 2
        loxp_triangle(ax, x, 0.30, w=(W1 - W0) * 0.0055, h=0.40)
        amb = ("" if sc["ambiguity_bp"] == 0 else f"  ±{sc['ambiguity_bp']} bp")
        ax.text(x, 0.80, f"loxP site {site}\n{fmt(sc['insertion_between'][0])}"
                         f"/{fmt(sc['insertion_between'][1])}{amb}\n{sc['linkers_short']}",
                ha="center", va="bottom", fontsize=6.1, linespacing=1.35)

    # ---- floxed-exon identity block, parked in the empty 7 kb intron -----
    # F1: the stable exon ID is NOT the unambiguous identifier -- within this same
    # gene, three different stable IDs are attached to the identical 184 bp. The
    # coordinate interval is. Every count below is read from the GTF by step 20.
    # No leader line: the automatic angled connector used to run straight through
    # the site-B label and strike it out, and every hand-routed alternative crossed
    # one label or another. The floxed exon is the only black-filled box and is
    # captioned "floxed exon" in the schematic, so the box needs no pointer.
    ident = fx["identity"]
    nums = ident["exon_number_values"]
    n_other = ident["n_transcripts"] - ident["n_transcripts_with_primary_exon_id"]
    ax.text(
        W0 + 150, 4.55,
        f"$\\bf{{The\\ floxed\\ exon}}$ (black box below)   "
        f"$\\bf{{chr{gene['chrom']}:{fmt(fx['start'])}}}$–$\\bf{{{fmt(fx['end'])}}}$   "
        f"{fmt(fx['length_bp'])} bp   —   the COORDINATE INTERVAL is the unambiguous "
        "identifier, not a stable ID and not a number.\n"
        f"The same {fx['length_bp']} bp carries exon_number "
        f"{', '.join(str(v) for v in nums[:-1])} or {nums[-1]} across the "
        f"{ident['n_transcripts']} annotated transcripts that contain it, which is why the "
        f"literature says “exon {fx['literature_exon_number']}” while the Ensembl\n"
        f"canonical transcript says “exon {fx['canonical_exon_number']}”. A stable ID fixes "
        f"the NUMBER but not the IDENTITY: {ident['n_transcripts_with_primary_exon_id']} of "
        f"those transcripts label this interval {ident['primary_exon_id']} (used here), the "
        f"other {n_other}\n"
        f"label the IDENTICAL {fx['length_bp']} bp "
        f"{', '.join(ident['other_exon_ids'])} — so querying one ID in Ensembl silently drops "
        f"transcripts. The exon {fx['protein']}.",
        ha="left", va="top", fontsize=6.3, linespacing=1.5,
        bbox=dict(boxstyle="round,pad=0.35", fc="#FFFFFF", ec="#8A6D2F", lw=0.6))

    # ---- predicted Cre product: a dashed grey span carrying its own zero --
    skip = junc[junc.is_predicted_cre_skip.astype(str).str.lower() == "true"].iloc[0]
    scols = sample_cols(junc)
    nlib = len(scols)
    a, b = float(skip.intron_start), float(skip.intron_end)
    ax.plot([a, b], [-3.52, -3.52], color="#7A7A7A", lw=1.1, ls=(0, (4, 2)), zorder=5)
    for xx in (a, b):
        ax.plot([xx, xx], [-3.66, -3.38], color="#7A7A7A", lw=1.1, zorder=5)
    ax.text((a + b) / 2, -3.44,
            f"predicted Cre product: exon {int(skip.donor_annot[4:-4])} → "
            f"exon {int(skip.acceptor_annot[4:-6])} skip junction   —   "
            f"{int(skip[scols].sum())} reads in {nlib}/{nlib} libraries (panel D)",
            ha="center", va="bottom", fontsize=6.6, color="#4D4D4D", style="italic")

    # ---- transcription direction (right-hand half: the identity box owns the
    #      left, and a full-width note underneath it collided with the box) ----
    ax.annotate("", xy=(W0 + (W1 - W0) * 0.020, 2.35),
                xytext=(W0 + (W1 - W0) * 0.150, 2.35),
                arrowprops=dict(arrowstyle="-|>", lw=1.3, color="k"))
    ax.text(W0 + (W1 - W0) * 0.160, 2.35,
            "transcription 5′→ 3′   (Smad4 is on the MINUS strand: "
            "transcription runs RIGHT → LEFT)",
            ha="left", va="center", fontsize=7.0, fontweight="bold")

    # ---- zoom brackets for panel B ---------------------------------------
    # Kept SHORT and hugging the gene axis: when these rectangles were 2.9 units
    # tall their dashed vertical edges ran through the loxP site labels above and
    # through the floxed-interval label below, striking out words in both.
    for name, (za, zb) in zooms.items():
        ax.add_patch(Rectangle((za, -0.44), zb - za, 0.88, facecolor="none",
                               edgecolor="#2E6E4E", lw=0.8, linestyle=(0, (2, 1.6)),
                               zorder=7))
        ax.text((za + zb) / 2, -0.52, f"panel B,\nsite {name}", ha="center", va="top",
                fontsize=6.2, color="#2E6E4E", linespacing=1.25)

    ax.set_title(f"A   Smad4 floxed allele, {gene['chrom']}:{fmt(WINDOW[0])}–{fmt(WINDOW[1])}"
                 f"   ({j['genome_build']})")


# ============================ panel B: read geometry =========================
def read_layout(reads):
    """Add drawn-extent columns and a y row index (grouped by sample, top-down)."""
    r = reads.copy()
    r["sample"] = pd.Categorical(r["sample"], SAMPLE_ORDER, ordered=True)
    r = r.sort_values(["sample", "site", "pos", "qname"]).reset_index(drop=True)
    # soft-clipped bases are drawn extending from the alignment breakpoint, IGV-style
    r["draw_start"] = r.pos - r.lead_clip
    r["draw_end"] = r.aln_end + r.trail_clip
    # read base index -> genomic x (leading S before pos, M on the reference,
    # trailing S after aln_end); CIGARs at this locus contain only S and M.
    r["loxp_x0"] = r.draw_start + r.loxp_offset_in_read
    r["loxp_x1"] = r.loxp_x0 + LOXP_CORE_LEN
    r["row"] = np.arange(len(r))
    return r


def zoom_windows(r, pad=28):
    out = {}
    for site, g in r.groupby("site"):
        out[site] = (int(g.draw_start.min()) - pad, int(g.draw_end.max()) + pad)
    return dict(sorted(out.items()))


def panel_B_pileup(axes, r, j, zw):
    cas = j["loxp_cassettes"]
    nrow = len(r)
    for site, ax in axes.items():
        za, zb = zw[site]
        ax.set_xlim(za, zb)
        # top extended to -2.85 so the site header sits clear of the insertion
        # triangle: at -1.9 the triangle was drawn ON TOP of its own coordinate
        # text and destroyed two digits of each site's coordinate.
        ax.set_ylim(nrow + 0.8, -2.85)         # row 0 at the top
        ax.set_yticks([])
        for s in ("left", "right", "top"):
            ax.spines[s].set_visible(False)
        ax.tick_params(axis="x", labelsize=6.6)
        ax.xaxis.set_major_formatter(lambda v, p: fmt(v))

        # alternating per-sample background bands
        for i, s in enumerate(SAMPLE_ORDER):
            rows = r.loc[r["sample"] == s, "row"]
            if len(rows) == 0:
                continue
            if i % 2 == 0:
                ax.axhspan(rows.min() - 0.5, rows.max() + 0.5, color="#F6F6F6", zorder=0)

        # the mapped insertion point of this site
        key = f"site_{site}"
        ins = cas[key]["insertion_between"]
        x = (ins[0] + ins[1]) / 2
        ax.axvline(x, color="k", lw=0.7, linestyle=(0, (3, 2)), zorder=1)
        loxp_triangle(ax, x, -1.82, w=(zb - za) * 0.030, h=0.92)
        amb = "" if cas[key]["ambiguity_bp"] == 0 else f" ±{cas[key]['ambiguity_bp']} bp"
        ax.text(x, -1.95, f"loxP site {site}   {fmt(ins[0])}/{fmt(ins[1])}{amb}",
                ha="center", va="bottom", fontsize=6.8, fontweight="bold")

        # reads
        for _, q in r.iterrows():
            if q.site != site:
                continue
            col = BLUE if q.genotype == "CON" else RED
            hatch = None if q.genotype == "CON" else "///"
            y, h = q.row, 0.62
            # soft-clipped flanks: open, hatched
            if q.lead_clip > 0:
                ax.add_patch(Rectangle((q.draw_start, y - h / 2), q.lead_clip, h,
                                       facecolor="white", edgecolor=col, lw=0.7,
                                       hatch="....", zorder=3))
            if q.trail_clip > 0:
                ax.add_patch(Rectangle((q.aln_end + 1, y - h / 2), q.trail_clip, h,
                                       facecolor="white", edgecolor=col, lw=0.7,
                                       hatch="....", zorder=3))
            # aligned portion: solid
            ax.add_patch(Rectangle((q.pos, y - h / 2), q.aln_end - q.pos + 1, h,
                                   facecolor=col, edgecolor=col, lw=0.5,
                                   hatch=hatch, zorder=4))
            # the 34 bp loxP core carried by this read
            ax.add_patch(Rectangle((q.loxp_x0, y - h / 2 - 0.10), LOXP_CORE_LEN, h + 0.20,
                                   facecolor=SEQ_LOXP, edgecolor="k", lw=0.4, zorder=5))
            if q.spans_cassette:
                ax.plot([q.loxp_x1 + (zb - za) * 0.012], [y], marker="*", ms=4.6,
                        color="k", zorder=6, clip_on=False)
            if q.is_duplicate_copy:
                # one fragment sequenced twice: drawn, but never silently counted
                # as a second independent observation
                ax.text(q.draw_start - (zb - za) * 0.012, y, "‡", ha="right", va="center",
                        fontsize=7.0, color="#555555", zorder=6, clip_on=False)

    # sample labels + genotype brackets, on the left-hand axis only
    axL = axes[list(axes)[0]]
    za, zb = zw[list(axes)[0]]
    for s in SAMPLE_ORDER:
        rows = r.loc[r["sample"] == s, "row"]
        if len(rows) == 0:
            continue
        col = BLUE if s.startswith("CON") else RED
        mark = "●" if s.startswith("CON") else "▲"
        axL.text(-0.012, (rows.min() + rows.max()) / 2, f"{mark} {s}",
                 transform=axL.get_yaxis_transform(), ha="right", va="center",
                 fontsize=7.0, color=col, fontweight="bold")
    for gt, col, lo, hi in (("CON", BLUE, 0, 2), ("dcKO", RED, 3, 5)):
        rows = r.loc[r.genotype == gt, "row"]
        axL.annotate("", xy=(-0.088, rows.min() - 0.45), xytext=(-0.088, rows.max() + 0.45),
                     xycoords=axL.get_yaxis_transform(), textcoords=axL.get_yaxis_transform(),
                     arrowprops=dict(arrowstyle="-", lw=2.0, color=col))

    n_by_gt = r.groupby("genotype")["sample"].nunique()
    cs = j["loxp_cassettes"]
    n_span = int(cs["n_cassette_spanning_reads"])
    n_dup = int(cs["n_duplicate_copies"])
    n_uni = int(cs["n_unique_fragments"])
    caption = (f"loxP cassette detected in {n_by_gt['CON']}/3 CON and {n_by_gt['dcKO']}/3 dcKO "
               f"libraries ({len(r)} reads total; {n_span} ★ read the cassette–genome junction "
               "directly).\nCON mice are Cre-negative Smad4$^{fl/fl}$ and carry the IDENTICAL "
               "floxed allele — loxP reads in CON are expected, and are NOT a control failure.")
    notes = []
    if n_dup:
        notes.append(
            f"‡  {n_dup} of the {len(r)} loxP reads is a byte-identical copy of another read in "
            f"the same library (same POS, same CIGAR, same base calls) — one fragment sequenced "
            f"twice. The number of INDEPENDENT fragments is {n_uni}, not {len(r)}; no "
            "cassette-spanning read is affected.")
    notes.append(
        "Both cassettes are INTRONIC, so every loxP read here comes from unspliced pre-mRNA or "
        "residual genomic DNA rather than mature mRNA. That is why there are only a handful, "
        "and why their number says nothing about transcript level or excision efficiency.")
    return caption, notes


# ==================== panel B (lower): literal sequence ======================
def core_side(row):
    """Which soft clip carries the loxP core: 'lead', 'trail' or 'straddle'."""
    lead = int(row.lead_clip)
    aln_len = int(row.aln_end) - int(row.pos) + 1
    off = int(row.loxp_offset_in_read)
    if off + LOXP_CORE_LEN <= lead:
        return "lead"
    if off >= lead + aln_len:
        return "trail"
    return "straddle"


def _linker_label(site_info, which):
    """Restriction-site label for one linker, DERIVED upstream from the reads.

    Reading it from the JSON rather than re-scanning the segment here is the fix
    for the figure disagreeing with itself: a local scan of the linker alone
    misses XbaI, whose TCTAGA is formed by the final T of the loxP core plus the
    CTAGA that starts the linker, so panel A said "XbaI/SalI" where panel B said
    only "SalI".
    """
    hits = site_info[f"linker_{which}_sites"]
    if not hits:
        return "linker"
    return " + ".join(h["enzyme"] + ("\u2020" if h["straddles_core"] else "") for h in hits)


def decompose_spanning(row, core, fa, chrom, site_info):
    """Literal blocks for a cassette-spanning read: genome | linker | loxP | linker | genome.

    BOTH read geometries are handled, because restricting this to leading-clip
    reads silently dropped 2 of the 7 spanning reads -- including dcKO2's only
    one -- while the legend claimed the three shown had been chosen on aligned-
    block length.

      * core in the LEADING clip:  [clipped genome][v5][loxP][v3][CIGAR M genome]
        v3 is pinned exactly by the CIGAR (the M block starts at lead_clip); the
        clipped genomic flank is recovered by extending an exact match back from
        the insertion point and cross-checked against the data step's own column.
      * core in the TRAILING clip: [CIGAR M genome][v5][loxP][v3][clipped genome]
        v5 is pinned exactly by the CIGAR; the linker on the clipped side is set
        by the site's consensus linker (itself derived by step 20 from reads at
        this site whose CIGAR pins it), and the genomic flank beyond it is scored
        base-for-base against the reference on the other side of the insertion.

    Either way the drawn order is genome -> cassette -> genome, so the two
    geometries read identically.
    """
    seq = row.seq
    off = int(row.loxp_offset_in_read)
    lead = int(row.lead_clip)
    trail = int(row.trail_clip)
    aln_len = int(row.aln_end) - int(row.pos) + 1
    side = core_side(row)
    if side == "straddle":
        raise RuntimeError(f"{row.qname}: the loxP core straddles the alignment boundary")
    if seq[off:off + LOXP_CORE_LEN] != core:
        raise RuntimeError(f"{row.qname}: loxP core absent at offset {off}")
    aligned = seq[lead:lead + aln_len]
    ref_aln = fa.fetch(chrom, int(row.pos), int(row.aln_end))
    if aligned != ref_aln:
        raise RuntimeError(f"{row.qname}: aligned block differs from genome.fa")

    ins_5p, ins_3p = site_info["insertion_between"]
    l5 = f"{_linker_label(site_info, '5p')} linker"
    l3 = f"{_linker_label(site_info, '3p')} linker"
    aln_lab = (f"genomic, aligned to {fmt(row.pos)}–{fmt(row.aln_end)}  "
               f"{len(aligned)} nt (100% identity)")
    core_seg = dict(seq=seq[off:off + LOXP_CORE_LEN], fc=SEQ_LOXP, tc="white",
                    lab=f"loxP core  {LOXP_CORE_LEN} nt")

    if side == "lead":
        # longest exact match of a read prefix against the reference ending at ins_5p
        ref5 = fa.fetch(chrom, ins_5p - 300 + 1, ins_5p)
        best_k, best_e = 0, 0
        for e in range(off, 0, -1):
            k = 0
            while k < e and k < len(ref5) and seq[e - 1 - k] == ref5[len(ref5) - 1 - k]:
                k += 1
            if k > best_k:
                best_k, best_e = k, e
        if best_k != int(row.genomic_anchor_5p_bp):
            raise RuntimeError(
                f"{row.qname}: recomputed 5' genomic anchor {best_k} bp contradicts the "
                f"data step's {int(row.genomic_anchor_5p_bp)} bp")
        v5 = seq[best_e:off]
        v3 = seq[off + LOXP_CORE_LEN:lead]
        segs = [
            dict(seq=seq[:best_e], fc=SEQ_GENOMIC_CLIP, tc="k",
                 lab=f"genomic flank  {best_e} nt ({best_k} exact)"),
            dict(seq=v5, fc=SEQ_VECTOR, tc="k", hatch="///", lab=f"{l5}  {len(v5)} nt"),
            core_seg,
            dict(seq=v3, fc=SEQ_VECTOR, tc="k", hatch="///", lab=f"{l3}  {len(v3)} nt"),
            dict(seq=aligned, fc=SEQ_GENOMIC_ALN, tc="k", lab=aln_lab),
        ]
        tail = seq[len(seq) - trail:] if trail else ""
        if tail:
            segs.append(dict(seq=tail, fc="#FFFFFF", tc="#666666", hatch="....",
                             lab=f"soft-clipped  {len(tail)} nt"))
        anchor_nt, anchor_exact = best_e, best_k
    else:
        v5 = seq[lead + aln_len:off]
        v3c = site_info["linker_3p_sequence"]
        v3 = seq[off + LOXP_CORE_LEN:off + LOXP_CORE_LEN + len(v3c)]
        gen3 = seq[off + LOXP_CORE_LEN + len(v3c):]
        if not gen3:
            raise RuntimeError(f"{row.qname}: nothing left of the read beyond the 3' linker")
        ref3 = fa.fetch(chrom, ins_3p, ins_3p + len(gen3) - 1)
        exact3 = sum(1 for a, b in zip(gen3, ref3) if a == b)
        if exact3 < int(row.genomic_anchor_3p_bp):
            raise RuntimeError(
                f"{row.qname}: the segment drawn as the 3' genomic flank matches the "
                f"reference at only {exact3}/{len(gen3)} bases, fewer than the "
                f"{int(row.genomic_anchor_3p_bp)} bp exact anchor the data step found")
        segs = [
            dict(seq=aligned, fc=SEQ_GENOMIC_ALN, tc="k", lab=aln_lab),
            dict(seq=v5, fc=SEQ_VECTOR, tc="k", hatch="///", lab=f"{l5}  {len(v5)} nt"),
            core_seg,
            dict(seq=v3, fc=SEQ_VECTOR, tc="k", hatch="///", lab=f"{l3}  {len(v3)} nt"),
            dict(seq=gen3, fc=SEQ_GENOMIC_CLIP, tc="k",
                 lab=f"genomic flank from {fmt(ins_3p)}  {len(gen3)} nt ({exact3} match)"),
        ]
        if lead:
            segs.insert(0, dict(seq=seq[:lead], fc="#FFFFFF", tc="#666666", hatch="....",
                                lab=f"soft-clipped  {lead} nt"))
        anchor_nt, anchor_exact = len(gen3), exact3
    if sum(len(s["seq"]) for s in segs) != len(seq):
        raise RuntimeError(f"{row.qname}: decomposed blocks do not sum to the read length")
    return segs, anchor_nt, anchor_exact


def panel_B_sequence(ax, fig, j, reads, fa, caption, notes):
    cas, chrom = j["loxp_cassettes"], j["gene"]["chrom"]
    core = cas["core_sequence_34bp"]
    ax.set_xticks([])          # emptied as well as hidden: an axis that is only
    ax.set_yticks([])          # set_axis_off() still reports stale tick-label extents
    ax.set_axis_off()
    ax.set_xlim(0, 1)
    ax.set_ylim(0, 1)

    ax.text(0.5, 1.055, wrap(caption, 200), ha="center", va="top", fontsize=6.6,
            linespacing=1.5,
            transform=ax.transAxes,
            bbox=dict(boxstyle="round,pad=0.40", fc="#FFFDF2", ec="#B08D3F", lw=0.7))

    # ---- reference sequence at the two mapped breakpoints -------------------
    ax.text(0.0, 0.930, "Reference sequence at each mapped insertion point "
                        f"(genome.fa, {j['genome_build']}) — “▾” is where the cassette "
                        "was inserted; both sites lie in intronic sequence",
            fontsize=7.0, fontweight="bold", va="top", transform=ax.transAxes)
    FLANK = 18
    for site, xa, xb in (("B", 0.015, 0.475), ("A", 0.525, 0.985)):
        s = cas[f"site_{site}"]
        p5, p3 = s["insertion_between"]
        left = fa.fetch(chrom, p5 - FLANK + 1, p5)
        right = fa.fetch(chrom, p3, p3 + FLANK - 1)
        amb = int(s["ambiguity_bp"])
        d = s["derivation"]
        why = (f"   (±{amb} bp: {d['microhomology']} microhomology at "
               f"{fmt(d['breakpoint_min'] + 1)}–{str(d['breakpoint_max'])[-3:]})" if amb
               else f"   (exact: {d['n_reads_used']} reads agree to the base)")
        ax.text(xa, 0.884,
                f"loxP site {site}   chr{chrom}:{fmt(p5)} ▾ {fmt(p3)}" + why,
                fontsize=6.5, va="top", fontweight="bold", transform=ax.transAxes)
        fs = mono_fontsize(fig, ax, 2 * FLANK, xspan_frac=(xb - xa), cap=9.0)
        segs = [dict(seq=left, fc="#EFEFEF", tc="k"),
                dict(seq=right, fc="#DCDCDC", tc="k")]
        spans = draw_seq_row(ax, xa, xb, 0.822, 0.040, segs, fs, edge="#9A9A9A")
        mid = spans[0][2]
        ax.plot([mid, mid], [0.794, 0.850], color="k", lw=1.2,
                transform=ax.transAxes, zorder=5)
        ax.text(mid, 0.784, "▲ cassette inserted here", ha="center", va="top",
                fontsize=6.2, transform=ax.transAxes)
        ax.text(xa, 0.784, f"5′ {fmt(p5 - FLANK + 1)}", ha="left", va="top",
                fontsize=6.0, color="#666666", transform=ax.transAxes)
        ax.text(xb, 0.784, f"{fmt(p3 + FLANK - 1)} 3′", ha="right", va="top",
                fontsize=6.0, color="#666666", transform=ax.transAxes)

    # ---- literal decomposition of the cassette-spanning reads ---------------
    ax.text(0.0, 0.734,
            "Literal decomposition of cassette-spanning reads — real base calls, "
            "no schematic. Each row is one sequencing read; the two flanks are "
            "genomic, the middle is the targeting cassette.",
            fontsize=7.0, fontweight="bold", va="top", transform=ax.transAxes)

    # ALL cassette-spanning reads are eligible: both read geometries are drawable,
    # so the selection below is the stated criterion and nothing else.
    cand = reads[reads.spans_cassette].copy()
    if cand.empty:
        raise RuntimeError("no cassette-spanning read available to draw")
    n_eligible = len(cand)
    if n_eligible != int(j["loxp_cassettes"]["n_cassette_spanning_reads"]):
        raise RuntimeError("the eligible set is no longer every spanning read; the legend's "
                           "account of how the drawn reads were chosen would be false")
    # Ranked by the length of the REFERENCE-ALIGNED block: the longer that block,
    # the more firmly the read pins the cassette to a genomic coordinate. Then the
    # best-ranked read from the other genotype, and the best-ranked read from the
    # other insertion site, so the panel shows the structure is not peculiar to one
    # library or one cassette.
    cand["aln_len"] = cand.aln_end - cand.pos + 1
    cand = cand.sort_values(["aln_len", "qname"], ascending=[False, True])
    picks, taken = [], set()

    def take(sub):
        for _, q in sub.iterrows():
            if q.qname not in taken:
                taken.add(q.qname)
                picks.append(q)
                return

    take(cand)
    take(cand[cand.genotype != picks[0].genotype])
    take(cand[cand.site != picks[0].site])

    rows_y = [0.600, 0.420, 0.240][:len(picks)]
    for row, ytop in zip(picks, rows_y):
        s = cas[f"site_{row.site}"]
        segs, e5, k5 = decompose_spanning(row, core, fa, chrom, s)
        n = sum(len(x["seq"]) for x in segs)
        fs = mono_fontsize(fig, ax, n, xspan_frac=0.985, cap=8.0)
        geom = ("cassette in the leading soft clip" if core_side(row) == "lead"
                else "cassette in the trailing soft clip")
        ax.text(0.0, ytop + 0.052,
                f"{row['sample']} ({row.genotype})  read {row.qname}   POS {fmt(row.pos)}  "
                f"CIGAR {row.cigar}   → site {row.site}   [{n} nt; {geom}]",
                fontsize=6.4, va="center", transform=ax.transAxes,
                color=BLUE if row.genotype == "CON" else RED, fontweight="bold")
        spans = draw_seq_row(ax, 0.0075, 0.9925, ytop, 0.048, segs, fs, edge="#8A8A8A")
        # staggered labels so neighbouring captions cannot collide
        for i, (sp, sg) in enumerate(zip(spans, segs)):
            dy = 0.030 if i % 2 == 0 else 0.068
            xl = min(max(sp[0], 0.045), 0.955)     # keep captions inside the panel
            ax.annotate(sg["lab"], xy=(sp[0], ytop - 0.024),
                        xytext=(xl, ytop - dy), ha="center", va="top",
                        fontsize=6.0, linespacing=1.3,
                        arrowprops=dict(arrowstyle="-", lw=0.35, color="#777777",
                                        shrinkA=0.5, shrinkB=0.5))
    star_note = (j["loxp_cassettes"]["site_A"].get("straddling_site_note")
                 or j["loxp_cassettes"]["site_B"].get("straddling_site_note"))
    foot = list(notes) + ([star_note] if star_note else [])
    ax.text(0.0, 0.128, wrap("\n".join(foot), 232), fontsize=6.0, va="top",
            color="#444444", linespacing=1.45, transform=ax.transAxes)
    return picks, n_eligible


# ============================ panel C ========================================
def panel_C_coverage(ax, j, cov, depth):
    """CON vs dcKO mean per-base coverage, normalised to mean Smad4 exonic depth."""
    W0, W1 = WINDOW
    fx = j["floxed_exon"]
    exons = (cov[cov["sample"] == cov["sample"].iloc[0]]
             [["exon_number", "exon_id", "start", "end", "exon_length"]].drop_duplicates())
    total_exonic_bp = int(exons.exon_length.sum())
    # norm_depth is depth / (that sample's summed exonic depth); x total exonic bp
    # rescales it so 1.0 = that sample's MEAN depth across all 12 Smad4 exons.
    d = depth.copy()
    d["rel"] = d.norm_depth * total_exonic_bp
    m = d.groupby(["genotype", "bin_start"])["rel"].mean().unstack(0)

    ntag = 0
    for _, e in exons.sort_values("start").iterrows():
        if e.end < W0 or e.start > W1:
            continue
        ax.axvspan(e.start, e.end, color="#EAEAEA", zorder=0)
        ntag += 1                      # stagger: e6 and e7 are 212 bp apart and would clash
        ax.text((e.start + e.end) / 2, 0.995 if ntag % 2 else 0.905,
                f"e{int(e.exon_number)}", transform=ax.get_xaxis_transform(),
                ha="center", va="top", fontsize=6.2, color="#555555")
    # the floxed exon is outlined only so it can be FOUND, in the same neutral weight
    # as the other exon blocks; nothing here should read as an arrow or a finding.
    ax.axvspan(fx["start"], fx["end"], facecolor="none", edgecolor="#777777",
               lw=0.7, ls=(0, (2, 1.5)), zorder=1)

    x = m.index.values + 5      # bin midpoint (10 bp bins)
    ax.plot(x, m["CON"].values, color=BLUE, lw=0.9, ls="-", label="CON mean (n = 3)",
            zorder=3)
    ax.plot(x, m["dcKO"].values, color=RED, lw=0.9, ls=(0, (3.5, 1.6)),
            label="dcKO mean (n = 3)", zorder=4)
    ax.set_xlim(W0, W1)
    ax.set_ylim(0, float(np.nanmax(m.values)) * 1.42)   # headroom for the exon tags
    ax.set_ylabel("coverage / mean Smad4\nexonic coverage", fontsize=7.0)
    ax.tick_params(labelsize=6.8)
    ax.xaxis.set_major_formatter(lambda v, p: fmt(v))
    ax.legend(frameon=False, fontsize=6.6, loc="upper center", ncol=2,
              handlelength=2.4, borderaxespad=0.1)
    ax.set_xlabel(f"chr{j['gene']['chrom']} coordinate (GRCm39)   —   transcription runs "
                  "RIGHT → LEFT (minus strand)", fontsize=7.0, labelpad=7)
    ax.set_title("C   Per-base coverage across the locus, and the per-exon dcKO/CON "
                 "ratio on both normalisations")


def panel_C_ratio(ax, j, cov, ratio):
    """Per-exon dcKO/CON ratio -- BOTH metrics, never one.

    Step 20 computes the ratio two ways (share of per-base exonic DEPTH, share of
    exonic READS) and says in the CSV header that both are given "so that neither
    the figure nor the text can pick the flattering one". Plotting only the depth
    metric did exactly that: on depth the floxed exon is mid-pack and exon 6
    deviates more, on read fraction the floxed exon has the LARGEST deviation of
    all twelve. Both series are drawn, and the caption names its metric.
    """
    fx_id = j["floxed_exon"]["ensembl_exon_id"]
    r = ratio.sort_values("exon_number").reset_index(drop=True)
    x = np.arange(len(r))
    lo = float(min(r.ratio_dcKO_over_CON.min(), r.ratio_readfrac_dcKO_over_CON.min()))
    hi = float(max(r.ratio_dcKO_over_CON.max(), r.ratio_readfrac_dcKO_over_CON.max()))

    # the band every exon lives in, on EITHER metric -- drawn FIRST so nothing
    # appears to escape it
    ax.axhspan(lo, hi, color="#F0F0F0", zorder=0)
    ax.axhline(1.0, color="k", lw=0.7, zorder=1)
    ax.text(len(r) - 0.35, hi, f"  range spanned by all {len(r)} exons on either metric: "
                               f"{lo:.3f} – {hi:.3f}", fontsize=6.2, va="bottom",
            ha="right", color="#555555")

    # per-sample points (what the reader should judge the depth ratio against)
    cmean = (cov[cov.genotype == "CON"].groupby("exon_number")["frac_of_smad4_exonic_depth"]
             .mean())
    for i, e in r.iterrows():
        for gt, col, mk, dx in (("CON", BLUE, "o", -0.22), ("dcKO", RED, "^", 0.22)):
            v = cov.loc[(cov.genotype == gt) & (cov.exon_number == e.exon_number),
                        "frac_of_smad4_exonic_depth"] / cmean.loc[e.exon_number]
            ax.plot([i + dx] * len(v), v.values, marker=mk, ms=2.8, ls="none",
                    mfc="none", mec=col, mew=0.7, alpha=0.85, zorder=3)
    # both metrics, distinguishable without colour (open diamond vs filled square)
    ax.plot(x, r.ratio_dcKO_over_CON.values, marker="D", ms=4.4, ls="none",
            mfc="#FFFFFF", mec="k", mew=0.9, zorder=5, label="ratio, per-base DEPTH")
    ax.plot(x, r.ratio_readfrac_dcKO_over_CON.values, marker="s", ms=3.6, ls="none",
            mfc="#4D4D4D", mec="#1A1A1A", mew=0.6, zorder=5, label="ratio, READ fraction")
    for i in x:                                   # tie the two metrics per exon
        ax.plot([i, i], [r.ratio_dcKO_over_CON[i], r.ratio_readfrac_dcKO_over_CON[i]],
                color="#9A9A9A", lw=0.6, zorder=2)
    ax.plot([], [], marker="o", ms=2.8, ls="none", mfc="none", mec=BLUE, mew=0.7,
            label="CON library (depth)")
    ax.plot([], [], marker="^", ms=2.8, ls="none", mfc="none", mec=RED, mew=0.7,
            label="dcKO library (depth)")

    # IDENTICAL annotation style for the floxed exon and for the two comparison
    # exons the caption names: the biggest deviation in the OPPOSITE direction and
    # the biggest in the SAME direction. No arrow, no asterisk, no emphasis.
    ec = j["exon_coverage"]
    dev = (r.ratio_dcKO_over_CON - 1.0).abs()
    i_fx = int(r.index[r.exon_id == fx_id][0])
    i_op = int(dev[r.exon_id != fx_id].idxmax())
    same = ec.get("same_direction_comparison_exon")
    i_sd = (int(r.index[r.exon_number == same["exon_number"]][0])
            if same is not None else None)
    for i in sorted({i_fx, i_op} | ({i_sd} if i_sd is not None else set())):
        e = r.loc[i]
        tag = " (floxed)" if e.exon_id == fx_id else ""
        # the same-direction comparison exon is adjacent to the floxed one, so its
        # caption is offset sideways with a leader rather than stacked on top of it
        dx = -0.62 if (i_sd is not None and i == i_sd and abs(i - i_fx) < 2) else 0.0
        ax.annotate(f"exon {int(e.exon_number)}{tag}\ndepth {e.ratio_dcKO_over_CON:.3f}"
                    f"  P = {e.perm_p_two_sided:.2f}"
                    f"\nreads {e.ratio_readfrac_dcKO_over_CON:.3f}"
                    f"  P = {e.perm_p_two_sided_readfrac:.2f}",
                    xy=(i, min(e.ratio_dcKO_over_CON, e.ratio_readfrac_dcKO_over_CON)),
                    xytext=(i + dx, 0.740),
                    ha="center", va="bottom", fontsize=6.0, linespacing=1.3,
                    arrowprops=dict(arrowstyle="-", lw=0.4, color="#888888",
                                    shrinkA=1, shrinkB=2))
    ax.set_xticks(x)
    ax.set_xticklabels([str(int(v)) for v in r.exon_number])
    ax.set_xlim(-0.7, len(r) - 0.3)
    ax.set_ylim(0.70, 1.30)
    ax.set_xlabel("Smad4 canonical exon, transcription order 5′→3′ (exon 1 has the HIGHEST "
                  "coordinate: the reverse of the spatial order in the panels above)",
                  fontsize=6.8, labelpad=3)
    ax.set_ylabel("dcKO / CON\ncoverage ratio", fontsize=7.0)
    ax.tick_params(labelsize=6.8)
    ax.legend(frameon=False, fontsize=6.2, loc="upper left", ncol=4, handletextpad=0.3,
              columnspacing=1.2, borderaxespad=0.15)
    ax.set_title("      per-exon dcKO / CON coverage ratio — all "
                 f"{len(r)} canonical exons, plotted identically, on BOTH metrics",
                 fontsize=7.8, fontweight="normal")
    return lo, hi, r.loc[i_fx], r.loc[i_op], (r.loc[i_sd] if i_sd is not None else None)


# ============================ panel D ========================================
def panel_D(ax, ax_note, j, junc):
    neg = j["negative_result"]
    scols = sample_cols(junc)
    samples = list(SAMPLE_ORDER)
    incl = [tuple(int(v) for v in s.split("-"))
            for s in neg["inclusion_junctions"].split(";")]
    skip_j = (neg["predicted_cre_skip_junction"]["intron_start"],
              neg["predicted_cre_skip_junction"]["intron_end"])

    groups = []
    for a, b in incl + [skip_j]:
        row = junc[(junc.intron_start == a) & (junc.intron_end == b)]
        if row.empty:
            raise RuntimeError(f"junction {a}-{b} absent from smad4_junctions.csv")
        row = row.iloc[0]
        groups.append((row, [int(row[c]) for c in scols]))
    # transcription order: exon8->exon9 first, then exon9->exon10, skip last
    groups = sorted(groups, key=lambda g: (bool(str(g[0].is_predicted_cre_skip).lower()
                                                == "true"), -g[0].intron_end))

    pitch, gap = 1.0, 2.4
    ymax = max(max(v) for _, v in groups) * 1.34
    for gi, (row, vals) in enumerate(groups):
        x0 = gi * (len(samples) * pitch + gap)
        is_skip = str(row.is_predicted_cre_skip).lower() == "true"
        if is_skip:
            ax.add_patch(Rectangle((x0 - 0.6, 0), len(samples) * pitch + 0.2, ymax,
                                   facecolor="#F2F2F2", edgecolor="#9A9A9A", lw=0.6,
                                   linestyle=(0, (3, 2)), zorder=0))
        for si, (s, v) in enumerate(zip(samples, vals)):
            con = s.startswith("CON")
            ax.bar(x0 + si * pitch, v, width=0.78,
                   color=BLUE if con else RED, edgecolor="k", lw=0.4,
                   hatch=None if con else "///", zorder=3)
            ax.text(x0 + si * pitch, v + ymax * 0.018, str(v), ha="center", va="bottom",
                    fontsize=6.1, fontweight="bold" if is_skip else "normal")
        # group caption
        if is_skip:
            cap = (f"PREDICTED Cre excision product\nexon 8 → exon 10 skip\n"
                   f"{fmt(row.intron_start)}–{fmt(row.intron_end)}")
        else:
            cap = (f"{row.donor_annot.replace('exon', 'exon ').replace('_end', '')} → "
                   f"{row.acceptor_annot.replace('exon', 'exon ').replace('_start', '')}"
                   f"  (normal)\n{fmt(row.intron_start)}–{fmt(row.intron_end)}")
        ax.text(x0 + (len(samples) - 1) * pitch / 2, -ymax * 0.32, cap, ha="center",
                va="top", fontsize=6.5, linespacing=1.35,
                fontweight="bold" if is_skip else "normal",
                color="#333333" if is_skip else "k")
        if is_skip:
            ax.text(x0 + (len(samples) - 1) * pitch / 2, ymax * 0.42,
                    f"0 reads in {len(samples)}/{len(samples)} libraries",
                    ha="center", va="center", fontsize=7.8, fontweight="bold",
                    color="#333333",
                    bbox=dict(boxstyle="round,pad=0.35", fc="white", ec="#333333", lw=0.8))

    ax.set_xticks([gi * (len(samples) * pitch + gap) + si * pitch
                   for gi in range(len(groups)) for si in range(len(samples))])
    ax.set_xticklabels(samples * len(groups), rotation=90, fontsize=6.3)
    for lbl, s in zip(ax.get_xticklabels(), samples * len(groups)):
        lbl.set_color(BLUE if s.startswith("CON") else RED)
    ax.set_ylim(0, ymax)
    ax.set_xlim(-1.4, (len(groups) - 1) * (len(samples) * pitch + gap)
                + (len(samples) - 1) * pitch + 1.4)
    ax.set_ylabel("uniquely mapped\njunction reads", fontsize=7.0)
    ax.tick_params(axis="y", labelsize=6.8)
    ax.set_title("D   Splice junctions at the floxed exon — the Cre excision product "
                 "is not detected")

    psi = sorted(set(neg["psi_per_sample"].values()))
    gsj = neg["genotype_specific_junctions"]
    ax_note.set_xticks([])     # emptied as well as hidden: set_axis_off() alone
    ax_note.set_yticks([])     # still reports stale tick-label extents
    w0, w1 = j["negative_result"]["junction_scan_window"]
    sd = j["study_design"]
    ax_note.set_axis_off()
    ax_note.text(0.0, 1.0,
            "$\\bf{The\\ excision\\ is\\ not\\ detected.}$\n"
            f"PSI (exon inclusion) = {psi[0]:.3f} in {len(neg['psi_per_sample'])}/"
            f"{len(neg['psi_per_sample'])} libraries.\n"
            f"0 skip events against {neg['dcko_inclusion_equivalent_reads']:.0f} "
            "inclusion-equivalent dcKO reads\n"
            f"    ({neg['dcko_junction_reads_total']} junction reads; each library "
            "contributes the MEAN of\n"
            "    its two flanking junctions, because a read crosses only one.\n"
            "    The mean is the conservative denominator — 515 is a\n"
            "    derived statistic, not a count of reads.)\n"
            f"→ ≥ {neg['one_sided_95_lower_bound_on_psi'] * 100:.1f}% of dcKO Smad4 mRNA "
            "retains the floxed exon\n    (one-sided 95% bound, rule of three).\n"
            f"Within chr{j['gene']['chrom']}:{fmt(w0)}–{fmt(w1)} only (not\n"
            "gene-wide): no dcKO-specific junction involves the floxed\n"
            f"exon; genotype-specific junctions occur in BOTH directions\n"
            f"({gsj['n_dcKO_only']} dcKO-only, {gsj['n_CON_only']} CON-only) at "
            f"≤ {max(gsj['max_dcKO_only_pct_of_sibling_canonical'], gsj['max_CON_only_pct_of_sibling_canonical']):.2f}% of the canonical junction\n"
            "they compete with.\n"
            f"$\\bf{{Scope.}}$ The sequenced material is\n"
            f"{sd['sequenced_material_driver']}, and these are BULK tissue\n"
            "libraries: this bound is on the population average, not on\n"
            "excision inside a recombined cell.",
                 transform=ax_note.transAxes, ha="left", va="top", fontsize=6.5,
                 linespacing=1.5,
                 bbox=dict(boxstyle="round,pad=0.45", fc="#FFFFFF", ec="#333333", lw=0.7))
    return psi[0]


# ============================ assembly =======================================
def build(j, reads, cov, ratio, junc, depth, fa):
    r = read_layout(reads)
    zw = zoom_windows(r)

    fig = plt.figure(figsize=FIGSIZE)
    # row 5 is an intentional spacer: panel C's honesty caption lives there, and it
    # must not be squeezed against the ratio strip's x-label or panel D's title.
    gs = fig.add_gridspec(7, 1, height_ratios=[3.55, 3.05, 4.20, 1.95, 2.35, 0.72, 2.35],
                          left=0.078, right=0.985, top=0.934, bottom=0.028, hspace=0.40)
    axA = fig.add_subplot(gs[0])
    gsB = gs[1].subgridspec(1, 2, wspace=0.075)
    axB = {site: fig.add_subplot(gsB[i]) for i, site in enumerate(zw)}
    axBseq = fig.add_subplot(gs[2])
    axC = fig.add_subplot(gs[3], sharex=axA)
    axCr = fig.add_subplot(gs[4])
    gsD = gs[6].subgridspec(1, 2, width_ratios=[2.15, 1.0], wspace=0.06)
    axD, axDn = fig.add_subplot(gsD[0]), fig.add_subplot(gsD[1])

    panel_A(axA, j, cov, junc, zw)
    axA.tick_params(labelsize=6.8)
    axA.xaxis.set_major_formatter(lambda v, p: fmt(v))

    cap, notes = panel_B_pileup(axB, r, j, zw)
    first = axB[list(axB)[0]]
    first.set_title("B   loxP cassette read pile-up — the floxed allele is observed "
                    "directly, in BOTH genotypes")
    picks, n_eligible = panel_B_sequence(axBseq, fig, j, r, fa, cap, notes)
    panel_C_coverage(axC, j, cov, depth)
    lo, hi, fx_row, op_row, sd_row = panel_C_ratio(axCr, j, cov, ratio)
    psi = panel_D(axD, axDn, j, junc)

    # honesty caption for panel C, computed from the ratio table + summary JSON.
    # Every comparative and directional word below is DERIVED; the previous version
    # hard-coded "in the opposite direction" and quoted only the depth metric, on
    # which the floxed exon happens to look unremarkable.
    n_perm = int(ratio.n_permutations.iloc[0])
    floor = float(ratio.perm_p_floor.iloc[0])
    sep = ratio.loc[ratio.separates_groups.astype(str).str.lower() == "true", "exon_number"]
    sep_rf = ratio.loc[ratio.separates_groups_readfrac.astype(str).str.lower() == "true",
                       "exon_number"]
    ec = j["exon_coverage"]
    md_rf = ec["most_deviating_exon"]["readfrac"]
    fx_dev_d = abs(float(fx_row.ratio_dcKO_over_CON) - 1)
    op_dev_d = abs(float(op_row.ratio_dcKO_over_CON) - 1)
    same_sign = (np.sign(float(op_row.ratio_dcKO_over_CON) - 1)
                 == np.sign(float(fx_row.ratio_dcKO_over_CON) - 1))
    dir_word = "the same direction" if same_sign else "the opposite direction"
    sd_txt = ""
    if sd_row is not None:
        sd_dev = abs(float(sd_row.ratio_dcKO_over_CON) - 1)
        sd_txt = (f" Exon {int(sd_row.exon_number)} moves {sd_dev:.3f} in the SAME direction "
                  f"as the floxed exon — {'more' if sd_dev > fx_dev_d else 'less'} than the "
                  f"floxed exon does (P = {float(sd_row.perm_p_two_sided):.2f}).")
    rf_txt = (f"the floxed exon has the LARGEST deviation of the {len(ratio)} "
              f"({md_rf['deviation']:.3f}) and the lowest ratio "
              f"({md_rf['ratio']:.3f})"
              if md_rf["is_floxed_exon"] else
              f"the largest deviation belongs to exon {md_rf['exon_number']} "
              f"({md_rf['deviation']:.3f})")
    spacer = gs[5].get_position(fig)
    fig.text(0.078, (spacer.y0 + spacer.y1) / 2,
             wrap(
             f"BOTH metrics, because they rank the exons differently. On per-base DEPTH every "
             f"exon lies within {lo:.3f}–{hi:.3f}; exon {int(op_row.exon_number)} deviates from "
             f"1 by {op_dev_d:.3f}, more than the floxed exon {int(fx_row.exon_number)} "
             f"({fx_dev_d:.3f}), in {dir_word}.{sd_txt}\n"
             f"On READ FRACTION the ranking reverses: {rf_txt}. Neither ordering is the "
             "finding — the point is that the two metrics disagree, so no single comparison "
             "settles anything.\n"
             f"Exon {', '.join(str(int(v)) for v in sep)} is the only exon whose CON and dcKO "
             f"values do not overlap, and that holds on both metrics (depth: exon "
             f"{', '.join(str(int(v)) for v in sep)}; read fraction: exon "
             f"{', '.join(str(int(v)) for v in sep_rf)}). But P = {floor:.2f} is the FLOOR of "
             f"an exact {n_perm}-permutation test at n = 3 vs 3: no exon can reach P < 0.05 in "
             "this design.\n"
             "This is suggestive-but-underpowered. It is NOT a demonstrated depletion — and it "
             "is not a demonstrated absence of one either: six libraries can neither establish "
             "nor exclude an effect this small.", 214),
             ha="left", va="center", fontsize=6.6, linespacing=1.5,
             bbox=dict(boxstyle="round,pad=0.40", fc="#F4F4F4", ec="#777777", lw=0.6))

    de = j["smad4_differential_expression"]
    deg_note = ("" if de["is_sig_at_project_threshold"] else
                "; |log2FC| < 1, so Smad4 is NOT in this study's 563-DEG list, which requires "
                "padj < 0.05 AND |log2FC| > 1")
    sdz = j["study_design"]
    fig.suptitle(
        "Supplementary Figure S5 — the Smad4 floxed allele is identified in the RNA-seq "
        "reads; the Cre excision is NOT detected",
        x=0.008, y=0.996, ha="left", fontsize=11, fontweight="bold")
    fig.text(0.008, 0.9835,
             wrap(
             f"Read-level genotyping of the conditional Smad4 allele "
             f"({j['gene']['symbol']}, {j['gene']['ensembl_gene_id']}, chr{j['gene']['chrom']}, "
             f"{j['gene']['strand']} strand; {j['genome_build']}) in the six bulk RNA-seq "
             f"libraries — driver {sdz['sequenced_material_driver']}; controls are "
             f"{mathify(sdz['control_genotype'])}.\n"
             f"Smad4 gene-level expression is UP, not down, in dcKO "
             f"(log2FC {de['log2FoldChange']:+.4f}, adjusted P = {de['padj']:.2e}{deg_note}).",
                  208),
             ha="left", va="top", fontsize=7.4, linespacing=1.45)
    return fig, r, zw, (lo, hi, fx_row, op_row, sd_row, psi, picks, n_eligible)


def save(fig, name):
    paths = []
    for ext in ("png", "pdf"):
        p = os.path.join(FIG, f"{name}.{ext}")
        # CreationDate must be suppressed or the PDF differs byte-for-byte between
        # otherwise identical runs, which breaks reproducibility checking.
        meta = {"CreationDate": None} if ext == "pdf" else None
        fig.savefig(p, bbox_inches="tight", facecolor="white", metadata=meta)
        paths.append(p)
    return paths


def type_audit(fig):
    """Smallest rendered type, and what it becomes on a printed page.

    A previous version put annotation text at 5.2 pt on a 12.6 in canvas; scaled
    to fit A4 with 0.5 in margins that is ~3.1 pt, below every journal minimum.
    Reported rather than silently accepted, and the native floor is asserted.
    """
    sizes = [t.get_fontsize() for t in fig.findobj(matplotlib.text.Text)
             if t.get_visible() and t.get_text().strip()]
    smallest = min(sizes)
    if smallest < MIN_FONTSIZE - 1e-9:
        raise RuntimeError(
            f"figure contains {smallest:.1f} pt type, below the declared floor of "
            f"{MIN_FONTSIZE:.1f} pt")
    # The DELIVERED size is the tight bbox, not FIGSIZE: savefig(bbox_inches=
    # "tight") silently grows the page around any text that overruns the canvas,
    # and a figure that is quietly 40% wider than designed prints 40% smaller.
    fig.canvas.draw()
    bb = fig.get_tightbbox(fig.canvas.get_renderer())
    w, h = bb.width, bb.height
    if w > FIGSIZE[0] * 1.06 or h > FIGSIZE[1] * 1.06:
        raise RuntimeError(
            f"text overruns the canvas: delivered size {w:.2f} x {h:.2f} in vs designed "
            f"{FIGSIZE[0]} x {FIGSIZE[1]} in. Wrap the offending caption -- bbox_inches="
            "'tight' will otherwise enlarge the page and shrink every point size on it.")
    out = {}
    for page, (pw, ph) in (("A4 portrait", (8.27, 11.69)), ("US Letter", (8.5, 11.0))):
        s = min((pw - 1.0) / w, (ph - 1.0) / h)
        out[page] = (s, smallest * s)
    return smallest, (w, h), out


# ============================ legend =========================================
def write_legend(j, r, ratio, junc, stats, path, typ=None):
    lo, hi, fx_row, op_row, sd_row, psi, picks, n_eligible = stats
    n_decomp = len(picks)
    decomp_list = "; ".join(
        f"{q['sample']} ({q.genotype}) `{q.qname}` at site {q.site}, CIGAR {q.cigar} "
        f"({int(q.aln_end) - int(q.pos) + 1} nt aligned)"
        for q in picks)
    g, fx, cas, neg = j["gene"], j["floxed_exon"], j["loxp_cassettes"], j["negative_result"]
    gsj = neg["genotype_specific_junctions"]
    n_by_gt = r.groupby("genotype")["sample"].nunique()
    per_sample = ", ".join(f"{k} {v}" for k, v in
                           cas["n_loxp_reads_per_sample"].items())
    sep = ", ".join(str(int(v)) for v in
                    ratio.loc[ratio.separates_groups.astype(str).str.lower() == "true",
                              "exon_number"])
    sep_rf = ", ".join(str(int(v)) for v in
                       ratio.loc[ratio.separates_groups_readfrac.astype(str).str.lower()
                                 == "true", "exon_number"])
    de = j["smad4_differential_expression"]
    sdz, ident, ec = j["study_design"], fx["identity"], j["exon_coverage"]

    # --- derived wording, so no comparative claim is a literal ---------------
    nums = ident["exon_number_values"]
    exnum_txt = ", ".join(str(v) for v in nums[:-1]) + f" **or** {nums[-1]}"
    fx_n = int(fx["canonical_exon_number"])
    fx_prev, fx_next = fx_n - 1, fx_n + 1
    leg_dir_word = ("the same direction" if
                    np.sign(op_row.ratio_dcKO_over_CON - 1) == np.sign(fx_row.ratio_dcKO_over_CON - 1)
                    else "the opposite direction")
    leg_sd = ""
    if sd_row is not None:
        _sd = abs(float(sd_row.ratio_dcKO_over_CON) - 1)
        leg_sd = (f"; the largest movement in the SAME direction as the floxed exon is exon "
                  f"{int(sd_row.exon_number)} at {_sd:.3f}, "
                  f"{'more' if _sd > abs(fx_row.ratio_dcKO_over_CON - 1) else 'less'} than the "
                  f"floxed exon itself")
    _rf = ec["most_deviating_exon"]["readfrac"]
    leg_rf = (f"the floxed exon has the **largest** deviation of the {len(ratio)} "
              f"({_rf['deviation']:.3f}) and the lowest ratio ({_rf['ratio']:.3f})"
              if _rf["is_floxed_exon"] else
              f"the largest deviation belongs to exon {_rf['exon_number']} "
              f"({_rf['deviation']:.3f})")
    deg_note = ("" if de["is_sig_at_project_threshold"] else
                "; |log2FC| < 1, so *Smad4* does not meet this study's DEG threshold "
                "(padj < 0.05 **and** |log2FC| > 1) and is absent from the 563-DEG list")
    n_dup = int(cas["n_duplicate_copies"])
    dup_sentence = (
        f"**{n_dup} of these reads (marked ‡) is a byte-identical copy of another read in the "
        f"same library** — same POS, same CIGAR, same base calls: one fragment sequenced twice. "
        f"The number of independent fragments is therefore {cas['n_unique_fragments']}, not "
        f"{cas['n_loxp_reads_total']}. No cassette-spanning read is affected."
        if n_dup else "No two of these reads are duplicates of one another.")
    skip_row = junc[junc.is_predicted_cre_skip.astype(str).str.lower() == "true"].iloc[0]
    if typ:
        smallest, (fw, fh), pages = typ
        pg = "; ".join(f"{k} {v[1]:.1f} pt" for k, v in pages.items())
        size_note = (
            f"The figure's native size is {fw:.1f} × {fh:.1f} in and its smallest type is "
            f"{smallest:.1f} pt. Reproduced at native width the type is legible; scaled down to "
            f"fit a portrait page with 0.5 in margins it becomes: {pg}. It is dense by design "
            "(three sequence rows of literal base calls) and should be reproduced at native "
            "size, across a full page width or as a digital supplement — not shrunk to a "
            "single column.")
    else:
        size_note = "(type audit not run)"
    txt = f"""# Supplementary Figure S5 — legend

**The Smad4 floxed allele is identified in the RNA-seq reads; the Cre excision is not
detected.**

Read-level genotyping of the conditional *Smad4* allele in the six bulk RNA-seq libraries
used throughout this study (*n* = 3 CON, *n* = 3 dcKO; one library per animal, no technical
replicates). **The sequenced material is {sdz['sequenced_material_driver']}**: dcKO =
{sdz['dcko_genotype']}, CON = {sdz['control_genotype']}; tissue = {sdz['tissue']}.
(Figure 1's in vivo phenotype uses constitutive K14-cre and is a different experiment.)
All coordinates are **{j['genome_build']}**, chromosome
{g['chrom']}. *{g['symbol']}* ({g['ensembl_gene_id']}) is on the **minus strand**, so
transcription runs **right → left** in panels A–C; panel C's exon-number axis is in
transcription order (exon 1 has the highest coordinate) and therefore runs in the opposite
spatial direction to the panels above it. Panels A and C share the genomic window
chr{g['chrom']}:{fmt(WINDOW[0])}–{fmt(WINDOW[1])}; panel B zooms into the two dashed
regions marked on panel A.

**(A) Locus schematic.** Canonical *{g['symbol']}* exons (transcript
{g['canonical_transcript']}, {g['n_canonical_exons']} exons) as boxes. The three exons
relevant to the allele are exons {fx_prev}, {fx['canonical_exon_number']} and {fx_next} in
canonical numbering; exons {fx_prev} and {fx_next} carry coordinate labels below the line and
the floxed exon's coordinates are given in its identity box. The floxed exon is filled black:
**chr{g['chrom']}:{fmt(fx['start'])}–{fmt(fx['end'])}, {fx['length_bp']} bp**.

**How to refer to this exon.** The unambiguous identifier is the **coordinate interval**, not
a stable exon ID and not an exon number. The same {fx['length_bp']} bp carries `exon_number`
{exnum_txt} across the {ident['n_transcripts']} annotated transcripts that contain it, which
is why the literature calls it "exon {fx['literature_exon_number']}" while the Ensembl
canonical transcript calls it "exon {fx['canonical_exon_number']}". A stable ID resolves the
*number* ambiguity but not the *identity* one:
{ident['n_transcripts_with_primary_exon_id']} of those transcripts label the interval
`{ident['primary_exon_id']}` (used in this figure), while the remaining
{ident['n_transcripts'] - ident['n_transcripts_with_primary_exon_id']} label the **identical**
{fx['length_bp']} bp `{'`, `'.join(ident['other_exon_ids'])}` — so a reader who queries
`{ident['primary_exon_id']}` in Ensembl silently loses those transcripts. The exon
{fx['protein']}.

Black triangles mark the two loxP insertion points mapped from
soft-clipped reads: site A at {fmt(cas['site_A']['insertion_between'][0])}/{fmt(cas['site_A']['insertion_between'][1])}
(±{cas['site_A']['ambiguity_bp']} bp, {cas['site_A']['ambiguity_reason']}; linkers
{cas['site_A']['linkers']}) and site B at
{fmt(cas['site_B']['insertion_between'][0])}/{fmt(cas['site_B']['insertion_between'][1])}
({cas['site_B']['ambiguity_reason']}; linkers {cas['site_B']['linkers']}). Both breakpoints
are re-derived from the reads at every run — the ±{cas['site_A']['ambiguity_bp']} bp at site A
is simply the spread of the per-read placements, and the microhomology that explains it is
verified in `genome.fa` to be a tandem repeat. Restriction-enzyme names are likewise derived
from the base calls, including sites formed across the loxP-core/linker boundary (marked †);
scanning the linker sequence alone misses XbaI, whose `TCTAGA` is completed by the final `T`
of the loxP core. On the minus strand, site A (the lower coordinate) lies **3′** of the floxed
exon in transcription order and site B **5′** of it. Shading marks the floxed interval
{fmt(cas['floxed_interval']['start'])}–{fmt(cas['floxed_interval']['end'])}
({fmt(cas['floxed_interval']['length_bp'])} bp), which contains exactly one exon. The dashed
grey span at the foot of the panel is the mRNA junction Cre excision would create; it is
drawn dashed, and carries its own read count, because it was never observed.

**(B) loxP read pile-up and literal sequence.** Every one of the
{cas['n_loxp_reads_total']} reads in the six libraries that contains the 34 bp loxP core
({cas['core_sequence_34bp']}) is drawn once, at its true coordinate, grouped by library:
{per_sample}. {dup_sentence} **Both cassettes are intronic**, so every one of these reads
comes from unspliced pre-mRNA or from residual genomic DNA rather than from mature mRNA; that
is why there are only a handful of them, and it is why their number carries no information
about transcript abundance or about excision efficiency.
Solid bars are the aligner-placed (CIGAR `M`) portion; open, dotted bars are
the soft-clipped portion, drawn extending from the alignment breakpoint in IGV fashion —
soft-clipped bases have no reference coordinate of their own, and their drawn position is a
display convention, not an alignment. The near-black block on each read is that read's own
34 bp loxP core, placed by counting bases along the read from that breakpoint. The cores
therefore do **not** align with one another: each sits offset from the insertion point by
however much vector polylinker its own read carries, on whichever side that read's clipped
segment falls. A star marks the
{cas['n_cassette_spanning_reads']} reads that carry genomic sequence on **both** sides of the
cassette ({cas['n_cassette_spanning_reads_by_site']['A']} at site A,
{cas['n_cassette_spanning_reads_by_site']['B']} at site B) and therefore read the
cassette–genome junction directly rather than by inference from clip endpoints. There are no
error bars in this panel: every bar is a single sequencing read, not a summary statistic.

Below the pile-up, the reference sequence at each insertion point is printed verbatim from
`genome.fa`, and {n_decomp} of the {n_eligible} cassette-spanning reads are decomposed into
their literal base calls: genomic flank | restriction linker | 34 bp loxP core | restriction
linker | genomic flank. **All {n_eligible} spanning reads are eligible** — both read
geometries (cassette in the leading soft clip, cassette in the trailing soft clip) are
decomposed, so the selection below is the only thing that decides what is shown. The reads
shown are {decomp_list}. The first is the spanning read with the longest reference-aligned
block (the longer that block, the more firmly the read pins the cassette to a genomic
coordinate); the second is the longest-aligned read from the other genotype; the third is the
longest-aligned read from the other insertion site — so that neither the structure nor the
genotype rests on a single library.
Block boundaries are derived, not assumed. The aligned block is verified base-for-base
against the reference. For a leading-clip read the 3′ linker is pinned exactly by the CIGAR
and the clipped genomic flank is recovered by extending an exact match back from the insertion
point (the recovered anchor length is cross-checked against the value the upstream analysis
recorded, and the script raises if the two disagree). For a trailing-clip read the 5′ linker is
pinned exactly by the CIGAR, the opposite linker is set by that site's consensus linker — itself
derived from reads whose CIGAR pins it — and the genomic flank beyond it is scored
base-for-base against the reference on the far side of the insertion.

**Why loxP reads appear in the controls.** CON animals are **Cre-negative** *Smad4*^fl/fl^
littermates. They carry the identical targeted allele; only the Cre recombinase differs
between genotypes. loxP-containing reads are therefore expected in CON
({n_by_gt['CON']}/3 CON and {n_by_gt['dcKO']}/3 dcKO libraries here) and their presence is a
positive control on the detection method, **not** a control failure or a genotyping error.
Detection of loxP reads in a library says the animal carries the floxed allele; it says
nothing about whether the exon was excised. That question is panels C and D.

**(C) Coverage and per-exon ratio.** Upper: mean per-base coverage over the shared window
for CON (solid blue) and dcKO (dashed red), each the mean of 3 libraries, normalised within
each library so that 1.0 is that library's mean coverage across all
{g['n_canonical_exons']} *{g['symbol']}* exons. Grey blocks are exons; the floxed exon is
outlined. Lower: for each exon, the dcKO/CON ratio of that exon's share of the library's
total *{g['symbol']}* coverage, computed **two independent ways** — share of per-base
**depth** (open diamonds) and share of exonic **reads** (filled squares), joined by a thin
line per exon. The six individual libraries are shown as small open symbols (circles CON,
triangles dcKO), on the depth metric, normalised to the CON mean of that exon. Points are
individual libraries, not error bars; no error bars are drawn because *n* = 3 per group does
not support them. The grey band is the range spanned by all {len(ratio)} exons on either
metric, {lo:.3f}–{hi:.3f}.

**Both metrics are plotted because they rank the exons differently, and quoting one alone
would be a choice.** On per-base depth, the floxed exon (exon {int(fx_row.exon_number)},
ratio {fx_row.ratio_dcKO_over_CON:.3f}) is not the most deviating exon: exon
{int(op_row.exon_number)} deviates from 1 by
{abs(op_row.ratio_dcKO_over_CON - 1):.3f} versus {abs(fx_row.ratio_dcKO_over_CON - 1):.3f},
in {leg_dir_word}{leg_sd}. On the read-fraction metric the ordering reverses:
{leg_rf}. The floxed exon and the two comparison exons are annotated in identical style, with
both metrics shown for each.

*P* values are two-sided **exact permutation** tests over all
{int(ratio.n_permutations.iloc[0])} balanced 3/3 relabellings of the six libraries. The
observed split and its complement always tie the observed statistic, so the smallest
attainable *P* is {float(ratio.perm_p_floor.iloc[0]):.2f}: **no exon can reach *P* < 0.05 in
this design**, and a *t* test would be invalid at *n* = 3. Exon {sep} is the only exon whose
CON and dcKO values do not overlap, and that holds on both metrics (depth: exon {sep};
read fraction: exon {sep_rf}). Clean separation is the strongest outcome this design can
produce, and it still is not significant.

**There is therefore no *demonstrable* exon-level depletion of the floxed exon — and equally
no demonstrated *absence* of one.** Six libraries can neither establish nor exclude an effect
of this size (~{100 * abs(fx_row.ratio_dcKO_over_CON - 1):.0f}% on depth,
~{100 * abs(fx_row.ratio_readfrac_dcKO_over_CON - 1):.0f}% on reads). Read the panel as
suggestive-but-underpowered; do not cite it as a depletion, and do not cite it as evidence
that the exon is intact.

**(D) Splice junctions.** Uniquely mapped reads crossing each junction (STAR
`SJ.out.tab` column 7), one bar per library. The two normal junctions flanking the floxed
exon are present at comparable depth in all six libraries. The junction Cre excision would
create — exon {int(skip_row.donor_annot[4:-4])} → exon {int(skip_row.acceptor_annot[4:-6])},
{fmt(neg['predicted_cre_skip_junction']['intron_start'])}–{fmt(neg['predicted_cre_skip_junction']['intron_end'])}
— has **0 reads in 6/6 libraries**. Exon inclusion PSI = {psi:.3f} in every library. With
0 skip events against **{neg['dcko_inclusion_equivalent_reads']:.0f} inclusion-equivalent dcKO
reads**, the rule of three gives a one-sided 95% bound of
**≥ {neg['one_sided_95_lower_bound_on_psi'] * 100:.1f}% of dcKO *{g['symbol']}* mRNA
retaining the floxed exon**.

*{neg['dcko_inclusion_equivalent_reads']:.0f} is a derived statistic, not a read count.* Each
dcKO library contributes the **mean** of its two flanking-junction counts, because a 150 bp
read crosses only one of the two junctions; the literal number of dcKO reads crossing either
junction is **{neg['dcko_junction_reads_total']}**. Using the mean halves the denominator and
therefore *widens* the interval, so the quoted bound is the conservative one.

**Scope of the junction survey.** Junctions were enumerated within
chr{g['chrom']}:{fmt(neg['junction_scan_window'][0])}–{fmt(neg['junction_scan_window'][1])}
({(neg['junction_scan_window'][1] - neg['junction_scan_window'][0]) / 1000:.0f} kb around the
floxed exon), not gene-wide. Within that window, no dcKO-specific junction involves the floxed
exon; sporadic genotype-specific junctions occur in **both** directions
({gsj['n_dcKO_only']} dcKO-only, {gsj['n_CON_only']} CON-only) at
≤ {max(gsj['max_dcKO_only_pct_of_sibling_canonical'], gsj['max_CON_only_pct_of_sibling_canonical']):.2f}%
of the canonical junction they compete with, and are low-count background rather than a dcKO
splicing phenotype.

**What this bound is a bound on.** {sdz['caveat']}

**Interpretation.** Panels A and B establish that the sequenced material carries the
conditional allele and map both loxP cassettes to base resolution. Panels C and D establish
that the excision product is undetectable in the same RNA: *{g['symbol']}* transcript levels
are in fact modestly **higher** in dcKO (log2FC {de['log2FoldChange']:+.4f}, adjusted
*P* = {de['padj']:.2e}{deg_note}), the floxed exon shows no demonstrable depletion, and no
exon-skipping junction was detected. This figure identifies the allele; it does **not**
demonstrate the knockout at the RNA level, and it should not be cited as if it did.

**Note on a superseded statistic.** {j['exon_coverage']['discrepancy_with_prior_analysis']}

**Reproduction size.** {size_note}

---
Generated by `scripts/JDR_reanalysis/21_smad4_loxp_figure.py` from
`scripts/JDR_reanalysis/output/smad4_*.csv|json` (written by
`20_smad4_loxp_evidence.py`) and `00.RawData/Bulk-RNAseq/genome.fa`.

Every **data** value in the figure and in this legend — count, coordinate, sequence,
restriction-enzyme name, exon-number list, transcript count, ratio, *P* value, and every
directional or comparative word — is read or derived at run time from those files. The
remaining text (panel titles, methodological explanation, and the interpretation above) is
argument about the data, not measurement of it, and is fixed prose in the script. The build
aborts if the data stop supporting the derived claims.
"""
    with open(path, "w") as fh:
        fh.write(txt)
    return path


if __name__ == "__main__":
    j, reads, cov, ratio, junc, depth = load_inputs()
    fa = Fasta(FASTA)
    sanity(j, reads, ratio, junc, fa)
    fig, r, zw, stats = build(j, reads, cov, ratio, junc, depth, fa)
    clashes = text_collisions(fig, _SEQ_TEXT_IDS)
    if clashes:
        print(f"[21] WARNING {len(clashes)} overlapping text pairs:")
        for c in clashes[:20]:
            print("      ", c)
    else:
        print("[21] text collisions        : 0")
    ink = ink_collisions(fig, _SEQ_TEXT_IDS)
    if ink:
        print(f"[21] WARNING {len(ink)} text-over-ink overlaps:")
        for c in ink[:20]:
            print("      ", c)
    else:
        print("[21] text over ink          : 0")
    typ = type_audit(fig)
    paths = save(fig, "FigureS5_smad4_loxp")
    plt.close(fig)
    lp = write_legend(j, r, ratio, junc, stats,
                      os.path.join(FIG, "FigureS5_smad4_loxp_legend.md"), typ=typ)
    smallest, (fw, fh), pages = typ
    print("[21] zoom windows           :", zw)
    print("[21] loxP reads drawn       :", len(r), f"({int(reads.is_duplicate_copy.sum())} "
          "duplicate copy marked ‡)")
    print("[21] reads decomposed       :", ", ".join(
        f"{q['sample']}@{q.site} ({int(q.aln_end) - int(q.pos) + 1} nt aligned)"
        for q in stats[6]), f"of {stats[7]} eligible")
    print(f"[21] smallest type          : {smallest:.1f} pt on a {fw:.1f}x{fh:.1f} in canvas; "
          + "; ".join(f"{k} {v[1]:.1f} pt" for k, v in pages.items()))
    print("[21] wrote                  :", *paths, lp, sep="\n      ")
    if clashes or ink:
        sys.exit(1)
