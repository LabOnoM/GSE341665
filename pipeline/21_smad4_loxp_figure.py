#!/usr/bin/env python3
"""
21_smad4_loxp_figure.py  --  Smad4 dcKO re-analysis pipeline (step 21)

Render Supplementary Figure S5: what the RNA-seq reads show about the Smad4
conditional allele, and why they cannot decide whether it was excised.

THE ONE JOB OF THIS FIGURE
    Demonstrate the LIMIT OF THE ASSAY.
The allele is identified; the excision product is undetectable; and neither
coverage nor junctions can distinguish "little recombination" from "recombination
followed by decay of the frameshifted transcript". A reader must not come away
believing the knockout was demonstrated at the RNA level, nor that the exon was
shown to be intact.

The figure shows DATA. The argument lives in the legend and in Supplementary
Methods B7 -- not in prose boxes drawn inside the axes.

Every DATA value -- number, coordinate, count, sequence, restriction-site name,
direction word -- is read or derived at run time from the outputs of
20_smad4_loxp_evidence.py and from genome.fa. If an input is missing the script
RAISES rather than substituting a value (project rule -- this repository once
shipped a figure built on fabricated fallback values).

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
FIGSIZE = (11.0, 11.0)

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
                f"summary JSON has no floxed_exon['{key}'] -- rerun step 20. The stable exon "
                "ID the figure uses to find the floxed exon must stay cross-checked.")
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
    ax.set_ylim(-3.05, 1.90)
    ax.set_yticks([])
    for sp in ("left", "right", "top"):
        ax.spines[sp].set_visible(False)
    ax.spines["bottom"].set_position(("data", -3.05))

    # ---- floxed interval -------------------------------------------------
    # Placed BELOW the panel-B zoom brackets: when the bracket rectangles were
    # tall enough to reach this label, their dashed edges struck through it.
    fi = cas["floxed_interval"]
    ax.axvspan(fi["start"], fi["end"], color="#F5E7CC", zorder=0)
    ax.annotate("", xy=(fi["start"], -1.42), xytext=(fi["end"], -1.42),
                arrowprops=dict(arrowstyle="<->", lw=0.8, color="#8A6D2F"))
    ax.text((fi["start"] + fi["end"]) / 2, -1.58,
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
        ax.annotate(f"exon {n}   {fmt(e.start)}–{fmt(e.end)}   ({fmt(e.exon_length)} bp)",
                    xy=((e.start + e.end) / 2, -0.28),
                    xytext=((e.start + e.end) / 2, -2.05), ha="center", va="top",
                    fontsize=6.3, linespacing=1.3,
                    arrowprops=dict(arrowstyle="-", lw=0.4, color="#666666",
                                    shrinkA=1, shrinkB=1))

    # ---- loxP triangles + their labels ----------------------------------
    # Coordinates only. The restriction linkers are drawn where they are DATA --
    # labelled on the decomposed read in panel B -- not narrated here.
    for site, key in (("A", "site_A"), ("B", "site_B")):
        sc = cas[key]
        x = (sc["insertion_between"][0] + sc["insertion_between"][1]) / 2
        loxp_triangle(ax, x, 0.30, w=(W1 - W0) * 0.0055, h=0.40)
        amb = ("" if sc["ambiguity_bp"] == 0 else f"  ±{sc['ambiguity_bp']} bp")
        ax.text(x, 0.86, f"loxP site {site}\n{fmt(sc['insertion_between'][0])}"
                         f"/{fmt(sc['insertion_between'][1])}{amb}",
                ha="center", va="bottom", fontsize=6.1, linespacing=1.35)

    # ---- predicted Cre product: a dashed grey span carrying its own zero --
    skip = junc[junc.is_predicted_cre_skip.astype(str).str.lower() == "true"].iloc[0]
    scols = sample_cols(junc)
    nlib = len(scols)
    a, b = float(skip.intron_start), float(skip.intron_end)
    ax.plot([a, b], [-2.85, -2.85], color="#7A7A7A", lw=1.1, ls=(0, (4, 2)), zorder=5)
    for xx in (a, b):
        ax.plot([xx, xx], [-2.97, -2.73], color="#7A7A7A", lw=1.1, zorder=5)
    ax.text((a + b) / 2, -2.79,
            f"exon {int(skip.donor_annot[4:-4])} → exon {int(skip.acceptor_annot[4:-6])} "
            f"junction Cre would create: {int(skip[scols].sum())} reads in {nlib}/{nlib}",
            ha="center", va="bottom", fontsize=6.6, color="#4D4D4D", style="italic")

    # ---- transcription direction ----------------------------------------
    ax.annotate("", xy=(W0 + (W1 - W0) * 0.865, 1.35),
                xytext=(W0 + (W1 - W0) * 0.995, 1.35),
                arrowprops=dict(arrowstyle="-|>", lw=1.3, color="k"))
    ax.text(W0 + (W1 - W0) * 0.855, 1.35, "transcription 5′→ 3′",
            ha="right", va="center", fontsize=7.0)

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

    ax.set_title(f"A   Smad4 locus, chr{gene['chrom']}:{fmt(WINDOW[0])}–{fmt(WINDOW[1])}, "
                 f"minus strand   —   floxed exon (black) = "
                 f"chr{gene['chrom']}:{fmt(fx['start'])}–{fmt(fx['end'])}, "
                 f"{fmt(fx['length_bp'])} bp")


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
                 fontsize=6.4, color=col, fontweight="bold")
    for gt, col, lo, hi in (("CON", BLUE, 0, 2), ("dcKO", RED, 3, 5)):
        rows = r.loc[r.genotype == gt, "row"]
        axL.annotate("", xy=(-0.088, rows.min() - 0.45), xytext=(-0.088, rows.max() + 0.45),
                     xycoords=axL.get_yaxis_transform(), textcoords=axL.get_yaxis_transform(),
                     arrowprops=dict(arrowstyle="-", lw=2.0, color=col))

    # ONE line of plain annotation plus the symbol key. No argument, no box:
    # the reasons loxP reads appear in CON, and why their number is uninformative,
    # are stated once in the legend.
    n_by_gt = r.groupby("genotype")["sample"].nunique()
    cs = j["loxp_cassettes"]
    n_span = int(cs["n_cassette_spanning_reads"])
    n_dup = int(cs["n_duplicate_copies"])
    dup_key = f";  ‡ = duplicate copy ({n_dup} of {len(r)})" if n_dup else ""
    caption = (f"All {len(r)} loxP-containing reads, {n_by_gt['CON']}/3 CON and "
               f"{n_by_gt['dcKO']}/3 dcKO libraries.   "
               f"★ = spans the cassette–genome junction ({n_span} reads){dup_key}.")
    return caption


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


def panel_B_sequence(ax, fig, j, reads, fa, caption):
    cas, chrom = j["loxp_cassettes"], j["gene"]["chrom"]
    core = cas["core_sequence_34bp"]
    ax.set_xticks([])          # emptied as well as hidden: an axis that is only
    ax.set_yticks([])          # set_axis_off() still reports stale tick-label extents
    ax.set_axis_off()
    ax.set_xlim(0, 1)
    ax.set_ylim(0, 1)

    ax.text(0.0, 1.055, caption, ha="left", va="top", fontsize=6.8,
            transform=ax.transAxes)

    # ---- reference sequence at the two mapped breakpoints -------------------
    ax.text(0.0, 0.955, "Reference sequence at each mapped insertion point "
                        "(genome.fa); both sites are intronic",
            fontsize=7.0, fontweight="bold", va="top", transform=ax.transAxes)
    FLANK = 18
    for site, xa, xb in (("A", 0.015, 0.475), ("B", 0.525, 0.985)):
        s = cas[f"site_{site}"]
        p5, p3 = s["insertion_between"]
        left = fa.fetch(chrom, p5 - FLANK + 1, p5)
        right = fa.fetch(chrom, p3, p3 + FLANK - 1)
        amb = int(s["ambiguity_bp"])
        why = f"   (±{amb} bp)" if amb else "   (exact)"
        ax.text(xa, 0.868,
                f"loxP site {site}   chr{chrom}:{fmt(p5)} ▾ {fmt(p3)}" + why,
                fontsize=6.5, va="top", fontweight="bold", transform=ax.transAxes)
        fs = mono_fontsize(fig, ax, 2 * FLANK, xspan_frac=(xb - xa), cap=9.0)
        segs = [dict(seq=left, fc="#EFEFEF", tc="k"),
                dict(seq=right, fc="#DCDCDC", tc="k")]
        spans = draw_seq_row(ax, xa, xb, 0.735, 0.072, segs, fs, edge="#9A9A9A")
        mid = spans[0][2]
        ax.plot([mid, mid], [0.685, 0.785], color="k", lw=1.2,
                transform=ax.transAxes, zorder=5)
        ax.text(mid, 0.673, "▲ cassette inserted here", ha="center", va="top",
                fontsize=6.2, transform=ax.transAxes)
        ax.text(xa, 0.673, f"5′ {fmt(p5 - FLANK + 1)}", ha="left", va="top",
                fontsize=6.0, color="#666666", transform=ax.transAxes)
        ax.text(xb, 0.673, f"{fmt(p3 + FLANK - 1)} 3′", ha="right", va="top",
                fontsize=6.0, color="#666666", transform=ax.transAxes)

    # ---- literal decomposition of ONE cassette-spanning read ----------------
    ax.text(0.0, 0.560,
            "One cassette-spanning read, decomposed into its literal base calls",
            fontsize=7.0, fontweight="bold", va="top", transform=ax.transAxes)

    # ALL cassette-spanning reads are eligible (both read geometries are drawable);
    # the single read shown is the one with the longest reference-aligned block, i.e.
    # the one that pins the cassette to a genomic coordinate most firmly. The count
    # check stays as a data-integrity guard: the drawn read must come from the set
    # step 20 calls cassette-spanning.
    cand = reads[reads.spans_cassette].copy()
    if cand.empty:
        raise RuntimeError("no cassette-spanning read available to draw")
    n_eligible = len(cand)
    if n_eligible != int(j["loxp_cassettes"]["n_cassette_spanning_reads"]):
        raise RuntimeError("the spanning-read table and the summary JSON disagree on how "
                           "many reads span the cassette")
    cand["aln_len"] = cand.aln_end - cand.pos + 1
    cand = cand.sort_values(["aln_len", "qname"], ascending=[False, True])
    picks = [cand.iloc[0]]

    for row, ytop in zip(picks, [0.320]):
        s = cas[f"site_{row.site}"]
        segs, e5, k5 = decompose_spanning(row, core, fa, chrom, s)
        n = sum(len(x["seq"]) for x in segs)
        fs = mono_fontsize(fig, ax, n, xspan_frac=0.985, cap=8.0)
        ax.text(0.0, ytop + 0.105,
                f"{row['sample']} ({row.genotype})  read {row.qname}   POS {fmt(row.pos)}  "
                f"CIGAR {row.cigar}   → site {row.site}   [{n} nt]",
                fontsize=6.4, va="center", transform=ax.transAxes,
                color=BLUE if row.genotype == "CON" else RED, fontweight="bold")
        spans = draw_seq_row(ax, 0.0075, 0.9925, ytop, 0.085, segs, fs, edge="#8A8A8A")
        # staggered labels so neighbouring captions cannot collide
        for i, (sp, sg) in enumerate(zip(spans, segs)):
            dy = 0.055 if i % 2 == 0 else 0.135
            xl = min(max(sp[0], 0.045), 0.955)     # keep captions inside the panel
            ax.annotate(sg["lab"], xy=(sp[0], ytop - 0.045),
                        xytext=(xl, ytop - dy), ha="center", va="top",
                        fontsize=6.0, linespacing=1.3,
                        arrowprops=dict(arrowstyle="-", lw=0.35, color="#777777",
                                        shrinkA=0.5, shrinkB=0.5))
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
    ax.set_ylabel("coverage / mean\nSmad4 exonic", fontsize=7.0)
    ax.tick_params(labelsize=6.8)
    ax.xaxis.set_major_formatter(lambda v, p: fmt(v))
    ax.legend(frameon=False, fontsize=6.6, loc="upper center", ncol=2,
              handlelength=2.4, borderaxespad=0.1)
    ax.set_xlabel(f"chr{j['gene']['chrom']} coordinate ({j['genome_build'].split()[0]})",
                  fontsize=7.0, labelpad=2)
    ax.set_title("C   Per-base coverage across the locus, and the per-exon dcKO/CON "
                 "ratio on both normalisations")


def panel_C_ratio(ax, j, cov, ratio):
    """Per-exon dcKO/CON ratio -- BOTH metrics, never one.

    Step 20 computes the ratio two ways (share of per-base exonic DEPTH, share of
    exonic READS) and says in the CSV header that both are given "so that neither
    the figure nor the text can pick the flattering one". Dropping one would be
    cherry-picking, and the two metrics rank the exons differently -- which is
    itself part of why coverage cannot settle the question.
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

    # IDENTICAL annotation style for the floxed exon and for ONE comparison exon:
    # the exon that deviates most on depth. No arrow, no asterisk, no emphasis.
    dev = (r.ratio_dcKO_over_CON - 1.0).abs()
    i_fx = int(r.index[r.exon_id == fx_id][0])
    i_op = int(dev[r.exon_id != fx_id].idxmax())
    for i in sorted({i_fx, i_op}):
        e = r.loc[i]
        tag = " (floxed)" if e.exon_id == fx_id else ""
        # one line, parked in the empty band below every plotted library point
        ax.annotate(f"exon {int(e.exon_number)}{tag}   depth {e.ratio_dcKO_over_CON:.3f}   "
                    f"reads {e.ratio_readfrac_dcKO_over_CON:.3f}",
                    xy=(i, min(e.ratio_dcKO_over_CON, e.ratio_readfrac_dcKO_over_CON)),
                    xytext=(i, 0.720),
                    ha="center", va="bottom", fontsize=6.0, linespacing=1.3,
                    arrowprops=dict(arrowstyle="-", lw=0.4, color="#888888",
                                    shrinkA=1, shrinkB=2))
    ax.set_xticks(x)
    ax.set_xticklabels([str(int(v)) for v in r.exon_number])
    ax.set_xlim(-0.7, len(r) - 0.3)
    ax.set_ylim(0.70, 1.30)
    n_perm = int(ratio.n_permutations.iloc[0])
    floor = float(ratio.perm_p_floor.iloc[0])
    ax.set_xlabel(f"Smad4 canonical exon (5′→3′)   —   P = {floor:.2f} is the floor of the "
                  f"exact {n_perm}-permutation test at n = 3 vs 3; no exon can reach P < 0.05",
                  fontsize=6.8, labelpad=2)
    ax.set_ylabel("dcKO / CON\ncoverage ratio", fontsize=7.0)
    ax.tick_params(labelsize=6.8)
    ax.legend(frameon=False, fontsize=6.2, loc="upper left", ncol=4, handletextpad=0.3,
              columnspacing=1.2, borderaxespad=0.15)
    ax.set_title("      per-exon dcKO / CON coverage ratio — all "
                 f"{len(r)} canonical exons, plotted identically, on BOTH metrics",
                 fontsize=7.8, fontweight="normal")
    return lo, hi, r.loc[i_fx], r.loc[i_op]


# ============================ panel D ========================================
def panel_D(ax, j, junc):
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
        ax.text(x0 + (len(samples) - 1) * pitch / 2, -ymax * 0.52, cap, ha="center",
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
    ax.set_title("D   Splice junctions at the floxed exon — the predicted Cre excision "
                 "junction has no supporting reads")
    # PSI, the rule-of-three bound and the survey-scope caveats live in the legend.
    psi = sorted(set(neg["psi_per_sample"].values()))
    return psi[0]


# ============================ assembly =======================================
def build(j, reads, cov, ratio, junc, depth, fa):
    r = read_layout(reads)
    zw = zoom_windows(r)

    fig = plt.figure(figsize=FIGSIZE)
    gs = fig.add_gridspec(6, 1, height_ratios=[4.40, 3.00, 3.40, 1.90, 2.20, 2.50],
                          left=0.098, right=0.972, top=0.912, bottom=0.055, hspace=0.50)
    axA = fig.add_subplot(gs[0])
    gsB = gs[1].subgridspec(1, 2, wspace=0.16)
    axB = {site: fig.add_subplot(gsB[i]) for i, site in enumerate(zw)}
    axBseq = fig.add_subplot(gs[2])
    axC = fig.add_subplot(gs[3], sharex=axA)
    axCr = fig.add_subplot(gs[4])
    axD = fig.add_subplot(gs[5])

    panel_A(axA, j, cov, junc, zw)
    axA.tick_params(labelsize=6.8)
    axA.xaxis.set_major_formatter(lambda v, p: fmt(v))

    cap = panel_B_pileup(axB, r, j, zw)
    first = axB[list(axB)[0]]
    first.set_title("B   loxP cassette read pile-up, both genotypes")
    picks, n_eligible = panel_B_sequence(axBseq, fig, j, r, fa, cap)
    panel_C_coverage(axC, j, cov, depth)
    lo, hi, fx_row, op_row = panel_C_ratio(axCr, j, cov, ratio)
    psi = panel_D(axD, j, junc)

    de = j["smad4_differential_expression"]
    sdz = j["study_design"]
    fig.suptitle(
        "Supplementary Figure S5 — RNA sequencing identifies the floxed Smad4 allele "
        "but cannot determine whether it was excised",
        x=0.008, y=0.998, ha="left", fontsize=11, fontweight="bold")
    fig.text(0.008, 0.9755,
             wrap(
             f"{j['gene']['symbol']} ({j['gene']['ensembl_gene_id']}), "
             f"chr{j['gene']['chrom']}, {j['gene']['strand']} strand, {j['genome_build']}. "
             f"Six bulk tooth-germ libraries, n = 3 dcKO vs 3 CON; driver "
             f"{sdz['sequenced_material_driver']}; controls "
             f"{mathify(sdz['control_genotype'])}. "
             f"Smad4 gene-level log2FC {de['log2FoldChange']:+.4f}, adjusted "
             f"P = {de['padj']:.2e}.", 165),
             ha="left", va="top", fontsize=7.4, linespacing=1.45)
    return fig, r, zw, (lo, hi, fx_row, op_row, psi, picks, n_eligible)


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
    if w < FIGSIZE[0] * 0.80 or h < FIGSIZE[1] * 0.80:
        raise RuntimeError(
            f"delivered size {w:.2f} x {h:.2f} in is far smaller than the designed "
            f"{FIGSIZE[0]} x {FIGSIZE[1]} in -- FIGSIZE no longer describes the figure, so "
            "every reproduction-size statement made about it would be wrong.")
    out = {}
    for page, (pw, ph) in (("A4 portrait", (8.27, 11.69)), ("US Letter", (8.5, 11.0))):
        s = min((pw - 1.0) / w, (ph - 1.0) / h)
        out[page] = (s, smallest * s)
    return smallest, (w, h), out


# ============================ legend =========================================
def write_legend(j, r, ratio, junc, stats, path):
    """Write the Figure S5 legend: 300-400 words, one paragraph per panel.

    Structure is fixed prose; every number, coordinate, sequence, ratio and
    P value in it is read or derived at run time from step 20's outputs. Anything
    needed for reproducibility but not for reading the figure -- the exon-numbering
    and stable-ID ambiguity, the breakpoint and linker derivations, the read-selection
    rule, the superseded prior statistic -- lives in Supplementary Methods B7 and is
    deliberately not duplicated here.
    """
    lo, hi, fx_row, op_row, psi, picks, n_eligible = stats
    g, fx, cas, neg = j["gene"], j["floxed_exon"], j["loxp_cassettes"], j["negative_result"]
    sdz = j["study_design"]
    n_by_gt = r.groupby("genotype")["sample"].nunique()

    # --- the drawn cassette-spanning read, named from the data ----------------
    q = picks[0]
    read_txt = f"{q['sample']}, `{q.cigar}` at {fmt(q.pos)}"

    # --- flanking-junction read range, straight out of the junction table -----
    scols = sample_cols(junc)
    incl = [tuple(int(v) for v in s.split("-")) for s in neg["inclusion_junctions"].split(";")]
    incl_vals = []
    for a, b in incl:
        row = junc[(junc.intron_start == a) & (junc.intron_end == b)]
        if row.empty:
            raise RuntimeError(f"junction {a}-{b} absent from smad4_junctions.csv")
        incl_vals.extend(int(v) for v in row.iloc[0][scols])
    j_lo, j_hi = min(incl_vals), max(incl_vals)

    skip_row = junc[junc.is_predicted_cre_skip.astype(str).str.lower() == "true"].iloc[0]
    skip_don = int(skip_row.donor_annot[4:-4])
    skip_acc = int(skip_row.acceptor_annot[4:-6])
    n_lib = len(scols)
    n_dup = int(cas["n_duplicate_copies"])
    dup_txt = (f"; \u2021 = a duplicate read, so {cas['n_unique_fragments']} independent "
               "fragments" if n_dup else "")
    _site = cas[f"site_{q.site}"]
    dagger = (" (\u2020 = a restriction site spanning the core/linker boundary)"
              if any(h["straddles_core"] for h in
                     _site["linker_5p_sites"] + _site["linker_3p_sites"]) else "")
    ambA = (f" (\u00b1{cas['site_A']['ambiguity_bp']} bp)"
            if cas["site_A"]["ambiguity_bp"] else " (exact)")
    ambB = (f" (\u00b1{cas['site_B']['ambiguity_bp']} bp)"
            if cas["site_B"]["ambiguity_bp"] else " (exact)")

    txt = f"""# Supplementary Figure S5 \u2014 legend

**RNA sequencing identifies the floxed *{g['symbol']}* allele but cannot determine whether it
was excised.** It shows what the reads establish and why they cannot settle excision or its
efficiency.

Six bulk tooth-germ libraries: *n* = 3 dcKO ({sdz['sequenced_material_driver']}) versus 3
{sdz['control_genotype']} controls. {j['genome_build']}, chr{g['chrom']}; *{g['symbol']}* is on
the **minus strand**, so transcription runs right \u2192 left in A\u2013C.

**(A)** Locus schematic. The floxed exon (black) is
**chr{g['chrom']}:{fmt(fx['start'])}\u2013{fmt(fx['end'])}, {fx['length_bp']} bp**; triangles are
the loxP insertion points, site A
{fmt(cas['site_A']['insertion_between'][0])}/{fmt(cas['site_A']['insertion_between'][1])}{ambA} and
site B {fmt(cas['site_B']['insertion_between'][0])}/{fmt(cas['site_B']['insertion_between'][1])}{ambB};
shading is the {fmt(cas['floxed_interval']['length_bp'])} bp floxed interval. Dashed = the exon
{skip_don} \u2192 exon {skip_acc} junction Cre would create, never observed.

**(B)** All {cas['n_loxp_reads_total']} reads carrying the 34 bp loxP core, in
{n_by_gt['CON']}/3 CON and {n_by_gt['dcKO']}/3 dcKO libraries (solid = aligned, open =
soft-clipped, black = core{dup_txt}). loxP reads are expected in Cre-negative CON, which carries
the identical allele; both cassettes are intronic, so read number is uninformative
about excision. Below: the reference sequence at each insertion point, and one cassette-spanning
read ({read_txt}) in literal base calls{dagger}.

**(C)** Per-base coverage across the locus, and for each of the {len(ratio)} canonical exons the
dcKO/CON ratio on both normalisations \u2014 per-base **depth** and exonic **reads**, with each
library plotted. All exons lie within {lo:.3f}\u2013{hi:.3f}; the floxed exon is
{fx_row.ratio_dcKO_over_CON:.3f} on depth (*P* = {fx_row.perm_p_two_sided:.2f}) and
{fx_row.ratio_readfrac_dcKO_over_CON:.3f} on reads
(*P* = {fx_row.perm_p_two_sided_readfrac:.2f}), and the metrics rank the exons differently.

**(D)** Uniquely mapped reads on the two flanking junctions and on the predicted excision
junction
(chr{g['chrom']}:{fmt(neg['predicted_cre_skip_junction']['intron_start'])}\u2013{fmt(neg['predicted_cre_skip_junction']['intron_end'])}),
one bar per library. Per-exon *P* values are in `output/smad4_exon_ratio.csv`. The Cre excision
is not detected: **0 reads, {n_lib}/{n_lib}
libraries**.

**The limit.** The allele is present; the excision product is undetectable \u2014 0 skip reads
in {n_lib}/{n_lib} libraries against {j_lo}\u2013{j_hi} reads on each flanking junction, PSI =
{psi:.3f} everywhere, one-sided 95% bound (rule of three)
\u2265 {neg['one_sided_95_lower_bound_on_psi'] * 100:.1f}% of dcKO *{g['symbol']}* mRNA
retaining the exon. Coverage cannot reach significance either:
*P* = {float(ratio.perm_p_floor.iloc[0]):.2f} is the floor of the exact
{int(ratio.n_permutations.iloc[0])}-permutation test at *n* = 3 vs 3. Both coverage and
junctions are measured on **surviving mRNA**, so neither can distinguish little recombination
from recombination followed by decay of the frameshifted transcript. Deciding excision and its
efficiency requires a DNA- or protein-level assay: genomic PCR across the floxed interval, an
amplicon spanning the excision junction, or SMAD4 immunostaining.

---
Generated by `21_smad4_loxp_figure.py` from the step-20 outputs and `genome.fa`; all values are
derived at run time. Details: Supplementary Methods B7.
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
                      os.path.join(FIG, "FigureS5_smad4_loxp_legend.md"))
    smallest, (fw, fh), pages = typ
    print("[21] zoom windows           :", zw)
    print("[21] loxP reads drawn       :", len(r), f"({int(reads.is_duplicate_copy.sum())} "
          "duplicate copy marked ‡)")
    print("[21] reads decomposed       :", ", ".join(
        f"{q['sample']}@{q.site} ({int(q.aln_end) - int(q.pos) + 1} nt aligned)"
        for q in stats[5]), f"of {stats[6]} eligible")
    print("[21] legend words           :", len(open(lp).read().split()))
    print(f"[21] smallest type          : {smallest:.1f} pt on a {fw:.1f}x{fh:.1f} in canvas; "
          + "; ".join(f"{k} {v[1]:.1f} pt" for k, v in pages.items()))
    print("[21] wrote                  :", *paths, lp, sep="\n      ")
    if clashes or ink:
        sys.exit(1)
