#!/usr/bin/env python3
"""
22_breeding_scheme.py -- render the mouse-line breeding and experimental-allocation
scheme as Supplementary Figure S6.

WHY THIS FIGURE EXISTS
----------------------
This study uses TWO Cre drivers on the same floxed allele, and the two drivers feed
different figures. That split has already been mis-recorded once (the RNA-sequenced
material was attributed to K14-cre in earlier drafts and corrected on 2026-07-27), and
the Figure 1 artwork still carries a column header that contradicts its own legend.
The scheme below draws the allocation stated in the Methods paragraph "Animals and
conditional deletion" so that the assignment is checkable at a glance rather than
reconstructed from prose.

SOURCES -- every fact drawn here is quoted or paraphrased from one of these files. Each
is a file in this repository; nothing is drawn on an undocumented recollection.

  [M] doc/Manuscript_JDR_Smad4_Identity_Metabolism.md
      "Animals and conditional deletion"  -> driver-to-figure allocation; Cre-negative
        littermate controls; 4-OHT 90 umol/L selected from a 50-300 umol/L titration;
        floxed-exon coordinates.
      "Phenotyping"  -> E16.5 / E18.5 / P0 dissection, Trowell culture, renal capsule
        recovered at 9 or 14 days, Faxitron MX-20, the four pathway inhibitors.
      "RNA sequencing and differential expression"  -> n = 3 vs 3, pooled germs, both
        sexes.
      Figure legends 1-3  -> which panels of Figures 2 and 3 carry which comparison.
  [B] 00.Manuscript/Submission_01_BBRC_2025-12/BBRC 投稿材料-2025.12/Manuscript-BBRC.docx
      (the original submission; the only place in this repository where the wet-lab
      detail is written down. Extract it with `unzip -p ... word/document.xml`.)
      2.1 Mouse lines  -> "K14-cre (JAX: 016230), K14-creER (JAX: 005107), and Smad4f/f
        (JAX: 017462) mice were maintained under specific pathogen-free (SPF)
        conditions. Genotyping was performed according to protocols provided by The
        Jackson Laboratory."
      2.2 Molar dissection and Trowell culture  -> "Since K14-cre; Smad4f/f embryos are
        embryonically lethal ..."; E16.5 K14-cre and P0 K14-creER dissection; DMEM +
        10% FBS (Lonsera S711-001S); the four inhibitors.
      2.4 Sub-renal implantation  -> "The 8-10-week-old WT male mice were anesthetized";
        no host strain is named, so none is drawn.
      2.5  -> Faxitron MX-20.
      Results  -> "Since K14-cre; Smad4f/f mice die shortly after birth [12], we
        transplanted molar germs from E18.5 ..."; 9 days (constitutive) and 14 days
        (inducible); the ten-point 50-300 umol/L titration that SELECTED 90 umol/L
        "for ... subsequent enamel development experiments"; and "we treated the dental
        epithelium of first molar germs from wild-type mice with various pathway
        inhibitors".
  [S5] scripts/JDR_reanalysis/21_smad4_loxp_figure.py and step 20  -> the floxed-exon
      coordinates and the frameshift, which are cross-referenced here, not repeated.

TWO POINTS ARE NOT RESOLVED BY THE SOURCES, and are marked on the panel itself with a
dagger rather than silently decided (a panel gets cropped away from its legend):
  dagger  the Figure 1 artwork column header says K14-Cre-ER while its own legend, the
          Results and the Methods all say K14-cre;
  ddagger the material for the inhibitor panels is [M] K14-creER germs but [B] wild-type
          germs.
Both are drawn as [M] states them, because [M] is the manuscript being submitted, and
both are flagged. Neither is asserted to be settled.

Group sizes are deliberately absent: they are unresolved [[n]] placeholders in every
affected figure legend and must be supplied by the authors. The panel points the reader
at the figure legends, which is where the field exists, not at the Methods, which only
says that group sizes are stated per experiment.

STRUCTURE -- one panel, one column per arm. The two arms are told apart by four cues that
all survive a greyscale conversion: POSITION (left/right); a cluster caption at the TOP of
each arm naming the driver and saying in words whether tamoxifen is used; the driver named
again inside every terminal figure box, so a box reproduced on its own still identifies
its mouse line; and a double-ruled tamoxifen node that exists in the inducible arm only.
Fill tone differs by about a dozen grey levels as a fifth, redundant cue -- it is not
load-bearing, and colour alone carries nothing.

The whole panel is width-bound: the type size everything ships at is EMBED_IN / panel
width, so a label a few characters too long, a cross node dragged sideways by its parents,
or a footer line wider than the drawing all cost every other label its legibility. Hence
MAX_LINE, the edge weights on the two backbones, and the footer carrying only the two
unresolved-point glosses -- housing, genotyping, floxed-exon coordinates and the E16.5
dissection live in the legend, where they cost nothing.

TWO ORDERING FACTS THE FIRST DRAFT OF THIS FIGURE GOT WRONG, now asserted rather than
proof-read: the titration that DEFINED 90 umol/L cannot descend from the 90 umol/L node,
and the inhibitor experiment is not a tamoxifen experiment at all. Both are enforced as
graph properties (TAM_NOT_DOWNSTREAM), so re-wiring either one aborts the build.

AUDITS (the script aborts rather than emitting a wrong figure):
  * graph reachability -- each terminal figure node must be reachable from the mouse
    lines that the Methods assign to it, and NOT reachable from the other driver
  * tamoxifen exposure as a graph property -- what must and must not sit downstream of
    the 90 umol/L node
  * tamoxifen may be named only in the inducible arm (or in an explicit negation)
  * every terminal figure box must name its driver inside the box
  * every node connected; every marker on the panel must be glossed in the footer, and
    every gloss must belong to a marker that is actually drawn
  * label lines are length-capped, and after rendering EVERY font size the .dot actually
    emits -- not just the largest -- is checked against the size it will reach on the
    page, together with the height that page would need

Output: figures/FigureS6_breeding_scheme.{png,pdf}
"""
import os
import re
import subprocess
from collections import deque

ROOT = os.path.dirname(os.path.dirname(os.path.dirname(os.path.abspath(__file__))))
FIG = os.path.join(ROOT, "scripts", "JDR_reanalysis", "figures")
os.makedirs(FIG, exist_ok=True)

NAME = "FigureS6_breeding_scheme"
EMBED_IN = 6.5          # width, in inches, at which 15_ embeds this figure
MIN_PT = 7.0            # refuse to ship ANY type smaller than this once scaled
MAX_PAGE_IN = 8.6       # a portrait text block
MAX_LINE = 38           # characters per node label line
MAX_FOOT = 112          # characters per footer line. Deliberately WIDER than a node
                        # line: the footer sits below the whole panel, so it does not set
                        # the panel width, and every extra footer line is 0.2 in of
                        # height the whole panel then pays for in scale.

# One point size for the body, one a step down for genotype superscripts and the footer.
# Both are checked against MIN_PT after rendering, so neither can quietly go sub-legible.
BODY_PT, SUP_PT, FOOT_PT, CLUSTER_PT = 13, 12, 12, 13

RED, INK = "#B2182B", "#222222"
# Constitutive arm is drawn cooler AND darker, inducible warmer AND lighter, so the two
# arms remain distinguishable after a greyscale conversion (~12 grey levels apart).
COOL_FILL, COOL_BG = "#DCE6F4", "#EBF1FA"      # Rec.601 luminance 229 / 240
WARM_FILL, WARM_BG = "#FAF0E4", "#FFFBF6"      # Rec.601 luminance 242 / 252

DAGGER, DDAGGER = "†", "‡"

# --------------------------------------------------------------------------- markup
# Labels are written as lists of lines in a two-character markup so that the plain
# text stays extractable for the audits below: **bold**, _italic_, ^{superscript}.

def _plain(line: str) -> str:
    line = re.sub(r"\*\*(.+?)\*\*", r"\1", line)
    line = re.sub(r"_(.+?)_", r"\1", line)
    return re.sub(r"\^\{(.+?)\}", r"\1", line)


def _html(line: str) -> str:
    line = line.replace("&", "&amp;").replace("<", "&lt;").replace(">", "&gt;")
    line = re.sub(r"\^\{(.+?)\}",
                  rf'<SUP><FONT POINT-SIZE="{SUP_PT}">\1</FONT></SUP>', line)
    line = re.sub(r"\*\*(.+?)\*\*", r"<B>\1</B>", line)
    return re.sub(r"_(.+?)_", r"<I>\1</I>", line)


def label(lines) -> str:
    return "<" + "<BR/>".join(_html(x) for x in lines) + ">"


# --------------------------------------------------------------------------- nodes
# (id, arm, lines, extra graphviz attributes)
#   arm "line" = parental mouse line (drawn outside both clusters)
#   arm "c"    = constitutive K14-cre arm      arm "i" = inducible K14-creER arm
LINE_BOX = 'shape=box, fillcolor="#EFEFEF", penwidth=1.4'
NOTE_BOX = ('shape=note, fillcolor="#F7F7F7", color="#AAAAAA", fontcolor="#3A3A3A", '
            f'style="filled,dashed", fontsize={FOOT_PT}')
TERMINAL = f'shape=box, fillcolor="white", color="{INK}", penwidth=2.2'
CROSS = 'shape=hexagon, fillcolor="white", margin="0.02,0.02"'

NODES = [
    # --- parental lines -----------------------------------------------------
    ("LINE_CRE", "line",
     ["**K14-cre**", "constitutive epithelial Cre", "JAX stock 016230"], LINE_BOX),
    ("LINE_FL", "line",
     ["**_Smad4_^{fl/fl}**", "floxed allele", "JAX stock 017462"], LINE_BOX),
    ("LINE_CREER", "line",
     ["**K14-creER**", "tamoxifen-inducible epithelial Cre", "JAX stock 005107"],
     LINE_BOX),

    # --- constitutive arm ---------------------------------------------------
    ("CROSS_C", "c", ["**K14-cre  ×  _Smad4_^{fl/fl}**"], CROSS),
    ("DCKO_C", "c",
     ["**dcKO**", "K14-cre;", "_Smad4_^{fl/fl}"], f'penwidth=1.8, color="{INK}"'),
    # The sibling pair is the widest rank in each arm, and the panel's type size is set
    # by its width, so these two boxes are broken over more lines than reads naturally.
    ("CON_C", "c",
     ["**CON**", "Cre-negative", "_Smad4_^{fl/fl}", "littermate"], ""),
    # tethered to the dcKO only: the note explains why the MUTANT cannot be studied
    # intact. The control is viable, so it must not be drawn as a parent of this note.
    ("WHY_C", "c",
     ["K14-cre; _Smad4_^{fl/fl} does not survive",
      "to enamel maturation — recorded as",
      "embryonic lethality and as death",
      "shortly after birth — so the phenotype",
      "cannot be seen in the intact mutant"], NOTE_BOX),
    ("PROC_C", "c",
     ["Tooth germs dissected at **E18.5** → renal",
      "capsule of 8–10-week-old wild-type",
      "male hosts → recovered at **9 days**"], ""),
    ("ASSAY_C", "c",
     ["X-ray (Faxitron MX-20), H&E histology,",
      "ameloblastin immunohistochemistry"], ""),
    ("FIG1", "c",
     [f"**MAIN FIGURE 1** {DAGGER}",
      "enamel phenotype of the graft",
      "material: **K14-cre**; no tamoxifen"], TERMINAL),

    # --- inducible arm ------------------------------------------------------
    ("CROSS_I", "i", ["**K14-creER  ×  _Smad4_^{fl/fl}**"], CROSS),
    ("DCKO_I", "i",
     ["**dcKO**", "K14-creER;", "_Smad4_^{fl/fl}"], f'penwidth=1.8, color="{INK}"'),
    ("CON_I", "i",
     ["**CON**", "Cre-negative", "_Smad4_^{fl/fl}", "littermate"], ""),
    ("DISS_I", "i",
     ["Mandibular first-molar germs dissected",
      "at **P0** → Trowell-type organ culture"], ""),
    # The titration is what DEFINED 90 umol/L, so it branches from the culture step in
    # parallel with the tamoxifen node and feeds INTO it. Drawing it downstream of the
    # 90 umol/L node would say the germs were dosed and then dose-ranged.
    ("FIG_S2", "i",
     ["Tamoxifen dose series — ten doses over",
      "50–300 µmol/L (50, 90, 120 shown)",
      "→ renal capsule → contact radiography",
      "→ **SUPPL. FIGURE S2** (K14-creER)"], TERMINAL),
    ("FIG3CE", "i",
     ["Cultured germs, each inhibitor versus",
      "vehicle: dorsomorphin, SB525334,",
      "SB203580, U0126 → RT-qPCR and",
      "immunoblot of enamel epithelium",
      f"→ **MAIN FIGURE 3C–E** (K14-creER {DDAGGER})"], TERMINAL),
    ("TAM", "i",
     ["**4-hydroxytamoxifen, 90 µmol/L, 24 h**",
      "added to the culture medium; this is",
      "the dose set by the titration above.",
      "**dcKO and CON germs are treated alike**,",
      "so tamoxifen is not a variable here"],
     f'fillcolor="#F7E3E3", color="{RED}", penwidth=1.8, peripheries=2, '
     'style="filled"'),
    ("FIG_S1", "i",
     ["Renal-capsule transplantation,",
      "recovered at **14 days** → H&E histology",
      "→ **SUPPL. FIGURE S1** (K14-creER)"], TERMINAL),
    ("BR_C", "i",
     ["Culture continued → enamel-organ",
      "epithelium peeled from mesenchyme"], ""),
    ("FIG23", "i",
     ["RT-qPCR and immunoblot, dcKO vs CON",
      "→ **MAIN FIGURE 2** and **FIGURE 3A, B**",
      "(K14-creER)"], TERMINAL),
    ("FIG456", "i",
     ["Total RNA → bulk RNA sequencing; n = 3",
      "dcKO vs 3 CON, each library a pool of",
      "germs from several pups, both sexes",
      "→ **MAIN FIGURES 4, 5, 6** (K14-creER)",
      "→ **SUPPL. FIGURES S3 and S5**"], TERMINAL),
]

TETHER = 'style=dashed, color="#AAAAAA", arrowsize=0.5'
SELECTS = (f'style=dashed, color="{RED}", arrowsize=0.6, fontsize={FOOT_PT}, '
           f'fontcolor="{RED}", fontname="Helvetica", label=" selects the dose"')

# (tail, head, attributes)
# weight=0 on the parental edges and weight>1 on the backbone of each arm: without
# this the cross node is dragged sideways by its parents and sits ~1 in off the column
# it feeds, which widens the whole panel and so shrinks every label on it.
EDGES = [
    ("LINE_CRE", "CROSS_C", "weight=0"),
    ("LINE_FL", "CROSS_C", "weight=0"),
    ("LINE_FL", "CROSS_I", "weight=0"),
    ("LINE_CREER", "CROSS_I", "weight=0"),

    ("CROSS_C", "DCKO_C", "weight=8"),
    ("CROSS_C", "CON_C", "weight=8"),
    ("DCKO_C", "WHY_C", TETHER),
    ("DCKO_C", "PROC_C", "weight=8"),
    ("CON_C", "PROC_C", "weight=8"),
    ("PROC_C", "ASSAY_C", "weight=12"),
    ("ASSAY_C", "FIG1", "weight=12"),

    ("CROSS_I", "DCKO_I", "weight=8"),
    ("CROSS_I", "CON_I", "weight=8"),
    ("DCKO_I", "DISS_I", "weight=8"),
    ("CON_I", "DISS_I", "weight=8"),
    ("DISS_I", "FIG_S2", "weight=12"),
    ("DISS_I", "FIG3CE", "weight=0"),
    ("DISS_I", "TAM", "weight=0"),
    ("FIG_S2", "TAM", SELECTS),
    ("TAM", "FIG_S1", "weight=12"),
    ("TAM", "BR_C", "weight=0"),
    ("BR_C", "FIG23", "weight=12"),
    ("BR_C", "FIG456", "weight=0"),
]

# Each arm is drawn as ONE column, and each destination carries the figure it feeds
# inside its own box. Setting the destinations side by side makes the panel ~27 in wide
# (3 pt type once scaled to page width); two sub-columns per arm still leaves it at
# 5-6 pt. The invisible edges below stack the destinations instead of ranking them
# abreast. WHY_C is stacked into the constitutive column rather than set beside it: the
# note is 48 characters wide, and placing it alongside the chain widened the whole panel
# by 4.5 in, which costs every other label ~2 pt. The constitutive arm is the short one,
# so it has the height to spare and none of the width.
SAME_RANK = [("DCKO_C", "CON_C"), ("DCKO_I", "CON_I")]
STACK = [("WHY_C", "PROC_C"),
         ("FIG_S2", "FIG3CE"), ("FIG3CE", "TAM"), ("FIG_S1", "BR_C"),
         ("FIG23", "FIG456")]

# The allocation the Methods paragraph states, written as an assertion rather than a
# comment: terminal node -> (lines it MUST descend from, lines it must NOT descend from)
# FIG3CE is included on the strength of [M]; [B] records wild-type germs for those
# panels, which is why that node carries the double dagger. The assertion tracks the
# manuscript, and the panel tells the reader the manuscript is not corroborated here.
ALLOCATION = {
    "FIG1":   (["LINE_CRE", "LINE_FL"],   ["LINE_CREER"]),
    "FIG_S1": (["LINE_CREER", "LINE_FL"], ["LINE_CRE"]),
    "FIG_S2": (["LINE_CREER", "LINE_FL"], ["LINE_CRE"]),
    "FIG3CE": (["LINE_CREER", "LINE_FL"], ["LINE_CRE"]),
    "FIG23":  (["LINE_CREER", "LINE_FL"], ["LINE_CRE"]),
    "FIG456": (["LINE_CREER", "LINE_FL"], ["LINE_CRE"]),
}

# Tamoxifen must not be reachable from anything the inhibitor panels descend from, and
# must be reachable by everything the manuscript describes as tamoxifen-induced.
TAM_DOWNSTREAM = ["FIG_S1", "FIG23", "FIG456"]
TAM_NOT_DOWNSTREAM = ["FIG1", "FIG_S2", "FIG3CE"]

CLUSTERS = [
    ("c", "Arm 1 — constitutive driver K14-cre\\nno tamoxifen at any step of this arm",
     COOL_FILL, COOL_BG),
    ("i", "Arm 2 — inducible driver K14-creER\\ndeletion is tamoxifen-dependent",
     WARM_FILL, WARM_BG),
]

# The footer is drawn as a plain-text NODE stacked under the constitutive arm, not as a
# graph label under the whole panel. The constitutive arm is much the shorter of the two,
# so the lower-left quadrant would otherwise be ~4 x 5 in of white space, and a graph
# label would add its own height on top of that. Putting the footer in the hole costs no
# height at all. The two markers are glossed here as well as in the legend, because a
# panel is often reproduced without its legend.
FOOTER_PARAS = [
    (f"{DAGGER} ", "The Figure 1 artwork is currently labelled K14-Cre-ER;Smad4 f/f, "
     "which disagrees with its own legend, with the Results and with the Methods. The "
     "scheme follows the text; that legend already flags the label for correction to "
     "K14-cre."),
    (f"{DDAGGER} ", "The material for the inhibitor panels is recorded as K14-creER "
     "germs in this manuscript's Methods and as germs from wild-type mice in the "
     "original submission. Drawn as this manuscript states it. Unresolved."),
    ("", "Both controls are Cre-negative littermates that still carry the floxed "
     "allele, which is why control libraries also contain loxP reads (Supplementary "
     "Figure S5). Group sizes are given in the individual figure legends. Housing, "
     "genotyping and the floxed-exon coordinates are given in the legend."),
]


def _footer_lines():
    """Wrap the footer to MAX_FOOT rather than hand-breaking it.

    The footer is kept to the two unresolved-point glosses and the one fact a reader
    needs to interpret the control libraries. Everything else that used to be here --
    housing, genotyping, the floxed-exon coordinates, the E16.5 dissection -- moved
    into the legend: on the panel it cost 1.3 in of height, which the whole panel then
    paid for twice over in scale.
    """
    import textwrap
    out = []
    for marker, para in FOOTER_PARAS:
        if not para:
            out.append("")
            continue
        wrapped = textwrap.wrap(marker + para, width=MAX_FOOT,
                                subsequent_indent="  " if marker else "")
        out.extend(wrapped)
    return out


FOOTER_LINES = _footer_lines()
FOOTER = "".join(x + "\\l" for x in FOOTER_LINES)


# --------------------------------------------------------------------------- audits
def audit() -> None:
    ids = [n[0] for n in NODES]
    assert len(ids) == len(set(ids)), "duplicate node id"
    arm = {n[0]: n[1] for n in NODES}

    for t, h, _ in EDGES:
        assert t in arm and h in arm, f"edge references an undeclared node: {t} -> {h}"
    for a, b in SAME_RANK + STACK:
        assert a in arm and b in arm, f"layout constraint on an undeclared node: {a}, {b}"

    # every node connected
    touched = {x for t, h, _ in EDGES for x in (t, h)}
    orphans = sorted(set(ids) - touched)
    assert not orphans, f"unconnected node(s): {orphans}"

    # line length caps -- the panel has to stay legible when scaled to page width
    for nid, _, lines, _ in NODES:
        for ln in lines:
            n = len(_plain(ln))
            assert n <= MAX_LINE, f"{nid}: label line of {n} chars exceeds {MAX_LINE}: {ln!r}"
    for ln in FOOTER_LINES:
        assert len(ln) <= MAX_FOOT, f"footer line of {len(ln)} chars exceeds {MAX_FOOT}: {ln!r}"

    # reachability: the driver -> figure allocation stated in the Methods
    adj = {}
    for t, h, att in EDGES:
        if "dir=none" in att:            # annotation tethers are not descent
            continue
        adj.setdefault(t, []).append(h)

    def reaches(src, dst):
        seen, q = {src}, deque([src])
        while q:
            cur = q.popleft()
            if cur == dst:
                return True
            for nxt in adj.get(cur, []):
                if nxt not in seen:
                    seen.add(nxt)
                    q.append(nxt)
        return False

    for term, (must, must_not) in ALLOCATION.items():
        for src in must:
            assert reaches(src, term), f"{term} is not reachable from {src}"
        for src in must_not:
            assert not reaches(src, term), (
                f"{term} is reachable from {src} — the diagram contradicts the Methods "
                "allocation of drivers to figures")

    # tamoxifen exposure, as a graph property rather than as prose. The titration and
    # the inhibitor panels must NOT sit downstream of the 90 umol/L node: the first
    # defines that dose and the second is not a tamoxifen experiment at all.
    for term in TAM_DOWNSTREAM:
        assert reaches("TAM", term), f"{term} should be downstream of tamoxifen"
    for term in TAM_NOT_DOWNSTREAM:
        assert not reaches("TAM", term), (
            f"{term} is drawn downstream of the 90 µmol/L tamoxifen node, which asserts "
            "an exposure the sources do not describe")

    # tamoxifen belongs to the inducible arm only (or to an explicit negation).
    # The constitutive arm carrying a tamoxifen step would be the same class of error
    # this figure exists to prevent, so it is asserted rather than merely proof-read.
    for nid, a, lines, _ in NODES:
        txt = " ".join(_plain(x) for x in lines)
        low = txt.lower()
        if "tamoxifen" in low:
            assert a == "i" or nid == "LINE_CREER" or "no tamoxifen" in low, (
                f"{nid} names tamoxifen but is neither in the inducible arm nor the "
                "K14-creER line box, and does not negate it")
    for key, title, _, _ in CLUSTERS:
        low = title.lower()
        if "tamoxifen" in low:
            assert key == "i" or "no tamoxifen" in low, (
                f"cluster {key} names tamoxifen without negating it")
    assert "no tamoxifen" in dict((k, t) for k, t, _, _ in CLUSTERS)["c"].lower(), (
        "the constitutive arm must state explicitly that no tamoxifen is used")

    # figure numbers must not stray across arms
    for nid, a, lines, _ in NODES:
        txt = " ".join(_plain(x) for x in lines)
        if "MAIN FIGURE 1" in txt:
            assert a == "c", f"{nid} names Figure 1 outside the constitutive arm"
        if re.search(r"FIGURES? S[123]|FIGURE 2|FIGURE 3|FIGURES 4", txt):
            assert a == "i", f"{nid} names an inducible-arm figure outside that arm"

    # every terminal figure box must name the material it came from, so that a box
    # reproduced on its own still identifies its mouse line
    for term in ALLOCATION:
        txt = " ".join(_plain(x) for x in
                       next(n[2] for n in NODES if n[0] == term))
        assert "K14-cre" in txt, f"{term} does not name its driver inside the box"

    # every marker drawn on the panel must be glossed in the footer, and every gloss
    # must belong to a marker that is actually drawn
    drawn = {m for m in (DAGGER, DDAGGER)
             for nid, _, lines, _ in NODES if m in " ".join(lines)}
    glossed = {m for m in (DAGGER, DDAGGER) if m in FOOTER}
    assert drawn == glossed == {DAGGER, DDAGGER}, (
        f"marker mismatch: drawn on the panel {sorted(drawn)}, glossed in the footer "
        f"{sorted(glossed)}")


# ----------------------------------------------------------------------------- dot
def build_dot() -> str:
    arm_of = {n[0]: n[1] for n in NODES}

    def emit(nid, lines, extra):
        att = f"label={label(lines)}"
        if extra:
            att += ", " + extra
        return f'  {nid} [{att}];\n'

    out = [
        'digraph breeding {\n'
        '  graph [rankdir=TB, splines=polyline, nodesep=0.16, ranksep=0.18,\n'
        '         fontname="Helvetica", bgcolor="white", labelloc=b, labeljust=l,\n'
        f'         fontsize={FOOT_PT}, fontcolor="#3A3A3A", label="{FOOTER}"];\n'
        f'  node  [shape=box, style="rounded,filled", fontname="Helvetica", '
        f'fontsize={BODY_PT},\n'
        '         margin="0.13,0.070", color="#666666", penwidth=1.0,\n'
        f'         fillcolor="{COOL_FILL}"];\n'
        f'  edge  [color="#555555", arrowsize=0.6, penwidth=1.0];\n'
    ]
    for nid, a, lines, extra in NODES:
        if a == "line":
            out.append(emit(nid, lines, extra))
    out.append("  {rank=same; LINE_CRE; LINE_FL; LINE_CREER;}\n")

    for key, title, fill, bg in CLUSTERS:
        out.append(f'  subgraph cluster_{key} {{\n'
                   f'    label="{title}"; labeljust=l; labelloc=t; '
                   f'fontsize={CLUSTER_PT};\n'
                   f'    style="rounded,filled"; color="#999999"; penwidth=1.4;\n'
                   f'    fillcolor="{bg}"; fontcolor="{INK}"; margin=8;\n'
                   f'    node [fillcolor="{fill}"];\n')
        for nid, a, lines, extra in NODES:
            if a == key:
                out.append("  " + emit(nid, lines, extra))
        # rank sets must be declared INSIDE the cluster; a rank set at graph level
        # silently evicts its members from the cluster (graphviz warns and the cluster
        # box then omits them)
        for x, y in SAME_RANK:
            if arm_of[x] == key and arm_of[y] == key:
                out.append(f"    {{rank=same; {x}; {y};}}\n")
        out.append("  }\n")

    for t, h, att in EDGES:
        out.append(f'  {t} -> {h}{(" [" + att + "]") if att else ""};\n')
    for a, b in STACK:
        out.append(f"  {a} -> {b} [style=invis, weight=6];\n")
    out.append("}\n")
    return "".join(out)


def render(name: str, dot: str) -> None:
    dp = os.path.join(FIG, f"{name}.dot")
    open(dp, "w").write(dot)
    for fmt, extra in (("png", ["-Gdpi=300"]), ("pdf", [])):
        subprocess.run(["dot", f"-T{fmt}", *extra, dp, "-o",
                        os.path.join(FIG, f"{name}.{fmt}")], check=True)
    from PIL import Image
    w, h = Image.open(os.path.join(FIG, f"{name}.png")).size
    win, hin = w / 300, h / 300
    scale = EMBED_IN / win
    page_h = hin * scale

    # Check EVERY font size the .dot actually emits, not just the largest one. The
    # previous version tested only the node body size and reported "8 pt" while the
    # footer shipped at 6.8 pt and the fl/fl superscripts at 5.5 pt.
    sizes = sorted({float(x) for x in
                    re.findall(r'fontsize\s*=\s*"?([0-9.]+)', dot)} |
                   {float(x) for x in
                    re.findall(r'POINT-SIZE\s*=\s*"([0-9.]+)"', dot)})
    scaled = {s: s * scale for s in sizes}
    worst = min(scaled.values())
    print(f"  {name}: {win:.1f} x {hin:.1f} in  ->  at {EMBED_IN} in wide, "
          f"height {page_h:.1f} in")
    print("    type on the page: " +
          ", ".join(f"{s:g} pt → {v:.1f} pt" for s, v in scaled.items()))
    assert worst >= MIN_PT, (
        f"the smallest type in this panel would reach {worst:.1f} pt on the page "
        f"(< {MIN_PT} pt). Shorten labels, raise that font size, or split the panel — "
        "do not ship type this small.")
    assert page_h <= MAX_PAGE_IN, (
        f"the scaled panel is {page_h:.1f} in tall, taller than a portrait text block "
        f"({MAX_PAGE_IN} in)")


if __name__ == "__main__":
    audit()
    render(NAME, build_dot())
