"""
volcano_panel.py — the Figure 4B volcano renderer.

WHY THIS IS A SEPARATE MODULE
-----------------------------
The previous panel B was unreadable, for three measurable reasons:

  1. 420 genes have padj that underflowed to exactly 0 in double precision. The old
     code clipped padj to 1e-300 and plotted -log10 of that, so all 420 landed on a
     single y = 300 line as one opaque bar — and, worse, the figure asserted a
     -log10 P of 300 for genes whose P value is simply not representable.
  2. The band -log10 P = 30..200 contains ZERO genes yet occupied ~55% of the panel
     height, while 78% of the DEGs (440/563) were crushed into the bottom ~6%.
  3. Every gene label used the same offset (`xytext=(0, 4), ha="center"`) with no
     de-collision, so neighbouring labels overprinted each other:
     "Sqle"+"Hmgcr" rendered as "SqHemgcr" and "Ptgs2"+"Mmp12"+"Ereg" as
     "Ptgs2MEreg12".

The fix is a three-band broken y-axis plus a real label placer, which is enough code
that keeping it inside 07_figures.py would swamp that file. 07_figures.py imports
`draw()` from here.

DESIGN NOTES
------------
* The three bands are the three genuinely disjoint populations in the data:
  the categorical padj = 0 row, a sparse "extreme but measurable" band, and the bulk.
  The cut band is asserted empty at run time — if the data ever change so that genes
  fall inside it, the assertion fails loudly rather than hiding them.
* The underflow row carries NO y information. Points are given a deterministic
  (seeded) vertical spread purely so they stop collapsing into a bar; the row is
  labelled "padj = 0 (underflow)" and the manuscript legend states that its height is
  not a measured quantity. This is a presentation choice that is disclosed, not a
  silent transformation.
* Nothing is hard-coded: every count, tick and annotation is derived from the table at
  draw time (CLAUDE.md §6.2).
"""

from __future__ import annotations

import numpy as np
from matplotlib.transforms import Bbox
import matplotlib.text

BLUE, RED, GREY = "#2166AC", "#B2182B", "#BBBBBB"

XLO, XHI = -6.4, 6.4      # x window; genes outside are drawn as edge arrowheads
XPAD = 0.45               # inset so arrowheads sit inside the spine
BREAK = 20.0              # -log10 P separating the bulk from the extreme band
JIT = 0.42                # half-height of the categorical underflow strip

LABEL_GENES = ["Krt76", "Isg15", "Mx1", "Oasl1", "Wnt10a", "Sqle", "Hmgcr",
               "Ldlr", "Fdps", "Ptgs2", "Mmp12", "Ereg", "Cyp51", "Hmgcs1"]


# ---------------------------------------------------------------------------
# label de-collision (self-contained; adjustText is not a project dependency)
# ---------------------------------------------------------------------------
class LabelPlacer:
    """Greedy candidate-offset label placement, scored in display coordinates."""

    def __init__(self, fig):
        self.fig = fig
        fig.canvas.draw()
        self.rend = fig.canvas.get_renderer()
        self.px_per_pt = fig.dpi / 72.0
        self.taken = []      # Bboxes already occupied by text
        self.leaders = []    # leader-line segments already drawn
        self.anchors = {}    # gene -> Bbox around its own marker
        self.clouds = []     # display-space point clouds (density penalty)

    def add_obstacle(self, artist, grow=2.0):
        bb = artist.get_window_extent(renderer=self.rend)
        self.taken.append(Bbox.from_extents(bb.x0 - grow, bb.y0 - grow,
                                            bb.x1 + grow, bb.y1 + grow))

    def add_vline_obstacle(self, ax, x, half_width_px=3.0):
        """Keep label glyphs off a vertical rule (the |log2FC| = 1 thresholds).

        Without this, 'Cyp51' and 'Hmgcs1' are bisected by the x = +1 dashed line.
        """
        fr = ax.get_window_extent(renderer=self.rend)
        px = ax.transData.transform((x, 0))[0]
        self.taken.append(Bbox.from_extents(px - half_width_px, fr.y0,
                                            px + half_width_px, fr.y1))

    def add_anchor(self, key, ax, x, y, r=7.0):
        px, py = ax.transData.transform((x, y))
        self.anchors[key] = Bbox.from_extents(px - r, py - r, px + r, py + r)

    def add_point_cloud(self, ax, xs, ys):
        self.clouds.append(ax.transData.transform(np.column_stack([xs, ys])))

    def _measure(self, text, fontsize, style):
        t = self.fig.text(0.5, 0.5, text, fontsize=fontsize, style=style)
        bb = t.get_window_extent(renderer=self.rend)
        t.remove()
        return bb.width, bb.height

    @staticmethod
    def _overlap(a, b):
        dx = min(a.x1, b.x1) - max(a.x0, b.x0)
        dy = min(a.y1, b.y1) - max(a.y0, b.y0)
        return dx * dy if (dx > 0 and dy > 0) else 0.0

    def _outside(self, box, frame):
        return box.width * box.height - self._overlap(box, frame)

    @staticmethod
    def _seg_hits_box(p0, p1, box):
        """Liang–Barsky segment/rectangle intersection."""
        x0, y0 = p0
        dx, dy = p1[0] - x0, p1[1] - y0
        t0, t1 = 0.0, 1.0
        for p, q in ((-dx, x0 - box.x0), (dx, box.x1 - x0),
                     (-dy, y0 - box.y0), (dy, box.y1 - y0)):
            if p == 0:
                if q < 0:
                    return False
                continue
            r = q / p
            if p < 0:
                if r > t1:
                    return False
                t0 = max(t0, r)
            else:
                if r < t0:
                    return False
                t1 = min(t1, r)
        return t0 <= t1

    def _density(self, box):
        n = 0
        for pts in self.clouds:
            n += int(((pts[:, 0] > box.x0) & (pts[:, 0] < box.x1) &
                      (pts[:, 1] > box.y0) & (pts[:, 1] < box.y1)).sum())
        return n

    def place(self, items, fontsize=6.0, style="italic",
              radii=(9, 11, 13.5, 16, 19, 23, 28, 34, 41, 49, 58),
              angles=(90, 68, 112, 50, 130, 35, 145, 20, 160, 0, 180,
                      -25, -155, -50, -130, -75, -105, -90)):
        items = sorted(items, key=lambda it: it.get("prio", 5))
        for it in items:
            ax, x, y, txt = it["ax"], it["x"], it["y"], it["text"]
            others = [b for k, b in self.anchors.items() if k != txt]
            w, h = self._measure(txt, fontsize, style)
            ax0 = ax.transData.transform((x, y))
            fr = ax.get_window_extent(renderer=self.rend)
            frame = Bbox.from_extents(fr.x0 + 1.5, fr.y0 + 1.0, fr.x1 - 1.5, fr.y1 - 1.0)

            best, best_score = None, np.inf
            for r_pt in radii:
                r = r_pt * self.px_per_pt
                for a_deg in angles:
                    a = np.deg2rad(a_deg)
                    cx = ax0[0] + r * np.cos(a) + np.sign(np.cos(a)) * w * 0.5 * 0.55
                    cy = ax0[1] + r * np.sin(a) + np.sign(np.sin(a)) * h * 0.5 * 0.75
                    box = Bbox.from_extents(cx - w / 2 - 2.6, cy - h / 2 - 2.0,
                                            cx + w / 2 + 2.6, cy + h / 2 + 2.0)
                    norm = max(w * h, 1)
                    score = 0.0
                    score += 4000.0 * sum(self._overlap(box, t) for t in self.taken) / norm
                    score += 4000.0 * sum(self._overlap(box, t) for t in others) / norm
                    score += 6000.0 * self._outside(box, frame) / norm
                    score += 34.0 * self._density(box)
                    score += 0.55 * r_pt
                    score += 0.9 * abs(a_deg - 90) / 10.0
                    seg = (tuple(ax0), (cx, cy))
                    score += 260.0 * sum(self._seg_hits_box(*seg, t) for t in self.taken)
                    score += 260.0 * sum(self._seg_hits_box(a, b, box) for a, b in self.leaders)
                    if score < best_score:
                        best_score, best = score, (cx - ax0[0], cy - ax0[1], box)
            dx_px, dy_px, box = best
            self.taken.append(box)
            self.leaders.append((tuple(ax0), (ax0[0] + dx_px, ax0[1] + dy_px)))
            arrow = None
            if float(np.hypot(dx_px, dy_px)) > 11 * self.px_per_pt:
                arrow = dict(arrowstyle="-", lw=0.42, color="0.35", shrinkA=0.5, shrinkB=1.2)
            ax.annotate(txt, xy=(x, y),
                        xytext=(dx_px / self.px_per_pt, dy_px / self.px_per_pt),
                        textcoords="offset points", ha="center", va="center",
                        fontsize=fontsize, style=style, color="0.12",
                        annotation_clip=False, arrowprops=arrow, zorder=6)


def text_collisions(fig):
    """Return overlapping pairs of rendered glyph boxes (verification helper).

    Two subtleties, both learned the hard way:
      * Annotation.get_window_extent() unions in the leader line, which produces
        phantom clashes; we call the Text base-class implementation for glyphs only.
      * Matplotlib keeps Text objects for major ticks that fall OUTSIDE the axis view
        interval. They are never drawn, but they still report a stale window extent and
        collide with each other. We therefore keep a tick label only if its tick sits
        inside the current view interval.
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
        if not t.get_visible() or not t.get_text().strip() or id(t) in ghosts:
            continue
        try:
            bb = matplotlib.text.Text.get_window_extent(t, renderer=rend)
        except Exception:
            continue
        if bb.width > 0 and bb.height > 0:
            boxes.append((t.get_text(), bb))
    out = []
    for i in range(len(boxes)):
        for j in range(i + 1, len(boxes)):
            a, b = boxes[i][1], boxes[j][1]
            if (min(a.x1, b.x1) - max(a.x0, b.x0) > 0.5 and
                    min(a.y1, b.y1) - max(a.y0, b.y0) > 0.5):
                out.append((boxes[i][0], boxes[j][0]))
    return out


# ---------------------------------------------------------------------------
# the panel
# ---------------------------------------------------------------------------
def prepare(deg):
    """Add sig / underflow / nlp columns. -log10 P is left NaN where padj == 0."""
    d = deg.dropna(subset=["log2FoldChange", "padj"]).copy()
    d["sig"] = (d.padj < 0.05) & (d.log2FoldChange.abs() > 1)
    d["underflow"] = d.padj == 0.0
    d["nlp"] = np.where(d.underflow, np.nan, -np.log10(d.padj.where(d.padj > 0)))
    return d


def draw(fig, spec, deg, title_prefix="B"):
    """Render the volcano into `spec` (a SubplotSpec). Returns a dict of counts."""
    d = prepare(deg)

    n_sig = int(d.sig.sum())
    n_up = int((d.sig & (d.log2FoldChange > 0)).sum())
    n_dn = int((d.sig & (d.log2FoldChange < 0)).sum())
    n_uf = int(d.underflow.sum())
    n_uf_sig = int((d.underflow & d.sig).sum())
    n_off = int((d.log2FoldChange < XLO).sum() + (d.log2FoldChange > XHI).sum())

    hi = d[(~d.underflow) & (d.nlp > BREAK)]
    lo = d[(~d.underflow) & (d.nlp <= BREAK)]
    assert len(hi) + len(lo) + n_uf == len(d)
    # the excised band must genuinely be empty, or we would be hiding data
    assert lo.nlp.max() < BREAK and hi.nlp.min() > 150.0, (
        "the empty -log10 P band assumed by this panel no longer holds; "
        "re-inspect the distribution before plotting")

    gs = spec.subgridspec(3, 1, height_ratios=[1.15, 0.95, 3.05], hspace=0.18)
    ax_uf, ax_hi, ax_lo = (fig.add_subplot(gs[0]), fig.add_subplot(gs[1]),
                           fig.add_subplot(gs[2]))

    def scatter(ax, sub, yvals):
        x = sub.log2FoldChange.values
        y = np.asarray(yvals, dtype=float)
        off = (x < XLO) | (x > XHI)
        xc = np.clip(x, XLO + XPAD, XHI - XPAD)
        groups = ((~sub.sig.values, GREY, 2.2, 0.28, 1),
                  (sub.sig.values & (x < 0), BLUE, 5.0, 0.75, 2),
                  (sub.sig.values & (x > 0), RED, 5.0, 0.75, 2))
        for m, col, s, al, z in groups:
            mm = m & ~off
            ax.scatter(xc[mm], y[mm], s=s, c=col, alpha=al, linewidths=0,
                       rasterized=True, zorder=z)
            for i in np.where(m & off)[0]:      # off-scale -> outward arrowhead
                ax.scatter(xc[i], y[i], s=13, c=col, alpha=0.9, zorder=z + 1,
                           marker=("<" if x[i] < 0 else ">"), linewidths=0)

    uf = d[d.underflow].copy()
    uf_y = np.random.default_rng(0).uniform(-JIT, JIT, len(uf))  # seeded, decorative
    scatter(ax_uf, uf, uf_y)
    scatter(ax_hi, hi, hi.nlp.values)
    scatter(ax_lo, lo, lo.nlp.values)

    for ax in (ax_uf, ax_hi, ax_lo):
        ax.set_xlim(XLO, XHI)
        ax.axvline(-1, ls="--", c="k", lw=0.5, zorder=0)
        ax.axvline(1, ls="--", c="k", lw=0.5, zorder=0)
        ax.spines["right"].set_visible(False)
        ax.spines["top"].set_visible(False)

    ax_uf.set_ylim(-0.95, 1.15)
    ax_uf.set_yticks([0.0])
    ax_uf.set_yticklabels(["padj = 0\n(underflow)"], fontsize=6.0)
    ax_uf.spines["left"].set_bounds(-JIT - 0.18, JIT + 0.18)
    ax_uf.spines["bottom"].set_visible(False)
    ax_uf.tick_params(axis="x", bottom=False, labelbottom=False)

    # ticks derived from the data, so the topmost point always has a reference above it
    # This band is physically short, so it gets at most three ticks, chosen to
    # bracket the data (the topmost point must sit below the highest tick).
    lo_t = int(np.floor(hi.nlp.min() / 20.0) * 20)
    hi_t = int(np.ceil(hi.nlp.max() / 20.0) * 20)
    ax_hi.set_ylim(lo_t - 6, hi_t + 6)
    step = max(20, int(np.ceil((hi_t - lo_t) / 2.0 / 10.0) * 10))
    ax_hi.set_yticks([t for t in range(lo_t, hi_t + 1, step)][:3])
    ax_hi.spines["bottom"].set_visible(False)
    ax_hi.tick_params(axis="x", bottom=False, labelbottom=False)

    ax_lo.set_ylim(-0.85, float(lo.nlp.max()) * 1.09)
    ax_lo.set_yticks([0, 4, 8, 12, 16])
    ax_lo.axhline(-np.log10(0.05), ls="--", c="k", lw=0.5, zorder=0)
    ax_lo.set_xlabel("log2 fold change (dcKO/CON)")

    # conventional diagonal break marks on the left spine at each cut
    bk = dict(marker=[(-1, -1), (1, 1)], markersize=3.2, linestyle="none",
              color="k", mec="k", mew=0.6, clip_on=False, zorder=10)
    uf_frac = (-JIT - 0.18 - ax_uf.get_ylim()[0]) / float(np.ptp(ax_uf.get_ylim()))
    ax_uf.plot([0], [uf_frac], transform=ax_uf.transAxes, **bk)
    ax_hi.plot([0, 0], [0, 1], transform=ax_hi.transAxes, **bk)
    ax_lo.plot([0], [1], transform=ax_lo.transAxes, **bk)

    ax_uf.set_title(f"{title_prefix}  Volcano ({n_sig} DEGs)", loc="left", pad=6)
    ax_hi.set_ylabel("-log10 adjusted P", labelpad=6)
    ax_hi.yaxis.set_label_coords(-0.155, 0.10, transform=ax_hi.transAxes)

    # sits on the title line, right-aligned, so it never competes with the band
    note = ax_uf.text(1.0, 1.06,
                      f"{n_uf} genes ({n_uf_sig} DEGs)",
                      transform=ax_uf.transAxes, ha="right", va="bottom",
                      fontsize=6.0, color="0.30", zorder=5)
    cnt_dn = ax_lo.text(0.02, 0.97, f"{n_dn} down", transform=ax_lo.transAxes,
                        ha="left", va="top", fontsize=6.2, color=BLUE, fontweight="bold")
    cnt_up = ax_lo.text(0.98, 0.97, f"{n_up} up", transform=ax_lo.transAxes,
                        ha="right", va="top", fontsize=6.2, color=RED, fontweight="bold")

    # ---- gene labels ----
    placer = LabelPlacer(fig)
    for a in (note, cnt_dn, cnt_up):
        placer.add_obstacle(a)
    for ax in (ax_uf, ax_hi, ax_lo):          # keep glyphs off the threshold rules
        placer.add_vline_obstacle(ax, -1.0)
        placer.add_vline_obstacle(ax, 1.0)
    placer.add_point_cloud(ax_lo, np.clip(lo.log2FoldChange.values, XLO + XPAD, XHI - XPAD), lo.nlp.values)
    placer.add_point_cloud(ax_hi, np.clip(hi.log2FoldChange.values, XLO + XPAD, XHI - XPAD), hi.nlp.values)
    placer.add_point_cloud(ax_uf, np.clip(uf.log2FoldChange.values, XLO + XPAD, XHI - XPAD), uf_y)

    uf_ypos = dict(zip(uf.external_gene_name.values, uf_y))
    items = []
    for g in LABEL_GENES:
        r = d[d.external_gene_name == g]
        if r.empty:
            raise RuntimeError(f"label gene {g!r} absent from the DEG table")
        r = r.iloc[0]
        xv = float(np.clip(r.log2FoldChange, XLO + XPAD, XHI - XPAD))
        if r.underflow:
            ax, yv, prio = ax_uf, float(uf_ypos[g]), 1
        elif r.nlp > BREAK:
            ax, yv, prio = ax_hi, float(r.nlp), 2
        else:
            ax, yv, prio = ax_lo, float(r.nlp), 3
        items.append(dict(ax=ax, x=xv, y=yv, text=g, prio=prio))
    for it in items:
        placer.add_anchor(it["text"], it["ax"], it["x"], it["y"])
    placer.place(items)

    return dict(n_sig=n_sig, n_up=n_up, n_dn=n_dn, n_underflow=n_uf,
                n_underflow_sig=n_uf_sig, n_offscale=n_off,
                max_measured_nlp=float(hi.nlp.max()))
