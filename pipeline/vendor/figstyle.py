"""Publication figure styling + plot archetype helpers.

Part of the `scientific-figure-master` skill — a venue-agnostic engine for
building publication figures from raw data + a human's intent (an instruction
or a manuscript caption/legend). The DEFAULT look is a neutral publication
style that essentially every journal accepts: sans-serif type at final print
size, honest axes, top/right spines off, colorblind-safe series, declared
error bars. Named presets (`prism`, `nature`) are opt-in deltas, chosen from
the target venue or the caption — never assumed. Validated by reproducing
figures across several papers, but none of that is baked into the defaults.
"""
from __future__ import annotations

import numpy as np
import matplotlib as mpl
import matplotlib.pyplot as plt
from matplotlib.patches import Polygon

MM = 1 / 25.4  # mm → inch

OKABE_ITO = ['#E69F00', '#56B4E9', '#009E73', '#F0E442',
             '#0072B2', '#D55E00', '#CC79A7', '#000000']

JOURNAL_MM = {  # (single column, double column) widths in mm
    'nature': (89, 183), 'science': (55, 175), 'cell': (85, 178),
    'peerj': (85, 178), 'plos': (83, 173), 'biofactors': (85, 178),
    'ieee': (88, 181), 'elsevier': (90, 190), 'default': (85, 178),
}

# Style presets: DEFAULTS for the archetype helpers, not hard rules. A preset
# is a small set of kwargs the recipes spread into e.g. barsdots(**bar_kwargs).
# 'neutral' is the venue-agnostic default; 'prism' emulates the GraphPad look
# many bio labs submit; 'nature' is a thinner, lowercase-panel variant.
STYLE_PRESETS = {
    'neutral': {
        'panel_case': 'upper', 'axes_linewidth': 0.8,
        'bar': dict(colors=None, edgecolor='black', errorbar='sd',
                    dot_size=9, open_dots=False),
        'dot': dict(open_markers=False),
    },
    'prism': {
        'panel_case': 'upper', 'axes_linewidth': 0.9,
        'bar': dict(colors=['white'], edgecolor='black', errorbar='sd',
                    dot_size=9, open_dots=True),
        'dot': dict(open_markers=True),
    },
    'nature': {
        'panel_case': 'lower', 'axes_linewidth': 0.6,
        'bar': dict(colors=None, edgecolor='none', errorbar='sem',
                    dot_size=7, open_dots=False),
        'dot': dict(open_markers=False),
    },
}
_ACTIVE_STYLE = 'neutral'


def setup(base_pt: float = 8, palette=None, style: str = 'neutral') -> dict:
    """Apply the publication rcParams for `style` (default 'neutral').
    Returns the preset dict so callers can spread its per-archetype defaults,
    e.g. `st = fs.setup(style='prism'); fs.barsdots(ax, g, **st['bar'])`."""
    global _ACTIVE_STYLE
    _ACTIVE_STYLE = style if style in STYLE_PRESETS else 'neutral'
    preset = STYLE_PRESETS[_ACTIVE_STYLE]
    mpl.rcParams.update({
        'font.family': 'sans-serif',
        'font.sans-serif': ['Arial', 'Helvetica', 'Liberation Sans',
                            'Nimbus Sans', 'DejaVu Sans'],
        'font.size': base_pt,
        'axes.labelsize': base_pt,
        'axes.titlesize': base_pt + 1,
        'xtick.labelsize': base_pt - 1,
        'ytick.labelsize': base_pt - 1,
        'legend.fontsize': base_pt - 1,
        'axes.linewidth': preset['axes_linewidth'],
        'xtick.major.width': preset['axes_linewidth'],
        'ytick.major.width': preset['axes_linewidth'],
        'xtick.major.size': 2.5, 'ytick.major.size': 2.5,
        'xtick.direction': 'out', 'ytick.direction': 'out',
        'axes.spines.top': False, 'axes.spines.right': False,
        'axes.unicode_minus': False,
        'legend.frameon': False,
        'savefig.dpi': 300,
        'pdf.fonttype': 42, 'ps.fonttype': 42,  # editable text in vectors
    })
    mpl.rcParams['axes.prop_cycle'] = plt.cycler(color=palette or OKABE_ITO)
    return preset


def style_kwargs(kind: str = 'bar', style: str | None = None) -> dict:
    """Default kwargs for an archetype under the active (or named) preset."""
    return dict(STYLE_PRESETS.get(style or _ACTIVE_STYLE, STYLE_PRESETS['neutral']).get(kind, {}))


def fig_mm(w_mm: float, h_mm: float, **kw):
    """Figure sized in mm of final print size."""
    return plt.subplots(figsize=(w_mm * MM, h_mm * MM), **kw)


def stars(p: float, ns: str = '') -> str:
    """p/FDR value → Prism-style significance string."""
    if p < 1e-4:
        return '****'
    if p < 1e-3:
        return '***'
    if p < 1e-2:
        return '**'
    if p < 0.05:
        return '*'
    return ns


def barsdots(ax, values_by_group, order=None, colors=None, bar_width=0.62,
             errorbar='sd', dot_size=9, dot_color='black', dot_colors=None,
             open_dots=False, jitter=0.10, edgecolor='black', center='mean', seed=0):
    """Prism-style bar (mean or median) + SD whiskers + jittered replicate dots.

    values_by_group: dict {label: 1D array}. Returns {label: x_position}.
    dot_colors: per-group list/dict overriding dot_color. open_dots: white-face
    black-edge circles (the common Prism replicate-dot style)."""
    rng = np.random.default_rng(seed)
    order = list(order or values_by_group)
    xs = {}
    for i, g in enumerate(order):
        v = np.asarray(values_by_group[g], dtype=float)
        v = v[~np.isnan(v)]
        c = 'white' if colors is None else (colors[i % len(colors)] if not isinstance(colors, dict) else colors[g])
        h = np.mean(v) if center == 'mean' else np.median(v)
        ax.bar(i, h, bar_width, facecolor=c, edgecolor=edgecolor, linewidth=0.9, zorder=2)
        if errorbar and len(v) > 1:
            err = np.std(v, ddof=1) if errorbar == 'sd' else np.std(v, ddof=1) / np.sqrt(len(v))
            ax.errorbar(i, h, yerr=err, fmt='none', ecolor='black',
                        elinewidth=0.9, capsize=2.5, capthick=0.9, zorder=3)
        if dot_size:
            dc = (dot_color if dot_colors is None else
                  (dot_colors[g] if isinstance(dot_colors, dict) else dot_colors[i % len(dot_colors)]))
            dkw = (dict(facecolor='none', edgecolor=dc, linewidth=0.7) if open_dots
                   else dict(facecolor=dc, edgecolor='none'))
            ax.scatter(i + rng.uniform(-jitter, jitter, len(v)), v, s=dot_size,
                       zorder=4, **dkw)
        xs[g] = i
    ax.set_xticks(range(len(order)))
    ax.set_xticklabels(order)
    ax.set_ylim(bottom=0)
    return xs


def dotstrip(ax, values_by_group, order=None, center='mean', dot_size=12,
             dot_colors=None, open_markers=False, jitter=0.10, bar_halfwidth=0.22,
             errorbar=None, seed=0):
    """Prism 'scatter with mean/median crossbar' (no bar rectangle) — the
    biof.2131 Fig 2G/3B/4B dialect. values_by_group: {label: 1D array}.
    dot_colors: per-group list/dict or None (black). errorbar: None|'sd'|'sem'.
    Returns {label: x_position}."""
    rng = np.random.default_rng(seed)
    order = list(order or values_by_group)
    xs = {}
    for i, g in enumerate(order):
        v = np.asarray(values_by_group[g], dtype=float)
        v = v[~np.isnan(v)]
        c = ('black' if dot_colors is None else
             (dot_colors[g] if isinstance(dot_colors, dict) else dot_colors[i % len(dot_colors)]))
        kw = (dict(facecolor='none', edgecolor=c, linewidth=0.8) if open_markers
              else dict(facecolor=c, edgecolor='none'))
        ax.scatter(i + rng.uniform(-jitter, jitter, len(v)), v, s=dot_size,
                   zorder=3, **kw)
        h = np.mean(v) if center == 'mean' else np.median(v)
        ax.hlines(h, i - bar_halfwidth, i + bar_halfwidth, color='black',
                  lw=1.4, zorder=4)
        if errorbar and len(v) > 1:
            err = np.std(v, ddof=1)
            if errorbar == 'sem':
                err /= np.sqrt(len(v))
            ax.errorbar(i, h, yerr=err, fmt='none', ecolor='black',
                        elinewidth=0.9, capsize=2.5, capthick=0.9, zorder=4)
        xs[g] = i
    ax.set_xticks(range(len(order)))
    ax.set_xticklabels(order)
    return xs


def grouped_barsdots(ax, table, group_order, cond_order, colors=None,
                     bar_width=0.36, errorbar='sd', dot_size=8, dot_color='black',
                     open_dots=False, seed=0):
    """Two-factor bars: table[(group, cond)] = 1D array. Returns {(g,c): x}.
    open_dots: white-face black-edge replicate circles."""
    rng = np.random.default_rng(seed)
    colors = colors or ['white', '#9e9e9e']
    xs = {}
    for gi, g in enumerate(group_order):
        for ci, c in enumerate(cond_order):
            v = np.asarray(table[(g, c)], dtype=float)
            v = v[~np.isnan(v)]
            x = gi + (ci - (len(cond_order) - 1) / 2) * bar_width
            ax.bar(x, v.mean(), bar_width * 0.92, facecolor=colors[ci % len(colors)],
                   edgecolor='black', linewidth=0.9, zorder=2,
                   label=c if gi == 0 else None)
            if len(v) > 1:
                err = np.std(v, ddof=1) if errorbar == 'sd' else np.std(v, ddof=1) / np.sqrt(len(v))
                ax.errorbar(x, v.mean(), yerr=err, fmt='none', ecolor='black',
                            elinewidth=0.9, capsize=2, capthick=0.9, zorder=3)
            if dot_size:
                dkw = (dict(facecolor='none', edgecolor=dot_color, linewidth=0.6)
                       if open_dots else dict(facecolor=dot_color, edgecolor='none'))
                ax.scatter(x + rng.uniform(-bar_width / 4, bar_width / 4, len(v)), v,
                           s=dot_size, zorder=4, **dkw)
            xs[(g, c)] = x
    ax.set_xticks(range(len(group_order)))
    ax.set_xticklabels(group_order)
    ax.set_ylim(bottom=0)
    return xs


def sig_bracket(ax, x1, x2, y, text, tick_frac=0.02, lw=0.8, fontsize=7, color='black'):
    """Significance bracket between x1 and x2 at height y (data coords)."""
    yspan = np.diff(ax.get_ylim())[0]
    t = tick_frac * yspan
    ax.plot([x1, x1, x2, x2], [y - t, y, y, y - t], lw=lw, c=color,
            solid_capstyle='butt', clip_on=False)
    ax.text((x1 + x2) / 2, y + t * 0.3, text, ha='center', va='bottom',
            fontsize=fontsize, color=color)


def violin(ax, values_by_group, order=None, dots=True, dot_size=3,
           median_color='red', body_color='0.87', seed=0):
    """Prism-style violin: outline body, solid red median, dashed red quartiles,
    per-observation dots inside. (PeerJ 10244 Fig 4 dialect.)"""
    rng = np.random.default_rng(seed)
    order = list(order or values_by_group)
    data = [np.asarray(values_by_group[g], dtype=float) for g in order]
    data = [d[~np.isnan(d)] for d in data]
    parts = ax.violinplot(data, positions=range(len(order)), widths=0.8,
                          showextrema=False)
    for b in parts['bodies']:
        b.set_facecolor(body_color)
        b.set_edgecolor('black')
        b.set_linewidth(0.8)
        b.set_alpha(1.0)
    for i, d in enumerate(data):
        q1, med, q3 = np.percentile(d, [25, 50, 75])
        ax.hlines(med, i - 0.22, i + 0.22, color=median_color, lw=1.2, zorder=4)
        ax.hlines([q1, q3], i - 0.15, i + 0.15, color=median_color, lw=0.8,
                  linestyle=(0, (3, 2)), zorder=4)
        if dots:
            ax.scatter(i + rng.uniform(-0.16, 0.16, len(d)), d, s=dot_size,
                       facecolor='black', edgecolor='none', alpha=0.7, zorder=3)
    ax.set_xticks(range(len(order)))
    ax.set_xticklabels(order)


def corr_panel(ax, x, y, dot_size=12, ci=0.95, line_color='black',
               annotate=True, loc='auto', fontsize=7, r2_as_percent=False,
               show=('eq', 'r2', 'p')):
    """Scatter + linear regression + CI band + annotation (Wald test on slope).
    show: any of 'eq' (formula), 'r' (Pearson r), 'r2', 'p' — pick the dialect
    (biof.2131 uses eq+r2+p; PeerJ 10244 uses r+P). loc: 'auto' places the text
    opposite the slope, or pass an (x, y) axes-fraction tuple. Returns
    (slope, intercept, r2, p)."""
    from scipy import stats
    x = np.asarray(x, dtype=float); y = np.asarray(y, dtype=float)
    m = ~(np.isnan(x) | np.isnan(y)); x, y = x[m], y[m]
    res = stats.linregress(x, y)
    xx = np.linspace(x.min(), x.max(), 100)
    yy = res.slope * xx + res.intercept
    n = len(x)
    tval = stats.t.ppf(0.5 + ci / 2, n - 2)
    resid = y - (res.slope * x + res.intercept)
    s_err = np.sqrt(np.sum(resid ** 2) / (n - 2))
    ci_band = tval * s_err * np.sqrt(1 / n + (xx - x.mean()) ** 2 / np.sum((x - x.mean()) ** 2))
    ax.scatter(x, y, s=dot_size, facecolor='black', edgecolor='none', zorder=3)
    ax.plot(xx, yy, color=line_color, lw=1.0, zorder=2)
    ax.plot(xx, yy + ci_band, color=line_color, lw=0.7, ls='--', zorder=2)
    ax.plot(xx, yy - ci_band, color=line_color, lw=0.7, ls='--', zorder=2)
    r2 = res.rvalue ** 2
    if annotate:
        sign = '+' if res.intercept >= 0 else '−'
        r2txt = f'{r2 * 100:.2f}%' if r2_as_percent else f'{r2:.3f}'
        if loc == 'auto':
            pos, ha = ((0.03, 0.97), 'left') if res.slope >= 0 else ((0.97, 0.97), 'right')
        elif isinstance(loc, (tuple, list)):
            pos, ha = tuple(loc), 'left'
        else:
            pos, ha = (0.03, 0.97), 'left'
        lines = []
        if 'eq' in show:
            lines.append(f'y = {res.slope:.4g}x {sign} {abs(res.intercept):.4g}')
        if 'r' in show:
            lines.append(f'r = {res.rvalue:.3f}')
        if 'r2' in show:
            lines.append(f'$R^2$ = {r2txt}')
        if 'p' in show:
            lines.append(f'P = {res.pvalue:.3g}')
        ax.annotate('\n'.join(lines), xy=pos, xycoords='axes fraction',
                    va='top', ha=ha, fontsize=fontsize)
    return res.slope, res.intercept, r2, res.pvalue


def _spearman_exact_p(a, b):
    """Exact permutation two-sided p for Spearman's rho at small n — matches
    GraphPad Prism's small-sample p (scipy's t-approximation differs). Falls
    back to scipy's asymptotic p for n>10 (exact enumeration too large)."""
    from scipy import stats
    from itertools import permutations
    a = np.asarray(a, float); b = np.asarray(b, float)
    n = len(a)
    rho = stats.spearmanr(a, b).correlation
    if n > 10:
        return rho, stats.spearmanr(a, b).pvalue
    ra = stats.rankdata(a)
    rb = stats.rankdata(b)
    obs = abs(np.corrcoef(ra, rb)[0, 1])
    perms = list(permutations(range(n)))
    cnt = sum(abs(np.corrcoef(ra, rb[list(p)])[0, 1]) >= obs - 1e-12 for p in perms)
    return rho, cnt / len(perms)


def _fmt_p(p, sig=4):
    """Prism-style p formatting: 4-decimal below 0.05, scientific with ×10 for
    tiny values (e.g. 0.0153, 4×10⁻⁶). Uses mathtext superscripts — Unicode
    superscript glyphs are missing from common fallback fonts (tofu boxes)."""
    if p >= 1e-3:
        return f'{p:.{sig}f}'
    mant, exp = f'{p:.0e}'.split('e')
    return f'{mant}×10$^{{{int(exp)}}}$'


def split_heatmap(ax, df, lower='spearman', upper='pearson', method=None,
                  cmap='RdBu_r', p_fontsize=5.5, label_fontsize=7, vmin=-1,
                  vmax=1, show_diag_line=True, diag='color', exact_spearman=False,
                  p_fmt=None, grid=False, highlight=None, annotate=True):
    """Correlation matrix, lower-left triangle = `lower` method, upper-right =
    `upper`; `method='pearson'` sets both (single-method heatmap). Significant
    cells annotated with bold white exact p-values (biof.1774 Fig 3A /
    biof.2131 Fig 6D dialect).

    diag: 'color' (r=1 cells), 'blank' (white), 'cross' (white with ×).
    exact_spearman: use permutation p for Spearman cells at small n (Prism).
    p_fmt: callable(p)->str for the annotation (default `_fmt_p`); None text
      when annotate=False.
    grid: draw white cell separators.
    highlight: (i0, i1, j0, j1) inclusive index box drawn with a green dashed
      border (biof.1774 green box dialect).
    df: samples × variables. Returns the image (for colorbar creation)."""
    if method:
        lower = upper = method
    from scipy import stats
    cols = list(df.columns)
    k = len(cols)
    r = np.zeros((k, k)); p = np.ones((k, k))
    for i in range(k):
        for j in range(k):
            a = df[cols[i]].astype(float); b = df[cols[j]].astype(float)
            m = ~(a.isna() | b.isna())
            if m.sum() < 3 or i == j:
                r[i, j] = 1.0 if i == j else np.nan
                continue
            meth = lower if i > j else upper
            if meth == 'spearman' and exact_spearman:
                r[i, j], p[i, j] = _spearman_exact_p(a[m], b[m])
            else:
                f = stats.spearmanr if meth == 'spearman' else stats.pearsonr
                r[i, j], p[i, j] = f(a[m], b[m])
    rr = r.copy()
    if diag != 'color':
        np.fill_diagonal(rr, np.nan)
    im = ax.imshow(rr, cmap=cmap, vmin=vmin, vmax=vmax)
    if grid:
        ax.set_xticks(np.arange(-.5, k), minor=True)
        ax.set_yticks(np.arange(-.5, k), minor=True)
        ax.grid(which='minor', color='white', lw=0.6)
    if diag == 'cross':
        for i in range(k):
            ax.text(i, i, '×', ha='center', va='center', fontsize=p_fontsize + 3,
                    color='0.4')
    fmt = p_fmt or _fmt_p
    if annotate:
        for i in range(k):
            for j in range(k):
                if i != j and p[i, j] < 0.05:
                    ax.text(j, i, fmt(p[i, j]), ha='center', va='center',
                            color='white', fontsize=p_fontsize, fontweight='bold')
    if show_diag_line:
        ax.plot([-0.5, k - 0.5], [-0.5, k - 0.5], color='black', lw=0.8)
    if highlight is not None:
        i0, i1, j0, j1 = highlight
        from matplotlib.patches import Rectangle
        ax.add_patch(Rectangle((j0 - 0.5, i0 - 0.5), j1 - j0 + 1, i1 - i0 + 1,
                               fill=False, edgecolor='#2ca02c', lw=1.2,
                               linestyle=(0, (4, 2)), zorder=5))
    ax.set_xticks(range(k)); ax.set_xticklabels(cols, rotation=45,
                                                ha='right', fontsize=label_fontsize)
    ax.set_yticks(range(k)); ax.set_yticklabels(cols, fontsize=label_fontsize)
    ax.tick_params(length=0)
    return im


def radar(ax, values_by_group, labels, colors=None, normalize='minmax',
          stat='median', alpha=0.25, label_fontsize=6, show_rticks=False,
          rtick_fontsize=5):
    """Radar/spider plot: per-variable min-max normalization over ALL pooled
    observations, then per-group medians (biof.2131 Fig 3D-F dialect). ax must
    be polar. values_by_group: {group: 2D array (observations × len(labels))}.
    show_rticks=True keeps the 0-1 radial ring labels visible."""
    groups = list(values_by_group)
    raw = {g: np.asarray(values_by_group[g], dtype=float) for g in groups}
    if normalize == 'minmax':
        # normalize per variable over ALL observations pooled (not over the
        # group aggregates — that degenerates to 0/1 polygons with 2 groups)
        allobs = np.vstack(list(raw.values()))
        lo, hi = np.nanmin(allobs, axis=0), np.nanmax(allobs, axis=0)
        span = np.where(hi - lo == 0, 1, hi - lo)
        raw = {g: (v - lo) / span for g, v in raw.items()}
    agg = {g: (np.nanmedian if stat == 'median' else np.nanmean)(v, axis=0)
           for g, v in raw.items()}
    M = np.vstack([agg[g] for g in groups])
    ang = np.linspace(0, 2 * np.pi, len(labels), endpoint=False)
    ang_c = np.concatenate([ang, ang[:1]])
    colors = colors or OKABE_ITO
    for i, g in enumerate(groups):
        v = np.concatenate([M[i], M[i][:1]])
        ax.plot(ang_c, v, color=colors[i], lw=1.2, label=g)
        ax.fill(ang_c, v, color=colors[i], alpha=alpha)
    ax.set_xticks(ang)
    ax.set_xticklabels(labels, fontsize=label_fontsize)
    if show_rticks:
        ax.set_rlim(0, 1)
        ax.set_rticks([0, 0.25, 0.5, 0.75, 1.0])
        ax.tick_params(axis='y', labelsize=rtick_fontsize)
    else:
        ax.set_yticklabels([])
    ax.tick_params(pad=-2)


def seahorse(ax, t, mean_by_group, err_by_group=None, colors=None,
             markers=None, sig=None, phase_lines=None, phase_labels=None,
             ylabel='OCR (pmol/min)', ms=3, lw=1.0):
    """Seahorse metabolic-flux time course (mito/glyco stress test) — mean ± SD
    per group over assay time, with optional injection phase separators and
    per-timepoint significance stars (eLife 83345 Fig 7E-H dialect).
    mean_by_group/err_by_group: {label: 1D over timepoints}. sig: {t_index: str}
    stars placed above the higher series. phase_lines: x positions of injections;
    phase_labels: list of (x, text) captions above the axis. Returns ax."""
    colors = colors or ['#c0392b', '#2a4d9b'] + OKABE_ITO
    markers = markers or ['s', 'o', '^', 'D']
    t = np.asarray(t, float)
    labels = list(mean_by_group)
    for i, g in enumerate(labels):
        m = np.asarray(mean_by_group[g], float)
        e = None if err_by_group is None else np.asarray(err_by_group[g], float)
        ax.errorbar(t, m, yerr=e, color=colors[i % len(colors)],
                    marker=markers[i % len(markers)], ms=ms, lw=lw,
                    elinewidth=0.7, capsize=1.5, label=g)
    for x in (phase_lines or []):
        ax.axvline(x, color='0.6', lw=0.7, ls='--', zorder=0)
    for x, txt in (phase_labels or []):
        ax.text(x, 1.02, txt, transform=ax.get_xaxis_transform(), ha='center',
                va='bottom', fontsize=6, color='0.35')
    if sig:
        top = np.vstack([np.asarray(mean_by_group[g], float) +
                         (0 if err_by_group is None else np.asarray(err_by_group[g], float))
                         for g in labels]).max(axis=0)
        for ti, s in sig.items():
            ax.text(t[ti], top[ti] * 1.04, s, ha='center', va='bottom', fontsize=7)
    ax.set_xlabel('Time (minutes)')
    ax.set_ylabel(ylabel)
    return ax


def timeseries_fdr(ax, t, mean_by_group, err_by_group=None, fdr=None,
                   colors=None, fdr_color='0.45', fdr_thresh=0.05,
                   ylabel_right='FDR', marker='o'):
    """Time series per group (mean ± err) + right-axis FDR values with a gray
    dashed threshold line (biof.2131 Fig 5 dialect). Returns the right axis."""
    colors = colors or OKABE_ITO
    for i, (g, m) in enumerate(mean_by_group.items()):
        e = None if err_by_group is None else err_by_group.get(g)
        ax.errorbar(t, m, yerr=e, color=colors[i], marker=marker, ms=3,
                    lw=1.0, elinewidth=0.7, capsize=1.5, label=g)
    ax2 = None
    if fdr is not None:
        ax2 = ax.twinx()
        ax2.plot(t, fdr, color=fdr_color, lw=0.9, ls=':', marker='.', ms=3)
        ax2.axhline(fdr_thresh, color='0.6', lw=0.8, ls='--')
        ax2.set_ylabel(ylabel_right, color=fdr_color)
        ax2.tick_params(axis='y', labelcolor=fdr_color)
        ax2.spines['top'].set_visible(False)
    return ax2


def volcano(ax, log2fc, padj, fc_thresh=1.0, p_thresh=0.05, dot_size=4,
            up_color='#D55E00', down_color='#0072B2', ns_color='0.75',
            genes=None, label_top=0, label_fontsize=5, label_by='padj'):
    """Volcano plot: x=log2FC, y=−log10(adjusted p), threshold dashed lines,
    colored up/down/ns. If `genes` (name array aligned to points) and
    `label_top>0`, label the top-N significant genes per direction with a
    greedy vertical declutter so labels don't overlap. label_by: 'padj'
    (most significant) or 'fc' (largest |log2FC|). Returns (n_up, n_down)."""
    log2fc = np.asarray(log2fc, float)
    padj = np.asarray(padj, float)
    y = -np.log10(np.clip(padj, 1e-300, None))
    up = (log2fc > fc_thresh) & (padj < p_thresh)
    dn = (log2fc < -fc_thresh) & (padj < p_thresh)
    ns = ~(up | dn)
    for m, c in ((ns, ns_color), (up, up_color), (dn, down_color)):
        ax.scatter(log2fc[m], y[m], s=dot_size, c=c, edgecolors='none',
                   rasterized=True)
    ax.axhline(-np.log10(p_thresh), ls='--', lw=0.7, c='0.4')
    for v in (-fc_thresh, fc_thresh):
        ax.axvline(v, ls='--', lw=0.7, c='0.4')
    ax.set_xlabel(r'log$_2$ fold change')
    ax.set_ylabel(r'−log$_{10}$ FDR')
    if genes is not None and label_top:
        genes = np.asarray(genes, dtype=object)
        rank = padj if label_by == 'padj' else -np.abs(log2fc)
        y0, y1 = ax.get_ylim()
        x0, x1 = ax.get_xlim()
        for mask, side in ((up, 1), (dn, -1)):
            idx = np.where(mask)[0]
            idx = idx[np.argsort(rank[idx])][:label_top]
            if len(idx) == 0:
                continue
            # stack labels in evenly spaced slots down from the top on that side
            idx = sorted(idx, key=lambda k: -y[k])
            slot_x = (x1 - 0.02 * (x1 - x0)) if side > 0 else (x0 + 0.02 * (x1 - x0))
            top_slot = y1 - 0.04 * (y1 - y0)
            step = 0.06 * (y1 - y0)
            for rank_i, i in enumerate(idx):
                sy = top_slot - rank_i * step
                ax.annotate(str(genes[i]), xy=(log2fc[i], y[i]),
                            xytext=(slot_x, sy),
                            fontsize=label_fontsize, va='center',
                            ha='right' if side > 0 else 'left',
                            arrowprops=dict(arrowstyle='-', lw=0.3, color='0.5',
                                            shrinkA=0, shrinkB=1))
    return int(up.sum()), int(dn.sum())


def pca_scatter(ax, M, groups=None, colors=None, dot_size=14, zscore=True,
                flip=(1, 1), zero_lines='dotted', color_map=None):
    """PCA scores plot from samples×variables matrix; % explained variance in
    the axis labels (in parentheses). PC signs are arbitrary — `flip=(sx, sy)`
    (each ±1) orients PC1/PC2 to match a target. color_map: {group: color}
    overrides the palette. zero_lines: 'dotted'/'solid'/None axes crosshair.
    Returns (scores, explained_pct)."""
    M = np.asarray(M, float)
    X = (M - M.mean(0)) / M.std(0, ddof=1) if zscore else M - M.mean(0)
    U, S, _ = np.linalg.svd(X, full_matrices=False)
    scores = U * S
    scores[:, 0] *= flip[0]
    scores[:, 1] *= flip[1]
    expl = S ** 2 / (S ** 2).sum() * 100
    colors = colors or OKABE_ITO
    if groups is None:
        ax.scatter(scores[:, 0], scores[:, 1], s=dot_size, c='black')
    else:
        groups = np.asarray(groups)
        for i, g in enumerate(dict.fromkeys(groups)):
            m = groups == g
            c = (color_map or {}).get(g, colors[i % len(colors)])
            ax.scatter(scores[m, 0], scores[m, 1], s=dot_size, color=c, label=str(g))
        ax.legend(fontsize=6)
    ax.set_xlabel(f'PC1 ({expl[0]:.2f}%)')
    ax.set_ylabel(f'PC2 ({expl[1]:.2f}%)')
    if zero_lines:
        ls = ':' if zero_lines == 'dotted' else '-'
        ax.axhline(0, lw=0.6, c='0.3', ls=ls, zorder=0)
        ax.axvline(0, lw=0.6, c='0.3', ls=ls, zorder=0)
    return scores, expl


def line_series(ax, x, series, colors=None, markers=None, smooth='akima',
                open_markers=True, lw=1.1, ms=4, fill_ref=None):
    """Single-measurement time-course lines with open markers and a
    monotone/Akima connecting curve (Prism 'connecting line' dialect;
    biof.1774 Fig 2I-T). NEVER adds error bars — for single measurements.
    series: {label: 1D y}. smooth: 'akima'|'pchip'|None (straight segments).
    Returns the axes."""
    colors = colors or ['black'] + OKABE_ITO
    markers = markers or ['o', '^', 'D', 's', 'v', 'P']
    x = np.asarray(x, float)
    for i, (lab, y) in enumerate(series.items()):
        y = np.asarray(y, float)
        c = colors[i % len(colors)]
        if smooth in ('akima', 'pchip') and len(x) > 2:
            try:
                from scipy.interpolate import Akima1DInterpolator, PchipInterpolator
                interp = (Akima1DInterpolator if smooth == 'akima' else PchipInterpolator)(x, y)
                xs = np.linspace(x.min(), x.max(), 200)
                ax.plot(xs, interp(xs), color=c, lw=lw, zorder=2)
            except Exception:
                ax.plot(x, y, color=c, lw=lw, zorder=2)
        else:
            ax.plot(x, y, color=c, lw=lw, zorder=2)
        mkw = (dict(mfc='none', mec=c) if open_markers else dict(mfc=c, mec=c))
        ax.plot(x, y, marker=markers[i % len(markers)], ls='none', ms=ms,
                markeredgewidth=1.0, label=lab, zorder=3, **mkw)
    return ax


def group_span_label(ax, x_from, x_to, text, y_axfrac=-0.34, text_gap=0.10,
                     lw=1.2, fontsize=8, fontweight='bold'):
    """Horizontal line under the tick labels spanning categories x_from..x_to
    (data coords) with a bold group name beneath (drug-span dialect)."""
    import matplotlib.transforms as mtrans
    trans = mtrans.blended_transform_factory(ax.transData, ax.transAxes)
    ax.plot([x_from - 0.3, x_to + 0.3], [y_axfrac, y_axfrac], transform=trans,
            color='black', lw=lw, clip_on=False, solid_capstyle='butt')
    ax.text((x_from + x_to) / 2, y_axfrac - text_gap, text, transform=trans,
            ha='center', va='top', fontsize=fontsize, fontweight=fontweight)


def broken_yaxis(fig, gridpos_or_ax, break_lo, break_hi, height_ratio=(1, 2.6),
                 d_mm=1.2, plot_fn=None, top_ylim=None, bot_ylim=None):
    """Replace one axes slot with a Prism-style broken y-axis: returns
    (ax_top, ax_bottom). The SAME artists must appear on both segments —
    pass `plot_fn(ax)` and it is called once per segment (prevents the silent
    data-loss trap of drawing on only one). top_ylim=(lo,hi) sets the upper
    segment (default (break_hi, None)); bot_ylim the lower (default (0, break_lo)).
    Diagonal break marks drawn on the left spine."""
    if hasattr(gridpos_or_ax, 'get_subplotspec'):
        ss = gridpos_or_ax.get_subplotspec()
        gridpos_or_ax.remove()
    else:
        ss = gridpos_or_ax
    sub = ss.subgridspec(2, 1, height_ratios=height_ratio, hspace=0.08)
    ax_top = fig.add_subplot(sub[0])
    ax_bot = fig.add_subplot(sub[1], sharex=ax_top)
    if plot_fn is not None:
        plot_fn(ax_top)
        plot_fn(ax_bot)
    ax_top.spines['bottom'].set_visible(False)
    ax_bot.spines['top'].set_visible(False)
    ax_top.tick_params(bottom=False, labelbottom=False)
    d = 0.012
    for ax, ys in ((ax_top, 0), (ax_bot, 1)):
        kw = dict(transform=ax.transAxes, color='black', clip_on=False, lw=0.8)
        ax.plot((-d * 2, d * 2), (ys - d, ys + d), **kw)
    ax_top.set_ylim(*(top_ylim or (break_hi, None)))
    ax_bot.set_ylim(*(bot_ylim or (0, break_lo)))
    return ax_top, ax_bot


def add_panel_labels(fig, axes, labels=None, dx_mm=-6.0, dy_mm=2.0,
                     fontsize=10, uppercase=True):
    """Bold panel letters aligned in figure coordinates near each axes' top-left."""
    import string
    labels = labels or list(string.ascii_uppercase if uppercase else string.ascii_lowercase)
    fw, fh = fig.get_size_inches()
    for ax, lab in zip(axes, labels):
        bb = ax.get_position()
        fig.text(bb.x0 + dx_mm * MM / fw, bb.y1 + dy_mm * MM / fh, lab,
                 fontsize=fontsize, fontweight='bold', ha='left', va='bottom')


def save_fig(fig, stem, formats=('pdf', 'png'), dpi=300, transparent=False):
    """Export to vector + raster. Never JPEG for data figures."""
    paths = []
    for ext in formats:
        p = f'{stem}.{ext}'
        fig.savefig(p, dpi=dpi, bbox_inches='tight', pad_inches=0.02,
                    transparent=transparent)
        paths.append(p)
    return paths
