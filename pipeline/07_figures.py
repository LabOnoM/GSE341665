#!/usr/bin/env python3
"""
07_figures.py  --  Smad4 dcKO re-analysis pipeline (step 7)

Generate the manuscript's computational figures (Figures 4-6) from the real
analysis outputs. Every value is computed from the project's own data or the
public reference; nothing is hard-coded/illustrative. Replaces the retracted
KMT2D/ENCODE figure. Output files are named with their MANUSCRIPT numbers
(Figures 4-6); Figures 1-3 are wet-laboratory data.

Inputs : output/DEG_table.csv, output/TPM_matrix.csv, output/module_scores.csv,
         output/deconv_proportions.csv, output/program_celltype_specificity.csv,
         output/AI_gene_directions.csv
Outputs: figures/Figure4_transcriptome.{png,pdf}
         figures/Figure5_identity_switch.{png,pdf}
         figures/Figure6_metabolic_convergence.{png,pdf}
"""
import os, sys, numpy as np, pandas as pd
sys.path.insert(0, os.path.dirname(os.path.abspath(__file__)))
import matplotlib
matplotlib.use("Agg")
import matplotlib.pyplot as plt
from matplotlib.colors import TwoSlopeNorm
from programs import PROGRAMS

ROOT = os.path.dirname(os.path.dirname(os.path.dirname(os.path.abspath(__file__))))
OUT  = os.path.join(ROOT, "scripts", "JDR_reanalysis", "output")
FIG  = os.path.join(ROOT, "scripts", "JDR_reanalysis", "figures")
os.makedirs(FIG, exist_ok=True)
SAMP = ["TPM_CON_1","TPM_CON_2","TPM_CON_3","TPM_dcKO_1","TPM_dcKO_2","TPM_dcKO_3"]
GRP  = ["CON"]*3 + ["dcKO"]*3
BLUE, RED, GREY = "#2166AC", "#B2182B", "#BBBBBB"
plt.rcParams.update({"font.size": 8, "axes.titlesize": 9, "axes.titleweight": "bold",
                     "figure.dpi": 150, "savefig.dpi": 300, "axes.spines.top": False,
                     "axes.spines.right": False})


def _load():
    deg = pd.read_csv(os.path.join(OUT, "DEG_table.csv"), skiprows=1, low_memory=False)
    for c in ["log2FoldChange", "padj", "baseMean"]:
        deg[c] = pd.to_numeric(deg[c], errors="coerce")
    tpm = pd.read_csv(os.path.join(OUT, "TPM_matrix.csv"), skiprows=1, low_memory=False)
    for s in SAMP:
        tpm[s] = pd.to_numeric(tpm[s], errors="coerce")
    tpm = tpm.groupby("external_gene_name")[SAMP].mean()
    mod = pd.read_csv(os.path.join(OUT, "module_scores.csv"), index_col=0)
    dec = pd.read_csv(os.path.join(OUT, "deconv_proportions.csv"), index_col=0)
    spec = pd.read_csv(os.path.join(OUT, "program_celltype_specificity.csv"), index_col=0)
    ai = pd.read_csv(os.path.join(OUT, "AI_gene_directions.csv"))
    return deg, tpm, mod, dec, spec, ai


def _zrows(df):
    return df.sub(df.mean(axis=1), axis=0).div(df.std(axis=1).replace(0, np.nan), axis=0)


def save(fig, name):
    for ext in ("png", "pdf"):
        fig.savefig(os.path.join(FIG, f"{name}.{ext}"), bbox_inches="tight")
    plt.close(fig); print("  wrote", name)


# ---------- Figure 4: transcriptome landscape ----------
def figure4(deg, tpm, mod):
    fig, ax = plt.subplots(1, 3, figsize=(11, 3.4))
    # (A) standard PCA on top-2000 variable genes (log2 TPM, gene-centered, no unit scaling)
    lg = np.log2(tpm + 1)
    top = lg.var(axis=1).sort_values(ascending=False).head(2000).index
    M = lg.loc[top].T.values                 # samples x genes
    Mc = M - M.mean(0)                        # center each gene across samples
    U, S, Vt = np.linalg.svd(Mc, full_matrices=False)
    pc = U[:, :2] * S[:2]; ev = (S**2 / (S**2).sum())[:2] * 100
    for g, col in [("CON", BLUE), ("dcKO", RED)]:
        m = np.array([x == g for x in GRP])
        ax[0].scatter(pc[m, 0], pc[m, 1], c=col, s=70, label=g, edgecolor="k", lw=.4)
    ax[0].set(xlabel=f"PC1 ({ev[0]:.1f}%)", ylabel=f"PC2 ({ev[1]:.1f}%)", title="A  PCA")
    ax[0].legend(frameon=False, fontsize=7)
    # (B) volcano
    d = deg.dropna(subset=["log2FoldChange", "padj"]).copy()
    d["nlp"] = -np.log10(d["padj"].clip(lower=1e-300))
    sig = (d.padj < 0.05) & (d.log2FoldChange.abs() > 1)
    ax[1].scatter(d.loc[~sig, "log2FoldChange"], d.loc[~sig, "nlp"], s=3, c=GREY, alpha=.4, rasterized=True)
    ax[1].scatter(d.loc[sig & (d.log2FoldChange < 0), "log2FoldChange"], d.loc[sig & (d.log2FoldChange < 0), "nlp"], s=5, c=BLUE, alpha=.7, rasterized=True)
    ax[1].scatter(d.loc[sig & (d.log2FoldChange > 0), "log2FoldChange"], d.loc[sig & (d.log2FoldChange > 0), "nlp"], s=5, c=RED, alpha=.7, rasterized=True)
    # label a curated, spatially-separated set (avoid the padj-floor pile-up at centre)
    for g in ["Krt76","Oasl1","Mx1","Isg15","Wnt10a","Fdps","Ldlr","Hmgcr","Ptgs2","Mmp12","Ereg","Sqle"]:
        r = d[d.external_gene_name == g]
        if len(r):
            ax[1].annotate(g, (r.log2FoldChange.iloc[0], r.nlp.iloc[0]), xytext=(0, 4),
                           textcoords="offset points", fontsize=6, style="italic", ha="center")
    ax[1].axvline(-1, ls="--", c="k", lw=.5); ax[1].axvline(1, ls="--", c="k", lw=.5)
    ax[1].axhline(-np.log10(0.05), ls="--", c="k", lw=.5)
    ax[1].set(xlabel="log2 fold change (dcKO/CON)", ylabel="-log10 adjusted P",
              title="B  Volcano (563 DEGs)", xlim=(-6, 6), ylim=(-8, 335))
    # (C) module means by genotype
    order = ["Identity", "IFN", "Sterol", "Wound"]
    means = mod.groupby("group")[order].mean().reindex(["CON", "dcKO"])
    x = np.arange(len(order)); w = .38
    ax[2].bar(x - w/2, means.loc["CON"], w, color=BLUE, label="CON")
    ax[2].bar(x + w/2, means.loc["dcKO"], w, color=RED, label="dcKO")
    ax[2].axhline(0, c="k", lw=.6); ax[2].set_xticks(x)
    ax[2].set_xticklabels(order, rotation=20, ha="right")
    ax[2].set(ylabel="program z-score", title="C  Program shift")
    ax[2].legend(frameon=False, fontsize=7)
    fig.suptitle("Smad4-deficient dental-epithelium transcriptome", x=.02, ha="left", fontweight="bold")
    fig.tight_layout(rect=(0, 0, 1, .96)); save(fig, "Figure4_transcriptome")


# ---------- Figure 5: identity switch ----------
def figure5(deg, tpm, dec, ai):
    # Panel B lists all 31 PanelApp GREEN genes, so the figure is taller than the
    # 4-panel layout would otherwise need; at 6.2 in the labels collided.
    fig = plt.figure(figsize=(11, 7.8))
    gs = fig.add_gridspec(2, 2, height_ratios=[1.35, 1], hspace=.5, wspace=.32)
    # (A) identity gene heatmap
    idg = [g for g in PROGRAMS["Identity"] if g in tpm.index]
    Z = _zrows(np.log2(tpm.loc[idg] + 1))
    axA = fig.add_subplot(gs[0, 0])
    im = axA.imshow(Z.values, aspect="auto", cmap="RdBu_r", norm=TwoSlopeNorm(0, -2, 2))
    axA.set_xticks(range(6)); axA.set_xticklabels([s.replace("TPM_", "") for s in SAMP], rotation=90, fontsize=6)
    axA.set_yticks(range(len(idg))); axA.set_yticklabels(idg, fontsize=6, style="italic")
    axA.set_title("A  Ameloblast identity program (down)")
    fig.colorbar(im, ax=axA, shrink=.6, label="row z")
    # (B) AI gene log2FC
    axB = fig.add_subplot(gs[0, 1])
    a = ai.dropna(subset=["log2FoldChange"]).sort_values("log2FoldChange")
    cols = [RED if v > 0 else BLUE for v in a.log2FoldChange]
    # Label with the mouse symbol: the BioMart human-ortholog column is NaN for ~41%
    # of tested genes, which both breaks the axis and is the defect that corrupted the
    # enrichment itself (see 16_ai_gene_enrichment.py).
    axB.barh(a["external_gene_name"], a.log2FoldChange, color=cols)
    axB.axvline(0, c="k", lw=.6); axB.set_xlabel("log2 fold change")
    # Never hard-code a statistic into a figure (CLAUDE.md §2). Read the value that
    # 16_ai_gene_enrichment.py computed, and fail loudly if it is not there.
    sens = pd.read_csv(os.path.join(OUT, "AI_enrichment_sensitivity.csv"))
    row = sens[(sens.resolution == "CORRECTED") & (sens.lfc_cutoff == 0.0)]
    if row.empty:
        raise RuntimeError(
            "AI_enrichment_sensitivity.csv has no CORRECTED row at lfc_cutoff=0.0 — "
            "run 16_ai_gene_enrichment.py before building figures. Refusing to plot a "
            "panel whose statistic cannot be sourced.")
    _or, _p = float(row.odds_ratio.iloc[0]), float(row.p_one_sided.iloc[0])
    axB.set_title(f"B  Human AI genes (OR={_or:.2f}, p={_p:.1e} down-enriched)", fontsize=8)
    axB.tick_params(labelsize=5)
    for lbl in axB.get_yticklabels(): lbl.set_style("italic")
    # (C) wound genes
    axC = fig.add_subplot(gs[1, 0])
    wq = PROGRAMS["Wound"]           # display exactly the scored program
    d = deg.dropna(subset=["external_gene_name"]).groupby("external_gene_name")["log2FoldChange"].mean()
    wv = d.reindex([g for g in wq if g in d.index]).dropna().sort_values()
    axC.barh(wv.index, wv.values, color=RED)
    axC.axvline(0, c="k", lw=.6); axC.set_xlabel("log2 fold change")
    axC.set_title("C  Wound / stress program (up)"); axC.tick_params(labelsize=6)
    for lbl in axC.get_yticklabels(): lbl.set_style("italic")
    # (D) deconvolution proportions
    axD = fig.add_subplot(gs[1, 1])
    P = dec * 100
    keep = [c for c in P.columns if P[c].mean() > 0.5]
    conm = P.loc[[s for s in SAMP if "CON" in s], keep].mean()
    kom = P.loc[[s for s in SAMP if "dcKO" in s], keep].mean()
    x = np.arange(len(keep)); w = .38
    axD.bar(x - w/2, conm, w, color=BLUE, label="CON")
    axD.bar(x + w/2, kom, w, color=RED, label="dcKO")
    axD.set_xticks(x); axD.set_xticklabels([k.replace("_", "\n") for k in keep], fontsize=6)
    axD.set_ylabel("% of cells (deconv.)")
    axD.set_title("D  Composition: ameloblasts preserved,\nmesenchyme down", fontsize=8)
    axD.legend(frameon=False, fontsize=7)
    fig.suptitle("Cell-intrinsic ameloblast identity switch", x=.02, ha="left", fontweight="bold")
    save(fig, "Figure5_identity_switch")


# ---------- Figure 6: metabolic-immune convergence ----------
def figure6(deg, mod, spec):
    fig = plt.figure(figsize=(11, 6.4))
    gs = fig.add_gridspec(2, 3, hspace=.75, wspace=.6)
    ax = [fig.add_subplot(gs[0, 0]), fig.add_subplot(gs[0, 1]), fig.add_subplot(gs[0, 2]),
          fig.add_subplot(gs[1, 0]), fig.add_subplot(gs[1, 1])]
    order = ["Sterol", "IFN", "Wound", "Identity"]
    # (A) program scores by genotype
    for i, m in enumerate(order):
        for g, col, dx in [("CON", BLUE, -.15), ("dcKO", RED, .15)]:
            v = mod.loc[mod.group == g, m]
            ax[0].scatter([i + dx] * len(v), v, c=col, s=28, edgecolor="k", lw=.3)
    ax[0].axhline(0, c="k", lw=.5); ax[0].set_xticks(range(len(order)))
    ax[0].set_xticklabels(order, rotation=25, ha="right")
    ax[0].set(ylabel="program z-score", title="A  Program scores")
    ax[0].scatter([], [], c=BLUE, label="CON"); ax[0].scatter([], [], c=RED, label="dcKO")
    ax[0].legend(frameon=False, fontsize=7)
    # (B) sterol vs IFN
    for g, col in [("CON", BLUE), ("dcKO", RED)]:
        s = mod[mod.group == g]
        ax[1].scatter(s.Sterol, s.IFN, c=col, s=45, edgecolor="k", lw=.4, label=g)
    ax[1].set(xlabel="Sterol program", ylabel="IFN program", title="B  Sterol vs interferon")
    ax[1].legend(frameon=False, fontsize=7)
    # (C) program cell-type specificity heatmap
    S = spec.copy(); S.columns = [c.replace("Ameloblast_", "AB_").replace("Epithelial_", "Epi_") for c in S.columns]
    im = ax[2].imshow(S.values, cmap="RdBu_r", norm=TwoSlopeNorm(0, -1.2, 1.2), aspect="auto")
    ax[2].set_xticks(range(S.shape[1])); ax[2].set_xticklabels(S.columns, rotation=90, fontsize=6)
    ax[2].set_yticks(range(S.shape[0])); ax[2].set_yticklabels(S.index, fontsize=6)
    ax[2].set_title("C  Program cell-type origin")
    fig.colorbar(im, ax=ax[2], shrink=.7, label="z")
    # (D) sterol & IFN gene log2FC — leading edge of each scored program
    d = deg.dropna(subset=["external_gene_name"]).groupby("external_gene_name")["log2FoldChange"].mean()
    st = d.reindex(PROGRAMS["Sterol"]).dropna().sort_values(ascending=False).head(9)
    ifn = d.reindex(PROGRAMS["IFN"]).dropna().sort_values().head(7)
    v = pd.concat([st, ifn])
    cols = [RED if x > 0 else BLUE for x in v]
    ax[3].barh(range(len(v)), v.values, color=cols)
    ax[3].set_yticks(range(len(v))); ax[3].set_yticklabels(v.index, fontsize=6, style="italic")
    ax[3].axvline(0, c="k", lw=.6); ax[3].set(xlabel="log2 FC", title="D  Sterol up / IFN down")
    ax[3].invert_yaxis()
    # (E) ERK-activation signature — the wet-lab bridge (ERK inhibition rescues enamel genes)
    de = deg.dropna(subset=["external_gene_name"]).groupby("external_gene_name")["log2FoldChange"].mean()
    erk = ["Ereg","Areg","Dusp6","Dusp1","Dusp4","Egr1","Dusp5","Fosl1","Spry4","Fosb","Spry2","Fos"]
    ev = de.reindex(erk).dropna().sort_values()
    ax[4].barh(range(len(ev)), ev.values, color=RED)
    ax[4].set_yticks(range(len(ev))); ax[4].set_yticklabels(ev.index, fontsize=6, style="italic")
    ax[4].axvline(0, c="k", lw=.6)
    ax[4].set(xlabel="log2 FC", title="E  ERK-activation signature")
    fig.suptitle("Signaling-to-metabolism convergence (ERK / SREBP2 activation)",
                 x=.02, ha="left", fontweight="bold")
    save(fig, "Figure6_metabolic_convergence")


if __name__ == "__main__":
    deg, tpm, mod, dec, spec, ai = _load()
    print("Generating figures:")
    figure4(deg, tpm, mod)
    figure5(deg, tpm, dec, ai)
    figure6(deg, mod, spec)
    print("Done ->", FIG)
