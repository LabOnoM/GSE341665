#!/usr/bin/env python3
"""
02_deg_analysis.py  --  Smad4 dcKO JDR re-analysis pipeline (step 2/7)

Differential-expression summary, curated gene-program directions, and enrichment
of human amelogenesis-imperfecta (AI) genes among down- vs up-regulated genes.

DEG threshold (project-specified): padj < 0.05 AND |log2FoldChange| > 1.
  -> reproduces the 563 significant DEGs (209 up / 354 down) in the source table.
Note: core enamel-matrix genes (Amelx/Ambn/Enam/Mmp20) have |log2FC|~0.4 and do
NOT pass |log2FC|>1 despite padj 1e-66..1e-103; they are analysed as a gene set.

Inputs : output/DEG_table.csv
Outputs: output/DEG_significant.csv, output/program_directions.csv,
         output/AI_gene_directions.csv, output/DEG_summary.txt
"""
import os, sys, pandas as pd, numpy as np
sys.path.insert(0, os.path.dirname(os.path.abspath(__file__)))
from scipy import stats

ROOT = os.path.dirname(os.path.dirname(os.path.dirname(os.path.abspath(__file__))))
OUT  = os.path.join(ROOT, "scripts", "JDR_reanalysis", "output")

PADJ, LFC = 0.05, 1.0
# AI gene panel: imported, never re-typed (lessons-learned §4). Sourced from
# Genomics England PanelApp panel 269 (GREEN genes) with MGI orthologs — see
# ai_genes.py and 17_fetch_ai_panel.py. The former hand-typed list from
# scripts/00_Mendelian_AI_Genes.R is QUARANTINED: its provenance column was found
# to be largely fabricated on 2026-07-27 (18/20 PMIDs cite unrelated papers).
import ai_genes
from programs import PROGRAMS as _P, IDENTITY_MATRIX, IDENTITY_TFS
PROGRAMS = {
    "Sterol_biosynthesis_SREBP2": _P["Sterol"],
    "Interferon_ISG": _P["IFN"],
    "Wound_stress_AP1_EGFR": _P["Wound"],
    "Ameloblast_master_TFs": IDENTITY_TFS,
    "Enamel_matrix_proteases": IDENTITY_MATRIX,
}


def main():
    df = pd.read_csv(os.path.join(OUT, "DEG_table.csv"), skiprows=1, low_memory=False)
    for c in ["baseMean","log2FoldChange","pvalue","padj"]:
        df[c] = pd.to_numeric(df[c], errors="coerce")
    # NOTE: do NOT resolve AI genes through hsapiens_homolog_associated_gene_name.
    # It is NaN for ~41% of tested genes, and `.astype(str)` turns NaN into the string
    # "nan", so `.isin()` silently scores those genes as NON-risk instead of flagging
    # them. That mis-scored 7 of the old 20 panel genes. Match on MGI mouse symbols
    # via ai_genes.resolve_in_table(), which raises rather than truncating.
    df["hs"] = df["hsapiens_homolog_associated_gene_name"].astype(str).str.upper()
    lines = []
    def log(s):
        print(s); lines.append(s)

    sig = df[(df.padj < PADJ) & (df.log2FoldChange.abs() > LFC)].copy()
    log(f"DEGs at padj<{PADJ} & |log2FC|>{LFC}: {len(sig)}  "
        f"(up={int((sig.log2FoldChange>0).sum())}, down={int((sig.log2FoldChange<0).sum())})")
    sig.sort_values("log2FoldChange").to_csv(os.path.join(OUT, "DEG_significant.csv"), index=False)

    # program directions (strict |LFC|>1 and lenient padj<.05,|LFC|>.25)
    log("\nProgram directions [strict |LFC|>1 : up/down] [lenient padj<.05,|LFC|>.25 : up/down] meanLFC")
    prog_rows = []
    for name, genes in PROGRAMS.items():
        sub = df[df.external_gene_name.isin(genes)]
        su = int(((sub.padj<PADJ)&(sub.log2FoldChange> 1)).sum())
        sd = int(((sub.padj<PADJ)&(sub.log2FoldChange<-1)).sum())
        lu = int(((sub.padj<PADJ)&(sub.log2FoldChange> 0.25)).sum())
        ld = int(((sub.padj<PADJ)&(sub.log2FoldChange<-0.25)).sum())
        log(f"  {name:28s} [{su}/{sd}] [{lu}/{ld}] {sub.log2FoldChange.mean():+.2f}")
        prog_rows.append(dict(program=name, n_tested=len(sub), strict_up=su, strict_down=sd,
                              lenient_up=lu, lenient_down=ld, mean_log2FC=round(sub.log2FoldChange.mean(),3)))
    pd.DataFrame(prog_rows).to_csv(os.path.join(OUT, "program_directions.csv"), index=False)

    # AI-gene enrichment among DOWN vs UP (Fisher, one-sided greater) at 3 thresholds
    log("\nMendelian AI-gene enrichment among DOWN vs UP:")
    for th, lab in [(1.0,"|LFC|>1"), (0.25,"|LFC|>0.25"), (0.0,"padj<0.05 any LFC")]:
        t = df[df.padj.notna()].copy()
        t["s"] = (t.padj<PADJ) & (t.log2FoldChange.abs()>th)
        s = t[t.s]; s_dir = np.where(s.log2FoldChange>0,"UP","DOWN")
        s_risk = ai_genes.resolve_in_table(s, strict=False)
        a=int(((s_dir=="DOWN")& s_risk).sum()); b=int(((s_dir=="DOWN")&~s_risk).sum())
        c=int(((s_dir=="UP")  & s_risk).sum()); d=int(((s_dir=="UP")  &~s_risk).sum())
        OR,p = stats.fisher_exact([[a,b],[c,d]], alternative="greater")
        log(f"  {lab:18s} risk_down={a} nonrisk_down={b} risk_up={c} nonrisk_up={d}  OR={OR:.2f} p={p:.3g}")

    ai = (df[ai_genes.resolve_in_table(df, strict=False)][["hsapiens_homolog_associated_gene_name","external_gene_name",
          "log2FoldChange","padj","baseMean"]].dropna(subset=["padj"]).sort_values("log2FoldChange"))
    ai.to_csv(os.path.join(OUT, "AI_gene_directions.csv"), index=False)
    log(f"\nAI panel ({ai_genes.PANEL_SOURCE}) genes tested: {len(ai)}  "
        f"(sig down={int(((ai.padj<PADJ)&(ai.log2FoldChange<-0.25)).sum())}, "
        f"sig up={int(((ai.padj<PADJ)&(ai.log2FoldChange>0.25)).sum())})")
    # KMT2D sanity (must be non-significant to justify dropping the KMT2D figure)
    k = df[df.external_gene_name=="Kmt2d"]
    if len(k):
        log(f"KMT2D check: log2FC={k.log2FoldChange.iloc[0]:+.2f} padj={k.padj.iloc[0]:.2g} "
            f"(NOT a DEG -> the KMT2D/ENCODE figure is unsupported)")
    with open(os.path.join(OUT, "DEG_summary.txt"), "w") as fh:
        fh.write("\n".join(lines) + "\n")


if __name__ == "__main__":
    main()
