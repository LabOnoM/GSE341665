#!/usr/bin/env python3
"""
14_workflow_diagram.py -- render the analysis workflow as publication figures.

The Supplementary Methods carries the workflow as a Mermaid block, which is
readable in markdown viewers but is NOT rendered by pandoc when building the
.docx. This script renders the identical workflow with Graphviz so the Word
supplementary contains real figures rather than a code listing.

The workflow is drawn as TWO panels. A single combined diagram is ~16 inches
wide; scaled into a portrait Word page its type falls to ~4 pt, which is not
legible in print. Splitting keeps each panel near page width with ~8 pt type.

  Panel A : primary bulk RNA-sequencing pipeline (Part A of the Methods)
  Panel B : re-analysis performed for this manuscript (Part B)

Output: figures/FigureS4A_pipeline.{png,pdf}, figures/FigureS4B_reanalysis.{png,pdf}
"""
import os, subprocess

ROOT = os.path.dirname(os.path.dirname(os.path.dirname(os.path.abspath(__file__))))
FIG = os.path.join(ROOT, "scripts", "JDR_reanalysis", "figures")
os.makedirs(FIG, exist_ok=True)

RED, GREY = "#B2182B", "#8A8A8A"
HEAD = """
digraph g {
  graph [rankdir=TB, splines=polyline, nodesep=0.22, ranksep=0.30,
         fontname="Helvetica", bgcolor="white"];
  node  [shape=box, style="rounded,filled", fontname="Helvetica", fontsize=13,
         margin="0.12,0.08", color="#666666", penwidth=1.0, fillcolor="#F4F7FB"];
  edge  [color="#555555", arrowsize=0.6, penwidth=1.0];
"""

PANEL_A = HEAD + r"""
  FQ  [label="Paired-end FASTQ, 150 bp\n3 CON + 3 dcKO; 21.0–22.5 M read pairs/sample", fillcolor="#EFEFEF"];
  TG  [label="A3  FastQC → Trim Galore 0.6.10 / cutadapt 4.9\n(Phred 20, adapter auto-detected) → FastQC"];
  ST  [label="A4  STAR 2.7.10b two-pass alignment\nEnsembl 116 GRCm39 index (--sjdbOverhang 149)\n93.66–93.84 % uniquely mapped"];
  SQ  [label="A5  StringTie merge → Salmon 1.10.0 (ISR)\n→ tximport → gene-level counts\n(featureCounts -s 2 -p -Q 20 -O also computed)", fillcolor="#E8F0FB"];
  FLT [label="A6  Filter ≥10 reads in ≥3 samples → 15,803 genes\nRUVSeq: 1 factor of unwanted variation (W_1)"];
  LRT [label="DESeq2 likelihood-ratio test\n~W_1 + genotype   vs   ~W_1"];
  SIG [label="padj < 0.05  and  |log2FC| > 1\n563 DEGs: 209 up, 354 down",
       fillcolor="#FBE9E9", color="__RED__", penwidth=2.0];
  GO  [label="A7  GO + KEGG enrichment at uncorrected P < 0.05\n— qualitative corroboration only",
       fillcolor="#F7F7F7", color="__GREY__", style="rounded,filled,dashed", fontcolor="#555555"];
  UN  [label="A8  Run but NOT used here (reasons in Supplementary Methods):\nDEXSeq · IsoformSwitchAnalyzeR · LeafCutter · CIRCexplorer2 · WGCNA",
       fillcolor="#F2F2F2", color="#AAAAAA", style="rounded,filled,dashed", fontcolor="#555555"];
  OUT [label="→ continues in Panel B", shape=plaintext, fillcolor="white", fontcolor="__RED__"];

  FQ -> TG -> ST -> SQ -> FLT -> LRT -> SIG -> GO -> UN -> OUT;
}
""".replace("__RED__", RED).replace("__GREY__", GREY)

PANEL_B = HEAD + r"""
  IN  [label="563 DEGs + gene-level TPM matrix (Panel A)", fillcolor="#EFEFEF"];
  P1  [label="01  Extract DEG and TPM tables", fillcolor="#E8F0FB"];
  P2  [label="02  DEG summary; enrichment of human amelogenesis-\nimperfecta genes, PanelApp 269 v4.37 GREEN\n(OR 9.21, P = 1.7e-4)", fillcolor="#E8F0FB"];
  P3  [label="03  Program scores (combined gene-set z-score)", fillcolor="#E8F0FB"];
  PR  [label="programs.py — single source of truth\nIdentity 16 · Sterol 20 · Wound 14 · IFN 18",
       shape=note, fillcolor="#FBE9E9", color="__RED__", penwidth=1.8];
  D1  [label="04  Download GSE146855 (P7 incisor, 6,248 cells)", fillcolor="#E8F0FB"];
  D2  [label="05  Seurat QC → PCA (20 PCs) → clusters (res 0.6)\n→ specificity-based de-novo annotation\n→ cell-type CPM signature", fillcolor="#E8F0FB"];
  D3  [label="11  UMAP + marker dot plot (Suppl. Figure S3)", fillcolor="#E8F0FB"];
  D4  [label="06  NNLS deconvolution (542 marker genes)", fillcolor="#E8F0FB"];
  F   [label="07  Figures 4–6", fillcolor="#E8F0FB"];
  DOC [label="10  Assemble manuscript\n(auto-numbered citations + embedded figures)", fillcolor="#E8F0FB"];
  GEO [label="12 / 13  GEO submission package", fillcolor="#E8F0FB"];

  IN -> P1 -> P2 -> PR -> P3 -> D1 -> D2;
  D2 -> D3; D2 -> D4;
  D3 -> F [style=invis]; D4 -> F;
  F -> DOC -> GEO;
  {rank=same; D3; D4;}
}
""".replace("__RED__", RED)

def render(name, dot):
    dp = os.path.join(FIG, f"{name}.dot")
    open(dp, "w").write(dot)
    for fmt, extra in (("png", ["-Gdpi=300"]), ("pdf", [])):
        out = os.path.join(FIG, f"{name}.{fmt}")
        subprocess.run(["dot", f"-T{fmt}", *extra, dp, "-o", out], check=True)
    from PIL import Image
    w, h = Image.open(os.path.join(FIG, f"{name}.png")).size
    print(f"  {name}: {w/300:.1f} x {h/300:.1f} in  ->  at 6.5in wide, type ≈ "
          f"{13*6.5/(w/300):.1f} pt")


if __name__ == "__main__":
    render("FigureS4A_pipeline", PANEL_A)
    render("FigureS4B_reanalysis", PANEL_B)
    # remove the superseded single-panel version so it cannot be embedded by mistake
    for ext in ("png", "pdf", "dot"):
        old = os.path.join(FIG, f"FigureS4_analysis_workflow.{ext}")
        if os.path.exists(old):
            os.remove(old)
