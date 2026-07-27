# Software environment

The pipeline is deterministic (fixed seeds). Step 05 (`05_scref_annotate.R`) is the only
stochastic step and sets `set.seed(1)`. Versions below are those used for the reported
analysis; the manuscript Supplementary Methods is the authoritative record.

## Primary bulk RNA-seq processing (upstream; see the analysis report)

| Tool | Version | Role |
|---|---|---|
| Trim Galore / Cutadapt | 0.6.10 / 4.9 | adapter + quality trimming |
| STAR | 2.7.10b | two-pass alignment to GRCm39 (Ensembl 116, `--sjdbOverhang 149`) |
| StringTie | — | de novo transcriptome assembly + merge |
| SQANTI3 | — | assembly annotation |
| Salmon | 1.10.0 | quasi-mapping quantification (libType A → ISR) |
| tximport | — | transcript → gene summarisation |
| RUVSeq | — | one factor of unwanted variation (k = 1, covariate `W_1`) |
| DESeq2 | — | likelihood-ratio test (`~W_1 + genotype` vs `~W_1`) |
| DEXSeq | — | differential exon usage (exploratory) |
| IsoformSwitchAnalyzeR | — | isoform switching + consequences (exploratory) |
| CIRCexplorer2 / circRNAprofiler | — | circular-RNA detection (exploratory) |

## This repository's re-analysis pipeline

**Python** (`pandas`, `numpy`, `scipy`, `matplotlib`) — steps 01–03, 06, 07, 10, 12–17.
**R** (`Seurat` 5.3.1, `Matrix`) — steps 05, 11.
**Pandoc** — `.docx` assembly (steps 10, 15).
**Graphviz** — workflow diagrams (step 14; Pandoc cannot render the Mermaid source directly).

External resources retrieved at run time (the code fails loudly rather than using a stored fallback):
- Genomics England **PanelApp** panel 269 "Amelogenesis imperfecta" v4.37 (pinned) — `17_fetch_ai_panel.py`
- **MGI** `HOM_MouseHumanSequence.rpt` for mouse↔human orthologs — `17_fetch_ai_panel.py`
- Public single-cell reference **GSE146855** — `04_download_scref.sh`

> Exact package versions not recorded by the original pipeline are marked as unrecorded in
> the Supplementary Methods rather than guessed.
