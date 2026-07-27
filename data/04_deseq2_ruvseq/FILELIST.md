# File manifest — Batch correction and differential expression

> **These files are not stored in Git** (too large for a code repository).
> They are archived on **Zenodo** (DOI: `[[Zenodo DOI — insert on release]]`) and,
> where noted, are also browsable in the interactive analysis report at
> <http://smad4kotooth.bs-gou.com/>. Download them into **this folder** to reproduce
> the pipeline; paths in the code resolve relative to the repository root.

| Expected file | Description | Approx. size | Source |
|---|---|---|---|
| `Parameters.Rdata` | Recorded pipeline parameters (Par_ReN, Par_SaN, thresholds) | <1 MB | Zenodo |
| `DESeq2_Processing.Rdata` | RUVSeq (k=1) + DESeq2 LRT objects | ~10-50 MB | Zenodo |
| `DEG_table_full.xlsx` | All 15,803 tested genes (log2FC, padj, annotation) | ~5-15 MB | Zenodo — small CSV mirror in results/supplementary_data/ |

_This folder is intentionally otherwise empty; the manifest keeps it under version control._
