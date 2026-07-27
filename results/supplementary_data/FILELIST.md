# Processed data tables

Small, versionable outputs of the pipeline. Large primary data are on Zenodo (see `data/*/FILELIST.md`).

| File | Contents |
|---|---|
| `Supplementary_Data_S1_DEGs_563.csv` | All 563 DEGs (adjusted *P* < 0.05, \|log2FC\| > 1) with DESeq2 statistics |
| `Supplementary_Data_S2_AI_gene_directions.csv` | Fold changes / adjusted *P* for the 31 PanelApp GREEN AI genes |
| `Supplementary_Data_S3_deconvolution_proportions.csv` | Per-sample cell-type proportions (NNLS) |
| `Supplementary_Data_S4_program_celltype_specificity.csv` | Cell-type specificity of each gene program |
| `Supplementary_Data_S5_module_scores.csv` | Per-sample program scores (combined gene-set *z*-score) |
| `AI_panel_PanelApp_v4.37.csv` | The pinned PanelApp panel with MGI mouse orthologs |
| `AI_enrichment_sensitivity.csv` | Fisher 2×2 counts, odds ratios and *P* across thresholds |
| `AI_gene_provenance_TableS1_S2.md` | Supplementary Tables S1–S2: verified OMIM/PMID provenance |

**Note on program-score statistics.** Program scores are reported as *effect sizes* with an
**exact permutation** *P*. With three animals per genotype the permutation floor is *P* = 0.10;
a parametric *t* test is not valid on these scores (see `pipeline/03_module_scores.py`).
