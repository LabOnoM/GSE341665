# Manuscript reproduction map

| Figure, table or result | Retained code | Required source and purpose |
|---|---|---|
| Main Figures 4–7; full Data S3 | 01_extract_tables.py, 03_module_scores.py, 06_deconvolution.py, 07_figures.py; volcano_panel.py, vendor/figstyle.py | Two source workbooks; DEG/TPM and ordered text-gene index; seven-cell-type CPM signature. Preserves original gene-level model statistics and all displayed genes |
| Four scored programs; Additional file 5 | programs.py, 03_module_scores.py | The single definition of 16 identity, 20 sterol, 14 wound/stress and 18 interferon genes; full-precision TPM; descriptive mean gene-wise z scores (no score-inference route) |
| Reference fractions/specificity; Additional files 6–7 | 06_deconvolution.py | Captured GSE146855 CPM signature; full-precision bulk TPM; sorted markers before NNLS |
| AI enrichment, Additional file 3 and Methods B6 | 16_ai_gene_enrichment.py, ai_genes.py, ai_gene_provenance.py | Pinned PanelApp 269 v4.37 GREEN panel and gene-level statistics; three reported threshold sensitivity tests |
| Pinned-panel source reconstruction | 17_fetch_ai_panel.py | Original PanelApp v4.37 JSON and MGI report in the deposit; explicit live download also tested. `all` reconstructs the mapped panel before enrichment |
| Data S1/S2; Figure 3C/E statistical evidence | 24_september_inhibitor_audit.py | The declared inhibitor data sheets in Figure_3.prism; 195 values and 52 ANOVA/two-sided Dunnett comparisons. The author-provided wet-laboratory artwork is preserved |
| Additional file 2 workbook; files 2–7; Table S4 | export_additional_data.py, ai_gene_provenance.py | Latest regenerated CSV/payload; public openpyxl. Table S4 contains the 31 active PanelApp genes; export checks its human–mouse pairs against the version-pinned analysis panel |
| Supplementary Figure S6 / reference reannotation | 04_download_scref.sh, 05_scref_annotate.R, 11_scref_umap_dotplot.R | Original deposited GSE146855 raw archive; exact R/library lock; fresh labels/signatures checked against captured reference before downstream use |
| Supplementary Figure S5 | 14_workflow_diagram.py | Single composite A–B; recorded primary methods and downstream analyses; matplotlib |
| Supplementary Figure S3 rendering | 21_smad4_loxp_figure.py | Panels A–C from locus GTF/evidence CSVs, summary JSON and Ensembl-116 sequence; complete junction data retained |
| BAM-level S3 extraction | 20_smad4_loxp_evidence.py | Six deposited regional BAM/indexes and filtered original STAR junctions; exact locus sequence/GTF; pysam. Original full-BAM extraction independently verified |
| Supplementary Figure S1 | 22_breeding_scheme.py | Author-confirmed experiment allocation in current Methods; includes RNA-sequencing Figures 4–7 |
| Sox2/Lgr5 supporting rows in Data S3 | 07_figures.py | Six full auxiliary STAR count files; kept separate from Salmon/tximport model inputs |
| Primary program inference; Data S3 program tables | program_gene_set_tests.R, program_model_functions.R, program_results.py | Archived counts/TPM/effective length/W_1; TMM/voom; adjusted ROAST and supplementary models; membership generated from programs.py |
| Upstream model/assembly inspection | Data deposit only | Saved DESeq2 processing/results/parameters and merged transcriptome GTF. No independent raw-count-model refit is supplied or claimed |


`reproduce.py` verifies the 56 source-input hashes in `data/INPUT_MANIFEST.json` before running. The original Western blot image deck and Figure 2 Prism project are additional experimental evidence; Figure 3 Prism is also a runtime input. `data/RAW_EVIDENCE_MANIFEST.json` identifies all three original evidence files.

Numerical validation fixtures are versioned with the code: five Additional-file CSVs, nine source CSVs, two inhibitor CSVs and the reader-facing Data S3 workbook in `results/validation/`. The data archive contains the inputs needed to generate these outputs. The ordered text-gene index selects the Data S3 gene-value rows without retaining manuscript prose or paragraph positions.

The complete workflow passed 35 numerical/model/table checks and 20 figure-file checks. All 15,803 retained genes, four primary and 12 supplementary ROAST rows, 31 AI panel genes, 195 inhibitor measurements and 52 comparisons agree with the benchmarks. All 6,248 reference barcodes, cluster and cell-type labels agree, as do the five locus evidence tables. Validation includes the exact zero/missing probability masks and preserves tiny nonzero P values.

Supplementary Figure S3 contains A–C. The junction observations used previously in panel D remain available in locus data and are represented by the predicted-junction annotation in A. Supplementary Figure S5 combines the primary and downstream analyses into one A–B figure. Neither change alters a scientific calculation.

Primary gene-level DESeq2 results are retained; program limma-voom models and ROAST tests are fitted during reproduction. The raw FASTQ-to-DESeq2 fit is outside the executable package. Figures 1–3 and S2/S4/S7 are author-supplied experimental artwork. Full computational results are recorded in [reproduction-validation.json](reproduction-validation.json).
