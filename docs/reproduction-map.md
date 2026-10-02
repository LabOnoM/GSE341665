# Current manuscript reproduction map

Authority: BMC Word manuscript dated 2 October 2026 and current Additional files 1–9. `reproduce.py` stages the deposited data in a new working directory, verifies all 48 input/export hashes, and invokes the scripts below. No unused historical analysis is part of this map.

| Figure, table or result | Retained code | Required source and purpose |
|---|---|---|
| Main Figures 4–7; full Data S8 | 01_extract_tables.py, 03_module_scores.py, 06_deconvolution.py, 07_figures.py; volcano_panel.py, vendor/figstyle.py | Two source workbooks; DEG/TPM and captured claim-source JSON; seven-cell-type CPM signature. Preserves original gene-level model statistics and all displayed genes |
| Four scored programs; Additional file 7 | programs.py, 03_module_scores.py | The single definition of 16 identity, 20 sterol, 14 wound/stress and 18 interferon genes; full-precision TPM; exhaustive 20 label assignments |
| Reference fractions/specificity; Additional files 8–9 | 06_deconvolution.py | Captured GSE146855 CPM signature; full-precision bulk TPM; sorted markers before NNLS |
| AI enrichment, Additional file 5 and Methods B6 | 16_ai_gene_enrichment.py, ai_genes.py, ai_gene_provenance.py | Pinned PanelApp 269 v4.37 GREEN panel and gene-level statistics; three reported threshold sensitivity tests |
| Optional pinned-panel retrieval | 17_fetch_ai_panel.py | PanelApp v4.37 and MGI homology records; network access. Main mode uses the captured panel and does not refresh it |
| Additional files 2–3; Figure 3C/E statistical evidence | 24_september_inhibitor_audit.py | The two September normalized PZFX files; 195 values and 52 ANOVA/two-sided Dunnett comparisons. The author-provided wet-laboratory artwork is preserved |
| Additional file 4 workbook; files 2–9; Tables S1–S2 | export_additional_data.py, ai_gene_provenance.py | Latest regenerated CSV/payload; public openpyxl. Table S1's retired-list provenance is retained because it is explicitly part of current Additional file 1 |
| Supplementary Figure S3 / reference reannotation | 04_download_scref.sh, 05_scref_annotate.R, 11_scref_umap_dotplot.R | GEO GSE146855 raw public reference for optional R/Seurat rerun; annotated Seurat object retained for reference inspection |
| Supplementary Figure S4 | 14_workflow_diagram.py | Recorded primary methods in Additional file 1 and current retained-code map; public matplotlib renderer |
| Supplementary Figure S5 rendering | 21_smad4_loxp_figure.py | Archived locus GTF/evidence CSVs/summary JSON and exact Ensembl-116 locus sequence; missing evidence stops execution |
| Optional BAM-level S5 extraction | 20_smad4_loxp_evidence.py | Original six BAM/indexes, STAR junctions and exact genome reference via explicit CLI paths; samtools required |
| Supplementary Figure S6 | 22_breeding_scheme.py | Author-confirmed experiment allocation in current Methods; includes RNA-sequencing Figures 4–7 |
| Sox2/Lgr5 supporting rows in Data S8 | 07_figures.py | Six full auxiliary STAR count files; kept separate from Salmon/tximport model inputs |
| Upstream model/assembly inspection | Data deposit only | Saved DESeq2 processing/results/parameters and merged transcriptome GTF. No independent raw-count-model refit is supplied or claimed |

Main Figures 1–3 and supplementary S1, S2 and S7 are author-supplied experimental artwork, not bioinformatics outputs. No code or normalized values reconstruct the missing complete raw scans, raw Ct values or raw band intensities. Tables S3–S5 transcribe experimental records, with unresolved source checks explicitly retained in Additional file 1.

`results/additional_files/` contains the seven current Additional-file CSVs; `results/source_csv/` contains all six Data S8 CSVs. Publication PDFs and the Excel workbook belong to the journal/data package; both Python execution modes produce their stated figure outputs. No duplicate binary publication artwork is needed in this code repository.

The GEO TPM matrix is exported with four decimal places; the source workbooks/plotting inputs preserve full precision. `--geo-root` checks all 15,803 × six retained TPM entries after that declared rounding. Captured Data S8 paragraph numbers refer to the source manuscript used to create its index, not the later BMC numbering. The exact captured source is deposited as `manuscript_claim_source.json`.

## Validation scope

Main and supplement modes completed successfully on 2 October 2026, using the pinned Python environment. All 13 CSVs and six scientific workbook sheets were checked against the current submission exports: 15,803 gene rows, 157 displayed-gene rows, 102 claim-index rows, 563 strict DEGs, 31 AI panel genes, 195 inhibitor values and 52 comparisons. Zero and missing masks agree exactly. Nonzero adjusted P values use zero absolute tolerance and relative tolerance 1e-12. Descriptive NNLS t-test P values agree within 1e-10 relative tolerance, reflecting insignificant floating-point changes from deterministic marker ordering. The workbook's Read me sheet is exporter documentation, not a scientific data sheet.

Figure S3's full R/Seurat rerun and S5's BAM-level re-extraction were not performed in this check. Original primary-fit tool versions/commands that were not recorded remain unknown. The current executable package starts from archived results and does not perform a full FASTQ-to-DESeq2 refit.

The code repository is private. The local cleaned Zenodo data/metadata still require live-record synchronization, using a new version if the record is already published. HTTP 404 for unauthenticated private GitHub, a GEO browser challenge and a Zenodo 403 do not establish public data access. The interactive report is historical context and includes excluded exploratory analyses; it is not the versioned archive.
