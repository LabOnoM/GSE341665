# Current manuscript reproduction map

Authority: BMC Word manuscript dated 2 October 2026 and current Additional files 1–7. `reproduce.py` stages the deposited data in a new working directory, verifies all 70 input/export hashes, and invokes the scripts below. No unused historical analysis is part of this map.

| Figure, table or result | Retained code | Required source and purpose |
|---|---|---|
| Main Figures 4–7; full Data S3 | 01_extract_tables.py, 03_module_scores.py, 06_deconvolution.py, 07_figures.py; volcano_panel.py, vendor/figstyle.py | Two source workbooks; DEG/TPM and captured claim-source JSON; seven-cell-type CPM signature. Preserves original gene-level model statistics and all displayed genes |
| Four scored programs; Additional file 5 | programs.py, 03_module_scores.py | The single definition of 16 identity, 20 sterol, 14 wound/stress and 18 interferon genes; full-precision TPM; exhaustive 20 label assignments |
| Reference fractions/specificity; Additional files 6–7 | 06_deconvolution.py | Captured GSE146855 CPM signature; full-precision bulk TPM; sorted markers before NNLS |
| AI enrichment, Additional file 3 and Methods B6 | 16_ai_gene_enrichment.py, ai_genes.py, ai_gene_provenance.py | Pinned PanelApp 269 v4.37 GREEN panel and gene-level statistics; three reported threshold sensitivity tests |
| Pinned-panel source reconstruction | 17_fetch_ai_panel.py | Original PanelApp v4.37 JSON and MGI report in the deposit; explicit live download also tested. `all` reconstructs the mapped panel before enrichment |
| Zenodo Data S1/S2; Figure 3C/E statistical evidence | 24_september_inhibitor_audit.py | The declared inhibitor data sheets in Figure_3.prism; 195 values and 52 ANOVA/two-sided Dunnett comparisons. The author-provided wet-laboratory artwork is preserved |
| Additional file 2 workbook; files 2–7; Table S4 | export_additional_data.py, ai_gene_provenance.py | Latest regenerated CSV/payload; public openpyxl. Table S4 contains the 31 active PanelApp genes; export checks its human–mouse pairs against the version-pinned analysis panel |
| Supplementary Figure S6 / reference reannotation | 04_download_scref.sh, 05_scref_annotate.R, 11_scref_umap_dotplot.R | Original deposited GSE146855 raw archive; exact R/library lock; fresh labels/signatures checked against captured reference before downstream use |
| Supplementary Figure S5 | 14_workflow_diagram.py | Recorded primary methods in Additional file 1 and current retained-code map; public matplotlib renderer |
| Supplementary Figure S3 rendering | 21_smad4_loxp_figure.py | Archived locus GTF/evidence CSVs/summary JSON and exact Ensembl-116 locus sequence; missing evidence stops execution |
| BAM-level S3 extraction | 20_smad4_loxp_evidence.py | Six deposited regional BAM/indexes and filtered original STAR junctions; exact locus sequence/GTF; pysam. Original full-BAM extraction independently verified |
| Supplementary Figure S1 | 22_breeding_scheme.py | Author-confirmed experiment allocation in current Methods; includes RNA-sequencing Figures 4–7 |
| Sox2/Lgr5 supporting rows in Data S3 | 07_figures.py | Six full auxiliary STAR count files; kept separate from Salmon/tximport model inputs |
| Upstream model/assembly inspection | Data deposit only | Saved DESeq2 processing/results/parameters and merged transcriptome GTF. No independent raw-count-model refit is supplied or claimed |

Main Figures 1–3 and supplementary S2, S4 and S7 are author-supplied experimental artwork, not bioinformatics outputs. No code or normalized values reconstruct the missing complete raw scans, raw Ct values or raw band intensities. Tables S1–S3 transcribe experimental records, with unresolved source checks explicitly retained in Additional file 1.

`results/additional_files/` contains the five current Additional-file CSVs; `results/source_csv/` contains all six Data S3 CSVs. Publication PDFs and the Excel workbook belong to the journal/data package; both Python execution modes produce their stated figure outputs. No duplicate binary publication artwork is needed in this code repository.

The GEO TPM matrix is exported with four decimal places; the source workbooks/plotting inputs preserve full precision. `--geo-root` checks all 15,803 × six retained TPM entries after that declared rounding. Captured Data S3 paragraph numbers refer to the source manuscript used to create its index, not the later BMC numbering. The exact captured source is deposited as `manuscript_claim_source.json`.

## Validation scope

All retained scientific workflows completed on 2 October 2026 and were rerun with required main-figure SVG validation on 3 October 2026 from an independent repository/data snapshot, with no private project helper imports or data downloads. The original GSE146855 count archive was rebuilt under R 4.3.2 and the recorded Seurat/SeuratObject/Matrix/ggplot2/patchwork/uwot versions. All 6,248 ordered barcodes, cluster labels and cell-type labels match exactly; the 21,059 × seven signature values, 17-cluster marker scores and seven cell-type counts match the archived CSVs exactly. The first 20 PCs agree after their arbitrary signs are aligned (maximum difference 9.68e-12). PCA axes can change sign across builds; see the [official R documentation](https://stat.ethz.ch/R-manual/R-devel/library/stats/html/prcomp.html). UMAP/PDF appearance is not checked for binary/pixel identity.

The five locus evidence CSVs match archived evidence exactly when extracted independently from the original full BAMs and from deposited regional BAMs. Every one of the 69,573 reference-window bases was verified against the original full genome FASTA. Live pinned PanelApp/MGI retrieval reproduces the archived 31-row mapped panel exactly. The live original-reference archive matches the original local reference archive by SHA-256.

Final `--mode all` validation passed 28 scientific tables/workbook sheets and 22 figure-file checks (including the four main-figure SVGs), including all 13 published CSV exports and six scientific workbook sheets. The 15,803 gene rows, 157 displayed-gene rows, 102 claim-index rows, 563 strict DEGs, 31 AI panel genes, 195 inhibitor values and 52 comparisons agree. Zero and missing masks agree exactly. Nonzero adjusted P values use zero absolute tolerance and relative tolerance 1e-12; descriptive NNLS t-test P values allow relative tolerance 1e-10. The workbook Read me sheet is exporter documentation. Full run details are in `docs/reproduction-validation.json`; the local preparation record retains commands, logs and hashes.

The saved primary model/assembly remain inspectable, but the missing executable FASTQ-to-DESeq2 fit has not been recovered or rerun. This audit cannot establish missing raw wet-laboratory records. The code repository remains private and the revised data snapshot still requires live Zenodo synchronization. Public accessibility and submission are separate from these verified local computations.

`reproduce.py` orchestrates all script entry points; `validate.py` checks numerical benchmarks. `reference_download.py` is exercised through `04_download_scref.sh`; `check_reference_environment.R` checks the full tested R dependency lock before analysis. `environment/setup_reference.R` restores/checks the declared library. Helper modules `programs.py`, `ai_genes.py`, `ai_gene_provenance.py`, `volcano_panel.py` and `vendor/figstyle.py` are exercised through their owning workflows. Seven scientific-input/probability regression tests are in `tests/test_reproduction.py`; CI runs those tests, without restricted study data.

Figures 4–7 are exported as PDF, SVG and PNG from the same plotting objects. SVG labels remain editable text. Validation requires each SVG to be valid XML with a viewBox and to contain no external resource references.

Figure 4 contains panels A–C (PCA, volcano plot and program scores). The five sex-marker rows remain in Data S3, labelled `Data S3 only` in the Displayed genes sheet/CSV, with empty current `figure_panels` fields in the complete statistics and claim index. Their expression values and statistics are unchanged. See [the reporting revision](sex-marker-reporting.md).
