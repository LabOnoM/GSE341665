# Panel table consolidation on 3 October 2026

Supplementary Table S4 is the 31-GREEN-gene PanelApp panel 269 v4.37 used by the current analysis. It was formerly Table S5. The retired 20-gene audit table and its export code were removed from active outputs; the historical audit remains in Git/project history. Scientific inputs, panel membership and numerical results are unchanged.

The exporter verifies all human/mouse pairs against the rebuilt panel and writes `Supplementary_Table_S4.md`. Figure 5B rejects any gene-direction input whose mouse membership differs from `ai_genes.py`. Supplementary Methods B6 and Figure S5B use the retained table identifier. Other supplementary identifiers are unchanged.

Validation on 3 October 2026: seven unit tests passed; archived PanelApp/MGI reconstruction exactly reproduced 31 rows and 186 cells; incomplete and wrong-ortholog inputs were rejected; a fresh README `--mode manuscript` run passed 19 numerical table/sheet checks and 12 figure-file checks (PDF/PNG/SVG for Figures 4–7). The unchanged R/Seurat and locus-extraction workflows were not rerun for this scoped change. Prior full-workflow validation remains a separate dated record.
