# FAQ

**Can I reproduce the figures without downloading the large data?**
Partly. Figures 4–6 need `DEG_table.csv` and `TPM_matrix.csv`, which are produced by
step 01 from the GEO matrices. The deconvolution panels additionally need the annotated
single-cell reference (step 05, from the public GSE146855).

**Why does the pipeline stop with an error about PanelApp?**
By design. Step 17 pins panel version 4.37; if PanelApp has published a newer version the
gene set would change silently, so the script stops. Review the difference, then update
the pin and the Methods together.

**Why is `padj` zero for some genes in my copy of the table?**
You read it with the default float parser. See `known-issues.md` §1.

**Why is the interferon programme's cellular origin "undetermined"?**
The deconvolution estimates immune and endothelial content at 0.0% in *both* genotypes for
this dissected material, so it can neither support nor exclude a compositional explanation.
An earlier draft claimed a compositional origin; that claim was not supported by our own
deconvolution and was withdrawn.

**Is the cholesterol/SREBP2 finding causal?**
No. It is correlative. A statin or *Srebf2* perturbation with free-cholesterol rescue would
be needed to test causality, and no further wet-laboratory experiments are planned for this
study.

**Can I get the raw, unfiltered reads?**
They do not exist in our hands — the sequencing provider returned adapter- and
quality-filtered files only. See `data-access.md`.
