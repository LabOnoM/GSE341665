# Known issues and caveats

Read this before running the pipeline or interpreting the tables. Each entry is a real
problem that was found and corrected during pre-submission auditing; they are recorded
here rather than quietly fixed, because they change how you should read the outputs.

## 1. `DEG_table.csv` must be read with `float_precision="round_trip"`

Very small adjusted *P* values are written in **full decimal notation** — *Amelx* has 91
leading zeros — and pandas' default float parser silently flushes them to `0.0`.

```python
pd.read_csv(path, skiprows=1, low_memory=False, float_precision="round_trip")
```

With the default parser, 420 genes appear to have `padj == 0` when only **5** genuinely
do, and the −log10 *P* distribution is fabricated.

**What this did and did not affect.** Gene counts and the enrichment were **unaffected**:
every threshold in this study is `padj < 0.05`, and 0 < 0.05, so mis-parsed genes stayed
on the correct side of every cutoff. It corrupted *P*-value magnitudes and anything
plotted from them.

## 2. Programme scores are effect sizes, not significance tests

Each gene is *z*-scored across the same six samples before averaging, which collapses
within-group variance. A parametric *t* test on the result returns *P* values (10⁻⁶–10⁻⁸)
that measure the arithmetic of the scoring procedure, not biological replication.

With n = 3 vs 3 the **exact permutation floor is *P* = 0.10**, and all four programmes
attain it. Statistical support for the identity claim rests on the gene-level DESeq2
statistics, not on the programme scores.

## 3. The amelogenesis-imperfecta panel is external and version-pinned

Genes come from Genomics England PanelApp panel 269, **version 4.37**, GREEN-rated only.
`pipeline/17_fetch_ai_panel.py` pins that version and **fails if PanelApp has moved on**,
rather than silently accepting a different gene set. If it fails, review the diff before
updating the pin.

A hand-curated predecessor was retired after an audit found 18 of its 20 PMIDs cited
unrelated papers. Verified provenance for all 31 current genes is in
`results/supplementary_data/AI_gene_provenance_TableS1_S2.md`.

## 4. Two Cre drivers

The sequenced material is **K14-creER** (tamoxifen-inducible). The in vivo enamel
phenotype in Figure 1 uses the constitutive **K14-cre** line. Do not assume one genotype
throughout.

## 5. The pipeline fails loudly by design

External resources (PanelApp, MGI, GEO GSE146855) are fetched at run time and there is no
stored fallback. A hard-coded fallback table is what produced a fabricated figure earlier
in this project's history, so a failed fetch now stops the run instead.
