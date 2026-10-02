#!/usr/bin/env python3
"""
06_deconvolution.py  --  Smad4 dcKO JDR re-analysis pipeline (step 6/7)

Estimate cell-type composition of the 6 bulk RNA-seq samples against the
annotated single-cell signature (non-negative least squares over a cell-type-
specific marker panel), and localize each DEG program to cell types.

Purpose: test whether the bulk changes are cell-INTRINSIC (ameloblast
reprogramming) or COMPOSITIONAL (shift in cell-type proportions).

Inputs : output/GSE146855_signature_cpm.csv, output/TPM_matrix.csv
Outputs: output/deconv_proportions.csv, output/program_celltype_specificity.csv,
         output/deconv_tests.txt

Caveats (stated in the manuscript): single P7-incisor reference, de novo
annotation, TPM(bulk)/CPM(reference) unit mismatch -> treat proportions as
relative; the CON-vs-dcKO comparison is more robust than absolute values.
"""
import os, sys, pandas as pd, numpy as np
sys.path.insert(0, os.path.dirname(os.path.abspath(__file__)))
from scipy.optimize import nnls
from scipy import stats

ROOT = os.path.dirname(os.path.dirname(os.path.dirname(os.path.abspath(__file__))))
OUT  = os.path.join(ROOT, "scripts", "JDR_reanalysis", "output")
SAMP = ["TPM_CON_1","TPM_CON_2","TPM_CON_3","TPM_dcKO_1","TPM_dcKO_2","TPM_dcKO_3"]
GROUP = pd.Series(["CON"]*3 + ["dcKO"]*3, index=SAMP)

from programs import PROGRAMS


def main():
    sig = pd.read_csv(os.path.join(OUT, "GSE146855_signature_cpm.csv"), index_col=0, float_precision="round_trip")
    t = pd.read_csv(os.path.join(OUT, "TPM_matrix.csv"), skiprows=1, low_memory=False, float_precision="round_trip")
    for s in SAMP:
        t[s] = pd.to_numeric(t[s], errors="coerce")
    bulk = t.dropna(subset=["external_gene_name"]).groupby("external_gene_name")[SAMP].mean()

    # cell-type-specific marker panel from the signature
    S = sig[sig.sum(axis=1) > 0]
    spec = S.div(S.mean(axis=1) + 1e-9, axis=0)
    markers = set()
    for ct in S.columns:
        cand = [g for g in spec[ct].sort_values(ascending=False).index
                if S.loc[g, ct] > 20 and spec.loc[g, ct] > 2][:120]
        markers.update(cand)
    markers = sorted(g for g in markers if g in bulk.index)

    Sm = S.reindex(markers)[S.columns].values
    Bm = bulk.reindex(markers)[SAMP].values
    props = []
    for j in range(Bm.shape[1]):
        x, _ = nnls(Sm, Bm[:, j]); x = x / x.sum() if x.sum() > 0 else x
        props.append(x)
    P = pd.DataFrame(props, index=SAMP, columns=S.columns)
    P.to_csv(os.path.join(OUT, "deconv_proportions.csv"))

    lines = [f"NNLS deconvolution over {len(markers)} cell-type-specific marker genes.",
             "\nEstimated cell-type proportions (%):", (P*100).round(1).to_string(),
             "\nCON vs dcKO (t-test):"]
    for ct in S.columns:
        con = P.loc[GROUP=="CON", ct]*100; ko = P.loc[GROUP=="dcKO", ct]*100
        tt, pp = stats.ttest_ind(con, ko)
        lines.append(f"  {ct:22s} CON={con.mean():5.1f}%  dcKO={ko.mean():5.1f}%  "
                     f"delta={ko.mean()-con.mean():+5.1f}  p={pp:.3f}")

    # program cell-type specificity (z across cell types) -> where each program lives
    Z = sig.T.apply(lambda c: (c - c.mean())/(c.std()+1e-9)).T
    rows = {}
    for name, gl in PROGRAMS.items():
        gl = [g for g in gl if g in Z.index]
        rows[name] = Z.reindex(gl).mean(axis=0)
    spec_tab = pd.DataFrame(rows).T
    spec_tab.to_csv(os.path.join(OUT, "program_celltype_specificity.csv"))
    lines += ["\nGene-program cell-type specificity (mean z across program genes; +=enriched):",
              spec_tab.round(2).to_string(),
              "\nInterpretation: identity/secretory -> ameloblast (intrinsic; ameloblast "
              "fraction preserved). interferon -> immune/endothelial (its decrease may be "
              "cellular origin UNDETERMINED — immune/endothelial fit to 0.0% in BOTH\n"
              "  groups, so this model cannot attribute the interferon change)."]
    print("\n".join(lines))
    with open(os.path.join(OUT, "deconv_tests.txt"), "w") as fh:
        fh.write("\n".join(lines) + "\n")


if __name__ == "__main__":
    main()
