#!/usr/bin/env python3
"""
03_module_scores.py  --  Smad4 dcKO JDR re-analysis pipeline (step 3/7)

Per-sample gene-program "program scores" (combined gene-set z-score: the mean
gene-wise z-score of log2 TPM; Lee et al. 2008, PMID 18989396) for
the four programs, the CON-vs-dcKO group tests, and the sterol<->interferon
reciprocal correlation across samples.

Inputs : output/TPM_matrix.csv
Outputs: output/module_scores.csv, output/module_tests.txt
"""
import os, sys, pandas as pd, numpy as np
sys.path.insert(0, os.path.dirname(os.path.abspath(__file__)))
from scipy import stats

ROOT = os.path.dirname(os.path.dirname(os.path.dirname(os.path.abspath(__file__))))
OUT  = os.path.join(ROOT, "scripts", "JDR_reanalysis", "output")
SAMP = ["TPM_CON_1","TPM_CON_2","TPM_CON_3","TPM_dcKO_1","TPM_dcKO_2","TPM_dcKO_3"]
GROUP = ["CON"]*3 + ["dcKO"]*3

from programs import PROGRAMS as MODULES


def program_score(tpm, genes):
    sub = tpm.reindex([g for g in genes if g in tpm.index])
    lg = np.log2(sub + 1)
    z = lg.sub(lg.mean(axis=1), axis=0).div(lg.std(axis=1).replace(0, np.nan), axis=0)
    return z.mean(axis=0)


def main():
    t = pd.read_csv(os.path.join(OUT, "TPM_matrix.csv"), skiprows=1, low_memory=False)
    for s in SAMP:
        t[s] = pd.to_numeric(t[s], errors="coerce")
    tpm = t.groupby("external_gene_name")[SAMP].mean()

    scores = pd.DataFrame({m: program_score(tpm, g) for m, g in MODULES.items()})
    scores["group"] = GROUP
    scores.to_csv(os.path.join(OUT, "module_scores.csv"))

    lines = ["Per-sample program z-scores:", scores.round(3).to_string(), ""]
    for m in MODULES:
        con = scores.loc[scores.group=="CON", m]; ko = scores.loc[scores.group=="dcKO", m]
        tt, pp = stats.ttest_ind(con, ko)
        lines.append(f"  {m:9s} CON={con.mean():+.2f} dcKO={ko.mean():+.2f} t={tt:+.1f} p={pp:.2g}")
    r, p = stats.pearsonr(scores["Sterol"], scores["IFN"])
    lines.append(f"\nSterol vs IFN across samples: Pearson r={r:.2f} p={p:.2g} (n=6; "
                 f"reflects genotype separation, not within-group coupling)")
    print("\n".join(lines))
    with open(os.path.join(OUT, "module_tests.txt"), "w") as fh:
        fh.write("\n".join(lines) + "\n")


if __name__ == "__main__":
    main()
