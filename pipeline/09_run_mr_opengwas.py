#!/usr/bin/env python3
"""
09_run_mr_opengwas.py -- two-sample MR: serum cholesterol -> dental outcomes.

Self-contained (no TwoSampleMR needed): queries the OpenGWAS REST API directly
for LD-clumped instruments and outcome associations, harmonizes, and computes
IVW (random-effects), MR-Egger (slope + intercept pleiotropy test), and
weighted-median estimates.

Token: read from ~/.opengwas_token (chmod 600) or $OPENGWAS_JWT. Never printed.
Scope caveats (for the manuscript): exposure = SERUM cholesterol (not ameloblast
sterol synthesis); outcomes = MATURE-tooth disease (no enamel-formation GWAS
exists). Supportive plausibility only.

Outputs: output/MR_results.csv, output/MR_per_snp/<exp>__<out>.csv
"""
import os, json, time, urllib.request, urllib.error, numpy as np, pandas as pd
from scipy import stats

ROOT = os.path.dirname(os.path.dirname(os.path.dirname(os.path.abspath(__file__))))
OUT = os.path.join(ROOT, "scripts", "JDR_reanalysis", "output")
os.makedirs(os.path.join(OUT, "MR_per_snp"), exist_ok=True)
TOKEN = os.environ.get("OPENGWAS_JWT") or open(os.path.expanduser("~/.opengwas_token")).read().strip()
BASE = "https://api.opengwas.io/api"
HDR = {"Authorization": f"Bearer {TOKEN}", "Content-Type": "application/json", "User-Agent": "research-mr"}

EXPOSURES = {"LDL_cholesterol": "ieu-a-300", "total_cholesterol": "ieu-a-301", "HDL_cholesterol": "ieu-a-299"}
# Developmental dental outcomes preferred (enamel formation biology). Note: no
# well-powered GWAS of amelogenesis imperfecta / MIH / enamel defects exists; the
# only developmental dental outcome in OpenGWAS is hypodontia (tooth agenesis,
# 250 cases -> severely underpowered). Caries retained as a mature-disease contrast.
OUTCOMES = {"hypodontia_FinnGen": "finn-b-K11_HYPODONTIA",   # developmental (tooth agenesis), underpowered
            "dental_caries_FinnGen": "finn-b-K11_CARIES"}      # mature-disease contrast
PALINDROMIC = {frozenset("AT"), frozenset("CG")}


def api(ep, body):
    for attempt in range(4):
        try:
            req = urllib.request.Request(f"{BASE}/{ep}", data=json.dumps(body).encode(), headers=HDR, method="POST")
            with urllib.request.urlopen(req, timeout=120) as r:
                return json.loads(r.read())
        except urllib.error.HTTPError as e:
            if e.code == 429:  # throttled
                time.sleep(5 * (attempt + 1)); continue
            raise
    raise RuntimeError("API retries exhausted")


def get_associations(rsids, oid, chunk=60):
    # OpenGWAS caps N(id)*N(variant) <= 64 per call, so chunk the instrument list
    out = []
    for i in range(0, len(rsids), chunk):
        out += api("associations", {"variant": rsids[i:i + chunk], "id": [oid]})
    return out


def harmonize(inst, outc):
    o = {d["rsid"]: d for d in outc}
    rows = []
    for x in inst:
        y = o.get(x["rsid"])
        if not y:
            continue
        ax, nx = x["ea"].upper(), x["nea"].upper()
        ay, ny = y["ea"].upper(), y["nea"].upper()
        if frozenset(ax + nx) in PALINDROMIC:          # drop palindromic (strand-ambiguous)
            continue
        if {ax, nx} != {ay, ny}:                        # allele set mismatch -> drop
            continue
        by = float(y["beta"]) * (1 if ay == ax else -1)  # align outcome to exposure effect allele
        bx = float(x["beta"])
        if bx < 0:                                        # orient all exposure effects positive
            bx, by = -bx, -by
        rows.append(dict(rsid=x["rsid"], bx=bx, sex=float(x["se"]),
                         by=by, sey=float(y["se"])))
    return pd.DataFrame(rows)


def mr(df):
    bx, by, sey = df.bx.values, df.by.values, df.sey.values
    w = 1.0 / sey**2
    # IVW through origin (random-effects)
    b_ivw = np.sum(w * bx * by) / np.sum(w * bx**2)
    Q = np.sum(w * (by - b_ivw * bx)**2); dfree = len(df) - 1
    se_ivw = np.sqrt(1 / np.sum(w * bx**2)) * np.sqrt(max(1.0, Q / dfree))
    p_ivw = 2 * stats.norm.sf(abs(b_ivw / se_ivw))
    # MR-Egger (WLS of by on bx with intercept)
    X = np.column_stack([np.ones_like(bx), bx]); W = np.diag(w)
    XtWX = X.T @ W @ X; beta = np.linalg.solve(XtWX, X.T @ W @ by)
    resid = by - X @ beta; sigma2 = (resid @ (W @ resid)) / (len(df) - 2)
    cov = sigma2 * np.linalg.inv(XtWX); se = np.sqrt(np.diag(cov))
    egger_int, egger_slope = beta[0], beta[1]
    p_int = 2 * stats.t.sf(abs(egger_int / se[0]), len(df) - 2)
    p_slope = 2 * stats.t.sf(abs(egger_slope / se[1]), len(df) - 2)
    # weighted median (bootstrap SE)
    ratios = by / bx; wr = bx**2 / sey**2

    def wmed(r, ww):
        idx = np.argsort(r); r, ww = r[idx], ww[idx]
        cum = np.cumsum(ww) - 0.5 * ww; cum /= np.sum(ww)
        return np.interp(0.5, cum, r)
    b_wm = wmed(ratios, wr)
    rng = np.random.default_rng(1)
    boot = [wmed(rng.normal(by, sey) / bx, wr) for _ in range(1000)]
    se_wm = np.std(boot); p_wm = 2 * stats.norm.sf(abs(b_wm / se_wm))
    return dict(n_snp=len(df),
                IVW_beta=b_ivw, IVW_se=se_ivw, IVW_p=p_ivw,
                Egger_beta=egger_slope, Egger_se=se[1], Egger_p=p_slope,
                Egger_intercept=egger_int, Egger_intercept_p=p_int,
                WeightedMedian_beta=b_wm, WeightedMedian_se=se_wm, WeightedMedian_p=p_wm,
                Q_het=Q, Q_p=stats.chi2.sf(Q, dfree))


def main():
    results = []
    for en, eid in EXPOSURES.items():
        inst = api("tophits", {"id": [eid], "pval": 5e-8, "clump": 1, "r2": 0.001, "kb": 10000, "pop": "EUR"})
        rsids = [d["rsid"] for d in inst]
        for on, oid in OUTCOMES.items():
            outc = get_associations(rsids, oid)
            h = harmonize(inst, outc)
            if len(h) < 3:
                print(f"  {en} -> {on}: only {len(h)} usable SNPs, skipped"); continue
            h.to_csv(os.path.join(OUT, "MR_per_snp", f"{en}__{on}.csv"), index=False)
            r = mr(h); r.update(exposure=en, outcome=on)
            results.append(r)
            print(f"  {en:18s} -> {on:22s} nSNP={r['n_snp']:2d}  "
                  f"IVW beta={r['IVW_beta']:+.3f} (p={r['IVW_p']:.3f})  "
                  f"Egger_int p={r['Egger_intercept_p']:.2f}")
    df = pd.DataFrame(results)
    cols = ["exposure", "outcome", "n_snp", "IVW_beta", "IVW_se", "IVW_p",
            "WeightedMedian_beta", "WeightedMedian_p", "Egger_beta", "Egger_p",
            "Egger_intercept", "Egger_intercept_p", "Q_het", "Q_p"]
    df = df[cols]
    df.to_csv(os.path.join(OUT, "MR_results.csv"), index=False)
    print("\n=== MR results (beta = per-SD-higher cholesterol effect on log-odds/liability of outcome) ===")
    print(df.round(3).to_string(index=False))
    print("\nSaved -> output/MR_results.csv")


if __name__ == "__main__":
    main()
