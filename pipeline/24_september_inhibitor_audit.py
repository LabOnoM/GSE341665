#!/usr/bin/env python3
"""Audit the September 2026 author-supplied normalized inhibitor data.

These are normalized values, not raw Ct values or raw band intensities. Biological
independence is reported in the author's September Methods, not established by a
project file. Read declared data sheets in Figure_3.prism and independently
recalculate ordinary ANOVA / two-sided Dunnett results. Native analyses are
retained in the unchanged project, separately from this recalculation.
"""
from pathlib import Path
import csv
import hashlib
import json
import math
import io
import zipfile

import numpy as np
import scipy
from scipy.integrate import quad
from scipy.special import roots_hermitenorm, ndtr, gammaln
from scipy.stats import f_oneway, t as student_t

ROOT = Path(__file__).resolve().parents[2]
SOURCE = ROOT / '00.Manuscript/prism_projects'
OUT = ROOT / '00.Manuscript/Submission_02_JMM_2026-07/4_Journal_JMM/Supplementary'
PROJECT = 'Figure_3.prism'


def read_prism_tables(path):
    """Read current declared inhibitor sheets, excluding cached/undeclared tables."""
    tables = []
    with zipfile.ZipFile(path) as archive:
        document = json.loads(archive.read('document.json'))
        sheets = document['sheets']['data']
        if len(sheets) != len(set(sheets)):
            raise ValueError('Duplicate declared Prism data sheet')
        for uid in sheets:
            sheet = json.loads(archive.read(f'data/sheets/{uid}/sheet.json'))
            title = sheet['title']
            if not title.startswith(('qPCR-', 'WB-')):
                continue
            prefix, gene = title.split('-', 1)
            table = sheet['table']
            if table['format'] != 'column' or table['dataFormat'] != 'y_single':
                raise ValueError(f'Unsupported inhibitor table format: {title}')
            labels = [json.loads(archive.read(f'data/sets/{x}.json'))['title']
                      for x in table['dataSets']]
            rows = list(csv.reader(io.StringIO(archive.read(
                f"data/tables/{table['uid']}/data.csv").decode('utf-8'))))
            if len(labels) != 5 or labels[0] != 'Control' or len(set(labels)) != 5:
                raise ValueError(f'Expected control and four distinct treatments: {title}')
            if len(rows) != 3 or any(len(row) != 5 for row in rows):
                raise ValueError(f'Expected three samples in each of five groups: {title}')
            values = [np.array([float(row[column]) for row in rows])
                      for column in range(5)]
            if not all(np.isfinite(v).all() for v in values):
                raise ValueError(f'Nonfinite normalized values: {title}')
            tables.append({'assay': 'qPCR' if prefix == 'qPCR' else 'immunoblot',
                           'target': gene, 'labels': labels, 'values': values,
                           'source_table': title})
    keys = [(x['assay'], x['target']) for x in tables]
    if len(keys) != len(set(keys)) or sum(x[0] == 'qPCR' for x in keys) != 8 or sum(x[0] == 'immunoblot' for x in keys) != 5:
        raise ValueError('Expected eight qPCR targets and five proteins in Figure_3.prism')
    return tables


def dunnett_p(t_value, df, k=4):
    """Equal-n two-sided Dunnett probability via deterministic quadrature.

    Conditional on the shared control normal variate and pooled variance, the
    k treatment-control contrasts are independent. Integrating those quantities
    gives the maximum-absolute-t distribution (correlation 1/2). Substitution
    s=abs(t)*sqrt(V/df) avoids loss of very small tail probabilities at large t.
    """
    a = abs(float(t_value))
    if a == 0:
        return 1.0
    z, weights = roots_hermitenorm(100)
    weights = weights / math.sqrt(2 * math.pi)
    logc = math.log(2) + (df / 2) * math.log(df / 2) - gammaln(df / 2) - df * math.log(a)

    def integrand(s):
        if s == 0:
            return 0.0
        inside = np.clip(ndtr(z + math.sqrt(2) * s) - ndtr(z - math.sqrt(2) * s), 0, 1)
        tail = np.sum(weights * (1 - inside**k))
        density = math.exp(logc + (df - 1) * math.log(s) - df * s*s / (2 * a*a))
        return density * tail

    value, error = quad(integrand, 0, 20, epsabs=1e-13, epsrel=2e-9, limit=250)
    assert error < max(1e-10, value * 1e-6)
    return min(1.0, max(0.0, value))


def stars(p):
    return '****' if p < 0.0001 else '***' if p < 0.001 else '**' if p < 0.01 else '*' if p < 0.05 else 'ns'


def main():
    # The same integral with one contrast must reduce to a two-sided t test.
    for tv in (0.1, 1, 2, 4, 10, 25):
        assert np.isclose(dunnett_p(tv, 10, 1), 2 * student_t.sf(tv, 10), rtol=2e-6, atol=1e-12)
    raw, comparisons, summaries = [], [], []
    hashes = {PROJECT: hashlib.sha256((SOURCE/PROJECT).read_bytes()).hexdigest()}
    for table in read_prism_tables(SOURCE/PROJECT):
        assay, gene = table['assay'], table['target']
        labels, values = table['labels'], table['values']
        assert len(values) == 5 and all(len(v) == 3 for v in values)
        assert labels[0] == 'Control'
        anova = f_oneway(*values)
        df_error = sum(map(len, values)) - len(values)
        mse = sum(np.sum((v-v.mean())**2) for v in values) / df_error
        for label, v in zip(labels, values):
            for i, value in enumerate(v, 1):
                raw.append(dict(assay=assay, target=gene, source_group=label, sample_index=i,
                                normalized_value=float(value), source_file=PROJECT, source_table=table['source_table'],
                                replication_status='n=3 biological samples reported by September Methods; sample IDs not supplied'))
            summaries.append(dict(assay=assay,target=gene,group=label,n=len(v),mean=float(v.mean()),sd=float(v.std(ddof=1))))
        for label, v in zip(labels[1:], values[1:]):
            tv = (v.mean()-values[0].mean()) / math.sqrt(mse * (1/len(v)+1/len(values[0])))
            pv = dunnett_p(tv, df_error)
            comparisons.append(dict(assay=assay,target=gene,treatment=label,control='Control',n_treatment=len(v),
                n_control=len(values[0]),mean_treatment=float(v.mean()),sd_treatment=float(v.std(ddof=1)),
                mean_control=float(values[0].mean()),sd_control=float(values[0].std(ddof=1)),
                mean_ratio=float(v.mean()/values[0].mean()),anova_F=float(anova.statistic),anova_df_between=4,
                anova_df_error=df_error,anova_P=float(anova.pvalue),dunnett_t=float(tv),dunnett_adjusted_P=pv,
                significance=stars(pv),source_file=PROJECT,source_table=table['source_table']))
    assert len(raw) == 195 and len(comparisons) == 52
    OUT.mkdir(exist_ok=True)
    for filename, rows in [('Supplementary_Data_S6_inhibitor_normalized_values.csv',raw),
                           ('Supplementary_Data_S7_inhibitor_ANOVA_Dunnett.csv',comparisons)]:
        with (OUT/filename).open('w',newline='') as f:
            w=csv.DictWriter(f,fieldnames=list(rows[0]));w.writeheader();w.writerows(rows)
    audit={'source_sha256':hashes,'scipy_version':scipy.__version__,'method':'ordinary one-way ANOVA; two-sided Dunnett, four comparisons within target, no correction across targets; deterministic quadrature',
        'input_level':'author-supplied normalized values; raw Ct and densitometry not included',
        'biological_n':'three per group according to the September Methods; project files do not establish independence',
        'source_format':'Native Prism archive; only declared data sheets are read. Original native analyses remain in the unchanged project.',
        'summary':summaries,'comparisons':comparisons}
    (OUT/'Inhibitor_Statistics_Provenance.json').write_text(json.dumps(audit,indent=2)+'\n')
    for r in comparisons:
        print(f"{r['assay']:10} {r['target']:6} {r['treatment']:23} ratio={r['mean_ratio']:.3f} Dunnett P={r['dunnett_adjusted_P']:.7g} {r['significance']}")


if __name__ == '__main__':
    main()
