#!/usr/bin/env python3
"""Reproduce the latest manuscript outputs in an explicit disposable workspace.

This reproduces analyses/figures from archived processed inputs. It does not fit
the original FASTQ-to-DESeq2 model. Large inputs are supplied by the data deposit.
"""
import argparse
import hashlib
import json
import shutil
import subprocess
import sys
from pathlib import Path

REPO = Path(__file__).resolve().parent


def run(command, cwd):
    print('Running:', ' '.join(map(str, command)), flush=True)
    subprocess.run(list(map(str, command)), cwd=cwd, check=True)


def main():
    ap = argparse.ArgumentParser(description=__doc__)
    ap.add_argument('--data-root', type=Path, required=True, help='Unpacked Zenodo data root')
    ap.add_argument('--work-dir', type=Path, required=True, help='New empty working directory')
    ap.add_argument('--destination', type=Path, required=True, help='New output directory')
    ap.add_argument('--mode', choices=['manuscript', 'supplement', 'reference', 'locus'], default='manuscript')
    ap.add_argument('--geo-root', type=Path, help='Optional unpacked GEO processed matrices')
    ap.add_argument('--bam-dir', type=Path, help='BAM, indexes and STAR junction files for locus mode')
    ap.add_argument('--reference-dir', type=Path, help='genome.fa, .fai and genome.gtf for locus mode')
    args = ap.parse_args()
    data = args.data_root.resolve()
    manifest = json.loads((REPO/'data/INPUT_MANIFEST.json').read_text())
    for item in manifest:
        p = data/item['path']
        if not p.is_file(): raise FileNotFoundError(f'Required deposited input missing: {p}')
        if hashlib.sha256(p.read_bytes()).hexdigest() != item['sha256']:
            raise ValueError(f'Deposited input checksum mismatch: {p}')
    work = args.work_dir.resolve(); dest = args.destination.resolve()
    for p in [work, dest]:
        if p.exists() and any(p.iterdir()): raise ValueError(f'Refusing to overwrite nonempty directory: {p}')
        p.mkdir(parents=True, exist_ok=True)
    pipe = work/'scripts/JDR_reanalysis'
    shutil.copytree(REPO/'pipeline', pipe)
    source = data/'01_processed_data/reproduction_inputs'
    shutil.copytree(source/'output', pipe/'output')
    shutil.copytree(source/'data', pipe/'data')
    raw = work/'00.RawData/Bulk-RNAseq';raw.mkdir(parents=True)
    shutil.copytree(source/'workbooks', raw/'02.TableOutput')
    for p in (data/'05_Raw_data_for_Figures/auxiliary_alignment_counts').glob('*.STAR.counts'):
        shutil.copy2(p, raw/p.name)
    shutil.copytree(data/'05_Raw_data_for_Figures/inhibitor_normalized',work/'00.Manuscript/ZXY_20260919_updated')
    supplementary = work/'00.Manuscript/Submission_02_JMM_2026-07/4_Journal_JMM/Supplementary'
    supplementary.mkdir(parents=True)
    shutil.copy2(data/'06_single_cell_reference/GSE146855_seurat.rds',pipe/'output/GSE146855_seurat.rds')
    if args.geo_root:
        # Record and check every retained TPM row against the GEO matrix.
        import numpy as np
        import pandas as pd
        t = pd.read_csv(pipe/'output/TPM_matrix.csv', skiprows=1, float_precision='round_trip')
        g = pd.read_csv(args.geo_root/'GSE_processed_gene_TPM.txt', sep='\t',float_precision='round_trip')
        samples=['CON1','CON2','CON3','dcKO1','dcKO2','dcKO3']
        idcol=g.columns[0]
        expected=t.set_index('ensembl_gene_id')[[f'TPM_{s[:-1]}_{s[-1]}' for s in samples]]
        observed=g.set_index(idcol).loc[expected.index,samples]
        # The deposited GEO export uses %.4f, whereas the source workbooks and
        # figure CSV preserve full precision. Verify the declared conversion.
        rounded=expected.map(lambda value: float(format(value,'.4f')))
        if not np.allclose(rounded.values,observed.values,rtol=0,atol=1e-14):
            raise ValueError('Archived TPM values do not agree with GEO processed matrix')
    py=sys.executable
    def step(name): run([py,pipe/name],work)
    if args.mode=='manuscript':
        # Extract the actual source workbooks, recompute four program scores,
        # NNLS fractions, reference specificity and pinned-panel enrichment.
        for script in ['01_extract_tables.py','03_module_scores.py','06_deconvolution.py','16_ai_gene_enrichment.py']:
            step(script)
        run([py,pipe/'07_figures.py','--root',work,'--destination',dest,
             '--claim-source',source/'manuscript_claim_source.json','--scratch',work/'payload'],work)
        step('24_september_inhibitor_audit.py')
        run([py,pipe/'export_additional_data.py','--workspace',work,'--destination',dest,
             '--payload',work/'payload/workbook_payload.json'],work)
    elif args.mode=='supplement':
        for script in ['14_workflow_diagram.py','21_smad4_loxp_figure.py','22_breeding_scheme.py']:step(script)
        shutil.copytree(pipe/'figures',dest/'Figures')
    elif args.mode=='reference':
        run(['bash',pipe/'04_download_scref.sh'],work)
        for script in ['05_scref_annotate.R','11_scref_umap_dotplot.R']:run(['Rscript',pipe/script],work)
        shutil.copytree(pipe/'figures',dest/'Figures')
    elif args.mode=='locus':
        if args.bam_dir is None or args.reference_dir is None:raise ValueError('Locus mode requires --bam-dir and --reference-dir')
        for p in args.bam_dir.resolve().glob('*'):
            if p.name.endswith(('.bam','.bai','.SJ.out.tab')):(raw/p.name).symlink_to(p)
        for name in ['genome.fa','genome.fa.fai','genome.gtf']:(raw/name).symlink_to(args.reference_dir.resolve()/name)
        step('20_smad4_loxp_evidence.py');step('21_smad4_loxp_figure.py')
        shutil.copytree(pipe/'figures',dest/'Figures')
    (dest/'EXECUTION.json').write_text(json.dumps({'mode':args.mode,'python':sys.version,
        'repository_code_hashes':{str(p.relative_to(REPO)):hashlib.sha256(p.read_bytes()).hexdigest()
          for p in sorted((REPO/'pipeline').rglob('*')) if p.is_file() and '__pycache__' not in p.parts and not p.name.startswith('._')},
        'data_inputs_verified':len(manifest),'raw_count_model_refit':False},indent=2))


if __name__=='__main__': main()
