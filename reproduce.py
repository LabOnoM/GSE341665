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
import os
import time
from datetime import datetime, timezone
from importlib.metadata import distributions
from pathlib import Path

REPO = Path(__file__).resolve().parent
STEPS = []
LOG_DIR = None


def run(command, cwd):
    print('Running:', ' '.join(map(str, command)), flush=True)
    start=time.monotonic()
    record={'command':list(map(str,command)),'started_utc':datetime.now(timezone.utc).isoformat()}
    STEPS.append(record)
    logfile=LOG_DIR/f'{len(STEPS):02d}_{Path(command[1]).name}.log'
    with logfile.open('w') as log:
        process=subprocess.Popen(record['command'],cwd=cwd,stdout=subprocess.PIPE,
                                 stderr=subprocess.STDOUT,text=True)
        for line in process.stdout:
            print(line,end='',flush=True);log.write(line);log.flush()
        status=process.wait()
    record.update(exit_code=status,seconds=round(time.monotonic()-start,3),log=str(logfile.name))
    if status:raise subprocess.CalledProcessError(status,record['command'])


def main():
    ap = argparse.ArgumentParser(description=__doc__)
    ap.add_argument('--data-root', type=Path, required=True, help='Unpacked Zenodo data root')
    ap.add_argument('--work-dir', type=Path, required=True, help='New empty working directory')
    ap.add_argument('--destination', type=Path, required=True, help='New output directory')
    ap.add_argument('--mode', choices=['manuscript', 'supplement', 'reference', 'locus', 'panel', 'all', 'reference-figure'], default='manuscript')
    ap.add_argument('--rscript', default='Rscript', help='Rscript executable for the recorded R 4.3.2 environment')
    ap.add_argument('--reference-archive', type=Path, help='Explicit original GSE146855 tar; defaults to deposited archive')
    ap.add_argument('--geo-root', type=Path, help='Optional unpacked GEO processed matrices')
    ap.add_argument('--bam-dir', type=Path, help='BAM, indexes and STAR junction files for locus mode')
    args = ap.parse_args()
    data = args.data_root.resolve()
    work = args.work_dir.resolve(); dest = args.destination.resolve()
    def overlaps(a,b):return a==b or a in b.parents or b in a.parents
    for a,b in [(work,dest),(work,data),(dest,data),(work,REPO),(dest,REPO)]:
        if overlaps(a,b):raise ValueError(f'Input, repository, work and output directories must be separate: {a}; {b}')
    for p in [work,dest]:
        if p.exists() and (not p.is_dir() or any(p.iterdir())):
            raise ValueError(f'Refusing to overwrite nonempty directory: {p}')
    manifest = json.loads((REPO/'data/INPUT_MANIFEST.json').read_text())
    for item in manifest:
        p = data/item['path']
        if data not in p.resolve().parents:raise ValueError(f'Input manifest escapes data root: {p}')
        if not p.is_file(): raise FileNotFoundError(f'Required deposited input missing: {p}')
        with p.open('rb') as f:actual=hashlib.file_digest(f,'sha256').hexdigest()
        if p.stat().st_size != item['bytes'] or actual != item['sha256']:
            raise ValueError(f'Deposited input checksum mismatch: {p}')
    for p in [work, dest]:
        p.mkdir(parents=True, exist_ok=True)
    global LOG_DIR
    LOG_DIR=dest/'logs';LOG_DIR.mkdir()
    for variable in ['OPENBLAS_NUM_THREADS','OMP_NUM_THREADS','MKL_NUM_THREADS','VECLIB_MAXIMUM_THREADS']:
        os.environ[variable]='1'
    os.environ['MPLBACKEND']='Agg'
    os.environ['PYTHON']=sys.executable
    pipe = work/'scripts/JDR_reanalysis'
    shutil.copytree(REPO/'pipeline', pipe,ignore=shutil.ignore_patterns('._*','__pycache__'))
    shutil.copytree(REPO/'environment',pipe/'environment',ignore=shutil.ignore_patterns('._*','__pycache__'))
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
    if args.mode in ('reference','all'):
        run([args.rscript,pipe/'check_reference_environment.R'],work)
        archive=args.reference_archive or data/'06_single_cell_reference/GSE146855_RAW.tar'
        run(['bash',pipe/'04_download_scref.sh','--archive',archive],work)
        for script in ['05_scref_annotate.R','11_scref_umap_dotplot.R']:run([args.rscript,pipe/script],work)
        from validate import check_reference
        check_reference(data,work)
        shutil.copytree(pipe/'figures',dest/'Figures',dirs_exist_ok=True)
    if args.mode in ('panel','all'):
        run([py,pipe/'17_fetch_ai_panel.py','--sources-dir',data/'03_external_panel_sources'],work)
        from validate import compare_csv
        compare_csv(source/'output/AI_panel_PanelApp_v4.37.csv',pipe/'output/AI_panel_PanelApp_v4.37.csv')
    if args.mode in ('manuscript','all'):
        # Extract the actual source workbooks, recompute four program scores,
        # NNLS fractions, reference specificity and pinned-panel enrichment.
        for script in ['01_extract_tables.py','03_module_scores.py','06_deconvolution.py','16_ai_gene_enrichment.py']:
            step(script)
        run([py,pipe/'07_figures.py','--root',work,'--destination',dest,
             '--claim-source',source/'manuscript_claim_source.json','--scratch',work/'payload'],work)
        step('24_september_inhibitor_audit.py')
        run([py,pipe/'export_additional_data.py','--workspace',work,'--destination',dest,
             '--payload',work/'payload/workbook_payload.json'],work)
    if args.mode in ('locus','all'):
        bam_dir=(args.bam_dir or data/'08_locus_alignments').resolve()
        for sample in ['CON1','CON2','CON3','dcKO1','dcKO2','dcKO3']:
            for suffix in ['.STAR.Aligned.sortedByCoord.out.bam','.STAR.Aligned.sortedByCoord.out.bam.bai','.STAR.SJ.out.tab']:
                if not (bam_dir/(sample+suffix)).is_file():raise FileNotFoundError(f'Missing locus input: {bam_dir/(sample+suffix)}')
        for p in bam_dir.glob('*'):
            if p.name.endswith(('.bam','.bai','.SJ.out.tab')):(raw/p.name).symlink_to(p)
        step('20_smad4_loxp_evidence.py');step('21_smad4_loxp_figure.py')
        from validate import check_locus
        check_locus(data,work)
        shutil.copytree(pipe/'figures',dest/'Figures',dirs_exist_ok=True)
    if args.mode in ('supplement','all'):
        scripts=['14_workflow_diagram.py','22_breeding_scheme.py'] if args.mode=='all' else ['14_workflow_diagram.py','21_smad4_loxp_figure.py','22_breeding_scheme.py']
        for script in scripts:step(script)
        shutil.copytree(pipe/'figures',dest/'Figures',dirs_exist_ok=True)
    if args.mode=='reference-figure':
        run([args.rscript,pipe/'check_reference_environment.R'],work)
        run([args.rscript,pipe/'11_scref_umap_dotplot.R'],work)
        shutil.copytree(pipe/'figures',dest/'Figures')
    for name,patterns in [('Reference',['GSE146855*csv','GSE146855*umap.rds','reference*json','reference*txt','scref*csv','scref*txt']),
                          ('Panel',['AI_panel*']),('Locus',['smad4*'])]:
        mode={'Reference':['reference','reference-figure','all'],'Panel':['panel','all'],'Locus':['locus','all']}[name]
        if args.mode in mode:
            target=dest/name;target.mkdir()
            for pattern in patterns:
                for p in (pipe/'output').glob(pattern):shutil.copy2(p,target/p.name)
    run([py,REPO/'validate.py','--data-root',data,'--work-dir',work,'--destination',dest,'--mode',args.mode],work)
    (dest/'EXECUTION.json').write_text(json.dumps({'passed':True,'mode':args.mode,'python':sys.version,
        'steps':STEPS,'data_manifest_sha256':hashlib.sha256((REPO/'data/INPUT_MANIFEST.json').read_bytes()).hexdigest(),
        'python_packages':{d.metadata['Name']:d.version for d in distributions()},
        'environment_hashes':{str(p.relative_to(REPO)):hashlib.sha256(p.read_bytes()).hexdigest()
           for p in sorted((REPO/'environment').rglob('*')) if p.is_file() and not p.name.startswith('._')},
        'threads':1,'entrypoint_sha256':hashlib.sha256(Path(__file__).read_bytes()).hexdigest(),
        'repository_code_hashes':{str(p.relative_to(REPO)):hashlib.sha256(p.read_bytes()).hexdigest()
          for p in sorted((REPO/'pipeline').rglob('*')) if p.is_file() and '__pycache__' not in p.parts and not p.name.startswith('._')},
        'data_inputs_verified':len(manifest),'raw_count_model_refit':False},indent=2))


if __name__=='__main__': main()
