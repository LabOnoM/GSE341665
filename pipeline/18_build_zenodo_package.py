#!/usr/bin/env python3
"""
18_build_zenodo_package.py — assemble the Zenodo deposition folder.

WHAT GOES TO ZENODO, AND WHAT DOES NOT
--------------------------------------
The three archives have different jobs, and duplicating between them is a mistake:

  * GEO/SRA  — the sequence files and the processed count/TPM matrices. This is where a
               reader expects RNA-seq data, and where the journal will look. (The files
               are the provider's adapter/quality-filtered *_Clean_Data* reads — the
               rawest available; unfiltered output was never returned.)
  * GitHub   — the analysis code, live and versioned (LabOnoM/GSE341665).
  * Zenodo   — a PURE DATA record: the de-novo transcriptome assembly, the saved pipeline
               objects, and the exploratory analyses (circRNA, isoform switching,
               differential exon usage) that the manuscript does not interpret. It holds
               NO code. An earlier version copied the pipeline in here, which duplicated
               all 24 scripts byte-for-byte against GitHub and created two things to keep
               in step. The citable code snapshot is produced instead by Zenodo's GitHub
               integration: publishing a release makes Zenodo archive that tag and mint a
               SEPARATE DOI for it — two records, two DOIs, no manual copying. The
               integration only sees PUBLIC repositories, so the code DOI can be minted
               once the repository is made public at publication.

Deliberately NOT deposited, with reasons:
  * FASTQ (15.4 GB)              -> belongs in GEO/SRA; depositing twice is redundant.
                                    (These are the provider's adapter/quality-filtered
                                    *_Clean_Data* files - the rawest available.)
  * Data.tar (16 GB)             -> byte-identical duplicate of those FASTQ; deposit neither
                                    it nor the *_val_*.fq.gz post-trimming output.
  * STAR genome index (24.4 GB)  -> regenerable from the public reference in one command.
  * Ensembl genome.fa / genome.gtf (9.2 GB) -> public at Ensembl release 116; citing the
                                    release is more useful than re-hosting it.
  * Pfam-A.hmm (2.1 GB)          -> public database, versioned upstream.
  * STAR BAM (12.8 GB)           -> derivable from FASTQ + reference. Staged only if
                                    --with-bam is passed; see the note in the manifest.

That leaves a record of roughly 2 GB, comfortably inside Zenodo's 50 GB per-record limit,
containing only what cannot be reconstructed from public sources.

Usage
-----
    python3 scripts/JDR_reanalysis/18_build_zenodo_package.py            # build + stage small files
    python3 scripts/JDR_reanalysis/18_build_zenodo_package.py --stage-large   # also copy the big files
    python3 scripts/JDR_reanalysis/18_build_zenodo_package.py --with-bam --stage-large

Large files are NOT copied by default: they are listed in MANIFEST.md with their source
paths and checksums so nothing is silently duplicated on disk.
"""

from __future__ import annotations

import argparse
import hashlib
import json
import os
import shutil
import sys

ROOT = os.path.dirname(os.path.dirname(os.path.dirname(os.path.abspath(__file__))))
RAW = os.path.join(ROOT, "00.RawData", "Bulk-RNAseq")
OUT_PIPE = os.path.join(ROOT, "scripts", "JDR_reanalysis", "output")
SUP = os.path.join(ROOT, "00.Manuscript", "Submission_02_JMM_2026-07", "4_Journal_JMM", "Supplementary")
DEST = os.path.join(ROOT, "00.Manuscript", "Submission_02_JMM_2026-07", "3_Zenodo")

SMALL_LIMIT = 60 * 1024 * 1024      # copy anything under this automatically

# (destination subfolder, source path, one-line description)
LARGE_ITEMS = [
    ("02_denovo_assembly", os.path.join(RAW, "00.Downloadable_Content",
                                        "DeNovoAssemblyTranscriptome.gtf.gz"),
     "StringTie-merged, SQANTI3-annotated de-novo transcriptome. Not reproducible without "
     "this exact pipeline run; every quantification depends on it."),
    ("04_pipeline_objects", os.path.join(RAW, "00.Results", "DESeq2_Processing.Rdata"),
     "DESeqDataSet, RUVSeq covariate W_1 per sample, and the log2 TPM matrix."),
    ("04_pipeline_objects", os.path.join(RAW, "00.Results", "DESeq2_Results.Rdata"),
     "DESeq2 likelihood-ratio-test results for all 15,803 tested genes."),
    ("04_pipeline_objects", os.path.join(RAW, "00.Results", "Parameters.Rdata"),
     "Recorded pipeline parameters, including Par_ReN, Par_SaN and Par_UseFDR4FEA."),
    ("03_exploratory_analyses", os.path.join(RAW, "00.Results", "circRNA_Results.Rdata"),
     "circRNAprofiler output: 184 circRNAs detected, none differential. Not interpreted "
     "in the manuscript."),
    ("03_exploratory_analyses", os.path.join(RAW, "00.Results", "aSwitchList.Rdata"),
     "IsoformSwitchAnalyzeR switch list with functional consequences. Not interpreted "
     "in the manuscript."),
    ("03_exploratory_analyses", os.path.join(RAW, "00.Results", "DEXseqResOut.Rdata"),
     "DEXSeq differential exon usage. Not interpreted in the manuscript."),
]
BAM_ITEM = ("05_alignments", os.path.join(RAW, "00.Downloadable_Content", "BAM_Files.zip"),
            "STAR two-pass alignments for all six libraries. Derivable from the FASTQ in "
            "GEO plus the public reference; included only for convenience.")


def md5(path: str, chunk: int = 1 << 20) -> str:
    h = hashlib.md5()
    with open(path, "rb") as fh:
        for block in iter(lambda: fh.read(chunk), b""):
            h.update(block)
    return h.hexdigest()


def human(n: int) -> str:
    for unit in ("B", "KB", "MB", "GB"):
        if n < 1024 or unit == "GB":
            return f"{n:.1f} {unit}" if unit != "B" else f"{n} B"
        n /= 1024.0
    return f"{n:.1f} GB"


def main() -> None:
    ap = argparse.ArgumentParser()
    ap.add_argument("--stage-large", action="store_true",
                    help="copy the large files in as well (default: manifest only)")
    ap.add_argument("--with-bam", action="store_true",
                    help="also include the 12.8 GB STAR BAM archive")
    args = ap.parse_args()

    items = list(LARGE_ITEMS) + ([BAM_ITEM] if args.with_bam else [])

    for sub in ("01_processed_data", "02_denovo_assembly",
                "03_exploratory_analyses", "04_pipeline_objects"):
        os.makedirs(os.path.join(DEST, sub), exist_ok=True)
    if args.with_bam:
        os.makedirs(os.path.join(DEST, "05_alignments"), exist_ok=True)

    staged, pending = [], []

    # ---- 02 processed data ------------------------------------------------
    proc_dst = os.path.join(DEST, "01_processed_data")
    n_proc = 0
    for src_dir, names in ((SUP, [f for f in os.listdir(SUP) if f.startswith("Supplementary_Data_")]),
                           (OUT_PIPE, ["AI_panel_PanelApp_v4.37.csv",
                                       "AI_enrichment_sensitivity.csv",
                                       "AI_gene_panel_status.csv",
                                       "module_scores.csv", "module_tests.txt",
                                       "deconv_proportions.csv",
                                       "program_celltype_specificity.csv"])):
        for f in names:
            p = os.path.join(src_dir, f)
            if os.path.exists(p):
                shutil.copy2(p, os.path.join(proc_dst, f))
                n_proc += 1
    print(f"  01_processed_data  : {n_proc} tables")

    # ---- large items ------------------------------------------------------
    for sub, src, desc in items:
        if not os.path.exists(src):
            pending.append((sub, src, desc, None, None, "SOURCE NOT FOUND"))
            continue
        size = os.path.getsize(src)
        if args.stage_large or size <= SMALL_LIMIT:
            dst = os.path.join(DEST, sub, os.path.basename(src))
            if not os.path.exists(dst) or os.path.getsize(dst) != size:
                shutil.copy2(src, dst)
            staged.append((sub, src, desc, size, md5(dst), "staged"))
        else:
            pending.append((sub, src, desc, size, None, "NOT copied — run with --stage-large"))

    # ---- manifest ---------------------------------------------------------
    lines = ["# Zenodo deposition manifest", "",
             "Generated by `scripts/JDR_reanalysis/18_build_zenodo_package.py`. "
             "Checksums are MD5 of the staged copy.", ""]
    if staged:
        lines += ["## Staged in this folder", "",
                  "| Folder | File | Size | MD5 | Description |", "|---|---|---|---|---|"]
        for sub, src, desc, size, h, _ in staged:
            lines.append(f"| `{sub}` | `{os.path.basename(src)}` | {human(size)} | `{h}` | {desc} |")
        lines.append("")
    if pending:
        lines += ["## Still to be added before upload", "",
                  "These are large and were not copied automatically. Either re-run with "
                  "`--stage-large`, or copy them manually from the paths below.", "",
                  "| Folder | Source path | Size | Status | Description |", "|---|---|---|---|---|"]
        for sub, src, desc, size, _, status in pending:
            lines.append(f"| `{sub}` | `{os.path.relpath(src, ROOT)}` | "
                         f"{human(size) if size else '—'} | {status} | {desc} |")
        lines.append("")
    lines += ["## Deliberately NOT deposited", "",
              "| Item | Size | Why not |", "|---|---|---|",
              "| Raw FASTQ | 15.4 GB | Deposited in GEO/SRA; duplicating it here would split the citation. |",
              "| STAR genome index | 24.4 GB | Regenerable from the public reference in one command. |",
              "| `genome.fa`, `genome.gtf` | 9.2 GB | Public at Ensembl release 116 (GRCm39); cite the release. |",
              "| `Pfam-A.hmm` | 2.1 GB | Public, versioned upstream. |",
              "| STAR BAM archive | 12.8 GB | Derivable from FASTQ + reference. Pass `--with-bam` to include. |",
              ""]
    open(os.path.join(DEST, "MANIFEST.md"), "w").write("\n".join(lines) + "\n")

    print(f"  staged             : {len(staged)}   pending: {len(pending)}")
    print(f"  wrote              : {os.path.relpath(os.path.join(DEST,'MANIFEST.md'), ROOT)}")
    total = sum(s[3] for s in staged if s[3])
    print(f"  package size so far: {human(total)}")


if __name__ == "__main__":
    main()
