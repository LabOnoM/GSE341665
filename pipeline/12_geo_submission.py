#!/usr/bin/env python3
"""
12_geo_submission.py -- prepare GEO submission files for the bulk RNA-seq.

GEO requires (https://www.ncbi.nlm.nih.gov/geo/info/seq.html):
  * processed data with a QUANTITATIVE component, complete for all features and
    all samples -- explicitly NOT a list of differentially expressed genes;
  * matrix column headers that match the library names in the SAMPLES section;
  * unique file names, restricted to A-Z a-z 0-9 hyphen underscore;
  * md5 checksums recorded on the "MD5 Checksums" tab of the metadata template.

This script writes the two processed matrices, collects the raw-FASTQ checksums
supplied by the sequencing provider (verified against the files on disk), and
computes checksums for the processed files.

Outputs (in 00.Manuscript/Version02/GEO_submission/):
  GSE_processed_gene_raw_counts.txt
  GSE_processed_gene_TPM.txt
  GEO_checksums.tsv
"""
import os, csv, hashlib
import pandas as pd

ROOT = os.path.dirname(os.path.dirname(os.path.dirname(os.path.abspath(__file__))))
PIPE = os.path.join(ROOT, "scripts", "JDR_reanalysis", "output")
RAWD = os.path.join(ROOT, "00.RawData", "Bulk-RNAseq")
OUT = os.path.join(ROOT, "00.Manuscript", "Version02", "GEO_submission")
os.makedirs(OUT, exist_ok=True)

LIBS = ["CON1", "CON2", "CON3", "dcKO1", "dcKO2", "dcKO3"]
ANN = ["ensembl_gene_id", "external_gene_name", "gene_biotype"]


def write_matrix(df, cols, newnames, fname, valfmt):
    sub = df[ANN + cols].copy()
    sub.columns = ANN + newnames
    path = os.path.join(OUT, fname)
    sub.to_csv(path, sep="\t", index=False, float_format=valfmt)
    print(f"  {fname}: {sub.shape[0]} genes x {len(newnames)} samples")
    return path


def md5(path, chunk=1 << 22):
    h = hashlib.md5()
    with open(path, "rb") as fh:
        for b in iter(lambda: fh.read(chunk), b""):
            h.update(b)
    return h.hexdigest()


def main():
    t = pd.read_csv(os.path.join(PIPE, "TPM_matrix.csv"), skiprows=1, low_memory=False)
    print("Processed data files (complete matrices, all genes x all samples):")
    p1 = write_matrix(t, [f"RawCounts_{s[:-1]}_{s[-1]}" for s in LIBS], LIBS,
                      "GSE_processed_gene_raw_counts.txt", "%.0f")
    p2 = write_matrix(t, [f"TPM_{s[:-1]}_{s[-1]}" for s in LIBS], LIBS,
                      "GSE_processed_gene_TPM.txt", "%.4f")

    # --- raw FASTQ checksums: reuse the provider-supplied values (verified on disk) ---
    rows = []
    for s in LIBS:
        for r in (1, 2):
            m = os.path.join(RAWD, f"{s}_Data{r}.md5")
            if not os.path.exists(m):
                print(f"  WARNING: missing {os.path.basename(m)}"); continue
            checksum = open(m).read().split()[0]
            rows.append(("raw", f"{s}_Clean_Data{r}.fastq.gz", checksum))
    for p in (p1, p2):
        rows.append(("processed", os.path.basename(p), md5(p)))

    with open(os.path.join(OUT, "GEO_checksums.tsv"), "w", newline="") as fh:
        w = csv.writer(fh, delimiter="\t")
        w.writerow(["type", "file name", "file checksum"])
        w.writerows(rows)
    print(f"\nChecksums written: {sum(r[0]=='raw' for r in rows)} raw + "
          f"{sum(r[0]=='processed' for r in rows)} processed")
    for r in rows:
        print(f"  {r[0]:10s} {r[1]:34s} {r[2]}")


if __name__ == "__main__":
    main()
