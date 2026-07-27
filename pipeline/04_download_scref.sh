#!/usr/bin/env bash
# 04_download_scref.sh  --  Smad4 dcKO JDR re-analysis pipeline (step 4/7)
#
# Download the public single-cell RNA-seq reference used for cell-type
# deconvolution: GSE146855 (Chiba et al. 2020, Front Cell Dev Biol),
# P7 mouse incisor dental epithelium+mesenchyme, 10x Genomics, 6,260 cells.
# Idempotent: skips the download if the tar is already present.
#
# Output: output/scref/GSE146855/  (10x matrix: barcodes/features/matrix + h5)
set -euo pipefail
HERE="$(cd "$(dirname "${BASH_SOURCE[0]}")" && pwd)"
DEST="$HERE/output/scref/GSE146855"
URL="https://ftp.ncbi.nlm.nih.gov/geo/series/GSE146nnn/GSE146855/suppl/GSE146855_RAW.tar"
mkdir -p "$DEST"
if [ ! -f "$DEST/GSE146855_RAW.tar" ]; then
  echo "Downloading GSE146855_RAW.tar (~75 MB)..."
  curl -sSL -o "$DEST/GSE146855_RAW.tar" "$URL"
fi
tar -xf "$DEST/GSE146855_RAW.tar" -C "$DEST"
# standardise 10x filenames into $DEST/mtx (Read10X expects barcodes/features/matrix)
mkdir -p "$DEST/mtx"
cp "$DEST"/*barcodes.tsv.gz "$DEST/mtx/barcodes.tsv.gz"
cp "$DEST"/*genes.tsv.gz    "$DEST/mtx/features.tsv.gz"
cp "$DEST"/*matrix.mtx.gz   "$DEST/mtx/matrix.mtx.gz"
echo "Reference ready at $DEST/mtx"
