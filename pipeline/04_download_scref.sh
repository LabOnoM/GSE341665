#!/usr/bin/env bash
# 04_download_scref.sh  --  Smad4 dcKO JDR re-analysis pipeline (step 4/7)
#
# Download the public single-cell RNA-seq reference used for cell-type
# deconvolution: GSE146855 (Chiba et al. 2020, Front Cell Dev Biol),
# P7 mouse incisor dental epithelium+mesenchyme, 10x Genomics, 6,260 cells.
# Requires --archive /path/to/original.tar or --download; verifies the SHA-256.
#
# Output: output/scref/GSE146855/  (10x matrix: barcodes/features/matrix + h5)
set -euo pipefail
HERE="$(cd "$(dirname "${BASH_SOURCE[0]}")" && pwd)"
exec "${PYTHON:-python3}" "$HERE/reference_download.py" "$@"
