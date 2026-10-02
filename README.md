# GSE341665

This repository contains code and numerical exports used by the latest BMC Oral Health manuscript, **Epithelial Smad4 loss reprograms ameloblast identity and metabolism**, and its Additional files. The study reports reduced ameloblast-associated expression with wound/stress and sterol-biosynthesis responses and reduced interferon-associated expression. Bulk deconvolution does not establish cell autonomy; sterol-gene expression does not measure SREBP2 activity or cholesterol flux; selective inhibitor transcript changes do not establish general enamel-protein rescue.

The current repository is private during manuscript preparation. The GEO and Zenodo identifiers are retained from the author manuscript; deposit access and the upload of the cleaned data snapshot still need verification. No public release is claimed.

## Code and data

| Location | Role |
|---|---|
| [GitHub GSE341665](https://github.com/LabOnoM/GSE341665) | Current analysis and plotting code; latest numerical exports |
| [GEO GSE341665](https://www.ncbi.nlm.nih.gov/geo/query/acc.cgi?acc=GSE341665) | Provider-filtered sequence files and complete count/TPM matrices |
| [Zenodo data DOI](https://doi.org/10.5281/zenodo.21637766) | Assembly, gene-level model objects, exact reproduction inputs, annotated reference and normalized inhibitor sources |
| [Interactive primary-analysis report](http://smad4kotooth.bs-gou.com/) | Browsable historical primary-analysis context; includes exploratory analyses excluded from the current paper and deposits |
| [GEO GSE146855](https://www.ncbi.nlm.nih.gov/geo/query/acc.cgi?acc=GSE146855) | Public external single-cell reference |

The provider's `*_Clean_Data*.fastq.gz` files are adapter- and quality-filtered reads, the rawest files returned to the authors. They are not unfiltered instrument output. The website is a supplementary browsing resource, not a versioned code or data archive.

## Reproduce the current manuscript

Install Python dependencies from `requirements-reproduction.txt`. Download the data snapshot matching `data/INPUT_MANIFEST.json`. From the repository root:

```bash
python reproduce.py --data-root /path/to/unpacked/zenodo \
  --work-dir /path/to/new/scratch/workspace \
  --destination /path/to/new/results --mode manuscript
```

The command extracts the source workbooks, recomputes exact program tests, NNLS fractions/reference specificity and pinned-panel enrichment, regenerates Figures 4–7, exports all six Data S8 CSVs, recalculates the 195 inhibitor values and 52 ANOVA/Dunnett comparisons, and writes Additional files 2–9 including the Data S8 workbook. The workbook uses public openpyxl. Supply `--geo-root /path/to/GEO` to cross-check the entire retained TPM matrix against GEO.

```bash
python reproduce.py --data-root /path/to/unpacked/zenodo \
  --work-dir /path/to/new/scratch/supplement \
  --destination /path/to/new/supplement --mode supplement
```

This regenerates Figures S4–S6 from the archived evidence. `--mode reference` downloads GSE146855 and regenerates the Seurat analysis/Figure S3; it additionally requires R/Seurat and network access. `--mode locus --bam-dir ... --reference-dir ...` re-extracts Figure S5 evidence from aligned reads and requires samtools and the named BAM/reference inputs. See [the reproduction map](docs/reproduction-map.md) for every retained script/input/output and the validation limits.

These commands reproduce archived-result analyses and figures. They do not execute a complete FASTQ-to-DESeq2 refit. Exact package versions missing from the original primary analysis remain unrecorded. The original alignment/quantification/model-fitting commands are described in Additional file 1 and the linked primary report; a complete executable upstream workflow has not been recovered into this repository.

Read adjusted P values with `float_precision="round_trip"`. Five source zeros represent numerical underflow; 613 missing adjusted P values remain missing. Gene-level inference uses the archived RUV-adjusted count model; heatmaps/PCA use unadjusted TPM. Program scores use exact label permutations, with a minimum two-sided P of 0.10 for three libraries per genotype.

## Scope and citation

Unused MR/GWAS, isoform/circular-RNA/exon-usage outputs, historical main-figure renderers, manuscript builders and deposit/publishing scripts are excluded from the current repository. Table S1 retains the provenance audit of the superseded gene list because that table is explicitly part of current Additional file 1. Historical commits remain available through Git; they do not define the current manuscript version.

Publication PDFs and the Excel workbook remain in the journal/data package; this code repository versions their numerical CSV evidence and regenerates the computational figures/workbook. Use `CITATION.cff` and record the exact commit. Code remains under the repository's existing MIT license; deposited data retain their existing CC-BY-4.0 terms. See [data access](docs/data-access.md) and [known limitations](docs/known-issues.md).
