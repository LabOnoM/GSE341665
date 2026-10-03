<!-- Current methods: ROAST migration, 3 October 2026. No live Zenodo action performed. -->
# GSE341665

This repository contains code and numerical exports used by the latest BMC Oral Health manuscript, **Epithelial Smad4 loss reprograms ameloblast identity and metabolism**, and its Additional files. The study reports reduced ameloblast-associated expression with wound/stress and sterol-biosynthesis responses and reduced interferon-associated expression. Bulk deconvolution does not establish cell autonomy; sterol-gene expression does not measure SREBP2 activity or cholesterol flux; selective inhibitor transcript changes do not establish general enamel-protein rescue.

The current repository is private during manuscript preparation. The GEO and Zenodo identifiers are retained from the author manuscript; deposit access and the upload of the cleaned data snapshot still need verification. No public release is claimed.

## Code and data

| Location | Role |
|---|---|
| [GitHub GSE341665](https://github.com/LabOnoM/GSE341665) | Current analysis and plotting code; latest numerical exports |
| [GEO GSE341665](https://www.ncbi.nlm.nih.gov/geo/query/acc.cgi?acc=GSE341665) | Provider-filtered sequence files and complete count/TPM matrices |
| [Zenodo data DOI](https://doi.org/10.5281/zenodo.21637766) | Assembly, gene-level model objects, exact reproduction inputs, annotated reference and normalized inhibitor sources |
| [Interactive primary-analysis report](http://smad4kotooth.nexitia.ai/) | Browsable historical primary-analysis context; includes exploratory analyses excluded from the current paper and deposits |
| [GEO GSE146855](https://www.ncbi.nlm.nih.gov/geo/query/acc.cgi?acc=GSE146855) | Public external single-cell reference |

The provider's `*_Clean_Data*.fastq.gz` files are adapter- and quality-filtered reads, the rawest files returned to the authors. They are not unfiltered instrument output. The website is a supplementary browsing resource, not a versioned code or data archive.

## Run every retained workflow

Use Python 3.12 and the exact tested library. Use R 4.3.2 for manuscript program inference and the single-cell reference. On macOS/Linux, from the repository root:

```bash
python3.12 -m venv /path/to/isolated/python-env
/path/to/isolated/python-env/bin/python -m pip install -r environment/python-lock.txt
Rscript environment/setup_program.R --library /path/to/isolated/program-library --install
export PROGRAM_R_LIBRARY=/path/to/isolated/program-library
# For all/reference modes, additionally restore the reference library:
Rscript environment/setup_reference.R --library /path/to/isolated/R-library --install
export R_LIBS_USER=/path/to/isolated/R-library
```

The R installer prints its plan and writes only to the specified library. It restores 149 tested package versions, including the six originally recorded reference packages. CRAN source builds require a C/C++/Fortran toolchain and the system libraries requested by those packages. This tested dependency lock is separate from the incompletely recorded original primary-fit environment. The lightweight `requirements-reproduction.txt` pins the direct Python dependencies; `environment/python-lock.txt` also pins their tested transitive dependencies.

Obtain and unpack the **matching** data snapshot from the author or the matching version after manual Zenodo upload. Its declared inputs must match `data/INPUT_MANIFEST.json`; the local staged snapshot has not yet been synchronized to the live record. All computation uses those declared files after software installation:

```bash
/path/to/isolated/python-env/bin/python reproduce.py \
  --data-root /path/to/unpacked/zenodo \
  --work-dir /path/to/new/scratch/workspace \
  --destination /path/to/new/results --mode all
```

`--mode all` rebuilds the reference from the original GSE146855 count archive, regenerates Figure S6, rebuilds the AI panel from original pinned PanelApp/MGI records, extracts source workbooks, fits the W_1-adjusted ROAST program model and supplementary sensitivities, recomputes descriptive scores and NNLS fractions, renders Figures 4–7, exports Additional files 2–7 and Zenodo Data S1/S2 and all Data S3 CSVs/workbook, re-extracts Figure S3 evidence from six deposited regional BAMs, and renders Figures S1, S3 and S5. Missing/mismatched inputs, environments or scientific benchmarks stop execution. No data downloads or private project helpers are required. Working and output directories must be empty, separate from each other and outside both the repository and data directory. The caller's current directory is unrestricted. Use `--rscript /path/to/Rscript` if multiple R installations exist.

Outputs include `Figures/`, `Additional_files/`, `Inhibitor_data/`, `Supplementary/`, `Program_analysis/`, `Reference/`, `Panel/`, `Locus/`, individual command logs, `EXECUTION.json` and `VALIDATION.json`. The final validation compares current scientific tables/sheets (counts in `VALIDATION.json`), verifies all 15,803 retained genes and adjusted-P masks, and checks 22 figure files (PDF/PNG for all retained figures, plus SVG for main Figures 4–7). It does not claim pixel-identical layouts across operating systems.

| Mode | Purpose | Extra environment |
|---|---|---|
| `manuscript` | Figures 4–7, all current computational numerical exports/workbook | Python plus pinned program R library |
| `supplement` | Figures S1, S3 and S5 from archived locus evidence | Python |
| `reference` | Raw reference rebuild, signature, labels and Figure S6 | Recorded R library |
| `reference-figure` | Figure S6 from the captured annotated reference object | Recorded R library |
| `panel` | Rebuild the panel from explicit original source records | Python |
| `locus` | Re-extract Figure S3 evidence and render it | Python/pysam |
| `all` | Every retained scientific workflow in dependency order | Python plus program and reference R libraries |

Use the same data/work/output arguments with any mode. `--bam-dir /path/to/full/alignments` optionally tests the original six full BAMs, indexes and STAR junction files. The deposited regional BAMs preserve all original alignment records overlapping chr18:73,768,343–73,837,915; full and regional extraction were independently checked against the same archived evidence. These subsets support this locus analysis and do not reconstruct genome-wide alignments. The exact 69,573-base Ensembl-116 sequence and locus GTF are included; neither a complete genome FASTA nor a separate samtools executable is required. Pysam 0.24.1 executes the actual bundled samtools 1.24 subcommands with the original flags.

Supply `--geo-root /path/to/unpacked/GEO` to cross-check all 15,803 × six retained TPM values against `GSE_processed_gene_TPM.txt` after the recorded four-decimal GEO export rounding. The source workbooks retain full precision. The production program stage additionally checks all 94,818 archived count entries against GEO when this option is supplied.

## Settings and independent checks

The reference uses genes detected in at least three cells, >200 and <7,000 features/cell, <25% mitochondrial counts, LogNormalize at 10,000, 2,000 vst variable features, 30 PCs, the first 20 downstream PCs, k=20 Annoy neighbors with 50 trees, and Louvain clustering at resolution 0.6. Global seed is 1; PCA/UMAP seed is 42; clustering seed is 0. UMAP uses uwot/cosine, 30 neighbors and min.dist 0.3. Cross-cluster marker z-score means are rounded to two decimals before choosing the first maximum. The maturation marker set includes **Gm17660**. Additional file 1 has been corrected to state these previously omitted details.

```bash
/path/to/isolated/python-env/bin/python validate.py \
  --data-root /path/to/unpacked/zenodo \
  --work-dir /path/to/completed/workspace \
  --destination /path/to/completed/results --mode all
/path/to/isolated/python-env/bin/python -m unittest discover -s tests -v
```

Ten regression checks plus production R synthetic-model checks protect input/output separation, existing work, interrupted downloads, missing source records, and tiny/missing probability states. GitHub Actions runs those checks and syntax validation and exercises both program-model designs on synthetic counts; it does not obtain restricted study data or certify a complete scientific rerun.

To deliberately re-download the public reference or panel, use `PYTHON=/path/to/python bash pipeline/04_download_scref.sh --download --destination /new/reference` or `python pipeline/17_fetch_ai_panel.py --download --save-sources /new/sources --output-dir /new/panel`. Both live routes were executed in the audit. Downloads verify the fixed reference SHA-256; the API explicitly requests panel v4.37. MGI changes over time, so live retrieval records its source hash and is separate from frozen manuscript reproduction. Offline runs explicitly read archived original records and never silently substitute a hard-coded panel.

These commands reproduce the computational evidence in the current manuscript. They do not refit the original FASTQ-to-DESeq2 model. Its executable upstream commands and some original primary-fit package versions have not been recovered. Wet-laboratory Figures 1–3/S2/S4/S7 and Tables S1–S3 are experimental artwork/records; plotting code cannot recover missing raw scans, Ct values or band intensities. See [the complete dependency map](docs/reproduction-map.md), [validation](docs/reproduction-validation.json) and [remaining limitations](docs/known-issues.md).

Read adjusted P values with `float_precision="round_trip"`. Five source zeros represent numerical underflow; 613 missing adjusted P values remain missing. Nonzero probability checks use zero absolute tolerance. Gene-level inference uses the archived RUV-adjusted count model; heatmaps/PCA use unadjusted TPM. The six per-library program z scores remain descriptive. Primary inference uses ROAST (mean statistic; two-sided directional PValue; 99,999 rotations; seed 20261003; ordinary-P BH across four programs) after reconstruction of lengthScaledTPM-style counts within the 15,803 retained genes, TMM and limma-voom (~ W_1 + genotype; dcKO minus CON; three residual degrees of freedom). Genotype-only directional ROAST (seed 20261004) is a sensitivity analysis. Mixed-direction ROAST in both designs is secondary. Membership is generated from `pipeline/programs.py`. The old program-score permutation test was valid for its statistic; its P=0.10 floor does not apply to the new model. The separate locus/exon exact permutation analysis is unchanged.

`original_DESeq2_refit=false` and `program_limma_voom_fit=true` distinguish reuse of gene-level inference from the new program fit. Programs were curated during earlier analysis; these are exploratory associations from the same six libraries. Archived W_1 derivation, donor/litter mapping and very large weights for abundant genes remain limitations. Raw-P 0.00001 is the finite-rotation resolution floor.

The canonical four-row `Supplementary_Data_S3_program_tests.csv` includes raw/BH ROAST probabilities and descriptive score means, with stable program IDs. `Supplementary_Data_S3_program_test_sensitivity.csv` has 12 secondary rows with explicit model and adjustment-family labels. Data S6 remains the six-library score table. GEO checking is optional for both archived counts and rounded TPM. The code plus matching data snapshot runs without a sibling GEO directory.

## Scope and citation

Unused MR/GWAS, isoform/circular-RNA/exon-usage outputs, historical main-figure renderers, manuscript builders and deposit/publishing scripts are excluded from the current repository. Supplementary Table S4 documents the 31 GREEN-rated PanelApp genes actually used. The retired 20-gene provenance audit and its export code are excluded from current manuscript outputs; its historical record remains in Git/project history. Historical commits remain available through Git; they do not define the current manuscript version.

Publication PDFs and the Excel workbook remain in the journal/data package; this code repository versions their numerical CSV evidence and regenerates the computational figures/workbook. Use `CITATION.cff` and record the exact commit. Code remains under the repository's existing MIT license; deposited data retain their existing CC-BY-4.0 terms. See [data access](docs/data-access.md) and [known limitations](docs/known-issues.md).

Experimental source projects `Figure_2.prism` and `Figure_3.prism` are staged for the Zenodo record, together with the original Western blot image deck. Data S1/S2 remain open CSVs in `results/inhibitor_data/` and the data deposit; reproduction exports them to `Inhibitor_data/`. They are no longer separate journal Additional-file uploads. Additional files 2–7 retain Data S3, S4, S5, S6, S7 and S8, respectively. Figure_3.prism is the executable source for the independent inhibitor recalculation; Figure_2.prism supplies additional experimental evidence. Duplicate normalized PZFX files are excluded. The matching local package retains 101 scientific files: 73 declared inputs/validation benchmarks, two additional original evidence files and 26 program-model/QC outputs. Root documentation and manifests are counted separately. Recalculation reads the declared native data sheets, preserves all measured values, and exports native group spellings and source-sheet titles. Inhibitor_data/Inhibitor_Statistics_Provenance.json records the source hash and calculation method.

Figures 4–7 are exported as PDF, SVG and PNG from the same plotting objects. SVG labels remain editable text. Validation requires each SVG to be valid XML with a viewBox and to contain no external resource references.

Experimental animals were confirmed by the authors on 3 October 2026 to be male, including all six RNA-seq libraries. Xist/Y-linked expression is retained in Data S3, with a brief note in Additional file 1. Figure 4 contains panels A–C; these marker values have no current figure-panel assignment and are not an independent determination of animal sex. Supplementary identifiers follow first citation in the current manuscript; see [the renumbering map](docs/supplementary-numbering.md).
