# GSE341665

Analysis and figure-generation code for the epithelial Smad4 study and its Additional files.

## Code and data

| Resource | Contents |
|---|---|
| [GitHub](https://github.com/LabOnoM/GSE341665) | Analysis scripts, environment specifications and numerical validation benchmarks |
| [GEO GSE341665](https://www.ncbi.nlm.nih.gov/geo/query/acc.cgi?acc=GSE341665) | Provider-filtered RNA-sequencing reads and complete count and TPM matrices |
| [Zenodo data record](https://doi.org/10.5281/zenodo.21637766) | Transcriptome assembly, saved gene-level model objects, reproduction inputs and original experimental evidence |
| [GEO GSE146855](https://www.ncbi.nlm.nih.gov/geo/query/acc.cgi?acc=GSE146855) | External mouse tooth single-cell reference |

Use the data snapshot matching `data/INPUT_MANIFEST.json`. The revised snapshot is prepared for manual author upload; until that upload is confirmed, obtain it from the authors. Reviewer credentials are supplied separately through the journal. Repository access settings are unchanged.

For assay-specific reference genes and loading controls, see [experimental measurements](docs/experimental-data.md).

The data package contains 56 required input files and two additional original evidence files. Native Figure 2 and Figure 3 Prism projects and the original Western blot image deck are retained. The 17 validation benchmarks are versioned with this code in `results/`, rather than duplicated in the data package. Journal documents, generated figures and regenerated program-model outputs are supplied separately or produced by the commands below.

## Reproduce the analyses

Use Python 3.12 and R 4.3.2 with the supplied package specifications:

```bash
python3.12 -m venv /path/to/python-env
/path/to/python-env/bin/python -m pip install -r environment/python-lock.txt
Rscript environment/setup_program.R --library /path/to/program-library --install
export PROGRAM_R_LIBRARY=/path/to/program-library
# Required additionally for all/reference modes:
Rscript environment/setup_reference.R --library /path/to/reference-library --install
export R_LIBS_USER=/path/to/reference-library
```

The R installers write only to the specified libraries. Source builds require the relevant compiler toolchain and system libraries. Verify the complete ZIP or extracted data directory before running. This checks consecutive folders, every file checksum, all manifests, and compatibility with this code version:

```bash
python verify_data_package.py /path/to/Smad4_Reproduction_Inputs_20261004.zip
```

The data folders are numbered consecutively from `01_processed_data` to `07_locus_alignments`; experimental sources are under `05_experimental_data`. The data README, dictionary and sample metadata describe their contents. Run the complete workflow with separate, empty working and output directories:

```bash
/path/to/python-env/bin/python reproduce.py \
  --data-root /path/to/unpacked/data \
  --work-dir /path/to/new/workspace \
  --destination /path/to/new/results --mode all
```

All 56 declared input hashes are checked before computation. `all` rebuilds the single-cell reference and pinned AI panel, extracts source workbooks, fits the ROAST program models, calculates descriptive program scores and NNLS proportions, recalculates inhibitor statistics, extracts locus evidence, and generates the computational figures and Additional-file tables. It requires no data downloads after environment installation. Use `--rscript /path/to/Rscript` to select a particular R installation.

| Mode | Outputs |
|---|---|
| `manuscript` | Figures 4–7, Additional files 2–7, Data S1/S2 and source tables |
| `supplement` | Figures S1, S3 and S5 from captured locus inputs |
| `reference` | Reference rebuild, labels, signatures and Figure S6 |
| `reference-figure` | Figure S6 from the annotated reference object |
| `panel` | AI panel rebuilt from pinned PanelApp/MGI records |
| `locus` | Locus evidence extracted from regional BAMs and Figure S3 |
| `all` | All workflows in dependency order |

Figure S3 contains panels A–C; complete splice-junction evidence remains in the data inputs and locus outputs. Figure S5 is a single composite containing panels A and B, exported as `FigureS5_workflow.pdf`, `.svg` and `.png`. Main Figures 4–7 also have PDF, SVG and PNG exports. Author-supplied experimental artwork is provided with the journal submission.

The final validation compares 35 numerical tables, workbook sheets and model/reference checks, and checks 20 figure files. Zero and missing probability masks must match exactly; nonzero values use zero absolute tolerance. Scientific values in all 15,803 retained gene rows, program results, AI enrichment, deconvolution and inhibitor comparisons are preserved. Layout may vary slightly across operating systems.

```bash
/path/to/python-env/bin/python validate.py \
  --data-root /path/to/unpacked/data \
  --work-dir /path/to/completed/workspace \
  --destination /path/to/completed/results --mode all
/path/to/python-env/bin/python -m unittest discover -s tests -v
Rscript tests/test_program_models.R
```

Outputs include figures, numerical tables, model objects, command logs, `EXECUTION.json` and `VALIDATION.json`. Python regression checks and the R synthetic-model checks cover input separation, missing/corrupt evidence, probability handling and the production program models. GitHub Actions runs code checks without restricted study inputs.

## Analysis settings

Gene-level statistics are retained from the DESeq2 likelihood-ratio test comparing `~ W_1 + genotype` with `~ W_1`. The original FASTQ-to-DESeq2 fit is not rerun by this package. The saved model and transcriptome assembly are included for inspection.

Primary gene-program inference uses lengthScaledTPM-style counts reconstructed within the 15,803 retained genes, TMM normalization, limma-voom and ROAST (`~ W_1 + genotype`; dcKO minus CON; three residual degrees of freedom). The mean statistic, two-sided directional P value, 99,999 rotations, seed 20261003 and Benjamini–Hochberg adjustment across four programs are used. Genotype-only directional ROAST is a sensitivity analysis; mixed-direction ROAST in both designs is secondary. Per-library z scores are descriptive. Membership is defined in `pipeline/programs.py`.

`--geo-root /path/to/GEO` optionally checks retained TPM and count entries against GEO exports. Source workbooks preserve full precision; GEO TPM was exported to four decimal places. `--bam-dir /path/to/alignments` optionally substitutes full alignments for the six deposited regional BAMs. Regional BAMs retain the original records overlapping chr18:73,768,343–73,837,915 and support the locus analysis.

The [reproduction map](docs/reproduction-map.md) describes each dependency, and [Supplementary Methods](docs/supplementary_methods.md) describe the analyses. The package reproduces computational outputs from the supplied inputs; it does not reconstruct raw Ct values or band intensities from normalized measurements. See [data access](docs/data-access.md) and [validation results](docs/reproduction-validation.json).

## Citation and licensing

Use `CITATION.cff` and record the exact code commit and matching data version. Code is licensed under MIT; deposited data retain CC-BY-4.0 terms. Archived gene-level outputs that serve as inputs are retained, while regenerated publication outputs are omitted from the source-data package.
