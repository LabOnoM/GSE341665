# Inputs from the data deposit

Download and unpack the data record DOI [10.5281/zenodo.21637766](https://doi.org/10.5281/zenodo.21637766). The local staged package is being updated to the 2 October manuscript scope; its live record must be synchronized before readers can use this manifest. `INPUT_MANIFEST.json` names and hashes every file in the new data snapshot. Code is in this repository; the data record contains no code. The entry point verifies hashes and creates the expected project-shaped working paths automatically.

The six full auxiliary `.STAR.counts` files support the nonretained Sox2/Lgr5 evidence. They are not the Salmon/tximport gene-level model inputs. The manuscript claim JSON preserves the captured source used by the existing Data S3 paragraph index, not later author metadata or BMC paragraph numbering.

The 70 declared files now include the original GSE146855 reference archive, pinned PanelApp/MGI records and six regional BAM/index/junction inputs. All retained workflows run without data downloads; `--mode all` validates source-reference labels, signatures, locus evidence and manuscript exports. These regional reads preserve the original alignment records and support only the stated locus.

The computational manifest retains 70 inputs/exports, using Figure_3.prism in place of the two duplicate PZFX files. Data S1/S2 measurements/statistics are unchanged; source labels, source-file fields and sheet titles now identify the native project. `RAW_EVIDENCE_MANIFEST.json` in the staged data deposit retains the original blot deck and the two editable Prism projects (Figure_3.prism also appears in the executable manifest); absence from the executable manifest never makes an experimental source disposable.
