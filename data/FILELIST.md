# Inputs from the data deposit

The existing manuscript identifier is [10.5281/zenodo.21637766](https://doi.org/10.5281/zenodo.21637766). The ROAST snapshot is prepared locally for manual upload; obtain the matching snapshot from the authors until its uploaded version is confirmed. `INPUT_MANIFEST.json` names and hashes every file in the new data snapshot. Code is in this repository; the data record contains no code. The entry point verifies hashes and creates the expected project-shaped working paths automatically.

The six full auxiliary `.STAR.counts` files support the nonretained Sox2/Lgr5 evidence. They are not the Salmon/tximport gene-level model inputs. The manuscript claim JSON preserves the captured source used by the existing Data S3 paragraph index, not later author metadata or BMC paragraph numbering.

Every required runtime source is declared in INPUT_MANIFEST.json. The archived DESeq2_Processing.Rdata now supplies new program-model inputs (counts, full-precision TPM/effective lengths and W_1), without refitting or modifying original DESeq2 results. No cached fitted program model is a runtime input. Original GSE146855, pinned PanelApp/MGI records, regional locus alignments, the source workbooks and Figure_3.prism support the other retained workflows. Counts derive from the current manifest rather than a historical hard-coded target.

RAW_EVIDENCE_MANIFEST.json separately retains the original blot deck and Figure_2.prism. The author will upload the matching validated data snapshot manually; the existing DOI is not evidence that revised files are already deposited.
