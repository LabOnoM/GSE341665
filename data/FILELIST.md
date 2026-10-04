# Reproduction inputs

`INPUT_MANIFEST.json` identifies the 56 required source inputs and their SHA-256 hashes. `RAW_EVIDENCE_MANIFEST.json` identifies the original Western blot image deck and both Prism projects; Figure 3 Prism overlaps the computational manifest. The matching data package therefore contains 58 scientific files, plus documentation and manifests.

The package retains the assembly, saved gene-level model, full-precision source workbooks, captured gene/expression tables, single-cell reference and raw count archive, pinned PanelApp/MGI records, and regional locus evidence. Archived outputs serving as inputs remain necessary because the gene-level model is reused. Fitted ROAST objects and journal results are regenerated and are not included.

Validation benchmarks are under `results/` in this repository. The matching revised data snapshot is prepared for manual author upload to [Zenodo](https://doi.org/10.5281/zenodo.21637766); obtain it from the authors until that version is confirmed. Private reviewer access details are supplied separately.
