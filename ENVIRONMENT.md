# Software environment

`requirements-reproduction.txt` pins the environment used for the October Python reproduction check. Python 3.12.14 is used for that check. It is separate from the historical analysis environment in Additional file 1. No primary-model package version is inferred from this new environment.

The supplement mode uses the same public matplotlib environment; Graphviz is not required. Reference mode needs R 4.3.2, Seurat 5.3.1, SeuratObject 5.2.0, Matrix 1.6.5, uwot 0.2.3, ggplot2 4.0.0 and patchwork 1.3.2, as recorded in Additional file 1. Locus mode needs samtools and the exact alignment/reference inputs. These optional modes are documented separately from the Python archived-result check.

The primary model used STAR 2.7.10b, Salmon 1.10.0, Ensembl 116 GRCm39, Trim Galore 0.6.10 and Cutadapt 4.9; original versions unrecorded for other tools remain unrecorded. `07_figures.py` and `export_additional_data.py` require no proprietary office or plotting runtime. NNLS marker iteration is sorted; Seurat scripts use seed 1; exact permutation tests enumerate all allocations.
