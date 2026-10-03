# Tested software environment

Python 3.12.14 and `environment/python-lock.txt` govern the reproduction workflow. The new program stage requires R 4.3.2 and eight exact packages in `environment/program-packages.tsv`: limma 3.58.1, edgeR 4.0.16, statmod 1.5.0, jsonlite 2.0.0, digest 0.6.37, Rcpp 1.1.0, locfit 1.5-9.12 and lattice 0.21-9. `setup_program.R --library /isolated/library --install` restores that library; omit `--install` to check it. Set `PROGRAM_R_LIBRARY` to that directory. The workflow checks the lock before fitting.

`all`, `reference` and `reference-figure` additionally require the 149-package reference lock and setup procedure described in README. Reference packages use a separate `R_LIBS_USER`; the program stage prepends its own library only for program execution. `all` therefore supports both environments without replacing system libraries.

Program RNG is Mersenne-Twister/Inversion/Rejection, with seeds 20261003 and 20261004. ROAST uses exactly 99,999 rotations. Single-cell seeds and the separate exhaustive locus permutation analysis remain unchanged. No package installation occurs automatically during reproduction.

Original FASTQ-to-DESeq2 processing used STAR 2.7.10b, Salmon 1.10.0, Ensembl 116 GRCm39, Trim Galore 0.6.10 and Cutadapt 4.9. Unrecorded original versions remain unrecorded; the new program library does not reconstruct them. `original_DESeq2_refit=false`; `program_limma_voom_fit=true` in manuscript/all modes.
