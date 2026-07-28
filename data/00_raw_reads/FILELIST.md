# File manifest — Raw sequencing reads

> **These files are not stored in Git** (too large for a code repository).
> They are archived on **Zenodo** (DOI: `[[Zenodo DOI — insert on release]]`) and,
> where noted, are also browsable in the interactive analysis report at
> <http://smad4kotooth.bs-gou.com/>. Download them into **this folder** to reproduce
> the pipeline; paths in the code resolve relative to the repository root.


> **Processing level.** These files are the sequencing provider's adapter- and
> quality-filtered output ("clean data"), not unfiltered instrument reads — the latter were
> never returned to the authors. They are the input to Trim Galore and the rawest data
> available for this study, and are what GEO holds in the "raw" role. Do not use
> `Data.tar` (a byte-identical duplicate) or `*_val_*.fq.gz` (post-trimming output).
> All 12 verified against their MD5 checksums on 2026-07-28 (12/12 match).

| Expected file | Description | Approx. size | Source |
|---|---|---|---|
| `CON{1,2,3}_R{1,2}.fastq.gz` | Paired-end 150 bp reads, 3 control libraries | ~3-5 GB each | Zenodo / SRA on GEO release |
| `dcKO{1,2,3}_R{1,2}.fastq.gz` | Paired-end 150 bp reads, 3 Smad4-dcKO libraries | ~3-5 GB each | Zenodo / SRA on GEO release |
| `md5sums.txt` | Checksums for all FASTQ | <1 KB | Zenodo |

_This folder is intentionally otherwise empty; the manifest keeps it under version control._
