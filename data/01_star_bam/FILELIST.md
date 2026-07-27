# File manifest — STAR two-pass alignments

> **These files are not stored in Git** (too large for a code repository).
> They are archived on **Zenodo** (DOI: `[[Zenodo DOI — insert on release]]`) and,
> where noted, are also browsable in the interactive analysis report at
> <http://smad4kotooth.bs-gou.com/>. Download them into **this folder** to reproduce
> the pipeline; paths in the code resolve relative to the repository root.

| Expected file | Description | Approx. size | Source |
|---|---|---|---|
| `{CON,dcKO}{1,2,3}.Aligned.sortedByCoord.out.bam` | STAR GRCm39 two-pass BAM (assembly + circRNA input) | ~4-8 GB each | Zenodo (BAM_Files.zip) |
| `*.STAR.Log.final.out` | Per-sample mapping statistics (93.7-93.8% unique) | <10 KB | Zenodo |

_This folder is intentionally otherwise empty; the manifest keeps it under version control._
