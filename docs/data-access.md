# Where the data lives

Four destinations, no overlap. Going to the wrong one is the commonest reason people
cannot find a file.

| Archive | Holds | Identifier |
|---|---|---|
| **NCBI GEO** | the 12 paired-end sequence files, and the gene-level count and TPM matrices | `[[GSE — on release]]` |
| **Zenodo** | de-novo transcriptome assembly, DESeq2/RUVSeq objects, exploratory analyses, frozen code snapshot | `[[DOI — on release]]` |
| **This repository** | the analysis code and the small result tables behind the figures | you are here |
| **Analysis report** | the complete primary analysis, browsable, including the exploratory arms | <http://smad4kotooth.bs-gou.com/> |

## About the sequence files

They are named `*_Clean_Data[12].fastq.gz` and are the **sequencing provider's adapter-
and quality-filtered output**, not unfiltered instrument reads — those were never returned
to us and cannot be deposited. They are the input to the Trim Galore step and are the
rawest data available for this study.

This matters if you are comparing read counts against a study that deposits unfiltered
reads: ours are already filtered once before our pipeline touches them.

## What is deliberately not deposited anywhere

The STAR genome index, the Ensembl release 116 GRCm39 genome and GTF, and the Pfam
database. All are public and regenerable; citing the release is more useful than
re-hosting several gigabytes.

## Data under embargo

Until the manuscript is published, the GEO and Zenodo records are private. Reviewers and
editors receive access tokens through the journal. If you need access before publication,
open a **Data access request** issue.
