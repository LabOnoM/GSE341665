# Contributing & getting in touch

This repository accompanies a manuscript on epithelial *Smad4* loss in the dental
epithelium. We welcome questions, corrections and reuse.

## Where to post what

| You want to… | Use |
|---|---|
| Ask how something works, or discuss the biology | **Discussions** → Q&A |
| Report a bug, a wrong number, or a broken script | **Issues** → Bug report |
| Ask for data we have not deposited | **Issues** → Data access request |
| Suggest an analysis or an extension | **Discussions** → Ideas |

## Reporting a problem with a number or figure

This project has a documented history of errors that survived into drafts, so we take
these seriously and would rather hear about one than not. Please include:

1. the exact value as published, and where (figure panel, table row, or line of text);
2. what you believe it should be, and how you obtained that;
3. the script and command you ran, if you re-ran anything.

Every figure value is computed at draw time — nothing is hard-coded — so a discrepancy is
usually reproducible, and we would like to reproduce it.

## Reusing the code

The pipeline is MIT-licensed; the data and figures are CC-BY-4.0. Please cite the
manuscript and the Zenodo DOI (see `CITATION.cff`).

Two things to know before you run it:

- **Read `DEG_table.csv` with `float_precision="round_trip"`.** Very small adjusted *P*
  values are stored in full decimal notation and pandas' default parser silently flushes
  them to `0.0`. See `docs/known-issues.md`.
- **The pipeline fetches external resources at run time** (the pinned PanelApp panel, the
  MGI ortholog report, GEO GSE146855) and **fails loudly** rather than falling back to a
  stored copy. If a source is unreachable the run stops — that is deliberate.

## Scope

We can help with the code and the analyses in this repository. We cannot supply raw
material, and no further wet-laboratory experiments are planned for this study.
