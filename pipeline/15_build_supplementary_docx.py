#!/usr/bin/env python3
"""
15_build_supplementary_docx.py -- assemble the single supplementary Word file.

Journals expect one supplementary document. This merges:
  * Supplementary Methods (Parts A and B of Supplementary_Methods_Pipeline.md)
  * Supplementary Figures S1-S6, each with its image embedded and its legend
    (S5's legend is not retyped here -- it is read from the file that
    21_smad4_loxp_figure.py generates alongside the figure, so the two can
    never drift apart)
  * the index of Supplementary Data files

The Mermaid workflow block is NOT passed to pandoc (pandoc would emit it as raw
code). It is replaced by the Graphviz-rendered Figure S4 produced by
14_workflow_diagram.py; the Mermaid source stays in the markdown for readers of
the repository.

Output: 00.Manuscript/Submission_02_JMM_2026-07/4_Journal_JMM/Supplementary/Supplementary_Information.docx
"""
import filecmp, os, re, shutil, subprocess, sys

sys.path.insert(0, os.path.dirname(os.path.abspath(__file__)))
import ai_gene_provenance as prov

ROOT = os.path.dirname(os.path.dirname(os.path.dirname(os.path.abspath(__file__))))
SUP = os.path.join(ROOT, "00.Manuscript", "Submission_02_JMM_2026-07", "4_Journal_JMM", "Supplementary")
FIGSRC = os.path.join(ROOT, "scripts", "JDR_reanalysis", "figures")
BUILD = os.path.join(SUP, "_supplementary_build.md")
OUT = "Supplementary_Information.docx"

TITLE = """# Supplementary Information

**Epithelial Smad4 loss reprograms ameloblast identity and metabolism**

Xinyu Zheng, Dian Su, Fei Hu, Xianghui Liu, Jinshuo Wu, Zhiwei He, Xiaohua Xie

---

## Contents

- **Supplementary Methods** — Part A: primary bulk RNA-sequencing pipeline; Part B: re-analysis performed for this manuscript
- **Supplementary Figure S1** — the inducible K14-creER driver reproduces the enamel phenotype
- **Supplementary Figure S2** — tamoxifen dose-dependence of the inducible phenotype
- **Supplementary Figure S3** — single-cell reference used for deconvolution (UMAP and marker dot plot)
- **Supplementary Figure S4** — complete analysis workflow (A: primary pipeline; B: re-analysis)
- **Supplementary Figure S5** — read-level genotyping of the conditional *Smad4* allele; the Cre excision product is not detected
- **Supplementary Figure S6** — mouse lines, crosses, and the allocation of each experiment to a figure
- **Supplementary Table S1** — provenance audit of the retired hand-curated AI-gene list
- **Supplementary Table S2** — verified provenance of the 31 PanelApp genes used in the enrichment
- **Supplementary Data S1–S5** — machine-readable result tables

---

"""

FIGURES = [
    ("S1", "FigureS1_inducible_phenotype.png", "4.6in",
     "**Supplementary Figure S1. The inducible K14-creER driver reproduces the enamel phenotype.** "
     "Hematoxylin-eosin sections of control (Smad4^fl/fl^) and K14-creER;Smad4^fl/fl^ molars after "
     "subrenal-capsule transplantation, shown at increasing magnification (C, D and E, F are "
     "progressive magnifications of A, B). am, ameloblasts; e, enamel space; d, dentin; arrows "
     "indicate enamel matrix."),
    ("S2", "FigureS2_tamoxifen_titration.png", "4.8in",
     "**Supplementary Figure S2. Tamoxifen dose-dependence of the inducible phenotype.** "
     "Contact radiographs of K14-creER;Smad4^fl/fl^ and control tooth germs after culture in 50, 90 "
     "or 120 µmol/L 4-hydroxytamoxifen; 90 µmol/L was used for subsequent experiments."),
    ("S3", "FigureS3_scref_umap_dotplot.png", "6.4in",
     "**Supplementary Figure S3. The single-cell reference used for deconvolution.** "
     "(A) UMAP embedding of 6,248 quality-controlled cells from GSE146855 (postnatal day 7 mouse "
     "incisor), coloured by the de-novo cell-type annotation. (B) Dot plot of the canonical markers "
     "used to assign each cell type; dot size indicates the percentage of cells expressing the gene "
     "and colour indicates scaled mean expression. Enamel-matrix transcripts (*Amelx*, *Ambn*, "
     "*Enam*, *Mmp20*) are detected at low level in every cluster, reflecting ambient RNA from this "
     "highly secretory tissue — the reason cell types were assigned by cross-cluster specificity "
     "rather than absolute expression."),
    ("S4A", "FigureS4A_pipeline.png", "5.6in",
     "**Supplementary Figure S4. Complete analysis workflow. (A) Primary bulk RNA-sequencing "
     "pipeline.** Steps are labelled with the corresponding section of the Supplementary Methods. "
     "The red-outlined node is the differential-expression result carried into this manuscript. "
     "Grey dashed nodes are analyses that the primary pipeline executed but that are **not** used "
     "here; the reason for each exclusion is given in Supplementary Methods A8."),
    ("S4B", "FigureS4B_reanalysis.png", "6.4in",
     "**(B) Re-analysis performed for this manuscript.** Numbers refer to the numbered scripts in "
     "`scripts/JDR_reanalysis/`. `programs.py` (red) is not a processing step but the single "
     "definition file from which every analysis and figure imports the four curated gene programs, "
     "so that scoring and display can never diverge; it is drawn in sequence for readability. "
     "Steps 11 and 06 both derive from the annotated single-cell reference produced by step 05."),
    # S6 is placed after S5 in the document; the figure itself is generated by
    # 22_breeding_scheme.py, which asserts the driver-to-figure allocation drawn below.
    ("S6", "FigureS6_breeding_scheme.png", "6.5in",
     "**Supplementary Figure S6. Mouse lines, crosses, and the allocation of each experiment "
     "to a figure.** Three lines were used, maintained under specific-pathogen-free "
     "conditions and genotyped by PCR from tail DNA using The Jackson Laboratory published "
     "protocols: the constitutive epithelial driver K14-cre (Jackson Laboratory stock 016230), "
     "the tamoxifen-inducible driver K14-creER (stock 005107), and the floxed allele "
     "*Smad4*^fl/fl^ (stock 017462). Two crosses were made, K14-cre × *Smad4*^fl/fl^ and "
     "K14-creER × *Smad4*^fl/fl^, and in each the control (CON) is a Cre-negative *Smad4*^fl/fl^ "
     "littermate from the same litter — which is why control libraries also contain loxP reads "
     "(floxed exon: Supplementary Figure S5). The constitutive arm uses no tamoxifen and feeds "
     "Figure 1 alone. In the inducible arm 4-hydroxytamoxifen (90 µmol/L, 24 h) was added to the "
     "culture medium of explanted germs of *both* genotypes, so it is not a variable in the "
     "comparison; that arm feeds Supplementary Figures S1 and S2, Figures 2 and 3, and the RNA "
     "sequencing behind Figures 4–6 and Supplementary Figures S3 and S5. The E16.5 dissection "
     "in the Methods feeds no figure and is not drawn. Group sizes are given in the "
     "individual figure legends. "
     "*Two points marked on the panel are unresolved: † the current Figure 1 artwork is labelled "
     "K14-Cre^ER^;Smad4^f/f^, contradicting its own legend, the Results and the Methods; and "
     "‡ the inhibitor panels (Figure 3C–E) are placed in the K14-creER arm by the Methods but "
     "recorded as wild-type germs in the original submission. Both are drawn as the Methods "
     "state them; neither is settled here.*"),
]

DATA_TABLE = """
---

## Supplementary Data files

| File | Contents |
|---|---|
| `Supplementary_Data_S1_DEGs_563.csv` | All 563 differentially expressed genes (adjusted *P* < 0.05, absolute log2 fold change > 1) with DESeq2 statistics and human homolog annotation. |
| `Supplementary_Data_S2_AI_gene_directions.csv` | Log2 fold changes and adjusted *P* values for the 31 GREEN-rated genes of the Genomics England PanelApp "Amelogenesis imperfecta" panel (269 v4.37). |
| `Supplementary_Data_S3_deconvolution_proportions.csv` | Per-sample cell-type proportions estimated by non-negative least squares. |
| `Supplementary_Data_S4_program_celltype_specificity.csv` | Cell-type specificity of each gene program in the single-cell reference. |
| `Supplementary_Data_S5_module_scores.csv` | Per-sample program scores for the four gene programs. |

Complete gene-level count and transcripts-per-million matrices for all samples are deposited with the
sequencing data in the Gene Expression Omnibus (see Data Availability).
"""


# Panels this pipeline generates and therefore owns: staged copies of these are
# refreshed from scripts/JDR_reanalysis/figures/ whenever they differ, so a figure and
# its legend cannot drift apart inside the submission package.
#
# S4A/S4B were held out of this set until 2026-07-29 on the grounds that the difference
# "predates this script". Resolved: the generator is the source of truth, and the staged
# copies were stale in two ways that mattered. They carried the SUPERSEDED enrichment
# statistic "OR ~ 9, P = 0.010" - retracted on 2026-07-27 in favour of OR 9.21,
# P = 1.7e-4 - so the journal package held a withdrawn number; and they were the old
# horizontally-spreading layout (5111 x 1757 px) that renders at ~4 pt on a portrait
# page, which is why 14_workflow_diagram.py was rewritten to draw vertical chains.
# Wet-laboratory panels (S1, S2) stay out of this set permanently: nothing in this
# pipeline generates them, so they are reported on divergence, never overwritten.
GENERATED = {"FigureS6_breeding_scheme.png",
             "FigureS4A_pipeline.png", "FigureS4B_reanalysis.png"}

S5_PNG = "FigureS5_smad4_loxp.png"
S5_LEGEND = "FigureS5_smad4_loxp_legend.md"


def figure_S5():
    """Figure S5 plus the legend generated by step 21, never a retyped copy.

    Every number in that legend is derived at run time from the analysis
    outputs, so retyping it here would create a second, unversioned copy that
    could silently disagree with the figure. It is read verbatim instead, with
    only its H1 title and its build-provenance footer removed and its headings
    demoted to fit under "Supplementary Figures".
    """
    src = os.path.join(FIGSRC, S5_LEGEND)
    if not os.path.exists(src):
        raise FileNotFoundError(
            f"{src} is missing. Run scripts/JDR_reanalysis/21_smad4_loxp_figure.py; "
            "this script will not substitute a hand-written legend.")
    dst = os.path.join(SUP, S5_PNG)
    if not os.path.exists(dst) or (os.path.getmtime(os.path.join(FIGSRC, S5_PNG))
                                   > os.path.getmtime(dst)):
        shutil.copy2(os.path.join(FIGSRC, S5_PNG), dst)

    txt = open(src).read()
    txt = txt.split("\n", 1)[1] if txt.startswith("# ") else txt   # drop the H1
    txt = txt.split("\n---\n")[0]                                  # drop the footer
    # bold lead-in instead of a heading level that would enter the ToC
    txt = re.sub(r"^#+\s*", "", txt, flags=re.M)
    # The legend carries its own bold headline, so prepend only the figure number.
    # Prepending a second headline here would duplicate it, and would also freeze a
    # wording that the figure script is free to change.
    body = "**Supplementary Figure S5.** " + txt.strip()
    # Guard on MEANING, not on an exact phrase: the panel exists to show that RNA-seq
    # cannot settle excision, and that the excision product was not observed. Check both
    # are still stated, so a legend rewrite cannot silently drop the point of the figure.
    missing = [need for need, alts in (
        ("the limit it demonstrates", ("cannot determine", "cannot decide", "cannot settle")),
        ("the zero result",           ("0 skip reads", "0 reads", "not detected")),
    ) if not any(a in body for a in alts)]
    if missing:
        raise RuntimeError(
            "the generated Figure S5 legend no longer states " + " and ".join(missing) +
            "; refusing to embed it. Figure S5 exists to show that RNA-seq alone cannot "
            "determine whether the floxed exon was excised - if that is no longer in the "
            "legend, fix the legend rather than this guard.")
    return f"![]({S5_PNG}){{width=6.5in}}\n\n{body}\n\n"


# RT-qPCR primers, transcribed verbatim from the wet-laboratory record
# (BBRC submission package, "Supplementary Table 1"). Held here rather than in a
# separate document so the table and its exon annotation cannot drift apart.
QPCR_PRIMERS = [
    ("Amelx",  "ACGGCTGCCTTATCATGCTCT",   "ATGGGGACCTGGATTTTGTTT"),
    ("Enam",   "TGCAGAAATCCGACTTCTCCT",   "CATCTGGAATGGCATGGCA"),
    ("Ambn",   "GAATACTGAGAAAGTGCGTTGA",  "TGATCCTGTTCCTGTCCCTAG"),
    ("Klk4",   "CGGGAGTCTTGGTGCATCC",     "CTTGGGAGCCTTTCAGGTTATG"),
    ("Mmp20",  "CTGGGAAGTTGGACCAGAATAC",  "AGGCGGTAGTTAGCCACATCA"),
    ("Smad4",  "ACACCAACAAGTAACGATGCC",   "GCAAAGGTTTCACTTTCCCCA"),
    ("Odam",   "TCTCGCAGCGTCTATTGTCT",    "AGGAATCCAGGAGTTGAAGG"),
    ("Amtn",   "ATCAGCCCAGTCATTACCAAAG",  "AGGTCTGACCCCAGAGTGAG"),
    ("Gapdh",  "AGGTCGGTGTGAACGGATTTG",   "GGGGTCGTTGATGGCAACA"),
]


def table_S3() -> str:
    """Supplementary Table S3 — RT-qPCR primers.

    The Smad4 row carries an explicit note on where the amplicon sits relative to the
    floxed exon. That placement is the reason gene-level counts and this assay both
    report total transcript rather than deletion, so it is stated in the table itself
    rather than left for a reader to work out.
    """
    rows = "\n".join(
        f"| *{g}* | {f} | {r} |" for g, f, r in QPCR_PRIMERS)
    return (
        "\n\n**Supplementary Table S3. Primer sequences used for RT-qPCR.** "
        "Sequences are as used in the wet-laboratory experiments of Figures 2 and 3. "
        "The *Smad4* pair amplifies a 63-base product lying within exon 2 "
        "(chromosome 18: 73,810,876–73,810,938, GRCm39), which is retained by the "
        "conditional allele. It therefore quantifies total *Smad4* transcript and is "
        "not positioned to detect excision of the floxed exon "
        "(chromosome 18: 73,782,797–73,782,980); see Supplementary Figure S5. "
        "*Gapdh* served as the normalisation reference.\n\n"
        "| Gene | Forward primer (5′ → 3′) | Reverse primer (5′ → 3′) |\n"
        "|---|---|---|\n" + rows + "\n")


def main():
    src = open(os.path.join(SUP, "Supplementary_Methods_Pipeline.md")).read()

    # keep Parts A and B; drop the Mermaid workflow section and the trailing S3 legend
    body = src.split("## Part C")[0]
    body = body.split("# Supplementary Methods — Complete Bioinformatics Pipeline", 1)[-1]
    body = body.replace("*Epithelial Smad4 loss reprograms ameloblast identity and metabolism*", "", 1)
    body = "## Supplementary Methods\n" + body.lstrip("\n")
    # demote heading levels so everything nests under the single document title
    body = re.sub(r"^## Part ", "### Part ", body, flags=re.M)
    body = re.sub(r"^### (A\d|B\d)\.", r"#### \1.", body, flags=re.M)

    def emit(entry):
        _, fname, width, legend = entry
        dst, src = os.path.join(SUP, fname), os.path.join(FIGSRC, fname)
        # Panels listed in GENERATED are re-staged whenever the pipeline's own output is
        # newer, so a re-rendered figure cannot be silently left behind while its legend
        # updates around it. Everything else is staged only when the package has no copy:
        # the wet-laboratory panels (S1, S2) exist only in the package, and the staged
        # S4A/S4B currently differ from the pipeline's output for reasons that predate
        # this script and are not its to decide. Those two are reported, not overwritten.
        if not os.path.exists(dst):
            shutil.copy2(src, dst)
        elif os.path.exists(src) and not filecmp.cmp(src, dst, shallow=False):
            if fname in GENERATED:
                shutil.copy2(src, dst)
            else:
                print(f"  NOTE: {fname} in the submission package differs from "
                      f"{os.path.relpath(src, ROOT)} and was left as it is")
        return f"![]({fname}){{width={width}}}\n\n{legend}\n\n"

    # S5 sits between the hand-written legends because its own legend is generated by
    # step 21 and read from disk; anything numbered after S5 is emitted after that call.
    def number(entry) -> int:
        return int(re.match(r"S(\d+)", entry[0]).group(1))

    figs = "\n---\n\n## Supplementary Figures\n\n"
    figs += "".join(emit(e) for e in FIGURES if number(e) < 5)
    figs += figure_S5()
    figs += "".join(emit(e) for e in FIGURES if number(e) > 5)

    tables = (
        "\n---\n\n## Supplementary Tables\n\n"
        "**Supplementary Table S1. Provenance audit of the retired hand-curated "
        "amelogenesis-imperfecta gene list.** The preliminary analysis used a 20-gene list "
        "hand-curated with OMIM and PMID annotations. On pre-submission audit (27 July 2026) "
        "that annotation was found unreliable: " + prov.summary().split('PANEL')[0].strip().rstrip('.') +
        ". Every identifier below was checked against NCBI E-utilities (PMIDs) and mim2gene_medgen "
        "(MIM class). Because of this, the analysis was moved to the versioned external panel "
        "documented in Supplementary Table S2. Genes marked “on panel: no” are not "
        "GREEN-rated on that panel and are not used here.\n\n"
        + prov.audit_markdown() + prov.audit_footnote() +
        "\n\n**Supplementary Table S2. Verified provenance of the 31 genes used in the "
        "amelogenesis-imperfecta enrichment.** Genes are the GREEN-rated (high-evidence) members of "
        "the Genomics England PanelApp “Amelogenesis imperfecta” panel (269 v4.37; retrieved "
        "27 July 2026). Each OMIM phenotype MIM was confirmed to be a phenotype (not gene) entry "
        "mapping to the stated gene, and each primary PMID was confirmed to resolve to a real report "
        "of that gene’s enamel or amelogenesis-imperfecta phenotype. Mouse orthologs are as "
        "carried in the Ensembl release 116 annotation. Two entries are flagged: PEX26 has no "
        "enamel-specific publication (the cited paper is the peroxisome-biogenesis-disorder primary "
        "report), and PLXNB2 has no OMIM phenotype MIM assigned yet.\n\n"
        + prov.panel_markdown() + "\n" + prov.panel_footnote() + "\n"
        + table_S3()
    )
    open(BUILD, "w").write(TITLE + body.rstrip() + "\n" + figs + tables + DATA_TABLE)
    subprocess.run(["pandoc", os.path.basename(BUILD), "-o", OUT,
                    "--toc-depth=3"], cwd=SUP, check=True)
    os.remove(BUILD)
    print("wrote", os.path.join(SUP, OUT))


if __name__ == "__main__":
    main()
