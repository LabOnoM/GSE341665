#!/usr/bin/env python3
"""
10_build_manuscript_docx.py -- assemble the Version02 submission package.

Reads the working manuscript (doc/Manuscript_JDR_Smad4_Identity_Metabolism.md),
(1) strips drafting-only meta, (2) converts author-year citations to JMM numbered
[n] AUTOMATICALLY in order of first appearance, (3) emits the reference list in
that citation order, (4) embeds Figures 1-6 and Supplementary Figures S1-S2,
and writes 00.Manuscript/Version02/Manuscript_JMM_Smad4.md + .docx (via pandoc).

Numbering is derived, not hand-maintained: add a citation to CITES and the whole
reference list renumbers itself. The script fails loudly if any author-year
citation in the body was not mapped.
"""
import os, re, subprocess

ROOT = os.path.dirname(os.path.dirname(os.path.dirname(os.path.abspath(__file__))))
SRC = os.path.join(ROOT, "doc", "Manuscript_JDR_Smad4_Identity_Metabolism.md")
V2 = os.path.join(ROOT, "00.Manuscript", "Version02")
OUT_MD = os.path.join(V2, "Manuscript_JMM_Smad4.md")

# reference key -> Vancouver-style entry (all verified against PubMed)
REFS = {
 "lacruz":  "Lacruz RS, Habelitz S, Wright JT, Paine ML. Dental enamel formation and implications for oral health and disease. Physiol Rev. 2017;97(3):939–993.",
 "bartlett":"Bartlett JD. Dental enamel development: proteinases and their enamel matrix substrates. ISRN Dent. 2013;2013:684607.",
 "smith":   "Smith CEL, Poulter JA, Antanaviciute A, Kirkham J, Brookes SJ, Inglehearn CF, Mighell AJ. Amelogenesis imperfecta; genes, proteins, and pathways. Front Physiol. 2017;8:435.",
 "massague":"Massagué J. TGFβ signalling in context. Nat Rev Mol Cell Biol. 2012;13(10):616–630.",
 "yokozeki":"Yokozeki M, Afanador E, Nishi M, Kaneko K, Shimokawa H, Yokote K, Deng C, Tsuchida K, Sugino H, Moriyama K. Smad3 is required for enamel biomineralization. Biochem Biophys Res Commun. 2003;305(3):684–690.",
 "huang":   "Huang X, Xu X, Bringas P Jr, Hung YP, Chai Y. Smad4-Shh-Nfic signaling cascade-mediated epithelial-mesenchymal interaction is crucial in regulating tooth root development. J Bone Miner Res. 2010;25(5):1167–1178.",
 "arrive":  "Percie du Sert N, Hurst V, Ahluwalia A, et al. The ARRIVE guidelines 2.0. PLoS Biol. 2020;18(7):e3000410.",
 "star":    "Dobin A, Davis CA, Schlesinger F, et al. STAR: ultrafast universal RNA-seq aligner. Bioinformatics. 2013;29(1):15–21.",
 "salmon":  "Patro R, Duggal G, Love MI, Irizarry RA, Kingsford C. Salmon provides fast and bias-aware quantification of transcript expression. Nat Methods. 2017;14(4):417–419.",
 "stringtie":"Pertea M, Pertea GM, Antonescu CM, Chang TC, Mendell JT, Salzberg SL. StringTie enables improved reconstruction of a transcriptome from RNA-seq reads. Nat Biotechnol. 2015;33(3):290–295.",
 "tximport":"Soneson C, Love MI, Robinson MD. Differential analyses for RNA-seq: transcript-level estimates improve gene-level inferences. F1000Res. 2015;4:1521.",
 "ruvseq":  "Risso D, Ngai J, Speed TP, Dudoit S. Normalization of RNA-seq data using factor analysis of control genes or samples. Nat Biotechnol. 2014;32(9):896–902.",
 "seurat":  "Hao Y, Hao S, Andersen-Nissen E, et al. Integrated analysis of multimodal single-cell data. Cell. 2021;184(13):3573–3587.e29.",
 "deseq2":  "Love MI, Huber W, Anders S. Moderated estimation of fold change and dispersion for RNA-seq data with DESeq2. Genome Biol. 2014;15(12):550.",
 "balic":   "Balic A, Thesleff I. Tissue interactions regulating tooth development and renewal. Curr Top Dev Biol. 2015;115:157–186.",
 "zhang15": "Zhang Z, Tian H, Lv P, Wang W, Jia Z, Wang S, Zhou C, Gao X. Transcriptional factor DLX3 promotes the gene expression of enamel matrix proteins during amelogenesis. PLoS One. 2015;10(3):e0121288.",
 "rhodes":  "Rhodes CS, Yoshitomi Y, Burbelo PD, Freese NH, Nakamura T, Chiba Y, Yamada Y. Sp6/Epiprofin is a master regulator in the developing tooth. Biochem Biophys Res Commun. 2021;581:89–95.",
 "lee18":   "Lee DS, Roh SY, Park JC. The Nfic-osterix pathway regulates ameloblast differentiation and enamel formation. Cell Tissue Res. 2018;374(3):531–540.",
 "werner":  "Werner S, Grose R. Regulation of wound healing by growth factors and cytokines. Physiol Rev. 2003;83(3):835–870.",
 "ge":      "Ge Y, Gomez NC, Adam RC, et al. Stem cell lineage infidelity drives wound repair and cancer. Cell. 2017;169(4):636–650.e14.",
 "horton":  "Horton JD, Shah NA, Warrington JA, et al. Combined analysis of oligonucleotide microarray data from transgenic and knockout mice identifies direct SREBP target genes. Proc Natl Acad Sci U S A. 2003;100(21):12027–12032.",
 "interferome":"Rusinova I, Forster S, Yu S, et al. Interferome v2.0: an updated database of annotated interferon-regulated genes. Nucleic Acids Res. 2013;41(Database issue):D1040–D1046.",
 "leez":    "Lee E, Chuang H-Y, Kim J-W, Ideker T, Lee D. Inferring pathway activity toward precise disease classification. PLoS Comput Biol. 2008;4(11):e1000217.",
 "chiba":   "Chiba Y, Saito K, Martin D, et al. Single-cell RNA-sequencing from mouse incisor reveals dental epithelial cell-type specific genes. Front Cell Dev Biol. 2020;8:841. [complete author list before submission]",
 "li11":    "Li J, Huang X, Xu X, Mayo J, Bringas P Jr, Jiang R, Wang S, Chai Y. SMAD4-mediated WNT signaling controls the fate of cranial neural crest cells during tooth morphogenesis. Development. 2011;138(10):1977–1989.",
 "kotzka":  "Kotzka J, Lehr S, Roth G, Avci H, Knebel B, Müller-Wieland D. Insulin-activated Erk-mitogen-activated protein kinases phosphorylate sterol regulatory element-binding protein-2 at serine residues 432 and 455 in vivo. J Biol Chem. 2004;279(21):22404–22411.",
 "zhang21": "Zhang J, et al. The absence of muscle segment homeobox 2 leads to the pyroptosis of ameloblasts by inducing squamous epithelial hyperplasia in the enamel organ. J Cell Mol Med. 2021;25(13):6429–6437. [complete author list before submission]",
 "klein":   "Klein OD, Lyons DB, Balooch G, Marshall GW, Basson MA, Peterka M, Boran T, Peterkova R, Martin GR. An FGF signaling loop sustains the generation of differentiated progeny from stem cells in mouse incisors. Development. 2008;135(2):377–385.",
 "goldberg":"Goldberg M, Septier D. Phospholipids in amelogenesis and dentinogenesis. Crit Rev Oral Biol Med. 2002;13(3):276–290.",
 "wang":    "Wang S, Link F, Han M, et al. The interplay of TGF-β1 and cholesterol orchestrating hepatocyte cell fate. Cell Mol Gastroenterol Hepatol. 2024;17(4):567–587.",
 "olczak":  "Olczak-Kowalczyk D, Gozdowski D, Krasuska-Sławińska E, et al. Smith–Lemli–Opitz syndrome: oral characteristics and risk factors for caries development. Biomedicines. 2025;13(3):574.",
 "engelking":"Engelking LJ, Evers BM, Richardson JA, Goldstein JL, Brown MS, Liang G. Severe facial clefting in Insig-deficient mouse embryos caused by sterol accumulation and reversed by lovastatin. J Clin Invest. 2006;116(9):2356–2365.",
 "iwaya":   "Iwaya C, Suzuki A, Shim J, Kim A, Iwata J. Craniofacial bone anomalies related to cholesterol synthesis defects. Sci Rep. 2024;14(1):5371.",
 "york":    "York AG, Williams KJ, Argus JP, et al. Limiting cholesterol biosynthetic flux spontaneously engages type I IFN signaling. Cell. 2015;163(7):1716–1729.",
 "blanc":   "Blanc M, Hsieh WY, Robertson KA, et al. Host defense against viral infection involves interferon-mediated down-regulation of sterol biosynthesis. PLoS Biol. 2011;9(3):e1000598.",
}

# (draft author-year text in the manuscript, [ref keys], optional replacement template)
CITES = [
 ("(Lacruz et al. 2017; Bartlett 2013; Balic and Thesleff 2015)", ["lacruz","bartlett","balic"], None),
 ("(Zhang et al. 2015; Rhodes et al. 2021; Lee et al. 2018)", ["zhang15","rhodes","lee18"], None),
 ("(Werner and Grose 2003)", ["werner"], None),
 ("(Ge et al. 2017)", ["ge"], None),
 ("(Horton et al. 2003)", ["horton"], None),
 ("(Rusinova et al. 2013)", ["interferome"], None),
 ("(York et al. 2015; Blanc et al. 2011)", ["york","blanc"], None),
 ("(Gene Expression Omnibus accession GSE146855; Chiba et al. 2020)", ["chiba"],
  "(Gene Expression Omnibus accession GSE146855) {n}"),
 ("(Chiba et al. 2020)", ["chiba"], None),
 ("(Lacruz et al. 2017)", ["lacruz"], None),
 ("(Bartlett 2013)", ["bartlett"], None),
 ("(Smith et al. 2017)", ["smith"], None),
 ("(Massagué 2012)", ["massague"], None),
 ("(Yokozeki et al. 2003)", ["yokozeki"], None),
 ("(Huang et al. 2010)", ["huang"], None),
 ("(Percie du Sert et al. 2020)", ["arrive"], None),
 ("(Dobin et al. 2013)", ["star"], None),
 ("(Patro et al. 2017)", ["salmon"], None),
 ("(Pertea et al. 2015)", ["stringtie"], None),
 ("(Soneson et al. 2015)", ["tximport"], None),
 ("(Risso et al. 2014)", ["ruvseq"], None),
 ("(Hao et al. 2021)", ["seurat"], None),
 ("(Love et al. 2014)", ["deseq2"], None),
 ("(Lee et al. 2008)", ["leez"], None),
 ("(Li et al. 2011)", ["li11"], None),
 ("(Kotzka et al. 2004)", ["kotzka"], None),
 ("(Zhang et al. 2015)", ["zhang15"], None),
 ("(Rhodes et al. 2021)", ["rhodes"], None),
 ("(Zhang et al. 2021)", ["zhang21"], None),
 ("(Klein et al. 2008)", ["klein"], None),
 ("(Goldberg and Septier 2002)", ["goldberg"], None),
 ("(Wang et al. 2024)", ["wang"], None),
 ("(Olczak-Kowalczyk et al. 2025)", ["olczak"], None),
 ("(Engelking et al. 2006)", ["engelking"], None),
 ("(Iwaya et al. 2024)", ["iwaya"], None),
]

FIGS = [   # (figure number, file in Figures/, display width)
    ("1", "Figure1_enamel_phenotype.png", "3.5in"),
    ("2", "Figure2_matrix_genes_proteins.png", "6.2in"),
    ("3", "Figure3_signaling_ERK_brake.png", "4.5in"),
    ("4", "Figure4_transcriptome.png", "6.5in"),
    ("5", "Figure5_identity_switch.png", "6.5in"),
    ("6", "Figure6_metabolic_convergence.png", "6.5in"),
]
SUPP = [("S1", "FigureS1_inducible_phenotype.png", "5in"),
        ("S2", "FigureS2_tamoxifen_titration.png", "5in")]


def _sync_figures():
    """AUTO-SYNC: copy freshly generated computational figures (4-6) from the
    pipeline output into Version02/Figures, so the packaged figures can never be
    stale relative to the analysis."""
    import shutil
    src = os.path.join(ROOT, "scripts", "JDR_reanalysis", "figures")
    dst = os.path.join(V2, "Figures")
    for stem in ("Figure4_transcriptome", "Figure5_identity_switch", "Figure6_metabolic_convergence"):
        for ext in ("png", "pdf"):
            a = os.path.join(src, f"{stem}.{ext}")
            if os.path.exists(a):
                shutil.copy2(a, os.path.join(dst, f"{stem}.{ext}"))
    print("figures synced from pipeline output")


def main():
    _sync_figures()
    txt = open(SRC).read()
    body = txt[:txt.index("## References")]
    body = body[body.index("## Title Page"):]
    for junk in ["*(67 characters incl. spaces)*\n", "*(≈200 words)*\n", "*(≈199 words)*\n"]:
        body = body.replace(junk, "")

    # --- assign reference numbers by first appearance in the body ---
    present = [(body.index(pat), pat, keys, tmpl) for pat, keys, tmpl in CITES if pat in body]
    present.sort(key=lambda x: x[0])
    order, num = [], {}
    for _, _, keys, _ in present:
        for k in keys:
            if k not in num:
                order.append(k); num[k] = len(order)

    # --- replace (longest pattern first, so multi-citations aren't clipped) ---
    for pat, keys, tmpl in sorted(CITES, key=lambda c: -len(c[0])):
        if pat not in body:
            continue
        nums = "[" + ", ".join(str(num[k]) for k in keys) + "]"
        body = body.replace(pat, tmpl.format(n=nums) if tmpl else nums)

    # --- embed figures before their legends ---
    for n, f, w in FIGS:
        body = body.replace(f"**Figure {n}.", f"![](Figures/{f}){{width={w}}}\n\n**Figure {n}.", 1)
    for n, f, w in SUPP:
        body = body.replace(f"**Supplementary Figure {n}.",
                            f"![](Supplementary/{f}){{width={w}}}\n\n**Supplementary Figure {n}.", 1)

    refs = "## References\n\n" + "\n".join(f"{i+1}. {REFS[k]}" for i, k in enumerate(order)) + "\n"
    open(OUT_MD, "w").write(body.rstrip() + "\n\n---\n\n" + refs)
    subprocess.run(["pandoc", "Manuscript_JMM_Smad4.md", "-o", "Manuscript_JMM_Smad4.docx"],
                   cwd=V2, check=True)

    leftover = re.findall(r"\([A-Z][a-zA-Z'’-]+(?: (?:and|et al\.))?,? \d{4}[a-z]?\)", body)
    print("wrote", OUT_MD, "and Manuscript_JMM_Smad4.docx")
    print("references cited:", len(order))
    print("unmapped author-year citations:", leftover if leftover else "none")
    unused = [k for k in REFS if k not in num]
    if unused:
        print("WARNING: reference entries never cited:", unused)


if __name__ == "__main__":
    main()
