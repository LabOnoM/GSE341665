"""
ai_gene_provenance.py — verified provenance for the amelogenesis-imperfecta (AI) gene
panel, and the audit trail of the retired hand-curated list.

Two datasets:

  AUDIT  — the 20 genes of the retired hand-typed panel
           (scripts/_quarantine/00_Mendelian_AI_Genes.R). For each: the OMIM MIM it
           claimed and what that MIM actually is, the PMID it claimed and what that PMID
           actually cites, and the verified correct values. This documents the
           2026-07-27 integrity audit that motivated switching to a versioned external
           panel.

  PANEL  — the 31 GREEN-rated genes of the Genomics England PanelApp "Amelogenesis
           imperfecta" panel (269 v4.37) actually used in the enrichment, each with its
           verified OMIM phenotype MIM, inheritance, mouse ortholog, and a verified
           primary PMID.

EVERY identifier below was confirmed on 2026-07-27 against authoritative sources:
  * OMIM phenotype-vs-gene MIM class — NCBI mim2gene_medgen (each MIM confirmed
    type == "phenotype", mapped to the stated gene);
  * PMID existence and topic — NCBI E-utilities esummary (each PMID confirmed to
    resolve to a real paper on the gene's enamel/AI phenotype);
  * mouse ortholog — MGI / NCBI Gene, as carried in the Ensembl-116 annotation.

Nothing here is recalled or guessed. Where evidence is genuinely incomplete it is
flagged in `note` (PEX26 has no enamel-specific paper; PLXNB2 has no OMIM phenotype MIM
yet; SP6's human evidence is recent and limited).

This module is imported by 15_build_supplementary_docx.py to render Supplementary
Tables S1 (audit) and S2 (verified panel). It is the single source for those tables.
"""

from __future__ import annotations

from dataclasses import dataclass


# ---------------------------------------------------------------------------
# Supplementary Table S1 — audit of the retired 20-gene hand-curated list
# ---------------------------------------------------------------------------
@dataclass(frozen=True)
class AuditRow:
    gene: str
    claimed_mim: str
    claimed_mim_actual: str      # what that MIM number actually is
    correct_mim: str             # verified phenotype MIM (or note)
    claimed_pmid: str
    claimed_pmid_actual: str     # what that PMID actually cites
    verified_pmid: str           # verified correct primary (or note)
    in_panel: bool               # retained on PanelApp GREEN (used) vs excluded


AUDIT: tuple[AuditRow, ...] = (
    AuditRow("AMELX",  "301200", "phenotype (correct)",           "301200",
             "8651276",  "Apert syndrome (FGFR2)",                "1916828", True),
    AuditRow("AMBN",   "616270", "phenotype (correct)",           "616270",
             "25327895", "Atlantic salmon conservation genomics", "24858907", True),
    AuditRow("ENAM",   "104500", "phenotype (correct)",           "104500",
             "11992261", "Noonan syndrome (PTPN11)",              "11487571", True),
    AuditRow("AMTN",   "617607", "phenotype (correct)",           "617607",
             "28861858", "does not exist in PubMed",              "27412008", False),
    AuditRow("MMP20",  "612529", "phenotype (correct)",           "612529",
             "19187894", "chronic anterior uveitis",              "15744043", True),
    AuditRow("KLK4",   "204700", "phenotype (correct)",           "204700",
             "14756670", "familial hypercholesterolemia",         "15235027", True),
    AuditRow("FAM83H", "130900", "phenotype (correct)",           "130900",
             "18327258", "THM1/hedgehog (mouse cilia)",           "18252228", True),
    AuditRow("ACP4",   "617297", "phenotype (correct)",           "617297",
             "28111907", "palaeontology essay",                   "27843125", True),
    AuditRow("SLC4A4", "(none)", "no MIM given",                  "604278",
             "15452090", "metipranolol / retina",                 "35260236", False),
    AuditRow("STIM1",  "612559", "different gene (leukaemia susceptibility 5)", "612783",
             "25193870", "JAK3 / T-ALL",                          "19420366", True),
    AuditRow("ORAI1",  "615887", "different gene (SLC24A4 AI locus)",           "612782",
             "24349437", "cre recombinase (retina)",              "20004786", True),
    AuditRow("CFTR",   "219700", "phenotype (correct)",           "219700",
             "24799079", "knee arthroplasty alignment",           "— (excluded)", False),
    AuditRow("LAMC2",  "150292", "gene MIM (not phenotype)",      "226700",
             "10798435", "congenital heart surgery nomenclature", "— (excluded)", False),
    AuditRow("LAMB3",  "150310", "gene MIM (not phenotype)",      "104530",
             "20034080", "HHS gene-patent news item",             "23958762", True),
    AuditRow("LAMA3",  "600805", "gene MIM (not phenotype)",      "619786",
             "36326426", "correct (LAMA3-associated AI)",         "36326426", True),
    AuditRow("COL17A1","113811", "gene MIM (not phenotype)",      "619787",
             "16685652", "POLG2 / progressive ophthalmoplegia",   "37979963", True),
    AuditRow("DLX3",   "600525", "gene MIM (not phenotype)",      "190320",
             "12649808", "HOXD13 brachydactyly",                  "9467018", True),
    AuditRow("PITX2",  "137600", "wrong disease (ant. seg. dysgenesis 4)",     "180500",
             "10508521", "CRB1 / retinitis pigmentosa",           "— (excluded)", False),
    AuditRow("WNT10A", "150400", "wrong disease (tooth agenesis 4)",           "257980",
             "17952058", "melanopsin / circadian clock",          "— (excluded)", False),
    AuditRow("KMT2D",  "147920", "phenotype (correct)",           "147920",
             "40103013", "correct (KMT2D enamel development)",     "— (excluded)", False),
)


# ---------------------------------------------------------------------------
# Supplementary Table S2 — verified provenance of the 31 PanelApp GREEN genes
# ---------------------------------------------------------------------------
@dataclass(frozen=True)
class PanelRow:
    human: str
    mouse: str
    phenotype_mim: str      # "—" where no OMIM phenotype MIM exists yet
    phenotype: str
    inheritance: str        # AD | AR | XL
    pmid: str
    citation: str           # FirstAuthor, Journal Year
    note: str = ""


PANEL: tuple[PanelRow, ...] = (
    PanelRow("ACP4","Acp4","617297","Amelogenesis imperfecta, type IJ","AR",
             "27843125","Seymen, Am J Hum Genet 2016"),
    PanelRow("AMBN","Ambn","616270","Amelogenesis imperfecta, type IF (hypoplastic)","AR",
             "24858907","Poulter, Hum Mol Genet 2014"),
    PanelRow("AMELX","Amelx","301200","Amelogenesis imperfecta, type IE (X-linked)","XL",
             "1916828","Lagerström, Genomics 1991"),
    PanelRow("C4orf26","Odaph","614832","Amelogenesis imperfecta, type IIA4 (hypomaturation)","AR",
             "22901946","Parry, Am J Hum Genet 2012","gene renamed ODAPH"),
    PanelRow("CNNM4","Cnnm4","217080","Jalili syndrome (cone-rod dystrophy + AI)","AR",
             "19200525","Parry, Am J Hum Genet 2009","syndromic"),
    PanelRow("COL17A1","Col17a1","619787","Amelogenesis imperfecta; junctional EB spectrum","AD",
             "37979963","Hany, J Med Genet 2024","also causes junctional EB"),
    PanelRow("DLX3","Dlx3","190320","Tricho-dento-osseous syndrome","AD",
             "9467018","Price, Hum Mol Genet 1998","syndromic"),
    PanelRow("ENAM","Enam","104500","Amelogenesis imperfecta, type IB (hypoplastic)","AD",
             "11487571","Rajpar, Hum Mol Genet 2001"),
    PanelRow("FAM20A","Fam20a","204690","Amelogenesis imperfecta, type IG (enamel-renal syndrome)","AR",
             "21549343","O'Sullivan, Am J Hum Genet 2011"),
    PanelRow("FAM20C","Fam20c","259775","Raine syndrome (with enamel/dental defects)","AR",
             "17924334","Simpson, Am J Hum Genet 2007","syndromic; enamel defects in survivors"),
    PanelRow("FAM83H","Fam83h","130900","Amelogenesis imperfecta, type IIIA (hypocalcified)","AD",
             "18252228","Kim JW, Am J Hum Genet 2008","gene renamed SACK1H (2026)"),
    PanelRow("GPR68","Gpr68","617217","Amelogenesis imperfecta, type IIA6 (hypomaturation)","AR",
             "27693231","Parry, Am J Hum Genet 2016"),
    PanelRow("ITGB6","Itgb6","616221","Amelogenesis imperfecta, type IH (hypoplastic)","AR",
             "24305999","Wang SK, Hum Mol Genet 2014"),
    PanelRow("KLK4","Klk4","204700","Amelogenesis imperfecta, type IIA1 (hypomaturation)","AR",
             "15235027","Hart PS, J Med Genet 2004"),
    PanelRow("LAMA3","Lama3","619786","Junctional epidermolysis bullosa; AI-spectrum enamel defect","AR",
             "36326426","Wang SK, Oral Dis 2023","also causes junctional EB"),
    PanelRow("LAMB3","Lamb3","104530","Amelogenesis imperfecta, type IA (hypoplastic)","AD",
             "23958762","Kim JW, J Dent Res 2013","also causes junctional EB"),
    PanelRow("LTBP3","Ltbp3","601216","Dental anomalies and short stature (DASS)","AR",
             "19344874","Noor, Am J Hum Genet 2009","syndromic"),
    PanelRow("MMP20","Mmp20","612529","Amelogenesis imperfecta, type IIA2 (hypomaturation)","AR",
             "15744043","Kim JW, J Med Genet 2005"),
    PanelRow("ORAI1","Orai1","612782","CRAC channelopathy (immunodeficiency 9) with defective enamel","AR",
             "20004786","McCarl, J Allergy Clin Immunol 2009","syndromic (recessive LOF)"),
    PanelRow("PEX1","Pex1","234580","Heimler syndrome 1 (AI + hearing loss)","AR",
             "26387595","Ratbi, Am J Hum Genet 2015","syndromic"),
    PanelRow("PEX26","Pex26","614872","Peroxisome biogenesis disorder 7A (Zellweger spectrum)","AR",
             "12851857","Matsumoto, Am J Hum Genet 2003",
             "PBD primary report; no enamel-specific paper — enamel inferred from spectrum"),
    PanelRow("PEX6","Pex6","616617","Heimler syndrome 2 (AI + hearing loss)","AR",
             "26387595","Ratbi, Am J Hum Genet 2015","syndromic"),
    PanelRow("PLXNB2","Plxnb2","—","AI, hearing loss and intellectual disability","AR",
             "38458752","Smith CEL, J Med Genet 2024","recently reported; no OMIM phenotype MIM yet"),
    PanelRow("RELT","Relt","618386","Amelogenesis imperfecta, type IIIC","AR",
             "30506946","Kim JW, Clin Genet 2019"),
    PanelRow("ROGDI","Rogdi","226750","Kohlschütter-Tönz syndrome (epilepsy + AI)","AR",
             "22424600","Schossig, Am J Hum Genet 2012","syndromic"),
    PanelRow("SLC10A7","Slc10a7","618363","Short stature, AI and skeletal dysplasia (SSASKS)","AR",
             "30082715","Dubail, Nat Commun 2018","syndromic"),
    PanelRow("SLC13A5","Slc13a5","615905","Developmental and epileptic encephalopathy 25 with AI","AR",
             "27600704","Schossig, J Med Genet 2017","syndromic"),
    PanelRow("SLC24A4","Slc24a4","615887","Amelogenesis imperfecta, type IIA5 (hypomaturation)","AR",
             "23375655","Parry, Am J Hum Genet 2013"),
    PanelRow("SP6","Sp6","620104","Amelogenesis imperfecta, type IK (hypoplastic)","AD",
             "32167558","Smith CEL, Hum Mol Genet 2020","human evidence recent and limited"),
    PanelRow("STIM1","Stim1","612783","CRAC channelopathy (immunodeficiency 10) with enamel defect","AR",
             "19420366","Picard, N Engl J Med 2009","syndromic (recessive LOF)"),
    PanelRow("WDR72","Wdr72","613211","Amelogenesis imperfecta, type IIA3 (hypomaturation)","AR",
             "19853237","El-Sayed, Am J Hum Genet 2009"),
)


def _md_escape(s: str) -> str:
    return s.replace("|", "\\|")


def audit_markdown() -> str:
    """Render Supplementary Table S1 as a Markdown table."""
    h = ("| Gene | Claimed MIM | What that MIM is | Correct phenotype MIM | "
         "Claimed PMID | What that PMID cites | Verified PMID | On panel |\n"
         "|---|---|---|---|---|---|---|---|\n")
    rows = ""
    for r in AUDIT:
        rows += "| {} | {} | {} | {} | {} | {} | {} | {} |\n".format(
            r.gene, r.claimed_mim, _md_escape(r.claimed_mim_actual), r.correct_mim,
            r.claimed_pmid, _md_escape(r.claimed_pmid_actual), r.verified_pmid,
            "yes" if r.in_panel else "no")
    return h + rows


def panel_markdown() -> str:
    """Render Supplementary Table S2 as a Markdown table."""
    h = ("| Gene (human) | Mouse | OMIM phenotype MIM | Phenotype | Inheritance | "
         "Primary reference (PMID) | Notes |\n"
         "|---|---|---|---|---|---|---|\n")
    rows = ""
    for r in PANEL:
        cite = f"{r.citation} (PMID {r.pmid})"
        rows += "| {} | *{}* | {} | {} | {} | {} | {} |\n".format(
            r.human, r.mouse, r.phenotype_mim, _md_escape(r.phenotype),
            r.inheritance, cite, _md_escape(r.note))
    return h + rows


def panel_footnote() -> str:
    """Footnote defining every abbreviation that appears in Supplementary Table S2.

    Kept next to the data so the two cannot drift; the abbreviation list was derived by
    extracting every acronym actually present in the rendered table, not written from
    memory.
    """
    return """
**Table S2 footnote — abbreviations and conventions.**

*Columns.* **Gene (human)** — HGNC-approved human gene symbol, as listed by PanelApp.
**Mouse** — the MGI-approved mouse ortholog symbol (italic), as carried in the Ensembl
release 116 annotation used for the differential-expression analysis; this is the symbol
on which genes were matched. **OMIM phenotype MIM** — the Mendelian Inheritance in Man
(OMIM) identifier of the *phenotype* entry, not the gene entry; each was confirmed to be
of type "phenotype" and to map to the stated gene in NCBI `mim2gene_medgen`.
**Phenotype** — the OMIM phenotype designation. **Inheritance** — mode of inheritance for
that phenotype. **Primary reference (PMID)** — the report establishing the gene-phenotype
link, with its PubMed identifier (PMID); each was confirmed to exist and to concern the
gene's enamel or amelogenesis-imperfecta phenotype. **Notes** — caveats affecting how the
row should be read.

*Inheritance codes.* **AD**, autosomal dominant; **AR**, autosomal recessive; **XL**,
X-linked.

*Disease and phenotype abbreviations.* **AI**, amelogenesis imperfecta; **EB**,
epidermolysis bullosa (**junctional EB**, the junctional subtype); **CRAC**, calcium
release-activated calcium channel; **DASS**, dental anomalies and short stature;
**SSASKS**, short stature, amelogenesis imperfecta and skeletal dysplasia with scoliosis;
**PBD**, peroxisome biogenesis disorder; **LOF**, loss of function.

*Amelogenesis-imperfecta subtype nomenclature.* OMIM numbers AI subtypes by a Roman
numeral for the defect class followed by a letter (and, where several genes cause the same
class, a digit): **type I** = hypoplastic (too little enamel formed, of normal
mineralisation); **type II** = hypomaturation (full-thickness but incompletely matured,
soft enamel); **type III** = hypocalcified (full-thickness but poorly mineralised). Hence
type IA, IB, IE, IF, IG, IH, IJ and IK are hypoplastic forms; types IIA1-IIA6 are
hypomaturation forms; and types IIIA and IIIC are hypocalcified forms. Where a parenthetical
descriptor follows the subtype it restates this class for readability; rows without one
follow the same Roman-numeral rule.

*Gene-symbol changes.* Two genes are listed under the symbol PanelApp uses, with the
current HGNC symbol given in Notes: **C4orf26** is now **ODAPH** (odontogenesis-associated
phosphoprotein), and **FAM83H** was renamed **SACK1H** (scaffolding CK1 anchoring protein
H) in February 2026. The mouse column retains *Odaph* and *Fam83h* because those are the
symbols present in the Ensembl release 116 annotation on which matching was performed.

*Symbols.* An em dash (**—**) in the OMIM column means no phenotype MIM has been assigned
yet; this applies only to PLXNB2, whose AI phenotype was first reported in 2024.

*Note-column meanings.* **"syndromic"** — enamel involvement occurs as one feature of a
broader syndrome rather than as isolated AI. **"also causes junctional EB"** — the same
gene causes both an AI-spectrum enamel defect and junctional epidermolysis bullosa.
**"recessive LOF"** — the enamel phenotype belongs to the recessive loss-of-function
disorder, not to the dominant gain-of-function disorder caused by the same gene.
**"PBD primary report; no enamel-specific paper"** — for PEX26 the cited reference
establishes the peroxisome biogenesis disorder; enamel involvement is inferred from the
disorder spectrum rather than from a dedicated study. **"human evidence recent and
limited"** — for SP6 the human AI association rests on one dominant family with functional
corroboration plus one de novo case.

*Panel definition.* All 31 genes are the GREEN-rated ("high evidence", diagnostic-grade)
members of the Genomics England PanelApp panel 269, "Amelogenesis imperfecta", version
4.37, retrieved 27 July 2026. Amber- and red-rated genes on the same panel were excluded
and are named in Supplementary Methods B6.
"""


def summary() -> str:
    n_wrong_pmid = sum(1 for r in AUDIT if "correct" not in r.claimed_pmid_actual)
    n_wrong_mim = sum(1 for r in AUDIT if "correct" not in r.claimed_mim_actual)
    return (f"AUDIT: {len(AUDIT)} candidate genes; {n_wrong_pmid}/{len(AUDIT)} claimed PMIDs "
            f"cite unrelated work; {n_wrong_mim}/{len(AUDIT)} claimed MIMs wrong. "
            f"PANEL: {len(PANEL)} PanelApp GREEN genes, all identifiers verified.")


if __name__ == "__main__":
    print(summary())
    assert len({r.human for r in PANEL}) == 31, "panel must have 31 distinct genes"
    assert len({r.gene for r in AUDIT}) == 20, "audit must have 20 distinct genes"
    print("\n### Supplementary Table S1 — audit\n")
    print(audit_markdown())
    print("\n### Supplementary Table S2 — verified panel\n")
    print(panel_markdown())
