#!/usr/bin/env python3
"""
programs.py -- SINGLE SOURCE OF TRUTH for the four curated gene programs.

Imported by 02_deg_analysis.py, 03_module_scores.py, 06_deconvolution.py and
07_figures.py so that every analysis and every figure scores exactly the same
gene sets. (Previously each script carried its own copy, and the sterol and
interferon lists had silently drifted apart.)

Program membership follows published definitions, cited in the manuscript
Methods:
  IDENTITY - enamel matrix genes and maturation proteases, plus dental-epithelial
             lineage transcription factors.
             Lacruz 2017 (PMID 28468833); Bartlett 2013 (PMID 24159389);
             Balic & Thesleff 2015 (PMID 26589925); Dlx3: Zhang 2015
             (PMID 25815730); Sp6: Rhodes 2021 (PMID 34662808); Nfic-Sp7:
             Lee 2018 (PMID 30091046).
             NOTE: Krt76 is deliberately EXCLUDED. It is strongly down-regulated
             in our data but has no ameloblast literature (the frequently
             mis-cited Duverger 2014 paper concerns KRT75, not KRT76), so it is
             reported in the text as an observation, not as a program member.
  STEROL   - the cholesterol-biosynthetic (mevalonate-sterol) pathway. The nine
             genes marked below are listed as direct SREBP targets in Horton 2003
             Table 1 (PMID 14512514); Srebf2 is the transcription factor itself;
             the remainder are further enzymes of the same pathway.
  WOUND    - canonical epithelial wound/stress response: EGFR-family ligands,
             matrix metalloproteinases, prostaglandin synthase, AP-1 factors.
             Werner & Grose 2003 (PMID 12843410); Ge 2017 (PMID 28434617).
             CAVEAT: Mmp12 is classically a macrophage metalloelastase and may be
             contributed in part by non-epithelial cells.
  IFN      - interferon-stimulated genes as annotated in the Interferome database
             (Rusinova 2013, PMID 23203888).
"""

# direct SREBP targets per Horton 2003 Table 1 (subset of STEROL, used for Methods text)
STEROL_HORTON_DIRECT = [
    "Hmgcr", "Hmgcs1", "Mvd", "Fdps", "Sqle", "Msmo1", "Idi1", "Ldlr", "Insig1",
]

PROGRAMS = {
    "Identity": [
        # matrix + maturation proteases
        "Amelx", "Ambn", "Enam", "Amtn", "Mmp20", "Klk4", "Acp4", "Lamc2",
        # dental-epithelial lineage transcription factors
        "Dlx3", "Sp6", "Shh", "Gli1", "Wnt10a", "Trp63", "Nfic", "Sp7",
    ],
    "Sterol": STEROL_HORTON_DIRECT + [
        "Srebf2",                                            # the TF itself
        "Fdft1", "Lss", "Cyp51", "Dhcr7", "Dhcr24",          # further pathway enzymes
        "Nsdhl", "Mvk", "Pmvk", "Acat2", "Sc5d",
    ],
    "Wound": [
        "Mmp10", "Mmp12", "Mmp13", "Ereg", "Areg", "Ptgs2", "Il1a", "Cxcl2",
        "Csf3", "Fosl1", "Mt1", "Mt2", "Slc7a11", "Sprr1a",
    ],
    "IFN": [
        "Irf7", "Irf9", "Stat1", "Isg15", "Mx1", "Oasl1", "Oas1a", "Oas2", "Oas3",
        "Ifit1", "Ifit2", "Ifit3", "Ifi44", "Ifi27l2a", "Rsad2", "Usp18", "Bst2",
        "Ccl5",
    ],
}

# convenience subsets used for reporting in 02_deg_analysis.py
IDENTITY_MATRIX = PROGRAMS["Identity"][:8]
IDENTITY_TFS = PROGRAMS["Identity"][8:]

# display labels
LABELS = {
    "Identity": "Ameloblast identity",
    "Sterol": "Sterol biosynthesis (SREBP2)",
    "Wound": "Wound / stress",
    "IFN": "Interferon-stimulated",
}

if __name__ == "__main__":
    for k, v in PROGRAMS.items():
        print(f"{k:10s} n={len(v):2d}  {', '.join(v)}")
