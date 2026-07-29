#!/usr/bin/env python3
"""
13_fill_geo_template.py -- populate GEO's seq_template.xlsx for this study.

Fills the Metadata worksheet (STUDY / SAMPLES / PROTOCOLS / PAIRED-END) and the
MD5 Checksums worksheet, following https://www.ncbi.nlm.nih.gov/geo/info/seq.html.
Cell positions were read from the blank template; the workbook's data-validation
lists (instrument model, library strategy) are preserved by editing in place.

Fields the authors must supply are written as "[[AUTHOR: ...]]" so they are
impossible to miss, rather than being invented here.

Input : seq_template.xlsx (blank GEO template; 2_GEO/ or Submission_02_JMM_2026-07/ root)
Output: 00.Manuscript/Submission_02_JMM_2026-07/2_GEO/GEO_metadata_Smad4_dcKO.xlsx
"""
import os, csv, shutil
import openpyxl

ROOT = os.path.dirname(os.path.dirname(os.path.dirname(os.path.abspath(__file__))))
V2 = os.path.join(ROOT, "00.Manuscript", "Submission_02_JMM_2026-07")
GEO_DIR = os.path.join(V2, "2_GEO")
OUT = GEO_DIR
os.makedirs(OUT, exist_ok=True)
# The blank GEO template is not kept in the repository (it is a downloadable GEO artifact,
# not our data). Fetch a fresh copy from https://www.ncbi.nlm.nih.gov/geo/info/seq.html and
# drop it into 2_GEO/ (or the submission root) before running this. The built output
# GEO_metadata_Smad4_dcKO.xlsx already carries the filled result.
SRC = next((c for c in (os.path.join(OUT, "seq_template.xlsx"),
                        os.path.join(V2, "seq_template.xlsx"))
            if os.path.exists(c)), None)
if SRC is None:
    raise SystemExit(
        "seq_template.xlsx not found in 2_GEO/ or the submission root.\n"
        "It is deliberately not tracked in git. Download the blank GEO seq template from\n"
        "https://www.ncbi.nlm.nih.gov/geo/info/seq.html and place it in 2_GEO/ before running.")
DST = os.path.join(OUT, "GEO_metadata_Smad4_dcKO.xlsx")

TITLE = ("Epithelial Smad4 loss reprograms ameloblast identity and metabolism in the "
         "mouse dental epithelium")

SUMMARY = (
 "Amelogenesis imperfecta is a group of inherited enamel defects, yet the intracellular "
 "signalling that maintains the specialised secretory identity of the enamel-forming "
 "ameloblast is poorly defined. Smad4 is the common intracellular mediator of bone "
 "morphogenetic protein (BMP) and transforming growth factor beta (TGF-beta) signalling. "
 "To determine what transcriptional state replaces the ameloblast secretory programme when "
 "this signalling hub is removed, we deleted Smad4 selectively in the dental epithelium and "
 "profiled the transcriptome of dissected dental epithelium from conditional-knockout and "
 "littermate control tooth germs by bulk RNA sequencing (three biological replicates per "
 "genotype). Smad4 loss produced a coordinated collapse of the ameloblast identity programme "
 "(enamel matrix genes, maturation proteases and dental-epithelial lineage transcription "
 "factors), which was significantly enriched for orthologues of human amelogenesis-imperfecta "
 "genes, together with activation of a wound-associated stress programme and, unexpectedly, "
 "coordinated activation of the SREBP2 cholesterol-biosynthesis programme. Cell-type "
 "deconvolution against a public mouse tooth single-cell reference indicated that the "
 "ameloblast fraction is preserved in mutant tissue, identifying the change as cell-intrinsic "
 "reprogramming rather than a shift in tissue composition.")

DESIGN = (
 "Bulk RNA sequencing of dissected dental epithelium from mouse mandibular molar tooth germs. "
 "Two groups were compared: Smad4 conditional knockout (dcKO) and Cre-negative littermate "
 "control (CON), with three biological replicates per group (six libraries total). Libraries "
 "were sequenced as 150 bp paired-end reads on an Illumina NovaSeq 6000. The comparison of "
 "interest is dcKO versus CON.")

CONTRIBUTORS = ["Xinyu,Zheng", "Dian,Su", "Fei,Hu", "Xianghui,Liu",
                "Jinshuo,Wu", "Zhiwei,He", "Xiaohua,Xie"]

LIBS = [
 # (library name, replicate, group, genotype, title)
 ("CON1", 1, "CON", "Smad4^fl/fl (Cre-negative littermate control)"),
 ("CON2", 2, "CON", "Smad4^fl/fl (Cre-negative littermate control)"),
 ("CON3", 3, "CON", "Smad4^fl/fl (Cre-negative littermate control)"),
 ("dcKO1", 1, "dcKO", "K14-creER; Smad4^fl/fl (tamoxifen-inducible dental-epithelium conditional knockout)"),
 ("dcKO2", 2, "dcKO", "K14-creER; Smad4^fl/fl (tamoxifen-inducible dental-epithelium conditional knockout)"),
 ("dcKO3", 3, "dcKO", "K14-creER; Smad4^fl/fl (tamoxifen-inducible dental-epithelium conditional knockout)"),
]

PROC1, PROC2 = "GSE_processed_gene_raw_counts.txt", "GSE_processed_gene_TPM.txt"

# Sources for the wet-lab protocols below, so they can be re-checked rather than trusted:
#   Mouse lines, husbandry, dissection, Trowell culture   -> BBRC Methods 2.1, 2.2
#   Tamoxifen concentration and 24 h exposure             -> BBRC Results (90 umol/L, 24 h)
#                                                            and the 50-300 umol/L titration
#   RNA extraction and library QC thresholds              -> BBRC Methods 2.7, 2.9
#   Developmental stage, pooling and sex composition      -> DERIVED FROM THE DEPOSITED DATA
#                                                            (see the assertions in verify())
# Anything not recoverable from those sources is marked [[AUTHOR: ...]] and must be settled
# with the wet-lab collaborators before the record is released. Do not silently fill them in.
PROTOCOLS = {
 57: ("growth protocol",
      "Mouse lines. Smad4-floxed mice (Smad4^fl/fl; Jackson Laboratory stock 017462) were "
      "crossed with the tamoxifen-inducible epithelial driver K14-creER (Jackson Laboratory "
      "stock 005107) to generate K14-creER;Smad4^fl/fl dental-epithelium conditional knockouts "
      "(dcKO). Cre-negative Smad4^fl/fl littermates from the same litters served as controls "
      "(CON). Colonies were maintained under specific-pathogen-free (SPF) conditions and "
      "genotyped by PCR from tail DNA using the protocols published by The Jackson Laboratory "
      "for each allele. [[AUTHOR: confirm the housing detail. Standard SPF caging with a 12 h "
      "light/dark cycle and food and water ad libitum is ASSUMED here and is not documented in "
      "the source records.]] "
      "Tissue collection and organ culture. Pups were collected at postnatal day 0 (P0). Heads "
      "were dissected in ice-cold PBS and held in complete DMEM in 12-well plates while "
      "genotypes were confirmed. Mandibular first-molar germs were then micro-dissected and "
      "grown at the medium-gas interface in Trowell-type organ culture in DMEM supplemented "
      "with 10% fetal bovine serum (Lonsera, S711-001S) and antibiotics, at 37 C in a "
      "humidified 5% CO2 atmosphere. The enamel-organ epithelium was subsequently separated "
      "from the mesenchyme and stored at -80 C until RNA extraction. [[AUTHOR: state the total "
      "time in culture between explant and epithelial harvest. The written record covers only "
      "the 24 h induction window, so the interval between the end of induction and harvest is "
      "unknown.]] "
      "Developmental stage, determined from the deposited data. The sequenced epithelium is at "
      "the enamel-secretory stage and had not entered maturation: in the processed TPM matrix "
      "Amelx averages 8.5e4 TPM and Ambn 1.8e4 TPM, whereas the maturation-stage markers Klk4 "
      "and Slc24a4 average only 17 and 0.4 TPM. "
      "Pooling and sex, determined from the deposited data. Each library is a pool of tooth "
      "germs from more than one pup and contains both sexes: all six libraries co-express Xist "
      "(91-136 TPM) together with the Y-linked genes Ddx3y, Eif2s3y, Uty and Kdm5d (Eif2s3y "
      "52-81 TPM), a combination no single animal can produce. The male fraction is "
      "systematically higher in the dcKO pools than in the CON pools - Xist is lower and all "
      "four Y-linked genes are higher in every dcKO library (Xist Benjamini-Hochberg adjusted "
      "P = 3.5e-54) - so sex composition is an unbalanced covariate in this design and users of "
      "these data should account for it. [[AUTHOR: state how many pups were pooled per "
      "library.]]"),
 58: ("treatment protocol",
      "Smad4 deletion in the K14-creER line is tamoxifen-dependent, so the inducer was applied "
      "to the sequenced material. 4-Hydroxytamoxifen (Sigma-Aldrich, H7904) was added to the "
      "culture medium at 90 umol/L for 24 h. That concentration was selected from a ten-point "
      "titration (50, 60, 70, 80, 90, 100, 120, 150, 200 and 300 umol/L) scored on tooth-germ "
      "development: 50 umol/L produced only partial enamel hypoplasia, while 120 umol/L and "
      "above arrested germ development. CON (Cre-negative Smad4^fl/fl) germs were treated "
      "identically, so 4-hydroxytamoxifen exposure is common to both groups and is not a "
      "variable in the dcKO-versus-CON comparison. No other pharmacological agent was applied "
      "to the sequenced libraries; the pathway inhibitors (dorsomorphin, SB525334, U0126, "
      "SB203580) were used only in separate, non-sequenced experiments."),
 59: ("*extract protocol",
      "The enamel-organ epithelial layer was separated from the mesenchyme of the mandibular "
      "first-molar germs and stored at -80 C. Total RNA was isolated with TRIzol (ThermoFisher "
      "Scientific, 15596018), quantified on a NanoDrop spectrophotometer and assessed for "
      "integrity on an Agilent Bioanalyzer 2100. Samples were released for library "
      "construction only if RNA concentration exceeded 50 ng/uL, RIN exceeded 7.0 and total "
      "RNA exceeded 1 ug. [[AUTHOR: confirm whether a DNase I step was included, and supply "
      "the per-sample RIN values if the provider's QC report is still available.]]"),
 60: ("*library construction protocol",
      "Strand-specific (reverse-stranded, dUTP-type) poly(A) mRNA libraries were prepared and "
      "sequenced as 150 bp paired-end reads on an Illumina NovaSeq 6000 (LC-Bio Technology, "
      "Hangzhou, China). Strandedness was confirmed computationally: Salmon automatic library "
      "detection returned ISR for all six libraries. [[AUTHOR: state the library preparation kit "
      "name and catalogue number, and confirm poly(A) selection versus rRNA depletion.]]"),
}

DATA_PROC = [
 # State the processing level of the deposited files explicitly. They are named
 # *_Clean_Data[12].fastq.gz and are the provider's adapter- and quality-filtered reads,
 # not unfiltered instrument output; reads prior to that filtering were never returned.
 # A reader comparing read counts with studies that deposit unfiltered reads needs this.
 ("*data processing step",
  "Processing level of the deposited sequence files: the raw files submitted here "
  "(*_Clean_Data[12].fastq.gz) are the adapter- and quality-filtered reads returned by the "
  "sequencing provider (LC-Bio, Hangzhou, China). Unfiltered instrument output was not "
  "returned to the submitters and is therefore not available for deposit; these files are "
  "the rawest data available for this study."),
 ("data processing step",
  "Provider-supplied sequence files were further quality- and adapter-trimmed with "
  "Trim Galore v0.6.10 (cutadapt v4.9), Phred quality cutoff 20, automatic adapter detection."),
 ("data processing step",
  "Trimmed reads were aligned to the mouse genome with STAR v2.7.10b in two-pass mode "
  "(--twopassMode Basic --chimSegmentMin 10 --outFilterMultimapNmax 20 "
  "--outFilterMismatchNoverLmax 0.04 --alignIntronMin 70 --alignIntronMax 1000000 "
  "--outFilterType BySJout), against an index built with --sjdbOverhang 149. Uniquely mapped "
  "read pairs ranged from 93.66% to 93.84% across the six libraries."),
 ("data processing step",
  "Transcript-level abundances were estimated with Salmon v1.10.0 (--libType A, auto-detected "
  "as ISR) against an index built on a StringTie-merged transcriptome, and summarised to gene "
  "level with tximport. Gene-level counts and TPM were retained."),
 ("data processing step",
  "Genes with at least 10 reads in at least 3 samples were retained, giving 15,803 genes. One "
  "factor of unwanted variation (W_1) was estimated with RUVSeq and included as a covariate."),
 ("data processing step",
  "Differential expression was assessed with DESeq2 using a likelihood-ratio test of the full "
  "design ~W_1 + genotype against the reduced design ~W_1. Genes with Benjamini-Hochberg "
  "adjusted P < 0.05 and absolute log2 fold change > 1 were considered differentially expressed."),
]

GENOME = ("GRCm39 (Ensembl release 116; Mus_musculus.GRCm39.dna.primary_assembly.fa with "
          "Mus_musculus.GRCm39.116.gtf)")

PROC_FORMAT = [
 (68, "*processed data files format and content",
  f"{PROC1}: tab-delimited text file of gene-level raw counts for all 15,803 analysed genes "
  "(rows) across all six libraries (columns). Columns 1-3 are ensembl_gene_id, "
  "external_gene_name and gene_biotype; the remaining column headers match the library names "
  "(CON1, CON2, CON3, dcKO1, dcKO2, dcKO3)."),
 (69, "processed data files format and content",
  f"{PROC2}: tab-delimited text file of gene-level transcripts per million (TPM), same rows, "
  "columns and header convention as the raw-counts file."),
]


def verify():
    """Re-derive every data-derived number quoted in the growth protocol, and fail loudly.

    The growth protocol asserts a developmental stage and a sex composition. Both are
    inferences from the expression matrix rather than facts copied from a lab record, so
    they are recomputed on every build: if the upstream tables ever change, this raises
    instead of letting a stale claim ship to a public archive.
    """
    import pandas as pd

    out = os.path.join(os.path.dirname(os.path.abspath(__file__)), "output")
    tpm = pd.read_csv(os.path.join(out, "TPM_matrix.csv"), header=1,
                      float_precision="round_trip")
    deg = pd.read_csv(os.path.join(out, "DEG_table.csv"), header=1,
                      float_precision="round_trip")
    cols = [c for c in tpm.columns if c.startswith("TPM_")]

    def mean_tpm(gene):
        row = tpm[tpm.external_gene_name.astype(str) == gene]
        if not len(row):
            raise SystemExit(f"verify: {gene} absent from TPM_matrix.csv")
        return row.iloc[0][cols].astype(float).mean()

    def span(gene):
        row = tpm[tpm.external_gene_name.astype(str) == gene]
        if not len(row):
            raise SystemExit(f"verify: {gene} absent from TPM_matrix.csv")
        v = row.iloc[0][cols].astype(float)
        return v.min(), v.max()

    # -- stage: secretory, not maturation --------------------------------------------
    for gene, lo, hi in (("Amelx", 8.0e4, 9.0e4), ("Ambn", 1.6e4, 2.0e4),
                         ("Klk4", 15.0, 19.0), ("Slc24a4", 0.2, 0.6)):
        m = mean_tpm(gene)
        if not lo <= m <= hi:
            raise SystemExit(f"verify: mean {gene} TPM {m:.4g} outside the stated {lo:g}-{hi:g}; "
                             "the developmental-stage claim in the growth protocol is stale")

    # -- pooling and sex: Xist and Y genes co-expressed in every library ---------------
    xlo, xhi = span("Xist")
    elo, ehi = span("Eif2s3y")
    if not (85 <= xlo and xhi <= 140):
        raise SystemExit(f"verify: Xist range {xlo:.1f}-{xhi:.1f} TPM contradicts the stated 91-136")
    if not (50 <= elo and ehi <= 85):
        raise SystemExit(f"verify: Eif2s3y range {elo:.1f}-{ehi:.1f} TPM contradicts the stated 52-81")
    for g in ("Ddx3y", "Eif2s3y", "Uty", "Kdm5d"):
        if span(g)[0] <= 0:
            raise SystemExit(f"verify: {g} is not detected in every library; the "
                             "mixed-sex pooling claim does not hold")

    # -- the direction and strength of the sex imbalance -------------------------------
    xr = deg[deg.external_gene_name.astype(str) == "Xist"]
    if not len(xr):
        raise SystemExit("verify: Xist absent from DEG_table.csv")
    lfc, padj = float(xr.iloc[0]["log2FoldChange"]), float(xr.iloc[0]["padj"])
    if lfc >= 0:
        raise SystemExit(f"verify: Xist log2FC {lfc:+.3f} is not negative; the growth protocol "
                         "states the dcKO pools are the more male of the two groups")
    if padj > 1e-50:
        raise SystemExit(f"verify: Xist padj {padj:.3g} weaker than the stated 3.5e-54")
    ncon = tpm[[c for c in cols if "CON" in c]]
    for g in ("Ddx3y", "Eif2s3y", "Uty", "Kdm5d"):
        row = tpm[tpm.external_gene_name.astype(str) == g].iloc[0]
        con = row[[c for c in cols if "CON" in c]].astype(float).max()
        ko = row[[c for c in cols if "dcKO" in c]].astype(float).min()
        if ko <= con:
            raise SystemExit(f"verify: {g} does not separate the groups completely "
                             f"(max CON {con:.1f} >= min dcKO {ko:.1f})")
    print("  verify            : stage, pooling and sex-imbalance claims re-derived from the data")


def main():
    verify()
    shutil.copy2(SRC, DST)
    wb = openpyxl.load_workbook(DST)
    ws = wb["Metadata"]

    # ---------------- STUDY ----------------
    ws["B12"] = TITLE
    ws["B13"] = SUMMARY
    ws["B14"] = DESIGN
    for i, c in enumerate(CONTRIBUTORS):
        ws.cell(row=15 + i, column=2, value=c)

    # ---------------- SAMPLES (rows 39-44) ----------------
    # Two template characteristics columns are repurposed. 'cell line' and 'batch' are both
    # inapplicable here - this is primary tissue, and all six libraries came off one flowcell
    # lane - whereas developmental stage and sex are the two facts a reader of these data
    # actually needs, and sex is not balanced across the groups.
    ws["F38"] = "**developmental stage"
    ws["J38"] = "**sex"
    for i, (lib, rep, grp, geno) in enumerate(LIBS):
        r = 39 + i
        title = (f"Dental epithelium, {'control' if grp=='CON' else 'Smad4 conditional knockout'}"
                 f", biological replicate {rep}")
        vals = {
            "A": lib, "B": title, "C": "RNA-Seq", "D": "Mus musculus",
            "E": "mandibular molar tooth germ, dental epithelium",
            "F": "postnatal day 0 (P0) at dissection; enamel-secretory stage at harvest",
            "G": "dental epithelium (enamel organ)",
            "H": geno,
            "I": "4-hydroxytamoxifen 90 umol/L, 24 h (applied to both genotypes)",
            "J": "pooled male and female",
            "K": "total RNA", "L": "paired-end", "M": "Illumina NovaSeq 6000",
            "N": (f"Pool of mandibular first-molar tooth germs from several P0 pups of both "
                  f"sexes. Column '{lib}' in the processed data matrices"),
            "O": PROC1, "P": PROC2,
            "Q": f"{lib}_Clean_Data1.fastq.gz", "R": f"{lib}_Clean_Data2.fastq.gz",
        }
        for col, v in vals.items():
            ws[f"{col}{r}"] = v

    # ---------------- PROTOCOLS ----------------
    for row, (_, text) in PROTOCOLS.items():
        ws.cell(row=row, column=2, value=text)
    for i, (_, text) in enumerate(DATA_PROC):
        ws.cell(row=62 + i, column=2, value=text)
    ws["B67"] = GENOME
    for row, _, text in PROC_FORMAT:
        ws.cell(row=row, column=2, value=text)

    # ---------------- PAIRED-END (data starts row 77) ----------------
    for i, (lib, *_rest) in enumerate(LIBS):
        ws.cell(row=77 + i, column=1, value=f"{lib}_Clean_Data1.fastq.gz")
        ws.cell(row=77 + i, column=2, value=f"{lib}_Clean_Data2.fastq.gz")

    # ---------------- MD5 checksums ----------------
    md = wb["MD5 Checksums"]
    rows = list(csv.DictReader(open(os.path.join(OUT, "GEO_checksums.tsv")), delimiter="\t"))
    raw = [r for r in rows if r["type"] == "raw"]
    pro = [r for r in rows if r["type"] == "processed"]
    for i, r in enumerate(raw):                      # cols A,B from row 8
        md.cell(row=8 + i, column=1, value=r["file name"])
        md.cell(row=8 + i, column=2, value=r["file checksum"])
    for i, r in enumerate(pro):                      # cols F,G from row 8
        md.cell(row=8 + i, column=6, value=r["file name"])
        md.cell(row=8 + i, column=7, value=r["file checksum"])

    wb.save(DST)
    print("wrote", DST)
    print(f"  STUDY filled | {len(LIBS)} samples | protocols + {len(DATA_PROC)} processing steps")
    print(f"  MD5: {len(raw)} raw + {len(pro)} processed")


if __name__ == "__main__":
    main()
