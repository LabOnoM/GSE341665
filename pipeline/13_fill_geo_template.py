#!/usr/bin/env python3
"""
13_fill_geo_template.py -- populate GEO's seq_template.xlsx for this study.

Fills the Metadata worksheet (STUDY / SAMPLES / PROTOCOLS / PAIRED-END) and the
MD5 Checksums worksheet, following https://www.ncbi.nlm.nih.gov/geo/info/seq.html.
Cell positions were read from the blank template; the workbook's data-validation
lists (instrument model, library strategy) are preserved by editing in place.

Fields the authors must supply are written as "[[AUTHOR: ...]]" so they are
impossible to miss, rather than being invented here.

Input : seq_template.xlsx (blank GEO template; GEO_submission/ or Version02/ root)
Output: 00.Manuscript/Version02/GEO_submission/GEO_metadata_Smad4_dcKO.xlsx
"""
import os, csv, shutil
import openpyxl

ROOT = os.path.dirname(os.path.dirname(os.path.dirname(os.path.abspath(__file__))))
V2 = os.path.join(ROOT, "00.Manuscript", "Version02")
OUT = os.path.join(V2, "GEO_submission")
os.makedirs(OUT, exist_ok=True)
# the blank template may sit at the Version02 root or inside GEO_submission/
SRC = next((c for c in (os.path.join(OUT, "seq_template.xlsx"),
                        os.path.join(V2, "seq_template.xlsx"))
            if os.path.exists(c)), os.path.join(OUT, "seq_template.xlsx"))
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

PROTOCOLS = {
 57: ("growth protocol",
      "[[AUTHOR: mouse strain background, housing, age/stage at collection. Mandibular first "
      "molar tooth germs were dissected and, where applicable, grown in Trowell-type organ "
      "culture; please state medium, supplements and culture duration.]]"),
 58: ("treatment protocol",
      "No pharmacological treatment was applied to the sequenced libraries. (Tamoxifen "
      "induction and pathway-inhibitor treatments were used only for separate, non-sequenced "
      "experiments in the associated study.) [[AUTHOR: confirm.]]"),
 59: ("*extract protocol",
      "Dental epithelium was separated from the mesenchyme of mandibular molar tooth germs and "
      "total RNA was extracted. [[AUTHOR: state the RNA extraction kit/method, DNase treatment, "
      "and the RNA integrity (RIN) values or Bioanalyzer QC used to accept libraries.]]"),
 60: ("*library construction protocol",
      "Strand-specific (reverse-stranded, dUTP-type) poly(A) mRNA libraries were prepared and "
      "sequenced as 150 bp paired-end reads on an Illumina NovaSeq 6000 (LC-Bio Technology, "
      "Hangzhou, China). Strandedness was confirmed computationally: Salmon automatic library "
      "detection returned ISR for all six libraries. [[AUTHOR: state the library preparation kit "
      "name and catalogue number, and confirm poly(A) selection versus rRNA depletion.]]"),
}

DATA_PROC = [
 ("*data processing step",
  "Sequence files supplied by the sequencing provider were quality- and adapter-trimmed with "
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


def main():
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
    for i, (lib, rep, grp, geno) in enumerate(LIBS):
        r = 39 + i
        title = (f"Dental epithelium, {'control' if grp=='CON' else 'Smad4 conditional knockout'}"
                 f", biological replicate {rep}")
        vals = {
            "A": lib, "B": title, "C": "RNA-Seq", "D": "Mus musculus",
            "E": "mandibular molar tooth germ, dental epithelium",
            "F": "", "G": "dental epithelium (enamel organ)",
            "H": geno, "I": "none", "J": "",
            "K": "total RNA", "L": "paired-end", "M": "Illumina NovaSeq 6000",
            "N": f"Column '{lib}' in the processed data matrices",
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
