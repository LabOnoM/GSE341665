# Current supplementary methods figures and tables

Additional file 1

# Supplementary Information

Epithelial Smad4 loss reprograms ameloblast identity and metabolism

Ziyi Wang1, Xinyu Zheng2,3, Mitsuaki Ono4, Dian Su3, Fei Hu3, Xianghui Liu3,5, Zhiwei He6**, Xiaohua Xie3*

1 Department of Advanced International and Information Dentistry, Faculty of Medicine, Dentistry and Pharmaceutical Sciences, Okayama University, Okayama, Japan

2 Department of Tissue and Developmental Biology, Graduate School of Dentistry, The University of Osaka, Suita, Osaka, Japan

3 Department of Stomatology, Second Affiliated Hospital of Harbin Medical University, Harbin, China

4 Department of Oral Rehabilitation and Regenerative Medicine, Okayama University Graduate School of Medicine, Dentistry and Pharmaceutical Sciences, Okayama, Japan

5 Qingdao Stomatology Hospital, Qingdao, China

6 Department of Stomatology, Shenzhen University General Hospital, Shenzhen University, Shenzhen, Guangdong, China

*Correspondence: Xiaohua Xie (xiexiaohua_79@126.com); **Correspondence: Zhiwei He (3563570@qq.com).

## Contents

Supplementary Methods — Part A: primary bulk RNA-sequencing pipeline; Part B: re-analysis performed for this manuscript

Supplementary Figure S1 — mouse lines, crosses, and the allocation of each experiment to a figure

Supplementary Figure S2 — genotyping records and dental epithelial region

Supplementary Figure S3 — read-level genotyping of the conditional Smad4 allele; the Cre excision product is not detected

Supplementary Figure S4 — tamoxifen dose-dependence of the inducible phenotype

Supplementary Figure S5 — complete analysis workflow (A: primary pipeline; B: re-analysis)

Supplementary Figure S6 — single-cell reference used for deconvolution (UMAP and marker dot plot)

Supplementary Figure S7 — the inducible K14-creER driver reproduces the enamel phenotype

Supplementary Table S1 — qPCR and genotyping primers

Supplementary Table S2 — culture and inhibitor regimens

Supplementary Table S3 — antibody records and identity checks

Supplementary Table S4 — verified provenance of the 31 PanelApp genes used in the enrichment

Supplementary Data S1–S8 — machine-readable result tables

## Supplementary Methods

This document records the exact workflow, parameters, and software versions for all computational analyses, in two parts: Part A, the primary bulk RNA-sequencing pipeline; and Part B, the re-analysis performed for this manuscript (gene-set (program) scoring and cell-type deconvolution).

### Part A — Primary bulk RNA-sequencing pipeline

#### A1. Experimental design

| Field | Value |
| --- | --- |
| Samples | 6: CON1, CON2, CON3 (control) and dcKO1, dcKO2, dcKO3 (Smad4 conditional knockout) |
| Genotype | dcKO: K14-creER; Smad4fl/fl (tamoxifen-inducible epithelial deletion at the secretory stage). CON: Cre-negative Smad4fl/fl littermates. The constitutive K14-cre line was used for the in vivo enamel phenotype (Figure 1), not for the sequenced material. |
| Groups | 2 groups × 3 biological replicates |
| Sex | All six study libraries were derived from male experimental animals. Each biological replicate pooled tooth germs from more than one male pup; exact pup counts and animal-to-library assignments remain incompletely documented. |
| Tissue | Dental epithelium (mandibular molar tooth germ) |
| Layout | Paired-end, 150 bp |
| Library strandedness | Reverse-stranded; automatically detected by Salmon as ISR, and used as -s 2 in featureCounts |
| Contrast | dcKO vs CON |

Sex was confirmed by the experimentalist. Expression of Xist, Ddx3y, Eif2s3y, Uty and Kdm5d is retained in Additional file 2 (Data S3) as supporting sample information, rather than an independent determination of sex.

#### A2. Reference genome and annotation

| Field | Value |
| --- | --- |
| Assembly | Mus musculus GRCm39,; Mus_musculus.GRCm39.dna.primary_assembly.fa |
| Annotation | Mus_musculus.GRCm39.116.gtf |
| Source | Ensembl release 116 (ftp://ftp.ensembl.org/pub/release-116/) |
| Annotation database | org.Mm.eg.db; biomaRt dataset mmusculus_gene_ensembl; KEGG organism mmu; UniProt taxon 10090 |

The GTF was filtered before index building (genome_filtered.gtf); genes lacking a gene_name attribute retained their gene_id.

#### A3. Read quality control and trimming

FastQC on the provider-supplied reads, on trimmed reads, and on aligned BAM files.

Processing level of the deposited sequence files. The sequencing provider (LC-Bio) returned adapter- and quality-filtered files, named *_Clean_Data[12].fastq.gz. Reads prior to that filtering were not returned and cannot be deposited. These files are what is deposited in GEO and what Trim Galore took as input, so the deposit represents the rawest data available for this study rather than unfiltered instrument output. The distinction matters when comparing read counts with studies that deposit unfiltered reads.

Trim Galore v0.6.10 (wrapping cutadapt v4.9) with default settings: quality Phred cutoff 20, ASCII+33 encoding, automatic adapter detection.

#### A4. Alignment (STAR)

STAR version 2.7.10b.

Index construction:

STAR --runThreadN 32 --runMode genomeGenerate \
     --genomeDir ./GenomeIndex_ENSEMBL_STAR_ReadLength149/ \
     --genomeFastaFiles genome.fa --sjdbGTFfile genome_filtered.gtf \
     --sjdbOverhang 149 --limitGenomeGenerateRAM 101781612134

Alignment (per sample; two-pass):

STAR --runMode alignReads --chimSegmentMin 10 \
     --genomeDir ./GenomeIndex_ENSEMBL_STAR_ReadLength149 \
     --readFilesIn <sample>_Clean_Data1_val_1.fq.gz <sample>_Clean_Data2_val_2.fq.gz \
     --outFileNamePrefix ./<sample>.STAR. --alignEndsType Local \
     --outSAMtype BAM SortedByCoordinate --outFilterMultimapNmax 20 \
     --outFilterMismatchNmax 999 --outFilterMismatchNoverLmax 0.04 \
     --alignIntronMin 70 --alignIntronMax 1000000 --alignMatesGapMax 1000000 \
     --alignSJoverhangMin 8 --alignSJDBoverhangMin 1 \
     --outSAMstrandField intronMotif --outFilterType BySJout \
     --twopassMode Basic --outSAMmapqUnique 60 \
     --limitGenomeGenerateRAM 101781612134 --runThreadN 32 --readFilesCommand zcat

Mapping statistics (from *.STAR.Log.final.out):

| Sample | Input read pairs | Uniquely mapped | Uniquely mapped % | Multi-mapped % | Chimeric % |
| --- | --- | --- | --- | --- | --- |
| CON1 | 21,303,455 | 19,975,969 | 93.77 | 3.85 | 0.30 |
| CON2 | 21,235,660 | 19,900,363 | 93.71 | 3.80 | 0.29 |
| CON3 | 20,986,103 | 19,694,042 | 93.84 | 3.81 | 0.28 |
| dcKO1 | 22,453,825 | 21,044,062 | 93.72 | 3.79 | 0.30 |
| dcKO2 | 21,439,266 | 20,080,247 | 93.66 | 3.74 | 0.31 |
| dcKO3 | 22,455,897 | 21,031,771 | 93.66 | 3.77 | 0.28 |

Mapping was uniform across samples (uniquely mapped 93.66–93.84%), with no genotype-associated difference in alignment quality.

#### A5. Quantification

Two quantifications were produced:

featureCounts (gene-level counts from the STAR BAMs): featureCounts -T 32 -s 2 -p -Q 20 -O -a ./genome_filtered.gtf -o <sample>.STAR.counts (-s 2 reverse-stranded; -p paired-end; -Q 20 minimum mapping quality; -O assign reads overlapping multiple features.)

Salmon v1.10.0 selective alignment (--libType A, auto-detected as ISR in every sample) against an index built from a StringTie-merged transcriptome (a de-novo assembly merged with the reference annotation); library type auto-detected as ISR.

The differential-expression analysis used the Salmon quantifications imported with tximport (pipeline parameter CountsFrom = "Tximport"). Transcripts per million (TPM) reported in this study derive from that import.

#### A6. Differential expression (DESeq2 with RUV correction)

This is the step most often mis-stated; the parameters below are transcribed from the saved DESeqDataSet and the pipeline parameter object.

| Field | Value |
| --- | --- |
| Package | DESeq2 |
| Counts | Salmon → tximport |
| Unwanted-variation correction | RUVSeq, one factor (k = 1), stored as covariate W_1. The specific RUVSeq function (RUVg/RUVs/RUVr) and the control-gene set are not recorded in the pipeline outputs; k = 1 is recoverable from the single W_1 column. |
| Design formula | ~ W_1 + RUV_x (where RUV_x is genotype: CON, dcKO) |
| Test | Likelihood-ratio test (LRT), full ~ W_1 + RUV_x vs reduced ~ W_1 (Par_Test = "LRT") |
| Contrast | RUV_x, dcKO vs CON |
| Expression filter | Genes with at least 10 reads in at least 3 samples (Par_ReN = 10, Par_SaN = 3) |
| Genes tested | 15,803 |
| Significance | padj < 0.05 (Benjamini–Hochberg) and \|log2 fold change\| > 1 (Par_FDR = 0.05, Par_LFC = 1, Par_PorFRR = "padj") |
| Log2 fold-change shrinkage | Not applied (fold changes are maximum-likelihood estimates) |
| Independent filtering | Enabled; 613 genes (3.9%) have padj = NA and are treated as non-significant |
| Result | 563 differentially expressed genes: 209 up-regulated, 354 down-regulated |

Because the LRT compares the full and reduced models, the reported P values test the overall genotype effect after adjusting for W_1; the accompanying log2 fold changes are the maximum-likelihood estimates for RUV_x dcKO vs CON.

#### A7. Functional enrichment

Over-representation analysis (ORA) and GSEA against GO (biological process, molecular function, cellular component) and KEGG (organism mmu), using the org.Mm.eg.db annotation package in R.

Significance was assessed on the raw P value, not the false-discovery rate. The report states that enriched GO and KEGG terms, in both GSEA and ORA, “were defined as P-value < 0.05”, and the pipeline parameter Par_UseFDR4FEA = FALSE confirms that FDR filtering was disabled. Enrichment results are therefore not multiplicity-corrected and are used in this manuscript only as qualitative corroboration of programs that were independently identified from the ranked differential-expression results; no enrichment statistic is reported as a finding.

The enrichment software is not named anywhere in the pipeline outputs. [[Authors: confirm the enrichment package and version before submission.]]

Sample clustering/heatmaps used Pearson distance with ward.D2 linkage (Par_dist = "pearson", Par_hclust = "ward.D2").

#### A8. Additional analyses performed by the primary pipeline

These were run but are not used in the present manuscript; they are listed for completeness and transparency.

| Analysis | Tool (version) | Outcome — stated accurately |
| --- | --- | --- |
| Differential exon usage | DEXSeq v1.48.0 | 700 significant exon bins across 331 genes at FDR < 0.1 — but tested on a restricted subset of 987 genes, not transcriptome-wide. Not used. |
| Isoform switching | IsoformSwitchAnalyzeR v2.2.0 with CPC2, Pfam, SignalP 6.0, IUPred3 | 1,288 isoform switches in 1,098 genes; no alternative-splicing type and no functional consequence reached enrichment significance. Not used. |
| Intron excision | LeafCutter | Only intron clustering was run; no differential intron-excision test was performed. LeafCutter must not be cited as a completed differential analysis. |
| Circular RNA | CIRCexplorer2 v2.3.8 + circRNAprofiler | 184 circRNAs detected from 210 back-spliced junctions; 0 differentially expressed. The RBP-motif module failed (the ATtRACT motif database could not be downloaded) and produced no output. Not used. |
| Co-expression networks | WGCNA | 27 modules. The soft-thresholding power actually selected, the network type, and the TOM type are not recorded, so the network cannot be reproduced exactly. With n = 6 and a binary trait, module–trait correlations are point-biserial coefficients on six points and are not interpretable as discoveries. Not used. |

Statement on excluded analyses. None of the analyses in the table above is used in the main text, and each is excluded for a stated reason: the circular-RNA arm returned no differentially expressed circRNAs; the isoform-switching arm returned switches but no enriched splicing type or functional consequence; the exon-usage arm was run on a restricted gene subset; LeafCutter differential testing was never executed; and the WGCNA network cannot be reproduced because its key parameters were not recorded and is in any case underpowered at n = 6. An earlier draft of this project additionally explored a KMT2D/ENCODE epigenetic model; that analysis was discarded because Kmt2d is not differentially expressed in these data (log2 fold change −0.08, adjusted P = 0.36).

#### A9. Parameters that could not be recovered

For transparency, the following were not recorded anywhere in the pipeline outputs and are therefore not stated in this manuscript. They would be required for a bit-exact reproduction of the primary pipeline:

The RUVSeq function used (RUVg, RUVs or RUVr) and its control-gene set.

The WGCNA soft-thresholding power selected, network type, TOM type and correlation type.

GSEA parameters (ranking metric, permutation number, minimum and maximum gene-set size) and the ORA background/universe set.

The DEXSeq design formula, in particular whether the RUV covariate was carried into the exon-usage model.

The sequencing library-preparation kit.

A further inconsistency is noted: the gene-level DESeq2 model and the isoform-level analysis each adjusted for a single latent factor, but these were different vectors (they differ in sign for one sample), so the two analyses did not correct for the same unwanted variation. This does not affect the gene-level results reported here, which are the only ones used.

### Part B — Re-analysis performed for this manuscript

Code for the analyses, figures and numerical exports used by the current manuscript is maintained at https://github.com/LabOnoM/GSE341665. The repository README, reproduction map and input manifest specify the retained scripts, exact inputs and execution commands. The cleaned package excludes unused legacy analyses and document/deposit builders.

#### B1. Tested standalone reproduction environment

| Component | Version |
| --- | --- |
| Python | 3.12.14 |
| numpy / pandas / scipy / matplotlib | 2.3.5 / 2.2.3 / 1.18.1 / 3.11.2 |
| R | 4.3.2 |
| Seurat / SeuratObject | 5.3.1 / 5.2.0 |
| Matrix / uwot / ggplot2 / patchwork | 1.6.5 / 0.2.3 / 4.0.0 / 1.3.2 |
| limma / edgeR / statmod | 3.58.1 / 4.0.16 / 1.5.0 |
| Operating system | macOS, arm64; software locks and session information in the repository |

#### B2. Descriptive scoring and ROAST gene-program inference

The four curated gene sets (programs) are defined once in pipeline/programs.py in the GitHub repository and imported by program-scoring, deconvolution, and expression-figure scripts. Membership is tabulated below.

Two membership decisions are recorded explicitly. (i) Krt76 is among the most strongly reduced transcripts in the mutant epithelium (log2 fold change −4.64) but was excluded from the ameloblast-identity program: no published work characterizes Krt76 in ameloblasts, and the paper most often associated with hair keratins in enamel concerns KRT75, not KRT76. It is therefore reported in the main text as an observation rather than as a literature-defined program member. (ii) Mmp12 is retained in the wound/stress program because it behaves with the rest of that program, but it is classically a macrophage metalloelastase and may be contributed in part by non-epithelial cells; this caveat is stated in the main text.

Program membership. The four programs follow published definitions and are defined once in pipeline/programs.py:

| Program | n | Members | Basis |
| --- | --- | --- | --- |
| Ameloblast identity | 16 | Amelx, Ambn, Enam, Amtn, Mmp20, Klk4, Acp4, Lamc2 (matrix and maturation proteases); Dlx3, Sp6, Shh, Gli1, Wnt10a, Trp63, Nfic, Sp7 (lineage transcription factors) | Lacruz et al. 2017; Bartlett 2013; Balic and Thesleff 2015; Zhang et al. 2015; Rhodes et al. 2021; Lee et al. 2018 |
| Sterol biosynthesis | 20 | Hmgcr, Hmgcs1, Mvd, Fdps, Sqle, Msmo1, Idi1, Ldlr, Insig1 (nine experimentally defined direct SREBP targets), Srebf2, and Fdft1, Lss, Cyp51, Dhcr7, Dhcr24, Nsdhl, Mvk, Pmvk, Acat2, Sc5d | Horton et al. 2003 |
| Wound / stress | 14 | Mmp10, Mmp12, Mmp13, Ereg, Areg, Ptgs2, Il1a, Cxcl2, Csf3, Fosl1, Mt1, Mt2, Slc7a11, Sprr1a | Werner and Grose 2003; Ge et al. 2017 |
| Interferon-stimulated | 18 | Irf7, Irf9, Stat1, Isg15, Mx1, Oasl1, Oas1a, Oas2, Oas3, Ifit1, Ifit2, Ifit3, Ifi44, Ifi27l2a, Rsad2, Usp18, Bst2, Ccl5 | Interferome (Rusinova et al. 2013) |

For each sample, a program score was computed as the mean of gene-wise z-scores of log2(TPM + 1) across program members — the combined gene-set z-score of Lee et al. (2008). These six-library values remain descriptive; ROAST tests gene expression through the fitted model rather than these averaged scores.

Primary inference. ROAST (Wu et al., 2010; doi:10.1093/bioinformatics/btq401) tests differential expression within each program while preserving inter-gene correlation. Archived DESeq2_Processing.Rdata supplied counts, full-precision TPM, effective transcript lengths and W_1 for the unchanged 15,803-gene universe. All 68 members mapped uniquely (Identity 16, Sterol 20, Wound 14, IFN 18). Counts were reconstructed within this universe as TPM × mean effective length, rescaled to archived retained count totals; no second length correction was applied. TMM normalization (Robinson and Oshlack, 2010; doi:10.1186/gb-2010-11-3-r25) and voom precision weights (Law et al., 2014; doi:10.1186/gb-2014-15-2-r29) preceded limma fitting (Ritchie et al., 2015; doi:10.1093/nar/gkv007). The primary design was ~ W_1 + genotype, dcKO minus CON, with three residual degrees of freedom. limma::mroast used the mean statistic, two-sided directional PValue, 99,999 rotations, seed 20261003 and midp = FALSE; ordinary-P BH adjustment covered four programs. The score means in the table are descriptive; its ROAST probabilities arise from the gene-level model.

Sensitivities. Genotype-only ROAST used ~ genotype, four residual degrees of freedom and seed 20261004. CAMERA (Wu and Smyth, 2012; doi:10.1093/nar/gks461) addressed the competitive question against the other retained genes, with inter.gene.cor = NA, allow.neg.cor = FALSE and use.ranks = FALSE in both models. Fixed correlation 0.01 CAMERA and mixed-direction ROAST were explicitly secondary. Each method/model/endpoint used its own four-program BH family. Primary and 28 supplementary sensitivity rows, membership and sample metadata are in Additional file 2 (Data S3) and matching CSVs. All analyses use the same six biological libraries and are not independent validation.

QC and limits. Inputs, sample alignment, contrast sign, design rank, count totals, TMM factors, voom weights/residuals, per-gene contributions and independent BH calculations were audited. Optional GEO checking matched all 94,818 archived count entries. R 4.3.2, limma 3.58.1, edgeR 4.0.16 and statmod 1.5.0 were used; package/RNG locks and complete session information are supplied. Small replication, incomplete W_1 derivation/control genes and donor/litter mapping, very large weights among abundant genes and post hoc curation limit inference. The programs were frozen for these runs, not prospectively preregistered. At 99,999 rotations, raw P = 0.00001 is the reporting floor. The old exact score-permutation P = 0.10 was valid for that statistic; the categorical claim that z-scoring invalidates t tests is unsupported and is retired. The separate locus/exon permutation analysis in B7 and Figure S3 is unchanged.

| Program | n genes | CON | dcKO | ROAST BH-adjusted P |
| --- | --- | --- | --- | --- |
| Ameloblast identity | 16 | +0.90 | −0.90 | 2.67e-05 |
| Sterol biosynthesis | 20 | −0.91 | +0.91 | 2.00e-05 |
| Wound / stress | 14 | −0.91 | +0.91 | 2.00e-05 |
| Interferon-stimulated | 18 | +0.88 | −0.88 | 5.00e-05 |

#### B3. Single-cell reference and annotation

| Field | Value |
| --- | --- |
| Dataset | Gene Expression Omnibus GSE146855 (Chiba et al. 2020), postnatal day 7 mouse incisor, 10x Genomics |
| Cells after quality control | 6,248 (of 6,260 deposited) |
| QC filters | nFeature_RNA > 200 and < 7,000; mitochondrial content < 25%; genes detected in ≥ 3 cells |
| Normalisation | NormalizeData (log-normalisation, scale factor 10,000) |
| Variable features | 2,000 (FindVariableFeatures, vst) |
| Dimensionality reduction | RunPCA, 30 components computed; 20 PCs used downstream |
| Clustering | FindNeighbors(dims = 1:20), FindClusters(resolution = 0.6) → 17 clusters |
| Embedding | RunUMAP(dims = 1:20) (uwot, cosine metric) |
| Seed | Global seed 1; PCA/UMAP seed 42; clustering seed 0. |

Annotation strategy (important). Enamel-matrix transcripts (Amelx, Ambn, Enam) are so abundant in this tissue that they appear as ambient RNA in every droplet — visible in Supplementary Figure S6B as low-level detection across all cell types. Assigning clusters by absolute marker expression would therefore label every cluster an ameloblast. Each cluster was instead assigned to the cell type whose canonical marker set showed the highest cross-cluster specificity, computed as the mean z-score of marker expression across clusters. Marker means were rounded to two decimals before assigning the first maximum in marker-set order.

Marker sets used: secretory ameloblast (Amelx, Ambn, Enam, Mmp20); maturation ameloblast (Klk4, Odam, Amtn, Slc24a4, Gm17660); other dental epithelium (Shh, Sfrp5, Krt14, Krt5, Epcam, Thbd, Gjb3, Cldn10); mesenchyme (Col1a1, Dcn, Lum, Prrx1, Pdgfra, Msx1, Sp7, Dspp); immune (Ptprc, Lyz2, C1qa, Cd68, Fcer1g); endothelial (Pecam1, Cdh5, Egfl7, Emcn); erythroid (Hba-a1, Hbb-bs, Hba-a2).

Resulting reference composition (Supplementary Figure S6):

| Cell type | Cells | % |
| --- | --- | --- |
| Mesenchyme | 2,527 | 40.4 |
| Epithelial (inner enamel epithelium / stratum intermedium) | 1,361 | 21.8 |
| Ameloblast, maturation-stage | 1,346 | 21.5 |
| Ameloblast, secretory-stage | 506 | 8.1 |
| Endothelial | 290 | 4.6 |
| Immune | 131 | 2.1 |
| Erythroid | 87 | 1.4 |

#### B4. Cell-type deconvolution

| Field | Value |
| --- | --- |
| Reference signature | Pseudobulk counts per annotated cell type, converted to counts per million |
| Marker panel | Genes with reference expression > 20 CPM in a cell type and > 2-fold above the cross-cell-type mean; up to 120 per cell type; union retained and intersected with the bulk matrix (542 genes) |
| Method | Non-negative least squares (scipy.optimize.nnls; Lawson–Hanson), solved independently per bulk sample |
| Normalisation | Coefficients scaled to sum to 1 per sample |
| Comparison | Two-tailed t test, CON vs dcKO |
| Program localisation | Cell-type specificity of each program = mean across-cell-type z-score of reference expression of its member genes |

Result: ameloblast maturation-stage 56.9% → 61.3% (P = 0.006); secretory-stage 33.7% → 35.1% (P = 0.084); mesenchyme 9.4% → 3.6% (P = 0.002); other populations below the detection limit of the fit in this dissected-epithelium material.

Limitations, stated explicitly. (i) The reference is postnatal-day-7 incisor whereas the bulk material is molar-stage, so proportions are relative rather than absolute; the CON-versus-dcKO contrast is the robust readout. (ii) Bulk values are TPM and reference values are CPM, an acknowledged unit mismatch. (iii) The reference was annotated de novo rather than using author labels, which are not distributed with the GEO submission. (iv) A stage-matched molar atlas (for example GSE146123) would refine the immune and endothelial estimates, on which the interferon interpretation depends.

#### B5. Reproducibility

The standalone workflow starts from checksum-verified source workbooks and archived gene-level objects. It preserves original DESeq2 results, recomputes descriptive scores and fits the new limma-voom ROAST program model and supplementary sensitivities, then regenerates Figures 4–7, fractions/specificity, pinned-panel enrichment, Data S3 and normalized inhibitor ANOVA/Dunnett exports. Manuscript mode requires the pinned program R library; all mode also rebuilds the GSE146855 reference, re-extracts locus evidence and renders the supplementary diagrams. Use the README setup commands and reproduce.py --data-root /data --work-dir /new/work --destination /new/results --mode all --rscript /path/to/Rscript from outside the project. No sibling GEO directory, data downloads, cached program fit or private helpers are required. --geo-root is optional. original_DESeq2_refit = false; program_limma_voom_fit = true. Full-data validation and synthetic-data CI are distinct.

Code location: https://github.com/LabOnoM/GSE341665. Sequencing files and complete processed matrices: GEO GSE341665 (https://www.ncbi.nlm.nih.gov/geo/query/acc.cgi?acc=GSE341665). The data DOI is https://doi.org/10.5281/zenodo.21637766. The data package contains exact archived-result inputs, source workbooks, gene-level model objects, the annotated reference, auxiliary alignment counts, normalized inhibitor sources and the original Western blot image deck. The editable Figure_2.prism and Figure_3.prism projects (saved in Prism 10.1.1), Data S1 measurements and Data S2 statistical comparisons are retained in the locally validated upload package for manual upload by the authors. CSV exports preserve access without Prism; matching normalized values, ANOVA statistics and all 52 significance classifications were checked against Figure_3.prism. It contains no code. Code and study deposits are restricted pending publication; deposit credentials are supplied to the editor separately, and code access during review may be requested from the corresponding authors. The matching ROAST snapshot is local; the existing DOI does not verify deposition of these revised files. No Zenodo account access or remote operation was performed.

The interactive primary-analysis report at http://smad4kotooth.nexitia.ai/ provides historical context, including exploratory analyses excluded from the current manuscript and code/data packages. It predates ROAST adoption and complements archived inputs and versioned code. Reference re-annotation/Figure S6 requires the recorded R/Seurat environment and GSE146855; BAM-level re-extraction of Figure S3 uses the deposited regional alignments and exact reference-locus sequence with pysam 0.24.0 (bundled samtools 1.24). Archived-evidence Figure S3 rendering uses the deposited reference-locus sequence. The current data snapshot also includes the original reference archive and PanelApp/MGI source records, allowing the retained workflows to run without data downloads.

#### B6. Amelogenesis-imperfecta gene panel, ortholog mapping and enrichment

Panel source. Human amelogenesis-imperfecta (AI) genes were taken from the Genomics England PanelApp panel 269, “Amelogenesis imperfecta”, version 4.37 (panel version dated 24 June 2026; retrieved 27 July 2026 by 17_fetch_ai_panel.py, which explicitly requests the pinned version; offline reproduction reads the declared, checksum-verified original PanelApp/MGI records without an implicit fallback). Only genes rated GREEN (high, diagnostic-grade evidence) were used, giving 31 genes. An external, versioned, publicly retrievable panel was used in preference to a list curated in-house so that the exact gene set is reproducible by a reader. The 15 amber- and red-rated genes on the same panel were excluded and are named here rather than filtered silently: AIRE, ALPL, AMTN, CLDN16, CLDN19, DSPP, ITGB4, KCNJ1, LAMC2, NECTIN4, SMARCD2, SMOC2, TMEM165, TP63, TUFT1. Per-gene provenance (OMIM phenotype MIM, inheritance, primary reference, mouse ortholog) is given in Supplementary Table S4.

Ortholog mapping. Human genes were mapped to mouse orthologs using the Mouse Genome Informatics vertebrate homology report (HOM_MouseHumanSequence), accepting only unambiguous one-to-one homology classes. Three genes required explicit manual resolution for nomenclature reasons, each documented in the analysis code:

| Gene | Resolution | Reason |
| --- | --- | --- |
| AMELX | → Amelx | its MGI homology class also contains the Y-linked paralog AMELY, so the 1:1 filter rejects it |
| C4orf26 | → Odaph | PanelApp retains the retired symbol; HGNC renamed the gene ODAPH |
| FAM83H | → Fam83h | renamed SACK1H by HGNC in February 2026, but carried as Fam83h in the Ensembl release 116 annotation used here |

Matching was performed on mouse symbols. The human-ortholog annotation column of the differential-expression table was deliberately not used: it is empty for 40.9% of tested genes, and because Series.isin() treats a missing value as a non-match, using it silently scored seven genuine, strongly down-regulated AI genes as non-risk. The resolution step therefore raises an error if any panel gene cannot be located, rather than proceeding on a truncated panel.

Enrichment and threshold sensitivity. Over-representation of panel genes among down- versus up-regulated genes was assessed by a one-sided Fisher exact test. The test used adjusted P < 0.05 with no fold-change threshold, because the core enamel genes that carry the signal have |log2 fold change| ≈ 0.4 and are excluded by the fold-change cutoff used to define DEGs. The result is threshold-dependent and is reported as such:

| Threshold | Panel down | Non-panel down | Panel up | Non-panel up | Odds ratio | P (one-sided) |
| --- | --- | --- | --- | --- | --- | --- |
| adjusted P < 0.05, no fold-change cutoff | 18 | 2060 | 2 | 2107 | 9.21 | 1.7 × 10⁻⁴ |
| adjusted P < 0.05 and \|log2FC\| > 0.25 | 15 | 1778 | 2 | 1656 | 6.99 | 1.9 × 10⁻³ |
| adjusted P < 0.05 and \|log2FC\| > 1 (the DEG rule) | 0 | 354 | 1 | 208 | — | 1.00 |

At the |log2 fold change| > 1 cutoff used to define the DEG list, only one panel gene reaches significance and the test is uninformative. The two significant up-regulated panel genes are ODAPH (C4orf26; log2 fold change +1.94) and ORAI1 (+0.72).

#### B7. The conditional allele, resolved from the aligned reads

Because gene-level counts showed Smad4 modestly increased in the mutant libraries (log2 fold change +0.32, adjusted P = 2.7 × 10⁻³) — the opposite of the naive expectation for a knockout — the locus was examined directly in the alignments rather than inferred. Two retained scripts document this analysis: pipeline/20_smad4_loxp_evidence.py derives evidence tables from BAM and reference inputs; pipeline/21_smad4_loxp_figure.py renders Supplementary Figure S3 from those tables and an archived window of the Ensembl-116 reference sequence. Required-input checks stop execution when evidence is missing.

Allele structure. loxP-containing reads are present in all six libraries, carried in the soft-clipped portion of reads that straddle each insertion point (the loxP core is not in the reference genome, so it cannot align). Fifteen such reads were recovered by scanning both orientations of the 34-base loxP core across chr18:73,778,000–73,792,000. One CON1 read (POS 73,784,710, CIGAR 70S79M) carries genomic sequence, a linker, the complete loxP core, a second linker and genomic sequence resuming — fixing that breakpoint unambiguously. Both breakpoints were confirmed against genome.fa.

| Field | Value |
| --- | --- |
| Floxed exon | chr18:73,782,797–73,782,980 (184 bp), entirely coding |
| Exon identity | exon 9 of ENSMUST00000025393 (canonical); exon 8 in the conventional 11-exon numbering |
| Exon accession | ENSMUSE00000143011 in 19 of the 22 transcripts containing this interval; the same interval is labelled ENSMUSE00001829774 (2 transcripts) and ENSMUSE00001829814 (1). The unambiguous identifier is the coordinate interval, not the accession or the ordinal. |
| loxP site A | inserted chr18:73,782,242/73,782,243 (±3 bp; AAC microhomology), intron 9/10 |
| loxP site B | inserted chr18:73,784,709/73,784,710 (exact; one read spans both flanks), intron 8/9 |
| Floxed interval | 2,467 bp, containing exactly one exon |
| Protein consequence | residues 318–379, amino-terminal MH2 domain; 184 mod 3 = 1, so excision is frameshifting and introduces a premature termination codon in the following exon |

Note that loxP reads occur in both genotypes. This is expected: CON animals are Cre-negative Smad4fl/fl and carry the cassette; only Cre is absent.

The excision product is not detected. Excision would splice the flanking exons directly together, creating a junction (chr18:73,782,125–73,790,265) absent from the control genome. That junction has zero supporting reads in all six libraries, against 121–187 uniquely-mapped reads supporting each flanking junction per library (percent-spliced-in = 1.000 in every sample). By the rule of three, ≥99.4 % of dcKO Smad4 mRNA retains the exon. The result is not an aligner artefact: it reproduces in STAR’s SJ.out.tab, by direct parsing of N CIGAR operations, and by an alignment-free k-mer search, and an unbiased scan of every junction across chr18:73,780,000–73,792,000 found no dcKO-specific junction anywhere in the gene.

Exon-level coverage. Each exon’s share of that sample’s total Smad4 exonic signal was computed under two normalisations (per-base depth and read count), both reported so that neither the figure nor the text can select the more favourable one. The floxed exon carries the lowest ratio in the gene (0.92 depth-fraction, 0.889 read-count) and is the only exon whose CON and dcKO values separate completely under both. This is a weak, non-significant indication in the expected direction, not evidence: the effect is 8–11 %, it cannot exceed the exact permutation floor of P = 0.10 at three-versus-three, and exon 6 deviates as far in the opposite direction under the depth normalisation. Coverage cannot in principle settle this question — if the recombined transcript is degraded, the surviving mRNA is entirely wild-type and a genuinely deleted allele still reads near 1.0.

Predicted fate of the recombined transcript. The coding consequence of excision was computed directly from the reference sequence rather than assumed. The floxed exon contributes 184 nucleotides of coding sequence to the canonical transcript (spanning residues 318–379 of the 551-residue protein, of which 319–378 are encoded entirely within it); because 184 is not a multiple of three, its removal shifts the reading frame from the exon 8 | exon 10 junction onward. In the reconstructed mutant coding sequence that junction falls mid-codon (after coding nucleotide 952, within codon 318), and the first in-frame termination codon appears at codon 321 — only 8 nucleotides into exon 10. The predicted product is therefore a 320-residue truncation: residues 1–317 of wild-type SMAD4 followed by three frameshifted residues (Gly-Cys-Thr) and a stop, removing the entire MH2 domain required for Smad complex formation.

This places the premature termination codon 299 nucleotides upstream of the final exon–exon junction, far beyond the 50–55 nucleotide threshold at which a terminating ribosome fails to displace the downstream exon-junction complex. The recombined transcript is therefore a canonical substrate for nonsense-mediated decay and would be expected to be degraded rather than to accumulate. Consistent with this, Ensembl annotates Smad4 transcripts lacking this exon with the nonsense_mediated_decay biotype.

Transcription is not prematurely terminated. An alternative explanation for the absent junction — that the recombined allele terminates transcription at the loxP scar, so the downstream exons are never transcribed — was tested and rejected. Premature termination predicts progressive depletion of the exons downstream of the deletion, with the final exon most affected. The observed dcKO/CON coverage ratios instead recover across that interval (exon 10, 0.937; exon 11, 0.991; exon 12, 0.989), and exon 8, which lies upstream of the deletion and cannot be affected by any downstream terminator, is lower (0.910) than the final exon. Transcription therefore proceeds through the entire locus in both genotypes. This is expected on first principles as well: the loxP cassette is present in the control allele too, so were it a terminator the control libraries would show the same truncation, and they do not.

Interpretation and its limit. Two explanations remain compatible with these observations and cannot be separated by RNA sequencing: nonsense-mediated decay of the frameshifted recombined transcript, which the analysis above shows is expected to be efficient, or incomplete recombination in this short-cultured explant material. We note that decay alone is not a complete explanation, because degradation of a substantial recombined pool would be expected to reduce total Smad4 transcript, whereas it is modestly increased; reconciling efficient recombination with these data would require an additional, untested assumption of transcriptional compensation at the locus. These bulk RNA measurements do not establish the proportion of recombined cells or whether effect sizes are attenuated. Quantifying genomic recombination requires a DNA-level assay across the floxed interval. An excision-junction RT-qPCR assay would test the surviving recombined transcript, and SMAD4 immunostaining would assess protein depletion; neither is a direct measure of the cellular DNA recombination fraction. Note that the RT-qPCR amplicon used in this study lies in exon 2 (Supplementary Table S1) and is therefore also unable to report the deletion.

## Supplementary Figures

Supplementary Figure S1. Mouse models and allocation of experiments. All experimental tooth-germ donors were male. The constitutive K14-cre;Smad4fl/fl arm uses E16.5 germs and no 4-OHT and supplies Figure 1. The inducible K14-creER;Smad4fl/fl arm uses P0 germs; Cre-negative Smad4fl/fl littermates are the genotype controls. Both genotypes receive 90 µmol/L 4-OHT for 24 h in the subsequent experiments, while the preceding dose series has its own concentrations. Mechanical epithelial separation supplies the molecular assays (Figures 2–7); the inhibitor experiments are assigned to this same inducible arm following the confirmed September record. Transplants are recovered after 14 days in eight-week-old male recipients, n = 3 paired recipients per comparison. The five inhibitor-table groups do not identify their individual Cre genotypes, so the scheme does not assert a specific mutant-rescue contrast. Stock 005107 identifies K14-creER and 017462 the Smad4-floxed line; the recorded constitutive-driver stock 016230 conflicts with the catalog and remains unresolved. Supplementary Figure S6 is an external single-cell reference, not a new experiment in these mice.

Supplementary Figure S2. Genotyping records and identification of dental epithelium. (A) Smad4 allele-screening gel. M, DNA ladder; lanes 1–4, individual sample lanes; (+), positive control as labeled; WT, wild-type control. The displayed pattern is consistent with floxed-homozygous samples in lanes 1, 2, and 4, a heterozygous sample in lane 3, and the lower wild-type band in WT. This lane assignment is inferred from the displayed bands, not from an animal-ID key. The listed Smad4 primer pair has been reported to yield approximately 500-bp floxed and 436-bp wild-type products, consistent with this pattern. (B) Gel labeled K14-CreER in the source, with bands visible in lanes 1 and 4 and no corresponding visible bands in lanes 2 and 3. These are the intended positive and negative lanes; a definitive Cre call requires resolution of the primer-record discrepancy in Supplementary Table S1, where the supplied sequences match a 324-bp internal-control assay. (C) Representative tooth-germ photograph and the same image with the selected dental epithelial region outlined. Epithelium was separated mechanically under a stereomicroscope after 24-h culture. No spatial scale is supplied, and none has been inferred. These images document allele screening and tissue selection; they do not quantify recombination or epithelial purity.

Supplementary Figure S3. RNA sequencing identifies the floxed Smad4 allele but cannot determine whether it was excised. It shows what the reads establish and why they cannot settle excision or its efficiency.

Six bulk dental-epithelium libraries from male donors: n = 3 dcKO (K14-creER (tamoxifen-inducible)) versus 3 Cre-negative Smad4fl/fl littermate controls. GRCm39 / Ensembl release 116, chr18; Smad4 is on the minus strand, so transcription runs right → left in A–C.

(A) Locus schematic. The floxed exon (black) is chr18:73,782,797–73,782,980, 184 bp; triangles are the loxP insertion points, site A 73,782,242/73,782,243 (±3 bp) and site B 73,784,709/73,784,710 (exact); shading is the 2,467 bp floxed interval. Dashed = the exon 8 → exon 10 junction Cre would create, never observed.

(B) All 15 reads carrying the 34 bp loxP core, in 3/3 CON and 3/3 dcKO libraries (solid = aligned, open = soft-clipped, black = core; ‡ = a duplicate read, so 14 independent fragments). loxP reads are expected in Cre-negative CON, which carries the identical allele; both cassettes are intronic, so read number is uninformative about excision. Below: the reference sequence at each insertion point, and one cassette-spanning read (CON1, 70S79M at 73,784,710) in literal base calls.

(C) Per-base coverage across the locus, and for each of the 12 canonical exons the dcKO/CON ratio on both normalisations — per-base depth and exonic reads, with each library plotted. All exons lie within 0.889–1.094; the floxed exon is 0.924 on depth (P = 0.10) and 0.889 on reads (P = 0.10), and the metrics rank the exons differently.

(D) Uniquely mapped reads on the two flanking junctions and on the predicted excision junction (chr18:73,782,125–73,790,265), one bar per library. Per-exon P values are in output/smad4_exon_ratio.csv. The Cre excision is not detected: 0 reads, 6/6 libraries.

The limit. The allele is present; the excision product is undetectable — 0 skip reads in 6/6 libraries against 121–187 reads on each flanking junction, PSI = 1.000 everywhere, one-sided 95% bound (rule of three) ≥ 99.4% of dcKO Smad4 mRNA retaining the exon. Coverage cannot reach significance either: P = 0.10 is the floor of the exact 20-permutation test at n = 3 vs 3. Both coverage and junctions are measured on surviving mRNA, so neither can distinguish little recombination from recombination followed by decay of the frameshifted transcript. Genomic PCR across the floxed interval is needed to assess DNA recombination. An excision-junction RNA assay and SMAD4 immunostaining provide complementary transcript- and protein-level evidence.

Supplementary Figure S4. 4-OHT dose series in the inducible model. Contact radiographs of male P0 K14-creER;Smad4fl/fl germs (A, C, E) and Cre-negative Smad4fl/fl littermate controls (B, D, F), exposed to 50, 90, or 120 µmol/L 4-OHT for 24 h before subrenal-capsule transplantation. Both genotypes received the same concentration within each comparison. The source reports three paired recipients per comparison, with grafts recovered after 14 days; this allocation is applied to each illustrated dose comparison. The three displayed concentrations come from a ten-concentration screen spanning 50–300 µmol/L. The 90 µmol/L condition was selected for subsequent experiments. The 120 µmol/L mutant graft shows loss of recognizable tooth architecture, while a control tooth remains visible; this image alone does not establish a genotype-independent toxic effect.

Supplementary Figure S5. Complete analysis workflow. (A) Primary bulk RNA-sequencing pipeline. Steps are labelled with the corresponding section of the Supplementary Methods. The highlighted node is the differential-expression result carried into this manuscript. Historical exploratory analyses are excluded from this workflow; their prior execution and reasons for exclusion are documented in Supplementary Methods A8.

(B) Current manuscript re-analysis and reproduction. Numbers identify retained scripts in the GitHub pipeline/ directory. Step 01 extracts all retained gene statistics and TPM; 03 produces descriptive scores; program_gene_set_tests.R fits the adjusted ROAST model and supplementary sensitivities from archived full-precision inputs; 04/05 obtain and annotate GSE146855; 06 estimates fractions and specificity; 17/16 retrieve the pinned PanelApp panel and compute directional enrichment; and the current step 07 generates Figures 4–7 and the complete Data S3 source CSVs. pipeline/programs.py supplies the single program definition. export_additional_data.py writes the Additional-file data and Table S4. Step 24 independently recalculates the normalized inhibitor statistics described separately. The diagram excludes retired analyses and manuscript/deposit publishing scripts.

Supplementary Figure S6. The single-cell reference used for deconvolution. (A) UMAP embedding of 6,248 quality-controlled cells from GSE146855 (postnatal day 7 mouse incisor), coloured by the de-novo cell-type annotation. (B) Dot plot of the canonical markers used to assign each cell type; dot size indicates the percentage of cells expressing the gene and colour indicates scaled mean expression. Enamel-matrix transcripts (Amelx, Ambn, Enam, Mmp20) are detected at low level in every cluster, reflecting ambient RNA from this highly secretory tissue — the reason cell types were assigned by cross-cluster specificity rather than absolute expression.

Supplementary Figure S7. Enamel phenotype in the inducible model. Male P0 K14-creER;Smad4fl/fl and Cre-negative control germs were both exposed to 90 µmol/L 4-OHT for 24 h, transplanted beneath the renal capsule, and recovered after 14 days. Hematoxylin-eosin sections are shown at increasing magnification: C and D enlarge the boxed regions in A and B; E and F enlarge the boxed regions in C and D. n = 3 paired recipients, with one experimental and one control graft per recipient. am, ameloblasts; e, enamel space; d, dentin; arrows identify the epithelial/enamel region. Scale bars: A and B, 500 µm; C and D, 100 µm; E, 50 µm. F is the paired higher-magnification image, but no separate scale bar is visible in the supplied panel.

## Supplementary Tables

Supplementary Table S1. Primer sequences for qPCR and genotyping.

| Target | Application | Forward (5′–3′) | Reverse (5′–3′) |
| --- | --- | --- | --- |
| Amelx | RT-qPCR | ACGGCTGCCTTATCATGCTCT | ATGGGGACCTGGATTTTGTTT |
| Ambn | RT-qPCR | GAATACTGAGAAAGTGCGTTGA | TGATCCTGTTCCTGTCCCTAG |
| Mmp20 | RT-qPCR | CTGGGAAGTTGGACCAGAATAC | AGGCGGTAGTTAGCCACATCA |
| Klk4 | RT-qPCR | CGGGAGTCTTGGTGCATCC | CTTGGGAGCCTTTCAGGTTATG |
| Amtn | RT-qPCR | ATCAGCCCAGTCATTACCAAAG | AGGTCTGACCCCAGAGTGAG |
| Odam | RT-qPCR | GTCACATCCTCACCACAGCA | GAGTTTCTGGAGCTGTGCCT |
| Enam | RT-qPCR | TGCAGAAATCCGACTTCTCCT | CATCTGGAATGGCATGGCA |
| Tuft1 | RT-qPCR | CTGAGGTGGACCCTAGCATGA | CTCGCTGATAACGGCTGAGTG |
| Gapdh | RT-qPCR | TCCAGAACATCATCCCTGCCTCTA | ACAAAGTGGTCGTTGAGGGCAATG |
| Smad4 | RT-qPCR  | ACACCAACAAGTAACGATGCC | GCAAAGGTTTCACTTTCCCCA |
| Smad4 | Genotyping PCR | TAAGAGCCACAGGGTCAAGC | TTCCAGGAAAAACAGGGCTA |
| K14-Cre (source label; see note 2) | Genotyping PCR | CTAGGCCACAGAATTGAAAGATCT | GTAGGTGGAAATTCTAGCATCATCC |

Sequence-identity sources: JAX internal-control protocol, p. 4; published Smad4 primer protocol.

Supplementary Table S2. Culture and inhibitor regimens. The September Methods and confirmed author note govern the concentrations and exposure times.

| Reagent | Role | Final concentration | Duration | Supplier and catalog |
| --- | --- | --- | --- | --- |
| 4-Hydroxytamoxifen | Inducible-model treatment, both genotypes | 90 µmol/L | 24 h | Sigma-Aldrich H7904-5MG |
| SB203580 | p38α/β inhibitor | 10 µmol/L | 24 h | Cell Signaling Technology 5633 |
| Dorsomorphin | BMP type I receptor inhibitor | 10 µmol/L | 24 h | Sigma-Aldrich P5499 |
| SB525334 | TGF-β type I receptor/ALK5 inhibitor | 1 µmol/L | 24 h | Tocris 3211 |
| U0126 | MEK1/2 inhibitor | 10 µmol/L | 24 h | Sigma-Aldrich U120 |
| DMSO | Vehicle | Matched across groups | 24 h | Not specified |

The dose-series images show 50, 90, and 120 µmol/L 4-OHT, distinct from the subsequent 90 µmol/L working condition. Figure 1 uses E16.5 constitutive-model germs (Smad4-cKO) without 4-OHT; the remaining experimental material uses P0 inducible-model germs (Smad4-dcKO).

Supplementary Table S3. Primary antibodies.

| Antibody as supplied | Host | Dilution | Application | Supplier | Catalog |
| --- | --- | --- | --- | --- | --- |
| Phospho-p38 MAPK | Rabbit | 1:1,000 | Western blot | Cell Signaling Technology | #9212 |
| Phospho-SMAD1/5/9 | Rabbit | 1:1,000 | Western blot | Cell Signaling Technology | #13820 |
| Phospho-SMAD2/3 | Rabbit | 1:1,000 | Western blot | Thermo Fisher | #PA5-110155 |
| AMELX | Rabbit | 1:500 | Western blot | Thermo Fisher | #PA5-114845 |
| AMBN | Rabbit | 1:500 | Western blot | Bioss antibodies | #bs-12467R |
| MMP20 | Rabbit | 1:1,000 | Western blot | Proteintech | #55467-1-AP |
| KLK4 | Rabbit | 1:1,000 | Western blot | Gentaur | #8941 |
| β-Tubulin | Rabbit | 1:1,000 | Western blot | Thermo Fisher | #600-401-880 |
| GAPDH | Rabbit | 1:1,000 | Western blot | Cell Signaling Technology | #2118 |
| SMAD4 | Rabbit | 1:1,000 | Western blot | Proteintech | #51069-2-AP |

Catalog checks: CST 9212; Thermo Fisher 600-401-880.

Supplementary Table S4. Verified provenance of the 31 genes used in the amelogenesis-imperfecta enrichment. Genes are the GREEN-rated (high-evidence) members of the Genomics England PanelApp “Amelogenesis imperfecta” panel (269 v4.37; retrieved 27 July 2026). Each OMIM phenotype MIM was confirmed to be a phenotype (not gene) entry mapping to the stated gene, and each primary PMID was confirmed to resolve to a real report of that gene’s enamel or amelogenesis-imperfecta phenotype. Mouse orthologs are as carried in the Ensembl release 116 annotation. Two entries are flagged: PEX26 has no enamel-specific publication (the cited paper is the peroxisome-biogenesis-disorder primary report), and PLXNB2 has no OMIM phenotype MIM assigned yet.

| Gene (human) | Mouse | OMIM phenotype MIM | Phenotype | Inheritance | Primary reference (PMID) | Notes |
| --- | --- | --- | --- | --- | --- | --- |
| ACP4 | Acp4 | 617297 | Amelogenesis imperfecta, type IJ | AR | Seymen, Am J Hum Genet 2016 (PMID 27843125) |  |
| AMBN | Ambn | 616270 | Amelogenesis imperfecta, type IF (hypoplastic) | AR | Poulter, Hum Mol Genet 2014 (PMID 24858907) |  |
| AMELX | Amelx | 301200 | Amelogenesis imperfecta, type IE (X-linked) | XL | Lagerström, Genomics 1991 (PMID 1916828) |  |
| C4orf26 | Odaph | 614832 | Amelogenesis imperfecta, type IIA4 (hypomaturation) | AR | Parry, Am J Hum Genet 2012 (PMID 22901946) | gene renamed ODAPH |
| CNNM4 | Cnnm4 | 217080 | Jalili syndrome (cone-rod dystrophy + AI) | AR | Parry, Am J Hum Genet 2009 (PMID 19200525) | syndromic |
| COL17A1 | Col17a1 | 619787 | Amelogenesis imperfecta; junctional EB spectrum | AD | Hany, J Med Genet 2024 (PMID 37979963) | also causes junctional EB |
| DLX3 | Dlx3 | 190320 | Tricho-dento-osseous syndrome | AD | Price, Hum Mol Genet 1998 (PMID 9467018) | syndromic |
| ENAM | Enam | 104500 | Amelogenesis imperfecta, type IB (hypoplastic) | AD | Rajpar, Hum Mol Genet 2001 (PMID 11487571) |  |
| FAM20A | Fam20a | 204690 | Amelogenesis imperfecta, type IG (enamel-renal syndrome) | AR | O’Sullivan, Am J Hum Genet 2011 (PMID 21549343) |  |
| FAM20C | Fam20c | 259775 | Raine syndrome (with enamel/dental defects) | AR | Simpson, Am J Hum Genet 2007 (PMID 17924334) | syndromic; enamel defects in survivors |
| FAM83H | Fam83h | 130900 | Amelogenesis imperfecta, type IIIA (hypocalcified) | AD | Kim JW, Am J Hum Genet 2008 (PMID 18252228) | gene renamed SACK1H (2026) |
| GPR68 | Gpr68 | 617217 | Amelogenesis imperfecta, type IIA6 (hypomaturation) | AR | Parry, Am J Hum Genet 2016 (PMID 27693231) |  |
| ITGB6 | Itgb6 | 616221 | Amelogenesis imperfecta, type IH (hypoplastic) | AR | Wang SK, Hum Mol Genet 2014 (PMID 24305999) |  |
| KLK4 | Klk4 | 204700 | Amelogenesis imperfecta, type IIA1 (hypomaturation) | AR | Hart PS, J Med Genet 2004 (PMID 15235027) |  |
| LAMA3 | Lama3 | 619786 | Junctional epidermolysis bullosa; AI-spectrum enamel defect | AR | Wang SK, Oral Dis 2023 (PMID 36326426) | also causes junctional EB |
| LAMB3 | Lamb3 | 104530 | Amelogenesis imperfecta, type IA (hypoplastic) | AD | Kim JW, J Dent Res 2013 (PMID 23958762) | also causes junctional EB |
| LTBP3 | Ltbp3 | 601216 | Dental anomalies and short stature (DASS) | AR | Noor, Am J Hum Genet 2009 (PMID 19344874) | syndromic |
| MMP20 | Mmp20 | 612529 | Amelogenesis imperfecta, type IIA2 (hypomaturation) | AR | Kim JW, J Med Genet 2005 (PMID 15744043) |  |
| ORAI1 | Orai1 | 612782 | CRAC channelopathy (immunodeficiency 9) with defective enamel | AR | McCarl, J Allergy Clin Immunol 2009 (PMID 20004786) | syndromic (recessive LOF) |
| PEX1 | Pex1 | 234580 | Heimler syndrome 1 (AI + hearing loss) | AR | Ratbi, Am J Hum Genet 2015 (PMID 26387595) | syndromic |
| PEX26 | Pex26 | 614872 | Peroxisome biogenesis disorder 7A (Zellweger spectrum) | AR | Matsumoto, Am J Hum Genet 2003 (PMID 12851857) | PBD primary report; no enamel-specific paper — enamel inferred from spectrum |
| PEX6 | Pex6 | 616617 | Heimler syndrome 2 (AI + hearing loss) | AR | Ratbi, Am J Hum Genet 2015 (PMID 26387595) | syndromic |
| PLXNB2 | Plxnb2 | — | AI, hearing loss and intellectual disability | AR | Smith CEL, J Med Genet 2024 (PMID 38458752) | recently reported; no OMIM phenotype MIM yet |
| RELT | Relt | 618386 | Amelogenesis imperfecta, type IIIC | AR | Kim JW, Clin Genet 2019 (PMID 30506946) |  |
| ROGDI | Rogdi | 226750 | Kohlschütter-Tönz syndrome (epilepsy + AI) | AR | Schossig, Am J Hum Genet 2012 (PMID 22424600) | syndromic |
| SLC10A7 | Slc10a7 | 618363 | Short stature, AI and skeletal dysplasia (SSASKS) | AR | Dubail, Nat Commun 2018 (PMID 30082715) | syndromic |
| SLC13A5 | Slc13a5 | 615905 | Developmental and epileptic encephalopathy 25 with AI | AR | Schossig, J Med Genet 2017 (PMID 27600704) | syndromic |
| SLC24A4 | Slc24a4 | 615887 | Amelogenesis imperfecta, type IIA5 (hypomaturation) | AR | Parry, Am J Hum Genet 2013 (PMID 23375655) |  |
| SP6 | Sp6 | 620104 | Amelogenesis imperfecta, type IK (hypoplastic) | AD | Smith CEL, Hum Mol Genet 2020 (PMID 32167558) | human evidence recent and limited |
| STIM1 | Stim1 | 612783 | CRAC channelopathy (immunodeficiency 10) with enamel defect | AR | Picard, N Engl J Med 2009 (PMID 19420366) | syndromic (recessive LOF) |
| WDR72 | Wdr72 | 613211 | Amelogenesis imperfecta, type IIA3 (hypomaturation) | AR | El-Sayed, Am J Hum Genet 2009 (PMID 19853237) |  |

Table S4 footnote — abbreviations and conventions.

Columns. Gene (human) — HGNC-approved human gene symbol, as listed by PanelApp. Mouse — the MGI-approved mouse ortholog symbol (italic), as carried in the Ensembl release 116 annotation used for the differential-expression analysis; this is the symbol on which genes were matched. OMIM phenotype MIM — the Mendelian Inheritance in Man (OMIM) identifier of the phenotype entry, not the gene entry; each was confirmed to be of type “phenotype” and to map to the stated gene in NCBI mim2gene_medgen. Phenotype — the OMIM phenotype designation. Inheritance — mode of inheritance for that phenotype. Primary reference (PMID) — the report establishing the gene-phenotype link, with its PubMed identifier (PMID); each was confirmed to exist and to concern the gene’s enamel or amelogenesis-imperfecta phenotype. Notes — caveats affecting how the row should be read.

Inheritance codes. AD, autosomal dominant; AR, autosomal recessive; XL, X-linked.

Disease and phenotype abbreviations. AI, amelogenesis imperfecta; EB, epidermolysis bullosa (junctional EB, the junctional subtype); CRAC, calcium release-activated calcium channel; DASS, dental anomalies and short stature; SSASKS, short stature, amelogenesis imperfecta and skeletal dysplasia with scoliosis; PBD, peroxisome biogenesis disorder; LOF, loss of function.

Amelogenesis-imperfecta subtype nomenclature. OMIM numbers AI subtypes by a Roman numeral for the defect class followed by a letter (and, where several genes cause the same class, a digit): type I = hypoplastic (too little enamel formed, of normal mineralisation); type II = hypomaturation (full-thickness but incompletely matured, soft enamel); type III = hypocalcified (full-thickness but poorly mineralised). Hence type IA, IB, IE, IF, IG, IH, IJ and IK are hypoplastic forms; types IIA1-IIA6 are hypomaturation forms; and types IIIA and IIIC are hypocalcified forms. Where a parenthetical descriptor follows the subtype it restates this class for readability; rows without one follow the same Roman-numeral rule.

Gene-symbol changes. Two genes are listed under the symbol PanelApp uses, with the current HGNC symbol given in Notes: C4orf26 is now ODAPH (odontogenesis-associated phosphoprotein), and FAM83H was renamed SACK1H (scaffolding CK1 anchoring protein H) in February 2026. The mouse column retains Odaph and Fam83h because those are the symbols present in the Ensembl release 116 annotation on which matching was performed.

Symbols. An em dash (—) in the OMIM column means no phenotype MIM has been assigned yet; this applies only to PLXNB2, whose AI phenotype was first reported in 2024.

Note-column meanings. “syndromic” — enamel involvement occurs as one feature of a broader syndrome rather than as isolated AI. “also causes junctional EB” — the same gene causes both an AI-spectrum enamel defect and junctional epidermolysis bullosa. “recessive LOF” — the enamel phenotype belongs to the recessive loss-of-function disorder, not to the dominant gain-of-function disorder caused by the same gene. “PBD primary report; no enamel-specific paper” — for PEX26 the cited reference establishes the peroxisome biogenesis disorder; enamel involvement is inferred from the disorder spectrum rather than from a dedicated study. “human evidence recent and limited” — for SP6 the human AI association rests on one dominant family with functional corroboration plus one de novo case.

Panel definition. All 31 genes are the GREEN-rated (“high evidence”, diagnostic-grade) members of the Genomics England PanelApp panel 269, “Amelogenesis imperfecta”, version 4.37, retrieved 27 July 2026. Amber- and red-rated genes on the same panel were excluded and are named in Supplementary Methods B6.

## Supplementary data and repository files

| File | Contents |
| --- | --- |
| Data S1 (Zenodo CSV) | Supplementary_Data_S1_inhibitor_normalized_values.csv: 195 normalized inhibitor measurements, also present in Figure_3.prism. Sample indices are within-table indices, not animal IDs. |
| Data S2 (Zenodo CSV) | Supplementary_Data_S2_inhibitor_ANOVA_Dunnett.csv: 52 treatment-versus-vehicle comparisons with ANOVA statistics, group means/SDs and independently recalculated two-sided Dunnett P values. Native Prism analyses remain in Figure_3.prism. |
| Additional_file_2.xlsx (Data S3) | Complete retained-gene statistics, per-library TPM, primary ROAST tests, supplementary sensitivities and provenance and manuscript gene indexes; nine accompanying CSV tables are supplied in GitHub and Zenodo. |
| Additional_file_3.csv (Data S4) | Log2 fold changes and adjusted P values for the 31 GREEN-rated genes of the Genomics England PanelApp “Amelogenesis imperfecta” panel (269 v4.37). |
| Additional_file_4.csv (Data S5) | All 563 differentially expressed genes (adjusted P < 0.05, absolute log2 fold change > 1) with DESeq2 statistics and human homolog annotation. |
| Additional_file_5.csv (Data S6) | Per-sample program scores for the four gene programs. |
| Additional_file_6.csv (Data S7) | Per-sample cell-type proportions estimated by non-negative least squares. |
| Additional_file_7.csv (Data S8) | Cell-type specificity of each gene program in the single-cell reference. |

Complete gene-level count and transcripts-per-million matrices are associated with GEO GSE341665. The GEO TPM export is rounded to four decimal places; the archived source workbooks and plotting inputs retain full precision. The data-file map and input checksums are supplied in the GitHub repository.

Additional file 2 (Supplementary Data S3). Complete gene-level expression evidence for Figures 4–7 and the main text. The workbook contains all 15,803 retained gene rows with fitted DESeq2 log2 fold changes, raw and BH-adjusted P values, six-library TPM values, program membership and panel mapping; focused displayed-gene and text-claim indexes; primary ROAST tests, supplementary sensitivity/model/membership/sample information; and descriptive deconvolution comparisons. Auxiliary alignment feature counts for Sox2 and Lgr5 are kept in a separate sheet and are not the Salmon/tximport DESeq2 inputs. Zero adjusted P values result from numerical underflow; missing results remain missing. The nine accompanying machine-readable CSV tables are hosted in the GitHub repository (results/source_csv/) and included in the staged Zenodo data package. The original DESeq2 model is not refitted; the new program limma-voom model is fitted. The matching data snapshot is prepared locally for manual Zenodo upload.
