Additional file 1

## Supplementary Information

**Epithelial Smad4 loss reprograms ameloblast identity and metabolism**

Ziyi Wang1, Xinyu Zheng2,3, Mitsuaki Ono4, Dian Su3, Fei Hu3, Xianghui Liu3,5, Zhiwei He6**, Xiaohua Xie3*

1 Department of Advanced International and Information Dentistry, Faculty of Medicine, Dentistry and Pharmaceutical Sciences, Okayama University, Okayama, Japan

2 Department of Tissue and Developmental Biology, Graduate School of Dentistry, The University of Osaka, Suita, Osaka, Japan

3 Department of Stomatology, Second Affiliated Hospital of Harbin Medical University, Harbin, China

4 Department of Oral Rehabilitation and Regenerative Medicine, Okayama University Graduate School of Medicine, Dentistry and Pharmaceutical Sciences, Okayama, Japan

5 Qingdao Stomatology Hospital, Qingdao, China

6 Department of Stomatology, Shenzhen University General Hospital, Shenzhen University, Shenzhen, Guangdong, China

*Correspondence: Xiaohua Xie (xiexiaohua_79@126.com); **Correspondence: Zhiwei He (3563570@qq.com).

## Contents

**Supplementary Methods** — Part A: primary bulk RNA-sequencing pipeline; Part B: re-analysis performed for this manuscript

**Supplementary Figure S1** — mouse lines, crosses, and the allocation of each experiment to a figure

**Supplementary Figure S2** — genotyping records and dental epithelial region

**Supplementary Figure S3** — read-level evidence at the conditional Smad4 locus

**Supplementary Figure S4** — 4-OHT dose series in the inducible model

**Supplementary Figure S5** — analysis workflow (A: RNA sequencing; B: gene-program and reference-based analyses)

**Supplementary Figure S6** — single-cell reference used for deconvolution (UMAP and marker dot plot)

**Supplementary Figure S7** — the inducible K14-creER driver reproduces the enamel phenotype

**Supplementary Table S1** — qPCR and genotyping primers

**Supplementary Table S2** — culture and inhibitor regimens

**Supplementary Table S3** — **primary **antibod**i**es

**Supplementary Table S4** — **hum**an A**I** gene**s and mouse ortholog**s used in the enrichment** analysis**

**Supplementary Data S1–S8** — machine-readable result tables

## Supplementary Methods

This document records the workflow, parameters, and software versions for the computational analyses in two parts: **Part A**, the primary bulk RNA-sequencing pipeline; and **Part B**, gene-program scoring and ROAST inference, cell-type deconvolution, AI-gene enrichment, and Smad4 allele analysis.

Part A — Primary bulk RNA-sequencing pipeline

A1. Experimental design

| Field | Value |
| --- | --- |
| Samples | 6: CON1, CON2, CON3 (control) and dcKO1, dcKO2, dcKO3 (Smad4 conditional knockout) |
| Genotype | dcKO: K14-creER; Smad4fl/fl (tamoxifen-inducible epithelial deletion in P0 tooth-germ cultures). CON: Cre-negative Smad4fl/fl littermates. The constitutive K14-cre line was used for the in vivo enamel phenotype (Figure 1), not for the sequenced material. |
| Groups | 2 groups × 3 biological replicates |
| Sex | All six libraries were derived from male donors. Each biological replicate pooled tooth germs from more than one pup. |
| Tissue | Dental epithelium (mandibular molar tooth germ) |
| Layout | Paired-end, 150 bp |
| Library strandedness | Reverse-stranded; automatically detected by Salmon as ISR, and used as -s 2 in featureCounts |
| Contrast | dcKO vs CON |

All tooth-germ donors were male. Expression of Xist, Ddx3y, Eif2s3y, Uty and Kdm5d is provided in Additional file 2 (Data S3) as supporting sample information.

A2. Reference genome and annotation

| Field | Value |
| --- | --- |
| Assembly | Mus musculus GRCm39,Mus_musculus.GRCm39.dna.primary_assembly.fa |
| Annotation | Mus_musculus.GRCm39.116.gtf |
| Source | Ensembl release 116 (ftp://ftp.ensembl.org/pub/release-116/) |
| Annotation database | org.Mm.eg.db; biomaRt dataset mmusculus_gene_ensembl; KEGG organism mmu; UniProt taxon 10090 |

The GTF was filtered before index building (genome_filtered.gtf); genes lacking a gene_name attribute retained their gene_id.

A3. Read quality control and trimming

**FastQC** on the provider-supplied reads, on trimmed reads, and on aligned BAM files.

The sequencing provider (LC-Bio) **s**u**ppli**ed adapter- and quality-filtered **pa**i**r**ed**-end** **reads (***_Clean_Data[12].fastq.gz**)**. These files **we**re deposited in GEO and **used as input for further trimming **w**it**h Trim Galore.

**Trim Galore v0.6.10** (wrapping **cutadapt v4.9**) with default settings: quality Phred cutoff 20, ASCII+33 encoding, automatic adapter detection.

A4. Alignment (STAR)

**STAR version 2.7.10b.**

Index construction:

STAR --runThreadN 32 --runMode genomeGenerate \     --genomeDir ./GenomeIndex_ENSEMBL_STAR_ReadLength149/ \     --genomeFastaFiles genome.fa --sjdbGTFfile genome_filtered.gtf \     --sjdbOverhang 149 --limitGenomeGenerateRAM 101781612134

Alignment (per sample; two-pass):

STAR --runMode alignReads --chimSegmentMin 10 \     --genomeDir ./GenomeIndex_ENSEMBL_STAR_ReadLength149 \     --readFilesIn <sample>_Clean_Data1_val_1.fq.gz <sample>_Clean_Data2_val_2.fq.gz \     --outFileNamePrefix ./<sample>.STAR. --alignEndsType Local \     --outSAMtype BAM SortedByCoordinate --outFilterMultimapNmax 20 \     --outFilterMismatchNmax 999 --outFilterMismatchNoverLmax 0.04 \     --alignIntronMin 70 --alignIntronMax 1000000 --alignMatesGapMax 1000000 \     --alignSJoverhangMin 8 --alignSJDBoverhangMin 1 \     --outSAMstrandField intronMotif --outFilterType BySJout \     --twopassMode Basic --outSAMmapqUnique 60 \     --limitGenomeGenerateRAM 101781612134 --runThreadN 32 --readFilesCommand zcat

**Mapping statistics (from *****.STAR.Log.final.out****):**

| Sample | Input read pairs | Uniquely mapped | Uniquely mapped % | Multi-mapped % | Chimeric % |
| --- | --- | --- | --- | --- | --- |
| CON1 | 21,303,455 | 19,975,969 | 93.77 | 3.85 | 0.30 |
| CON2 | 21,235,660 | 19,900,363 | 93.71 | 3.80 | 0.29 |
| CON3 | 20,986,103 | 19,694,042 | 93.84 | 3.81 | 0.28 |
| dcKO1 | 22,453,825 | 21,044,062 | 93.72 | 3.79 | 0.30 |
| dcKO2 | 21,439,266 | 20,080,247 | 93.66 | 3.74 | 0.31 |
| dcKO3 | 22,455,897 | 21,031,771 | 93.66 | 3.77 | 0.28 |

Mapping was uniform across samples (uniquely mapped 93.66–93.84%), with no genotype-associated difference in alignment quality.

A5. Quantification

Two quantifications were produced:

**featureCounts** (gene-level counts from the STAR BAMs): featureCounts -T 32 -s 2 -p -Q 20 -O -a ./genome_filtered.gtf -o <sample>.STAR.counts (-s 2 reverse-stranded; -p paired-end; -Q 20 minimum mapping quality; -O assign reads overlapping multiple features.)

**Salmon v1.10.0** selective alignment (--libType A, auto-detected as **ISR** in every sample) against an index built from a **StringTie**-merged transcriptome (a de-novo assembly merged with the reference annotation); library type auto-detected as **ISR**.

**The differential-expression analysis used the Salmon quantifications imported with tximport** (pipeline parameter CountsFrom = "Tximport"). Transcripts per million (TPM) reported in this study derive from that import.

A6. Differential expression (DESeq2 with RUV correction)

Differential expression was analyzed with DESeq2 using the design and parameters below.

| Field | Value |
| --- | --- |
| Package | DESeq2 |
| Counts | Salmon → tximport |
| Unwanted-variation correction | RUVSeq, one factor (k = 1), included as the sample covariate W_1. Per-library values are provided in Additional file 2 (Data S3). |
| Design formula | ~ W_1 + RUV_x (where RUV_x is genotype: CON, dcKO) |
| Test | Likelihood-ratio test (LRT), full ~ W_1 + RUV_x vs reduced ~ W_1 (Par_Test = "LRT") |
| Contrast | RUV_x, dcKO vs CON |
| Expression filter | Genes with at least 10 reads in at least 3 samples (Par_ReN = 10, Par_SaN = 3) |
| Genes tested | 15,803 |
| Significance | padj < 0.05 (Benjamini–Hochberg) and \|log2 fold change\| > 1 (Par_FDR = 0.05, Par_LFC = 1, Par_PorFRR = "padj") |
| Log2 fold-change shrinkage | Not applied (fold changes are maximum-likelihood estimates) |
| Independent filtering | Enabled; 613 genes (3.9%) have padj = NA and are treated as non-significant |
| Result | 563 differentially expressed genes: 209 up-regulated, 354 down-regulated |

Because the LRT compares the full and reduced models, the reported *P* values test the overall genotype effect after adjusting for W_1; the accompanying log2 fold changes are the maximum-likelihood estimates for RUV_x dcKO vs CON.

A7. Functional enrichment

Over-representation analysis (ORA) and GSEA against **GO** (biological process, molecular function, cellular component) and **KEGG** (organism mmu), using the org.Mm.eg.db annotation package in R.

GO and KEGG enrichment used an unadjusted P-value threshold of 0.05. These analyses were exploratory; program-level inference used the ROAST tests described in B2.

Sample clustering/heatmaps used Pearson distance with **ward.D2** linkage (Par_dist = "pearson", Par_hclust = "ward.D2").

Part B — Gene-program and reference-based analyses

Analysis and figure-generation code, software versions, input-file checksums, and execution instructions are provided at https://github.com/LabOnoM/GSE341665.

B1. Software environment

| Component | Version |
| --- | --- |
| Python | 3.12.14 |
| numpy / pandas / scipy / matplotlib | 2.3.5 / 2.2.3 / 1.18.1 / 3.11.2 |
| R | 4.3.2 |
| Seurat / SeuratObject | 5.3.1 / 5.2.0 |
| Matrix / uwot / ggplot2 / patchwork | 1.6.5 / 0.2.3 / 4.0.0 / 1.3.2 |
| limma / edgeR / statmod | 3.58.1 / 4.0.16 / 1.5.0 |
| Operating system | macOS, arm64; software locks and session information in the repository |

B2. Descriptive scoring and ROAST gene-program inference

The four gene sets are defined in pipeline/programs.py and listed below.

*Krt76* was analyzed as an individual transcript and was not included in the ameloblast-identity gene set. *Mmp12* was included in the wound/stress gene set; its expression may include contributions from non-epithelial cells.

**Gene-set membership** and **supporting liter**a**tu**re.

| Program | n | Members | Basis |
| --- | --- | --- | --- |
| Ameloblast identity | 16 | Amelx, Ambn, Enam, Amtn, Mmp20, Klk4, Acp4, Lamc2 (enamel-associated proteins and proteases); Dlx3, Sp6, Shh, Gli1, Wnt10a, Trp63, Nfic, Sp7 (differentiation-associated transcription factors and signaling genes) | Lacruz et al. 2017; Bartlett 2013; Balic and Thesleff 2015; Zhang et al. 2015; Rhodes et al. 2021; Lee et al. 2018 |
| Sterol biosynthesis | 20 | Hmgcr, Hmgcs1, Mvd, Fdps, Sqle, Msmo1, Idi1, Ldlr, Insig1 (nine experimentally defined direct SREBP targets), Srebf2, and Fdft1, Lss, Cyp51, Dhcr7, Dhcr24, Nsdhl, Mvk, Pmvk, Acat2, Sc5d | Horton et al. 2003 |
| Wound / stress | 14 | Mmp10, Mmp12, Mmp13, Ereg, Areg, Ptgs2, Il1a, Cxcl2, Csf3, Fosl1, Mt1, Mt2, Slc7a11, Sprr1a | Werner and Grose 2003; Ge et al. 2017 |
| Interferon-stimulated | 18 | Irf7, Irf9, Stat1, Isg15, Mx1, Oasl1, Oas1a, Oas2, Oas3, Ifit1, Ifit2, Ifit3, Ifi44, Ifi27l2a, Rsad2, Usp18, Bst2, Ccl5 | Interferome (Rusinova et al. 2013) |

For each sample, a **program score** was computed as the mean of gene-wise *z*-scores of log2(TPM + 1) across program members — the combined gene-set *z*-score of Lee et al. (2008). These six-library values remain descriptive; ROAST tests gene expression through the fitted model rather than these averaged scores.

**Primary inference.** ROAST (Wu et al., 2010; doi:10.1093/bioinformatics/btq401) tested differential expression within each gene set while preserving inter-gene correlation. Counts, full-precision TPM, effective transcript lengths and W_1 were obtained from DESeq2_Processing.Rdata for the 15,803 retained genes. All 68 gene-set members mapped uniquely: ameloblast identity, 16; sterol biosynthesis, 20; wound/stress response, 14; and interferon response, 18. Counts were reconstructed as TPM × mean effective transcript length and rescaled to the retained count totals, without a second length correction. TMM normalization (Robinson and Oshlack, 2010; doi:10.1186/gb-2010-11-3-r25) and voom precision weights (Law et al., 2014; doi:10.1186/gb-2014-15-2-r29) preceded limma fitting (Ritchie et al., 2015; doi:10.1093/nar/gkv007). The primary design was ~ W_1 + genotype, with dcKO minus CON as the contrast and three residual degrees of freedom. limma::mroast used the mean statistic, a two-sided directional test, 99,999 rotations, seed 20261003 and midp = FALSE. Benjamini–Hochberg adjustment was applied across the four gene sets. The score means below are descriptive; ROAST probabilities were obtained from the fitted model.

**Sensitivity analyses.** Genotype-only ROAST used ~ genotype, four residual degrees of freedom and seed 20261004. Mixed-direction ROAST was performed for both designs. Each model and endpoint used a separate four-gene-set Benjamini–Hochberg adjustment. The four primary and 12 sensitivity results, gene-set membership and sample covariates are provided in Additional file 2 (Data S3). All analyses used the same six biological libraries.

Gene sets were curated after inspection of the differential-expression results, and program analyses were exploratory. Analyses used R 4.3.2, limma 3.58.1, edgeR 4.0.16 and statmod 1.5.0. At 99,999 rotations, the minimum reported raw P value was 0.00001. Per-sample W_1 values and gene-set membership are provided in Additional file 2 (Data S3).

| Program | n genes | CON | dcKO | ROAST BH-adjusted P |
| --- | --- | --- | --- | --- |
| Ameloblast identity | 16 | +0.90 | −0.90 | 2.67e-05 |
| Sterol biosynthesis | 20 | −0.91 | +0.91 | 2.00e-05 |
| Wound / stress | 14 | −0.91 | +0.91 | 2.00e-05 |
| Interferon-stimulated | 18 | +0.88 | −0.88 | 5.00e-05 |

B3. Single-cell reference and annotation

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

**Annotation strategy.** Enamel-matrix transcripts (*Amelx*, *Ambn*, *Enam*) were detected across reference cell types, consistent with ambient RNA contamination (Supplementary Figure S6B). Annotation was therefore based on relative marker expression across clusters. For each marker gene, mean expression was converted to z-scores across clusters; these z-scores were then averaged within each marker set. Marker-set scores were rounded to two decimals, and each cluster was assigned to the cell type with the highest score, using the first maximum in marker-set order to resolve ties. No ambient-RNA correction was applied.

Marker sets used: secretory ameloblast (*Amelx, Ambn, Enam, Mmp20*); maturation ameloblast (*Klk4, Odam, Amtn, Slc24a4, Gm17660*); other dental epithelium (*Shh, Sfrp5, Krt14, Krt5, Epcam, Thbd, Gjb3, Cldn10*); mesenchyme (*Col1a1, Dcn, Lum, Prrx1, Pdgfra, Msx1, Sp7, Dspp*); immune (*Ptprc, Lyz2, C1qa, Cd68, Fcer1g*); endothelial (*Pecam1, Cdh5, Egfl7, Emcn*); erythroid (*Hba-a1, Hbb-bs, Hba-a2*).

**Resulting reference composition** (Supplementary Figure S6):

| Cell type | Cells | % |
| --- | --- | --- |
| Mesenchyme | 2,527 | 40.4 |
| Epithelial (inner enamel epithelium / stratum intermedium) | 1,361 | 21.8 |
| Ameloblast, maturation-stage | 1,346 | 21.5 |
| Ameloblast, secretory-stage | 506 | 8.1 |
| Endothelial | 290 | 4.6 |
| Immune | 131 | 2.1 |
| Erythroid | 87 | 1.4 |

B4. Cell-type deconvolution

| Field | Value |
| --- | --- |
| Reference signature | Pseudobulk counts per annotated cell type, converted to counts per million |
| Marker panel | Genes with reference expression > 20 CPM in a cell type and > 2-fold above the cross-cell-type mean; up to 120 per cell type; union retained and intersected with the bulk matrix (542 genes) |
| Method | Non-negative least squares (scipy.optimize.nnls; Lawson–Hanson), solved independently per bulk sample |
| Normalisation | Coefficients scaled to sum to 1 per sample |
| Comparison | Exploratory two-tailed unpaired t tests, CON vs dcKO; no adjustment across cell types. |
| Program localisation | Each gene’s reference CPM values were converted to z-scores across the seven cell types, then averaged across the available genes in each program. |

Estimated representation (CON → dcKO): maturation-stage ameloblasts, 56.9% → 61.3% (*P* = 0.006); secretory-stage ameloblasts, 33.7% → 35.1% (*P* = 0.084); mesenchyme, 9.4% → 3.6% (*P* = 0.002). The remaining fitted coefficients were zero. These estimates describe relative contributions to the expression mixture; zero coefficients do not establish absence of the corresponding cells.

**Interpretation****.** The reference comprised postnatal-day-7 incisor cells, whereas the bulk samples were derived from P0 molar tooth-germ cultures. Reference profiles were normalized to CPM and bulk expression to TPM. Developmental-stage, tooth-type, annotation, and normalization differences limit interpretation of the fitted coefficients as absolute cell fractions. The reference was annotated de novo using the procedure in B3. Deconvolution describes possible compositional contributions but does not identify the cells responsible for a bulk expression change or establish cell autonomy.

B5. Reproducibility

The reproduction workflow uses saved gene-level analysis objects and source workbooks. It retains the original DESeq2 estimates, recomputes descriptive scores and cell-type estimates, fits the limma-voom ROAST models, and generates **Figures 4–7** and the associated numerical tables. The complete workflow also reconstructs the GSE146855 reference, calculates AI-gene enrichment from the version-pinned PanelApp and MGI records, analyzes regional alignments for **Figure S3**, and generates the supplementary workflow diagrams. Normalized inhibitor measurements are read from Figure_3.prism for ANOVA and Dunnett comparisons.

RNA-sequencing reads and complete count/TPM matrices are deposited under GEO accession GSE341665. The de novo transcriptome assembly and saved gene-level analysis objects are held in Zenodo (https://doi.org/10.5281/zenodo.21637766). Analysis code, execution instructions, and machine-readable results are maintained at https://github.com/LabOnoM/GSE341665. Additional analysis inputs and experimental source records are available from the corresponding authors. **Additional files 2–7** contain the computational results reported in the manuscript.

Reference annotation and **Figure S6** use the R/Seurat environment specified in B1. **Figure S3** uses regional alignments and the Ensembl release 116 reference-locus sequence. Input-file checksums and the software environment are specified in the repository.

B6. Amelogenesis-imperfecta gene panel, ortholog mapping and enrichment

**Panel source.** Human amelogenesis-imperfecta genes were obtained from Genomics England PanelApp panel 269, version 4.37 (retrieved 27 July 2026). The analysis included 31 GREEN-rated genes with high diagnostic-grade evidence. The 15 amber- and red-rated genes were excluded: *AIRE, ALPL, AMTN, CLDN16, CLDN19, DSPP, ITGB4, KCNJ1, LAMC2, NECTIN4, SMARCD2, SMOC2, TMEM165, TP63* and* TUFT1*. Human and mouse gene symbols, OMIM phenotype identifiers, inheritance patterns and references are provided in Supplementary Table S4.

**Ortholog mapping.** Human genes were mapped to mouse orthologs using the Mouse Genome Informatics vertebrate homology report (HOM_MouseHumanSequence), accepting only unambiguous one-to-one homology classes. Three genes required explicit manual resolution for nomenclature reasons, each documented in the analysis code:

| Gene | Resolution | Reason |
| --- | --- | --- |
| AMELX | → Amelx | its MGI homology class also contains the Y-linked paralog AMELY, so the 1:1 filter rejects it |
| C4orf26 | → Odaph | PanelApp retains the retired symbol; HGNC renamed the gene ODAPH |
| FAM83H | → Fam83h | renamed SACK1H by HGNC in February 2026, but carried as Fam83h in the Ensembl release 116 annotation used here |

Mouse orthologs were matched to the differential-expression results by m*o*use gene symbol. All 31 panel genes were represented after resolving the nomenclature exceptions listed above.

**Enrichment and threshold sensitivity.** Over-representation of panel genes among down- versus up-regulated genes was assessed by a one-sided Fisher exact test. The test used adjusted *P* < 0.05 with no fold-change threshold, because the core enamel genes that carry the signal have |log2 fold change| ≈ 0.4 and are excluded by the fold-change cutoff used to define DEGs. The result is threshold-dependent and is reported as such:

| Threshold | Panel down | Non-panel down | Panel up | Non-panel up | Odds ratio | P (one-sided) |
| --- | --- | --- | --- | --- | --- | --- |
| adjusted P < 0.05, no fold-change cutoff | 18 | 2060 | 2 | 2107 | 9.21 | 1.7 × 10⁻⁴ |
| adjusted P < 0.05 and \|log2FC\| > 0.25 | 15 | 1778 | 2 | 1656 | 6.99 | 1.9 × 10⁻³ |
| adjusted P < 0.05 and \|log2FC\| > 1 (the DEG rule) | 0 | 354 | 1 | 208 | — | 1.00 |

At the |log2 fold change| > 1 cutoff used to define the DEG list, only one panel gene reaches significance and the test is uninformative. The two significant up-regulated panel genes are *ODAPH* (C4orf26; log2 fold change +1.94) and *ORAI1* (+0.72).

B7. Read-level analysis of the conditional Smad4 allele

The Smad4 locus was examined in the six RNA-sequencing libraries to identify loxP-containing reads, exon coverage and splice junctions. Gene-level *Smad4* expression was modestly **increased** in dcKO samples (log2 fold change +0.32; adjusted *P* = 2.7 × 10⁻³). Read-level evidence was extracted from the alignments using pipeline/20_smad4_loxp_evidence.py and plotted with pipeline/21_smad4_loxp_figure.py (**Supplementary Figure S3**).

**Allele structure.** loxP-containing reads are present in all six libraries, carried in the *soft-clipped* portion of reads that straddle each insertion point (the loxP core is not in the reference genome, so it cannot align). Fifteen such reads were recovered by scanning both orientations of the 34-base loxP core across chr18:73,778,000–73,792,000. One CON1 read (POS 73,784,710, CIGAR 70S79M) carries genomic sequence, a linker, the complete loxP core, a second linker and genomic sequence resuming — fixing that breakpoint unambiguously. Both breakpoints were confirmed against genome.fa.

| Field | Value |
| --- | --- |
| Floxed exon | chr18:73,782,797–73,782,980 (184 bp), entirely coding |
| Exon identity | exon 9 of ENSMUST00000025393 (canonical); exon 8 in the conventional 11-exon numbering |
| Exon accession | ENSMUSE00000143011 in 19 of the 22 transcripts containing this interval; the same interval is labelled ENSMUSE00001829774 (2 transcripts) and ENSMUSE00001829814 (1). The unambiguous identifier is the coordinate interval, not the accession or the ordinal. |
| loxP site A | inserted chr18:73,782,242/73,782,243 (±3 bp; AAC microhomology), intron 9/10 |
| loxP site B | inserted chr18:73,784,709/73,784,710 (exact; one read spans both flanks), intron 8/9 |
| Floxed interval | 2,467 bp, containing exactly one exon |
| Protein consequence | residues 318–379, amino-terminal MH2 domain; 184 mod 3 = 1, so excision is frameshifting and introduces a premature termination codon in the following exon |

Both genotypes carried the floxed Smad4 allele; loxP-containing reads were therefore expected in Cre-negative controls as well as dcKO samples.

**Exon coverage.** Each exon’s contribution to total *Smad4* exonic signal was calculated using per-base depth and exonic read counts. For the floxed exon, the dcKO/CON ratios were 0.924 and 0.889, respectively (exact permutation *P* = 0.10 for each). With three libraries per genotype, the minimum attainable two-sided permutation *P* value was 0.10. Coverage describes the surviving RNA pool.

**Predicted fate of the recombined transcript.** The coding consequence of excision was computed directly from the reference sequence rather than assumed. The floxed exon contributes 184 nucleotides of coding sequence to the canonical transcript (spanning residues 318–379 of the 551-residue protein, of which 319–378 are encoded entirely within it); because 184 is not a multiple of three, its removal shifts the reading frame from the exon 8 | exon 10 junction onward. In the reconstructed mutant coding sequence that junction falls mid-codon (after coding nucleotide 952, within codon 318), and the first in-frame termination codon appears at codon 321 — only 8 nucleotides into exon 10. The predicted product is therefore a 320-residue truncation: residues 1–317 of wild-type SMAD4 followed by three frameshifted residues (Gly-Cys-Thr) and a stop, removing the entire MH2 domain required for Smad complex formation.

This places the premature termination codon 299 nucleotides upstream of the final exon–exon junction, far beyond the 50–55 nucleotide threshold at which a terminating ribosome fails to displace the downstream exon-junction complex. The recombined transcript is therefore a canonical substrate for nonsense-mediated decay and would be expected to be degraded rather than to accumulate. Consistent with this, Ensembl annotates *Smad4* transcripts lacking this exon with the nonsense_mediated_decay biotype.

## Supplementary Figures

**Supplementary Figure S1. Mouse models and allocation of experiments.**

All tooth-germ donors were male. The constitutive K14-cre;Smad4fl/fl model supplied E16.5 tooth germs for Figure 1 and received no 4-OHT. The inducible K14-creER;Smad4fl/fl model supplied P0 tooth germs, with Cre-negative Smad4fl/fl littermates as controls. Both genotypes received 90 µmol/L 4-OHT for 24 h, except in the dose-series experiment. Mechanically isolated dental epithelium was used for molecular analyses (Figures 2–7). All five inhibitor groups, including DMSO vehicle, used Smad4-dcKO tooth germs. For each transplantation comparison, three eight-week-old male recipients each received one experimental and one control germ in the same kidney. Grafts were recovered after 14 days and assessed qualitatively.

**Supplementary Figure S2. PCR images and isolation of dental epithelium.**

**(A)** Smad4 allele-screening PCR. M, DNA ladder; lanes 1–4, individual samples; (+), positive control; WT, wild-type control. The Smad4 primer pair produces an approximately 500-bp floxed product and a 436-bp wild-type product. **(B)** PCR image labeled K14-CreER, with products visible in lanes 1 and 4. **(C)** Representative tooth-germ photograph and the selected dental epithelial region, outlined in yellow. Dental epithelium was separated mechanically under a stereomicroscope after 24 h of culture. These images document allele screening and tissue selection.

**Supplementary Figure S3.** **Read-level evidence at the conditional *****Smad4***** locus.**

Six bulk dental-epithelium libraries from male donors: *n* = 3 dcKO (K14-creER (tamoxifen-inducible)) versus 3 Cre-negative Smad4fl/fl littermate controls. GRCm39 / Ensembl release 116, chr18; *Smad4* is on the **minus strand**, so transcription runs right → left in A–C. **(A)** Locus schematic. The floxed exon (black) is **chr18:73,782,797–73,782,980, 184 bp**; triangles are the loxP insertion points, site A 73,782,242/73,782,243 (±3 bp) and site B 73,784,709/73,784,710 (exact); shading is the 2,467 bp floxed interval. Dashed = the exon 8 → exon 10 junction Cre would create, never observed. **(B)** All 15 reads carrying the 34 bp loxP core, in 3/3 CON and 3/3 dcKO libraries (solid = aligned, open = soft-clipped, black = core; ‡ = a duplicate read, so 14 independent fragments). loxP reads are expected in Cre-negative CON, which carries the identical allele; both cassettes are intronic, so read number is uninformative about excision. Below: the reference sequence at each insertion point, and one cassette-spanning read (CON1, 70S79M at 73,784,710) in literal base calls. **(C)** Per-base coverage across the locus, and for each of the 12 canonical exons the dcKO/CON ratio on both normalisations — per-base depth and exonic reads, with each library plotted. All exons lie within 0.889–1.094; the floxed exon is 0.924 on depth (*P* = 0.10) and 0.889 on reads (*P* = 0.10), and the metrics rank the exons differently.

**Supplementary Figure S4. 4-OHT dose series in the inducible model.**

Contact radiographs of male P0 K14-creER;Smad4fl/fl tooth germs **(A, C, E)** and Cre-negative Smad4fl/fl controls **(B, D, F)** exposed to 50, 90, or 120 µmol/L 4-OHT for 24 h before subrenal-capsule transplantation. Both genotypes received the same concentration within each comparison. Three recipients per comparison each received one experimental and one control graft, and were recovered after 14 days. Representative images are shown without statistical testing. The displayed concentrations were selected from a screen spanning 50–300 µmol/L. The 90 µmol/L condition was used for subsequent experiments.

**Supplementary Figure S5. Analysis workflow. **

**(A)** Bulk RNA-sequencing analysis, from read processing to gene-level differential expression and functional enrichment. Differentially expressed genes were defined by Benjamini–Hochberg adjusted *P* < 0.05 and |log₂ fold change| > 1. **(B)** Gene-program and reference-based analyses using bulk expression data, mouse orthologs of human AI genes, and the GSE146855 single-cell reference. ROAST tests, descriptive program scores, AI-gene enrichment, and cell-type deconvolution provide the data for Figures 4–7 and the Additional-file tables.

**Supplementary Figure S6. The single-cell reference used for deconvolution.**

**(A)** UMAP embedding of 6,248 quality-controlled cells from GSE146855 (postnatal day 7 mouse incisor), coloured by the de-novo cell-type annotation. **(B)** Dot plot of the canonical markers used to assign each cell type; dot size indicates the percentage of cells expressing the gene and colour indicates scaled mean expression. Enamel-matrix transcripts (*Amelx*, *Ambn*, *Enam*, *Mmp20*) were detected across reference cell types, consistent with ambient RNA contamination. Annotation used relative marker expression across clusters as described in **Supplementary Methods B3**.

**Supplementary Figure S7. Enamel phenotype in the inducible model.**

Male P0 K14-creER;Smad4fl/fl and Cre-negative control tooth germs were exposed to 90 µmol/L 4-OHT for 24 h, transplanted beneath the renal capsule, and recovered after 14 days. Hematoxylin-eosin sections are shown at increasing magnification: **C** and **D** enlarge the boxed regions in **A** and **B**; **E** and **F** enlarge the boxed regions in **C** and **D**. Three recipients each received one experimental and one control graft; representative images are shown without statistical testing. am, ameloblasts; e, enamel space; d, dentin; arrows indicate the epithelial/enamel region. Scale bars: A and B, 500 µm; C and D, 100 µm; E, 50 µm.

## Supplementary Tables

**Supplementary Table S1. Primer sequences for qPCR and genotyping.**

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
| PCR internal control | Genotyping PCR | CTAGGCCACAGAATTGAAAGATCT | GTAGGTGGAAATTCTAGCATCATCC |

**Supplementary Table S2. Culture conditions and pathway inhibitors.**

| Reagent | Role | Final concentration | Duration | Supplier and catalog |
| --- | --- | --- | --- | --- |
| 4-Hydroxytamoxifen | Inducible-model treatment, both genotypes | 90 µmol/L | 24 h | Sigma-Aldrich H7904-5MG |
| SB203580 | p38α/β inhibitor | 10 µmol/L | 24 h | Cell Signaling Technology 5633 |
| Dorsomorphin | BMP type I receptor inhibitor | 10 µmol/L | 24 h | Sigma-Aldrich P5499 |
| SB525334 | TGF-β type I receptor/ALK5 inhibitor | 1 µmol/L | 24 h | Tocris 3211 |
| U0126 | MEK1/2 inhibitor | 10 µmol/L | 24 h | Sigma-Aldrich U120 |
| DMSO | Vehicle | Matched across groups | 24 h | Not specified |

The dose-series images show 50, 90, and 120 µmol/L 4-OHT, distinct from the subsequent 90 µmol/L working condition. **Figure 1** uses E16.5 constitutive-model germs (Smad4-cKO) without 4-OHT; the remaining experimental material uses P0 inducible-model germs (Smad4-dcKO). All five inhibitor groups, including the DMSO vehicle group, used Smad4-dcKO germs.

**Supplementary Table S3. Primary antibodies.**

| Antibody | Host | Dilution | Application | Supplier | Catalog |
| --- | --- | --- | --- | --- | --- |
| p38 immunoblot antibody | Rabbit | 1:1,000 | Western blot | Cell Signaling Technology | #9212 |
| Phospho-SMAD1/5/9 | Rabbit | 1:1,000 | Western blot | Cell Signaling Technology | #13820 |
| Phospho-SMAD2/3 | Rabbit | 1:1,000 | Western blot | Thermo Fisher | #PA5-110155 |
| AMELX | Rabbit | 1:500 | Western blot | Thermo Fisher | #PA5-114845 |
| AMBN | Rabbit | 1:500 | Western blot | Bioss antibodies | #bs-12467R |
| MMP20 | Rabbit | 1:1,000 | Western blot | Proteintech | #55467-1-AP |
| KLK4 | Rabbit | 1:1,000 | Western blot | Gentaur | #8941 |
| β-Tubulin | Rabbit | 1:1,000 | Western blot | Thermo Fisher | #600-401-880 |
| GAPDH | Rabbit | 1:1,000 | Western blot | Cell Signaling Technology | #2118 |
| SMAD4 | Rabbit | 1:1,000 | Western blot | Proteintech | #51069-2-AP |

**Supplementary Table S4. Human amelogenesis-imperfecta genes and mouse orthologs used in the enrichment analysis.** The 31 genes were GREEN-rated members of Genomics England PanelApp panel 269, version 4.37 (retrieved 27 July 2026). The table lists OMIM phenotype identifiers, inheritance patterns, primary references and mouse orthologs in Ensembl release 116. For PEX26, the cited reference describes the peroxisome biogenesis disorder; no enamel-specific study was identified. PLXNB2 has no assigned OMIM phenotype identifier.

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

**Table S4 footnote — abbreviations and conventions.**

***Columns******:*** **Gene (human)***,* human gene symbol*;* **Mouse***,* mouse ortholog symbol used for *matching to *the expression *d*at*a;* **OMIM phenotype MIM***,* identifier of the *associated phenotype*; **Phenotype***,* OMIM phenotype designation*;* **Inheritance***,* mode of inheritance*;* **Primary reference (PMID)***,* *sup*porting p*ub*li*ca*tion*;* **Notes***,* *additional* *inform*ati*o*n re*lev*a*nt to interpretation*.

***Inheritance codes****.* **AD**, autosomal dominant; **AR**, autosomal recessive; **XL**, X-linked.

***Disease and phenotype abbreviations****.* **AI**, amelogenesis imperfecta; **EB**, epidermolysis bullosa (**junctional EB**, the junctional subtype); **CRAC**, calcium release-activated calcium channel; **DASS**, dental anomalies and short stature; **SSASKS**, short stature, amelogenesis imperfecta and skeletal dysplasia with scoliosis; **PBD**, peroxisome biogenesis disorder; **LOF**, loss of function.

***Amelogenesis-imperfecta subtype nomenclature******.*** OMIM numbers AI subtypes by a Roman numeral for the defect class followed by a letter (and, where several genes cause the same class, a digit): **type I** = hypoplastic (too little enamel formed, of normal mineralisation); **type II** = hypomaturation (full-thickness but incompletely matured, soft enamel); **type III** = hypocalcified (full-thickness but poorly mineralised). Hence type IA, IB, IE, IF, IG, IH, IJ and IK are hypoplastic forms; types IIA1-IIA6 are hypomaturation forms; and types IIIA and IIIC are hypocalcified forms. Where a parenthetical descriptor follows the subtype it restates this class for readability; rows without one follow the same Roman-numeral rule.

***Gene-symbol changes****.* Two genes are listed under the symbol PanelApp uses, with the current HGNC symbol given in Notes: **C4orf26** is now **ODAPH** (odontogenesis-associated phosphoprotein), and **FAM83H** was renamed **SACK1H** (scaffolding CK1 anchoring protein H) in February 2026. The mouse column retains *Odaph* and *Fam83h* because those are the symbols present in the Ensembl release 116 annotation on which matching was performed.

***Symbols****.* An em dash (**—**) in the OMIM column means no phenotype MIM has been assigned yet; this applies only to PLXNB2, whose AI phenotype was first reported in 2024.

***Note-column meanings****.* **“syndromic”** — enamel involvement occurs as one feature of a broader syndrome rather than as isolated AI. **“also causes junctional EB”** — the same gene causes both an AI-spectrum enamel defect and junctional epidermolysis bullosa. **“recessive LOF”** — the enamel phenotype belongs to the recessive loss-of-function disorder, not to the dominant gain-of-function disorder caused by the same gene. **“PBD primary report; no enamel-specific paper”** — for PEX26 the cited reference establishes the peroxisome biogenesis disorder; enamel involvement is inferred from the disorder spectrum rather than from a dedicated study. **“human evidence recent and limited”** — for SP6 the human AI association rests on one dominant family with functional corroboration plus one de novo case.

***Panel definition.*** All 31 genes are the GREEN-rated (“high evidence”, diagnostic-grade) members of the Genomics England PanelApp panel 269, “Amelogenesis imperfecta”, version 4.37, retrieved 27 July 2026. Amber- and red-rated genes on the same panel were excluded and are named in **Supplementary Methods B6**.

## Supplementary data and repository files

| File | Contents |
| --- | --- |
| Data S1 (GitHub CSV) | Supplementary_Data_S1_inhibitor_normalized_values.csv: 195 normalized inhibitor measurements, also present in Figure_3.prism. Sample indices are within-table indices, not animal IDs. |
| Data S2 (GitHub CSV) | Supplementary_Data_S2_inhibitor_ANOVA_Dunnett.csv: 52 treatment-versus-vehicle comparisons with ANOVA statistics, group means/SDs and independently recalculated two-sided Dunnett P values. Native Prism analyses remain in Figure_3.prism. |
| Additional_file_2.xlsx (Data S3) | Gene-level statistics, per-library TPM, ROAST tests and sensitivity analyses, gene-set membership, sample covariates and gene values reported in the text and figures. |
| Additional_file_3.csv (Data S4) | Log2 fold changes and adjusted P values for the 31 GREEN-rated genes of the Genomics England PanelApp “Amelogenesis imperfecta” panel (269 v4.37). |
| Additional_file_4.csv (Data S5) | All 563 differentially expressed genes (adjusted P < 0.05, absolute log2 fold change > 1) with DESeq2 statistics and human homolog annotation. |
| Additional_file_5.csv (Data S6) | Per-sample program scores for the four gene programs. |
| Additional_file_6.csv (Data S7) | Per-sample cell-type proportions estimated by non-negative least squares. |
| Additional_file_7.csv (Data S8) | Cell-type specificity of each gene program in the single-cell reference. |

Complete gene-level count and transcripts-per-million matrices are associated with GEO GSE341665. The GEO TPM export is rounded to four decimal places; the archived source workbooks and plotting inputs retain full precision. The data-file map and input checksums are supplied in the GitHub repository.

**Additional file 2** (**Supplementary Data S3**). Gene-level expression and statistical results for Figures 4–7 and the main text. The workbook contains all 15,803 retained genes with DESeq2 log2 fold changes, raw and Benjamini–Hochberg-adjusted P values, TPM values for the six libraries, gene-set membership and figure assignments. Separate sheets provide genes reported in the text and figures, ROAST tests and sensitivity analyses, sample covariates, and deconvolution comparisons. Auxiliary alignment counts for Sox2 and Lgr5 are reported separately from the Salmon/tximport model inputs. Zero adjusted P values indicate numerical underflow; missing test results remain missing.
