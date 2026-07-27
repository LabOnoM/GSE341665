#!/usr/bin/env Rscript
# 05_scref_annotate.R  --  Smad4 dcKO JDR re-analysis pipeline (step 5/7)
#
# Cluster and annotate the GSE146855 reference (de novo, since GEO ships raw
# counts without author labels), then export a cell-type pseudobulk CPM
# signature for deconvolution.
#
# Annotation method: assign each Seurat cluster to the cell type whose canonical
# marker set has the highest CROSS-CLUSTER SPECIFICITY (z-scored mean expression).
# This is essential here because enamel-matrix transcripts (Amelx/Ambn/Enam) are
# so hyper-abundant they contaminate every droplet as ambient RNA; absolute
# expression would (wrongly) label every cluster an ameloblast.
#
# Inputs : output/scref/GSE146855/mtx (from 04_download_scref.sh)
# Outputs: output/GSE146855_signature_cpm.csv        (genes x 7 cell types, CPM)
#          output/GSE146855_celltype_specificity.csv (cluster x cell type z)
#          output/GSE146855_seurat.rds
# Runtime: ~2 min.  seed fixed for reproducibility.
suppressMessages(library(Seurat)); set.seed(1)
HERE <- dirname(sub("--file=", "", grep("--file=", commandArgs(FALSE), value = TRUE)))
if (length(HERE) == 0 || HERE == "") HERE <- "scripts/JDR_reanalysis"
OUT <- file.path(HERE, "output")

so <- CreateSeuratObject(Read10X(file.path(OUT, "scref/GSE146855/mtx"), gene.column = 2),
                         min.cells = 3, min.features = 200)
so[["percent.mt"]] <- PercentageFeatureSet(so, pattern = "^mt-")
so <- subset(so, subset = nFeature_RNA > 200 & nFeature_RNA < 7000 & percent.mt < 25)
so <- NormalizeData(so, verbose = FALSE)
so <- FindVariableFeatures(so, nfeatures = 2000, verbose = FALSE)
so <- ScaleData(so, verbose = FALSE)
so <- RunPCA(so, npcs = 30, verbose = FALSE)
so <- FindNeighbors(so, dims = 1:20, verbose = FALSE)
so <- FindClusters(so, resolution = 0.6, verbose = FALSE)
cat("cells:", ncol(so), " clusters:", length(levels(Idents(so))), "\n")

markers <- list(
  Ameloblast_secretory  = c("Amelx", "Ambn", "Enam", "Mmp20"),
  Ameloblast_maturation = c("Klk4", "Odam", "Amtn", "Slc24a4", "Gm17660"),
  Epithelial_IEE_SI     = c("Shh", "Sfrp5", "Krt14", "Krt5", "Epcam", "Thbd", "Gjb3", "Cldn10"),
  Mesenchyme            = c("Col1a1", "Dcn", "Lum", "Prrx1", "Pdgfra", "Msx1", "Sp7", "Dspp"),
  Immune                = c("Ptprc", "Lyz2", "C1qa", "Cd68", "Fcer1g"),
  Endothelial           = c("Pecam1", "Cdh5", "Egfl7", "Emcn"),
  Erythroid             = c("Hba-a1", "Hbb-bs", "Hba-a2"))
markers <- lapply(markers, function(g) intersect(g, rownames(so)))

avg <- AverageExpression(so, assays = "RNA", slot = "data", verbose = FALSE)$RNA
zt  <- t(scale(t(as.matrix(avg)))); zt[is.na(zt)] <- 0          # z across clusters
score <- sapply(markers, function(gs) if (length(gs)) colMeans(zt[gs, , drop = FALSE]) else rep(-Inf, ncol(zt)))
score <- as.data.frame(round(score, 2))
assign <- colnames(score)[max.col(as.matrix(score), ties.method = "first")]
names(assign) <- sub("^g", "", rownames(score))
so$celltype <- unname(assign[as.character(Idents(so))])
cat("\ncell-type counts:\n"); print(table(so$celltype))

rawc  <- GetAssayData(so, assay = "RNA", slot = "counts")
types <- sort(unique(so$celltype))
sig   <- sapply(types, function(ct) Matrix::rowSums(rawc[, so$celltype == ct, drop = FALSE]))
cpm   <- sweep(sig, 2, colSums(sig), "/") * 1e6
write.csv(round(cpm, 3), file.path(OUT, "GSE146855_signature_cpm.csv"))
write.csv(score,          file.path(OUT, "GSE146855_celltype_specificity.csv"))
saveRDS(so,               file.path(OUT, "GSE146855_seurat.rds"))
cat("signature genes:", nrow(cpm), " cell types:", ncol(cpm), "\nDONE\n")
