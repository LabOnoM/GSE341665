#!/usr/bin/env Rscript
# 11_scref_umap_dotplot.R -- Supplementary figures for the single-cell reference
#
# Loads the annotated GSE146855 Seurat object (from 05_scref_annotate.R), computes
# a UMAP embedding, and produces two supplementary panels documenting the reference
# used for deconvolution:
#   S6A  UMAP coloured by the de-novo cell-type annotation
#   S6B  Dot plot of the canonical markers used to assign each cell type
#
# Also writes the exact package versions used, for the supplementary methods.
#
# Output: figures/FigureS6_scref_umap_dotplot.{png,pdf}
#         output/scref_umap_versions.txt, output/scref_celltype_counts.csv
suppressMessages({library(Seurat); library(ggplot2); library(patchwork)})
set.seed(1)
HERE <- dirname(sub("--file=", "", grep("--file=", commandArgs(FALSE), value = TRUE)))
if (length(HERE) == 0 || HERE == "") HERE <- "scripts/JDR_reanalysis"
OUT <- file.path(HERE, "output"); FIG <- file.path(HERE, "figures")
dir.create(FIG, showWarnings = FALSE)

so <- readRDS(file.path(OUT, "GSE146855_seurat.rds"))
cat("cells:", ncol(so), " cell types:", length(unique(so$celltype)), "\n")

# UMAP on the same 20 PCs used for clustering in 05_scref_annotate.R
so <- RunUMAP(so, dims = 1:20, umap.method = "uwot", metric = "cosine",
              n.neighbors = 30, min.dist = 0.3, seed.use = 42, verbose = FALSE)

ct_order <- c("Ameloblast_secretory", "Ameloblast_maturation", "Epithelial_IEE_SI",
              "Mesenchyme", "Immune", "Endothelial", "Erythroid")
ct_order <- ct_order[ct_order %in% unique(so$celltype)]
so$celltype <- factor(so$celltype, levels = ct_order)
Idents(so) <- "celltype"

pal <- c(Ameloblast_secretory = "#B2182B", Ameloblast_maturation = "#E8783A",
         Epithelial_IEE_SI = "#F4A582", Mesenchyme = "#2166AC",
         Immune = "#4DAF4A", Endothelial = "#984EA3", Erythroid = "#999999")

# ---- S6A: UMAP ----
p1 <- DimPlot(so, reduction = "umap", cols = pal, pt.size = .3, label = TRUE,
              label.size = 2.6, repel = TRUE) +
  ggtitle("A  Single-cell reference (GSE146855, P7 incisor): de-novo cell-type annotation") +
  theme(legend.position = "right", legend.text = element_text(size = 6),
        plot.title = element_text(size = 9, face = "bold"),
        axis.title = element_text(size = 7), axis.text = element_text(size = 6))

# ---- S6B: marker dot plot (the markers used for annotation) ----
markers <- c("Amelx", "Ambn", "Enam", "Mmp20",          # secretory ameloblast
             "Klk4", "Odam", "Amtn", "Slc24a4",          # maturation ameloblast
             "Shh", "Sfrp5", "Krt14", "Epcam",           # other dental epithelium
             "Col1a1", "Dcn", "Lum", "Pdgfra",           # mesenchyme
             "Ptprc", "Lyz2", "C1qa",                    # immune
             "Pecam1", "Cdh5", "Egfl7",                  # endothelial
             "Hba-a1", "Hbb-bs")                         # erythroid
markers <- markers[markers %in% rownames(so)]
p2 <- DotPlot(so, features = markers, dot.scale = 3.4) +
  scale_colour_gradient2(low = "#2166AC", mid = "white", high = "#B2182B") +
  ggtitle("B  Canonical marker expression used for cell-type assignment") +
  theme(axis.text.x = element_text(angle = 90, hjust = 1, vjust = .5, size = 6, face = "italic"),
        axis.text.y = element_text(size = 6), axis.title = element_blank(),
        legend.text = element_text(size = 6), legend.title = element_text(size = 6),
        plot.title = element_text(size = 9, face = "bold"))

p <- p1 / p2 + plot_layout(heights = c(1.15, 1))
ggsave(file.path(FIG, "FigureS6_scref_umap_dotplot.png"), p, width = 9.5, height = 8.2, dpi = 300)
ggsave(file.path(FIG, "FigureS6_scref_umap_dotplot.pdf"), p, width = 9.5, height = 8.2)

# ---- provenance for the supplementary methods ----
tab <- as.data.frame(table(so$celltype)); colnames(tab) <- c("cell_type", "n_cells")
tab$percent <- round(100 * tab$n_cells / sum(tab$n_cells), 1)
write.csv(tab, file.path(OUT, "scref_celltype_counts.csv"), row.names = FALSE)
print(tab)

pk <- c("Seurat", "SeuratObject", "Matrix", "ggplot2", "patchwork", "uwot")
v <- sapply(pk, function(x) tryCatch(as.character(packageVersion(x)), error = function(e) NA))
writeLines(c(paste("R", getRversion()),
             paste(names(v), v, sep = " ")),
           file.path(OUT, "scref_umap_versions.txt"))
cat("\nversions:\n"); print(v)
cat("DONE\n")
saveRDS(so, file.path(OUT,"GSE146855_seurat_umap.rds"))
