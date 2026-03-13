# ==========================================
# Seurat scRNA-seq Analysis for T cells (Lymph Node)
# Input: E-HCAD-8.aggregated_filtered_counts.mtx (+ _cols, _rows)
# ==========================================

# -------------------------------
# Load required packages
# -------------------------------
if (!requireNamespace("Seurat", quietly = TRUE)) install.packages("Seurat")
if (!requireNamespace("Matrix", quietly = TRUE)) install.packages("Matrix")
if (!requireNamespace("readxl", quietly = TRUE)) install.packages("readxl")
if (!requireNamespace("DoubletFinder", quietly = TRUE)) install.packages("DoubletFinder")
if (!requireNamespace("biomaRt", quietly = TRUE)) install.packages("biomaRt")
if (!requireNamespace("ggsci", quietly = TRUE)) install.packages("ggsci")
if (!requireNamespace("RColorBrewer", quietly = TRUE)) install.packages("RColorBrewer")
if (!requireNamespace("patchwork", quietly = TRUE)) install.packages("patchwork")
if (!requireNamespace("plyr", quietly = TRUE)) install.packages("plyr")
if (!requireNamespace("dplyr", quietly = TRUE)) install.packages("dplyr")
if (!requireNamespace("tidyr", quietly = TRUE)) install.packages("tidyr")

library(Seurat)
library(Matrix)
library(DoubletFinder)
library(biomaRt)
library(ggplot2)
library(readxl)
library(ggsci)
library(RColorBrewer)
library(patchwork)
library(plyr)
library(dplyr)
library(tidyr)

# -------------------------------
# Load raw data
# -------------------------------
counts <- readMM("E-HCAD-8.aggregated_filtered_counts.mtx")
cell_info <- read.table("E-HCAD-8.aggregated_filtered_counts.mtx_cols", header = FALSE, stringsAsFactors = FALSE)
gene_info <- read.table("E-HCAD-8.aggregated_filtered_counts.mtx_rows", header = FALSE, stringsAsFactors = FALSE)

rownames(counts) <- gene_info$V1
colnames(counts) <- cell_info$V1

# -------------------------------
# Subset T cells from lymph nodes
# -------------------------------
meta <- read_excel("meta.xlsx")
Tcells <- meta[meta$cell_type == "T cell", ]
lymph <- Tcells[Tcells$part == "lymph node", ]
counts_lymph <- counts[, lymph$cellname]

# -------------------------------
# Create Seurat object and filter lowly expressed genes
# -------------------------------
seurat_obj <- CreateSeuratObject(counts = counts_lymph, min.features = 100)
counts_data <- GetAssayData(seurat_obj, slot = "counts")
nonzero_genes <- counts_data > 0
keep_genes <- Matrix::rowSums(nonzero_genes) >= 10
counts_filtered <- counts_data[keep_genes, ]
seurat_obj <- CreateSeuratObject(counts_filtered)

# -------------------------------
# Normalize, scale, PCA, UMAP
# -------------------------------
seurat_obj <- NormalizeData(seurat_obj, normalization.method = "LogNormalize")
seurat_obj <- FindVariableFeatures(seurat_obj, selection.method = "vst", nfeatures = 2000)
all_genes <- rownames(seurat_obj)
seurat_obj <- ScaleData(seurat_obj, features = all_genes)
seurat_obj <- RunPCA(seurat_obj, features = VariableFeatures(seurat_obj))
seurat_obj <- RunUMAP(seurat_obj, dims = 1:15)
seurat_obj <- FindNeighbors(seurat_obj, dims = 1:15)
seurat_obj <- FindClusters(seurat_obj, resolution = 0.5)

# -------------------------------
# Convert Ensembl IDs to gene symbols
# -------------------------------
mart <- useMart("ensembl", dataset = "hsapiens_gene_ensembl")
gene_conversion <- getBM(
  attributes = c("ensembl_gene_id", "hgnc_symbol"),
  filters = "ensembl_gene_id",
  values = rownames(seurat_obj),
  mart = mart
)
gene_conversion <- gene_conversion[gene_conversion$hgnc_symbol != "", ]
rownames(seurat_obj@assays$RNA@counts) <- gene_conversion$hgnc_symbol[
  match(rownames(seurat_obj@assays$RNA@counts), gene_conversion$ensembl_gene_id)
]
rownames(seurat_obj@assays$RNA@data) <- rownames(seurat_obj@assays$RNA@counts)
rownames(seurat_obj@assays$RNA@scale.data) <- rownames(seurat_obj@assays$RNA@counts)

# Ensure unique gene names
counts_final <- seurat_obj[["RNA"]]@counts
genes <- rownames(counts_final)
genes[is.na(genes) | genes == ""] <- paste0("Gene_", seq_len(sum(is.na(genes) | genes == "")))
rownames(counts_final) <- make.unique(genes)
tcell_final <- CreateSeuratObject(counts = counts_final)

# -------------------------------
# Normalize, variable features, scaling, PCA
# -------------------------------
tcell_final <- NormalizeData(tcell_final)
tcell_final <- FindVariableFeatures(tcell_final)
all_genes_final <- rownames(tcell_final)
tcell_final <- ScaleData(tcell_final, features = all_genes_final)
tcell_final <- RunPCA(tcell_final, npcs = 10)

# -------------------------------
# Doublet detection
# -------------------------------
sweep.res <- paramSweep_v3(tcell_final, PCs = 1:10, sct = FALSE)
sweep.stats <- summarizeSweep(sweep.res)
bcmvn <- find.pK(sweep.stats)
print(head(bcmvn[order(-bcmvn$BCmetric), ], 10))

optimal.pK <- 0.05
pN <- 0.25
doublet_rate <- 0.075
nExp <- round(ncol(tcell_final) * doublet_rate)

tcell_final <- doubletFinder_v3(
  tcell_final,
  PCs = 1:10,
  pN = pN,
  pK = optimal.pK,
  nExp = nExp,
  reuse.pANN = FALSE,
  sct = FALSE
)

df_col <- grep("DF.classifications", colnames(tcell_final@meta.data), value = TRUE)
prop_doublet <- sum(tcell_final@meta.data[[df_col]] == "Doublet") / ncol(tcell_final)

if (!"umap" %in% names(tcell_final@reductions)) {
  tcell_final <- RunUMAP(tcell_final, dims = 1:10)
}
DimPlot(tcell_final, group.by = df_col, pt.size = 0.8) +
  ggtitle(paste0("Doublet classification (nExp=", nExp, ")"))

# Keep only singlets
tcell_clean <- subset(tcell_final, subset = DF.classifications_0.25_0.05_1689 == "Singlet")

# Final processing: normalization, PCA, UMAP, clustering
tcell_clean <- NormalizeData(tcell_clean)
tcell_clean <- FindVariableFeatures(tcell_clean, selection.method = "vst", nfeatures = 2000)
all_genes_clean <- rownames(tcell_clean)
tcell_clean <- ScaleData(tcell_clean, features = all_genes_clean)
tcell_clean <- RunPCA(tcell_clean, features = VariableFeatures(tcell_clean))
tcell_clean <- RunUMAP(tcell_clean, dims = 1:20)
tcell_clean <- FindNeighbors(tcell_clean, dims = 1:20)
tcell_clean <- FindClusters(tcell_clean, resolution = 0.5)

# -------------------------------
# Rename clusters
# -------------------------------
cluster_markers <- FindAllMarkers(
  tcell_clean,
  only.pos = TRUE,
  min.pct = 0.25,
  logfc.threshold = 0.25
)

old_clusters <- 0:13
new_labels <- paste0("Cluster ", 1:14)
tcell_clean$seurat_clusters <- plyr::mapvalues(
  x = tcell_clean$seurat_clusters,
  from = old_clusters,
  to = new_labels
)
tcell_clean$cluster_number <- as.numeric(gsub("Cluster ", "", tcell_clean$seurat_clusters))

# -------------------------------
# FigureS2A: UMAP of all clusters
# -------------------------------
n_clusters <- length(unique(tcell_clean$cluster_number))
cols <- c(brewer.pal(8, "Set2"), brewer.pal(8, "Pastel1")[1:6])

figureS2A <- DimPlot(
  object = tcell_clean,
  group.by = "cluster_number",
  label = TRUE,
  label.size = 2.5,
  reduction = "umap",
  pt.size = 0.2,
  cols = cols
) +
  ggtitle("FigureS2A") +
  theme(
    plot.title = element_blank(),
    axis.title = element_text(family = "Helvetica", size = 8, face = "plain", color = "black"),
    axis.text = element_text(family = "Helvetica", size = 8, face = "plain", color = "black"),
    legend.position = "none",
    panel.background = element_rect(fill = NA, color = NA),
    panel.grid.major = element_blank(),
    panel.grid.minor = element_blank(),
    plot.background = element_rect(fill = NA, color = NA)
  ) +
  guides(color = guide_legend(override.aes = list(size = 1)))

ggsave("FigureS2A.png", plot = figureS2A, width = 2.5, height = 2, units = "in", dpi = 900, bg = "transparent")

# -------------------------------
# FigureS2B: FeaturePlots of CD4 and CD8A
# -------------------------------
figureS2B_CD4 <- FeaturePlot(
  object = tcell_clean,
  features = "CD4",
  reduction = "umap",
  pt.size = 0.0001, order = TRUE
) +
  theme_minimal() +
  theme(
    plot.title = element_text(family = "Helvetica", size = 8, face = "italic", color = "black", hjust = 0.5, vjust = 0.1),
    axis.text = element_text(family = "Helvetica", size = 8, color = "black"),
    legend.background = element_rect(fill = alpha('white',0.6), color = NA),
    legend.key.size = unit(0.3, "cm"),
    legend.position = c(1, 0.07),
    legend.justification = c("right", "bottom"),
    panel.background = element_rect(fill = NA, color = NA),
    plot.background = element_rect(fill = NA, color = NA)
  ) +
  guides(color = guide_colorbar(barwidth = 0.3, barheight = 2.5))
ggsave("FigureS2B_CD4.png", plot = figureS2B_CD4, width = 2.5, height = 2.2, units = "in", dpi = 900, bg = "transparent")

figureS2B_CD8 <- FeaturePlot(
  object = tcell_clean,
  features = "CD8A",
  reduction = "umap",
  pt.size = 0.0001, order = TRUE
) +
  theme_minimal() +
  theme(
    plot.title = element_text(family = "Helvetica", size = 8, face = "italic", color = "black", hjust = 0.5, vjust = 0.1),
    axis.text = element_text(family = "Helvetica", size = 8, color = "black"),
    legend.background = element_rect(fill = alpha('white',0.6), color = NA),
    legend.key.size = unit(0.3, "cm"),
    legend.position = c(1, 0.07),
    legend.justification = c("right", "bottom"),
    panel.background = element_rect(fill = NA, color = NA),
    plot.background = element_rect(fill = NA, color = NA)
  ) +
  guides(color = guide_colorbar(barwidth = 0.3, barheight = 2.5))
ggsave("FigureS2B_CD8A.png", plot = figureS2B_CD8, width = 2.5, height = 2.2, units = "in", dpi = 900, bg = "transparent")

# -------------------------------
# FigureS2C: DotPlot
# -------------------------------
figureS2C <- DotPlot(
  tcell_clean,
  features = c("CD4","CD8A","IGHG1","IGHG2","IGHG3","IGHG4",
               "IGLC2","IGHM","GZMA","GZMB","GZMK","GZMM",
               "GNLY","PRF1","TBX21","GATA3","RORC","PDCD1",
               "TIGIT","CXCL13","FOXP3","ZNF683","ITGAE",
               "IFNG","TNF","TGFB1","CXCR5"),
  group.by = "seurat_clusters",
  cols = c("blue","firebrick1"),
  scale = TRUE,
  dot.scale = 3
) + RotatedAxis() +
  theme(
    text = element_text(family = "Helvetica", size = 8, color = "black"),
    axis.text.x = element_text(angle = 45, hjust = 1, size = 8),
    axis.text.y = element_text(size = 8),
    axis.title = element_blank(),
    plot.title = element_text(size = 8),
    legend.title = element_text(family = "Helvetica", size = 8),
    legend.text = element_text(family = "Helvetica", size = 8)
  )
figureS2C <- figureS2C +
  theme(
    legend.key.size = unit(0.32, "cm"),
    legend.spacing.y = unit(0.1, "cm")
  )
ggsave("FigureS2C.png", plot = figureS2C, width = 6.5, height = 2.5, units = "in", dpi = 900, bg = "transparent")

# -------------------------------
# FigureS2D: Violin plots of IGH genes
# -------------------------------
genes <- c("IGHG1","IGHG2","IGHG3","IGHG4","IGLC2","IGHM")
cols <- c(brewer.pal(8, "Set2"), brewer.pal(8, "Pastel1")[1:6])

data_long <- tcell_clean@assays$RNA@data[genes, ] %>%
  as.data.frame() %>%
  tibble::rownames_to_column(var = "Gene") %>%
  pivot_longer(-Gene, names_to = "cell", values_to = "Expression")

meta <- tcell_clean@meta.data %>%
  tibble::rownames_to_column(var = "cell") %>%
  select(cell, cluster_number)

data_long <- data_long %>% left_join(meta, by = "cell")

vln_list <- lapply(genes, function(g){
  df <- subset(data_long, Gene == g)
  ggplot(df, aes(x = factor(cluster_number), y = Expression, fill = factor(cluster_number))) +
    geom_violin(scale = "width", adjust = 1.5, alpha = 1, color = "gray50", trim = TRUE) +
    geom_jitter(width = 0.15, size = 0.5, alpha = 1) +
    scale_x_discrete(expand = c(0,0)) +
    scale_fill_manual(values = cols) +
    coord_cartesian(ylim = c(0,6), clip = "off") +
    geom_hline(yintercept = 0, color = "black", size = 0.5) +
    geom_vline(xintercept = 0.5, color = "black", size = 0.5) +
    ggtitle(g) +
    xlab(NULL) + ylab(NULL) +
    theme(
      axis.text.x = element_text(family = "Helvetica", size = 8, color = "black"),
      axis.text.y = element_text(family = "Helvetica", size = 8, color = "black"),
      axis.title = element_text(family = "Helvetica", size = 8, face = "italic", color = "black"),
      strip.text = element_text(family = "Helvetica", size = 8, color = "black"),
      plot.title = element_text(family = "Helvetica", size = 8, face = "italic", color = "black", hjust = 0.5),
      plot.background = element_rect(fill = NA, color = NA),
      panel.background = element_rect(fill = NA, color = NA),
      panel.grid.major = element_blank(),
      panel.grid.minor = element_blank(),
      legend.position = "none"
    )
})

figureS2D <- wrap_plots(vln_list, ncol = 3) +
  plot_annotation(theme = theme(plot.background = element_rect(fill = NA, color = NA)))


ggsave("FigureS2D.png", plot = figureS2D, width = 6.5, height = 3.5, units = "in", dpi = 900, bg = "transparent")
