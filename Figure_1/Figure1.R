###############################################
# T Cell RNA-seq Analysis – Figures 1A–1E
# Author: Tuan Nguyen
###############################################

# ===============================
# 0. Load necessary packages
# ===============================
library(Seurat)
library(ggplot2)
library(DoubletFinder)
library(plyr)
library(patchwork)
library(RColorBrewer)
library(readxl)
library(Matrix)

# ===============================
# 1. Load raw UMI matrix and metadata
# ===============================
# Load UMI count matrix
umi_matrix <- read.delim("celseq_matrix_ru10_molecules.tsv", header = TRUE, sep = "\t")
rownames(umi_matrix) <- umi_matrix$gene
umi_matrix <- umi_matrix[, -1]

# Load metadata
metadata <- read_excel("metadata.xlsx")
tcell_metadata <- metadata[metadata$type == "T cell", ]

# ===============================
# 2. Quality control filtering
# ===============================
filter_mt <- tcell_metadata$percent_mt_molecules <= 0.25
filter_ngenes <- tcell_metadata$genes_detected >= 1000
filter_samples <- !tcell_metadata$sample %in% c("300-0153","300-0451","301-0121","301-0122","301-0134")
qc_filter <- filter_mt & filter_ngenes & filter_samples

tcell_metadata_filtered <- tcell_metadata[qc_filter, ]

# Subset UMI matrix for filtered cells
umi_tcell <- umi_matrix[, tcell_metadata_filtered$cell_name]
umi_tcell[is.na(umi_tcell)] <- 0

# ===============================
# 3. Create Seurat object and filter lowly expressed genes
# ===============================
seurat_obj <- CreateSeuratObject(counts = umi_tcell, min.features = 100)
counts <- GetAssayData(seurat_obj, slot = "counts")
keep_genes <- Matrix::rowSums(counts > 0) >= 10
counts_filtered <- counts[keep_genes, ]
seurat_filtered <- CreateSeuratObject(counts_filtered)

# ===============================
# 4. Normalization, variable feature selection, scaling, PCA
# ===============================
seurat_filtered <- NormalizeData(seurat_filtered, normalization.method = "LogNormalize")
seurat_filtered <- FindVariableFeatures(seurat_filtered, selection.method = "vst", nfeatures = 2000)
seurat_filtered <- ScaleData(seurat_filtered, features = rownames(seurat_filtered))
seurat_filtered <- RunPCA(seurat_filtered, features = VariableFeatures(seurat_filtered))

# Check elbow plot for PCs
ElbowPlot(seurat_filtered)
PCs_to_use <- 1:20

# ===============================
# 5. Doublet detection
# ===============================
sweep_res <- paramSweep_v3(seurat_filtered, PCs = PCs_to_use, sct = FALSE)
sweep_stats <- summarizeSweep(sweep_res)
bcmvn <- find.pK(sweep_stats)
print(head(bcmvn[order(-bcmvn$BCmetric), ], 10)) # Display top pK values

optimal_pK <- 0.01
pN <- 0.25
doublet_rate <- 0.05
nExp <- round(ncol(seurat_filtered) * doublet_rate)

seurat_filtered <- doubletFinder_v3(seurat_filtered,
                                    PCs = PCs_to_use,
                                    pN = pN,
                                    pK = optimal_pK,
                                    nExp = nExp,
                                    reuse.pANN = FALSE,
                                    sct = FALSE)

df_col <- grep("DF.classifications", colnames(seurat_filtered@meta.data), value = TRUE)
prop_doublet <- sum(seurat_filtered@meta.data[[df_col]] == "Doublet") / ncol(seurat_filtered)
message("Observed doublet fraction: ", round(prop_doublet, 3))

# Run UMAP if not already done
if (!"umap" %in% names(seurat_filtered@reductions)) {
  seurat_filtered <- RunUMAP(seurat_filtered, dims = PCs_to_use)
}

# ===============================
# 6. Keep only singlets (safe method)
# ===============================
# Extract cells classified as Singlet
singlet_cells <- rownames(seurat_filtered@meta.data)[seurat_filtered@meta.data[[df_class_col]] == "Singlet"]

# Subset Seurat object by these cells
seurat_singlets <- subset(seurat_filtered, cells = singlet_cells)

# ===============================
# 7. Re-normalize, scale, PCA, UMAP, clustering
# ===============================
seurat_singlets <- NormalizeData(seurat_singlets)
seurat_singlets <- FindVariableFeatures(seurat_singlets, selection.method = "vst", nfeatures = 2000)
seurat_singlets <- ScaleData(seurat_singlets)
seurat_singlets <- RunPCA(seurat_singlets, features = VariableFeatures(seurat_singlets))
seurat_singlets <- RunUMAP(seurat_singlets, dims = 1:20)
seurat_singlets <- FindNeighbors(seurat_singlets, dims = 1:20)
seurat_singlets <- FindClusters(seurat_singlets, resolution = 0.5)

# ===============================
# 8. Re-label clusters
# ===============================
old_clusters <- 0:(length(unique(seurat_singlets$seurat_clusters)) - 1)
new_labels <- paste0("Cluster ", seq_along(old_clusters))
seurat_singlets$seurat_clusters <- plyr::mapvalues(
  x = seurat_singlets$seurat_clusters,
  from = old_clusters,
  to = new_labels
)

# ===============================
# 9. Add disease metadata
# ===============================
rownames(tcell_metadata_filtered) <- tcell_metadata_filtered$cell_name
cells_in_seurat <- colnames(seurat_singlets)
metadata_subset <- tcell_metadata_filtered[cells_in_seurat, , drop = FALSE]
seurat_singlets@meta.data$disease <- metadata_subset$disease

# ===============================
# 10. Figures
# ===============================

# Fig1A: UMAP colored by cluster
cols_clusters <- c(brewer.pal(8, "Set2"), brewer.pal(8, "Pastel1")[1:3])
fig1a <- DimPlot(seurat_singlets, group.by = "seurat_clusters",
                 label = TRUE, label.size = 2.5,
                 pt.size = 0.2, cols = cols_clusters) +
  theme_classic(base_size = 8)
ggsave("Fig1a.png", fig1a, width = 3.5, height = 2, dpi = 900, bg = "transparent")

# Fig1B: UMAP colored by disease
fig1b <- DimPlot(seurat_singlets, group.by = "disease",
                 reduction = "umap", pt.size = 0.2,
                 cols = c("OA" = "black", "RA" = "lightgray")) +
  theme_classic(base_size = 8)
ggsave("Fig1b.png", fig1b, width = 3.5, height = 2, dpi = 900)

# Fig1C: FeaturePlot CD4/CD8A
fig1c_CD4 <- FeaturePlot(seurat_singlets, features = "CD4", reduction = "umap", pt.size = 0.2)
fig1c_CD8 <- FeaturePlot(seurat_singlets, features = "CD8A", reduction = "umap", pt.size = 0.2)
ggsave("Fig1c_CD4.png", fig1c_CD4, width = 3.5, height = 2, dpi = 900, bg = "transparent")
ggsave("Fig1c_CD8.png", fig1c_CD8, width = 3.5, height = 2, dpi = 900, bg = "transparent")

# Fig1D: DotPlot of marker genes
marker_genes <- c("CD4","CD8A","IGHG1","IGHG2","IGHG3","IGHG4",
                  "IGLC2","IGHM","GZMA","GZMB","GZMK","GZMM",
                  "GNLY","PRF1","TBX21","GATA3","RORC","PDCD1",
                  "TIGIT","CXCL13","FOXP3","ZNF683","ITGAE",
                  "IFNG","TNF","TGFB1","CXCR5")
fig1d <- DotPlot(seurat_singlets, features = marker_genes,
                 group.by = "seurat_clusters", cols = c("blue","firebrick1"), dot.scale = 3) +
  RotatedAxis()
ggsave("Fig1d.png", fig1d, width = 6.5, height = 2.5, dpi = 900, bg = "transparent")

# Fig1E: Violin plot of immunoglobulin genes
seurat_singlets$cluster_number <- as.numeric(seurat_singlets$seurat_clusters)
genes_igh <- c("IGHG1","IGHG2","IGHG3","IGHG4","IGLC2","IGHM")
vln_list <- VlnPlot(seurat_singlets, features = genes_igh, group.by = "cluster_number", pt.size = 0.2, combine = FALSE)
fig1e <- wrap_plots(vln_list, ncol = 3)
ggsave("Fig1e.png", fig1e, width = 6.5, height = 3.5, dpi = 900, bg = "transparent")

# ===============================
# 11. Save session info for reproducibility
# ===============================
writeLines(capture.output(sessionInfo()), "sessionInfo.txt")
