# ============================================================
# B cell filtering and comparison with Cluster 7 T cells
# ============================================================

library(Seurat)
library(dplyr)
library(reshape2)
library(ggplot2)

# ============================================================
# 1. B CELL FILTERING
# ============================================================

# Subset metadata for B cells
bcell_metadata <- metadata[metadata$type == "B cell", ]

# QC filtering
filter_mt_b <- bcell_metadata$percent_mt_molecules <= 0.25
filter_ngenes_b <- bcell_metadata$genes_detected >= 1000
filter_samples_b <- !bcell_metadata$sample %in% c(
  "300-0153","300-0451","301-0121","301-0122","301-0134"
)

# Combine filters
filter_final_b <- filter_mt_b & filter_ngenes_b & filter_samples_b

# Subset metadata for passing B cells
bcell_metadata_filtered <- bcell_metadata[filter_final_b, ]

# Subset UMI matrix for filtered B cells
umi_bcell <- umi_matrix[, bcell_metadata_filtered$cell_name, drop = FALSE]


# ============================================================
# 2. EXTRACT CLUSTER 7 TCELLS
# ============================================================

cluster7_cells <- colnames(seurat_singlets)[
  Idents(seurat_singlets) == 6
]

cluster7_cells_use <- intersect(
  cluster7_cells,
  colnames(umi_tcell)
)

umi_cluster7 <- umi_tcell[, cluster7_cells_use, drop = FALSE]


# ============================================================
# 3. MERGE B CELLS AND CLUSTER 7 CELLS
# ============================================================

umi_merge <- cbind(
  umi_bcell,
  umi_cluster7
)

umi_merge <- as.matrix(umi_merge)
umi_merge[is.na(umi_merge)] <- 0
mode(umi_merge) <- "numeric"


# ============================================================
# 4. CREATE METADATA
# ============================================================

cell_group <- data.frame(
  cell_name = colnames(umi_merge),
  group = c(
    rep("B cell", ncol(umi_bcell)),
    rep("Cluster 7", ncol(umi_cluster7))
  ),
  stringsAsFactors = FALSE
)

cell_type <- rep("Other", ncol(umi_merge))
names(cell_type) <- colnames(umi_merge)

cell_type[colnames(umi_merge) %in% colnames(umi_bcell)] <- "B_cell"
cell_type[colnames(umi_merge) %in% cluster7_cells_use] <- "Cluster7_Tcell"

meta.data <- data.frame(
  cell_type = cell_type,
  group = cell_group$group,
  row.names = colnames(umi_merge)
)


# ============================================================
# 5. CREATE SEURAT OBJECT
# ============================================================

seu <- CreateSeuratObject(
  counts = umi_merge,
  meta.data = meta.data
)

seu@meta.data$cell_type <- as.factor(seu@meta.data$cell_type)

Idents(seu) <- seu@meta.data$cell_type

seu <- NormalizeData(seu)


# ============================================================
# 6. EXTRACT GENE EXPRESSION
# ============================================================

genes_of_interest <- c(
  "MS4A1",   # CD20
  "CD79A",
  "CD79B",
  "CD19",
  "SDC1",    # CD138
  "IGHG1",
  "IGHG2",
  "IGHG3",
  "IGHG4",
  "IGHM",
  "IGLC2"
)

df <- FetchData(
  seu,
  vars = c(genes_of_interest, "group"),
  slot = "data"
)

df$Cell <- rownames(df)


# ============================================================
# 7. CONVERT TO LONG FORMAT
# ============================================================

df_long <- melt(
  df,
  id.vars = c("group", "Cell"),
  variable.name = "Gene",
  value.name = "Expression"
)


# Rename gene symbols for plotting
df_long$Gene <- dplyr::recode(
  df_long$Gene,
  "MS4A1" = "CD20",
  "SDC1"  = "CD138"
)

# Gene order for plotting
gene_order <- c(
  "CD19",
  "CD20",
  "CD79A",
  "CD79B",
  "CD138",
  "IGHG1",
  "IGHG2",
  "IGHG3",
  "IGHG4",
  "IGHM",
  "IGLC2"
)

df_long$Gene <- factor(df_long$Gene, levels = gene_order)


# ============================================================
# 8. PLOT EXPRESSION (ALL CELLS)
# ============================================================

new_cols <- c(
  "B cell" = "#66C2A5",
  "Cluster 7" = "burlywood2"
)

figs1a <- ggplot(df_long, aes(x = Gene, y = Expression, fill = group)) +
  
  geom_boxplot(
    color = "black",
    width = 0.6,
    position = position_dodge(width = 0.9),
    outlier.shape = NA,
    linewidth = 0.5
  ) +
  
  stat_summary(
    fun.min = min,
    fun.max = max,
    geom = "errorbar",
    width = 0.25,
    linewidth = 0.5,
    color = "black",
    position = position_dodge(width = 0.9)
  ) +
  
  scale_fill_manual(values = new_cols) +
  
  coord_cartesian(ylim = c(0, 8)) +
  
  scale_y_continuous(
    expand = expansion(mult = c(0.05, 0.05))
  ) +
  
  scale_x_discrete(
    expand = expansion(add = 0.5)
  ) +
  
  theme_classic() +
  theme(
    axis.text.x = element_text(angle = 45, hjust = 1, size = 8, family = "Helvetica"),
    axis.text.y = element_text(size = 8, family = "Helvetica"),
    legend.position = "right"
  ) +
  
  xlab(NULL) +
  ylab(NULL)
ggsave(
  "figs1a.png",
  plot = figs1a,
  width = 6,
  height = 2,
  units = "in",
  dpi = 900,
  bg = "transparent"
)

# ============================================================
# 9. REMOVE CLUSTER 7 CELLS EXPRESSING B-CELL MARKERS
# ============================================================

genes_check <- c(
  "CD19",
  "CD20",
  "CD79A",
  "CD79B",
  "CD138"
)

cells_remove <- df_long %>%
  filter(
    group == "Cluster 7",
    Gene %in% genes_check,
    Expression > 1
  ) %>%
  distinct(Cell) %>%
  pull(Cell)

df_filter <- df_long %>%
  filter(!Cell %in% cells_remove)


# ============================================================
# 10. PLOT AFTER FILTERING
# ============================================================

figs1b <- ggplot(df_filter, aes(x = Gene, y = Expression, fill = group)) +
  
  geom_boxplot(
    color = "black",
    width = 0.6,
    position = position_dodge(width = 0.9),
    outlier.shape = NA,
    linewidth = 0.5
  ) +
  
  stat_summary(
    fun.min = min,
    fun.max = max,
    geom = "errorbar",
    width = 0.25,
    linewidth = 0.5,
    color = "black",
    position = position_dodge(width = 0.9)
  ) +
  
  scale_fill_manual(values = new_cols) +
  
  coord_cartesian(ylim = c(0, 8)) +
  
  scale_y_continuous(
    expand = expansion(mult = c(0.05, 0.05))
  ) +
  
  scale_x_discrete(
    expand = expansion(add = 0.5)
  ) +
  
  theme_classic() +
  theme(
    axis.text.x = element_text(angle = 45, hjust = 1, size = 8, family = "Helvetica"),
    axis.text.y = element_text(size = 8, family = "Helvetica"),
    legend.position = "right"
  ) +
  
  xlab(NULL) +
  ylab(NULL)
ggsave(
  "figs1b.png",
  plot = figs1b,
  width = 6,
  height = 2,
  units = "in",
  dpi = 900,
  bg = "transparent"
)

