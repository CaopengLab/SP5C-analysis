library(Seurat)
library(dplyr)
library(ggplot2)
library(patchwork)

# ============================================================
# 0. Parameters and output directory


out_dir <- paste0(
  "/public/home/wutong/Project/PJ6_Brain_Clab/",
  "Re/step4_Neu_recluster_SP5C_region"
)

dir.create(
  out_dir,
  recursive = TRUE,
  showWarnings = FALSE
)

# Analysis parameters
n_hvg <- 2000
n_pcs <- 30
use_dims <- 1:20
neighbor_k <- 20
resolution_use <- 0.1
random_seed <- 1234


# ============================================================
# 1. Extract neurons


neu <- subset(
  integrated,
  subset = celltype_manual %in% c("Exc Neu")
)

DefaultAssay(neu) <- "RNA"

cat("神经元数量：", ncol(neu), "\n")
print(table(neu$celltype_manual))

# Automatically identify the sample metadata column
if ("sample" %in% colnames(neu[[]])) {
  sample_col <- "sample"
} else if ("orig.ident" %in% colnames(neu[[]])) {
  sample_col <- "orig.ident"
} else {
  stop("metadata中没有找到sample或orig.ident")
}

cat("使用的样本信息列：", sample_col, "\n")


# ============================================================
# 2. Normalize the neuronal subset and identify highly variable genes


neu <- NormalizeData(
  neu,
  normalization.method = "LogNormalize",
  scale.factor = 10000,
  verbose = FALSE
)

neu <- FindVariableFeatures(
  neu,
  selection.method = "vst",
  nfeatures = n_hvg,
  verbose = FALSE
)

cat(
  "神经元特异HVG数量：",
  length(VariableFeatures(neu)),
  "\n"
)


# ============================================================
# 4. Scale data and run PCA

neu <- ScaleData(
  neu,
  features = VariableFeatures(neu),
  verbose = FALSE
)

neu <- RunPCA(
  neu,
  features = VariableFeatures(neu),
  npcs = n_pcs,
  seed.use = random_seed,
  verbose = FALSE
)

p_elbow <- ElbowPlot(
  neu,
  ndims = n_pcs
) +
  ggtitle("Neuronal PCA")

ggsave(
  filename = file.path(
    out_dir,
    "SP5C_Neu_recluster_ElbowPlot.pdf"
  ),
  plot = p_elbow,
  width = 6,
  height = 5
)


# NC paper：SNN K = 20
neu <- FindNeighbors(
  neu,
  reduction = "pca",
  dims = use_dims,
  k.param = neighbor_k,
  verbose = FALSE
)

# Louvain clustering
neu <- FindClusters(
  neu,
  resolution = 0.5,
  algorithm = 1,
  random.seed = random_seed,
  verbose = FALSE
)

Idents(neu) <- "seurat_clusters"

cat("各cluster细胞数量：\n")
print(sort(table(Idents(neu)), decreasing = TRUE))

neu <- RunUMAP(
  neu,
  reduction = "pca",
  dims = 1:10,
  n.neighbors = 150,
  min.dist = 0.5,
  spread = 0.5,
  metric = "euclidean",
  seed.use = 001,
  verbose = FALSE
)

DimPlot(
  neu,
  reduction = "umap",
  group.by = "CellTypeN",
  label = TRUE,
  label.size = 5,
  repel = TRUE,
  pt.size = 0.7,
  shuffle = TRUE,
  seed = random_seed,
  raster = FALSE
) 
saveRDS(neu2,"/public/home/wutong/Project/PJ6_Brain_Clab/Re/step4_Neu_recluster_SP5C_region/Exc_Nue.rds")


# Merge subclusters with the same top DEGs
# cluster_to_CellTypeN <- c(
#   "12" = "N1",
#   "16" = "N1",
#   
#   "2"  = "N2",
#   "7"  = "N2",
#   "18" = "N2",
#   
#   "0"  = "N3",
#   "6"  = "N3",
#   "11" = "N3",
#   "19" = "N3",
#   
#   "13" = "N4",
#   "20" = "N4",
#   
#   "1"  = "N5",
#   "3"  = "N6",
#   "4"  = "N7",
#   "8"  = "N8",
#   "9"  = "N9",
#   "14" = "N10",
#   "17" = "N11",
#   
#   "5"  = "N12",  # previously N8
#   "10" = "N13",  # previously N11
#   "15" = "N14"   # previously N13
# )

cluster_to_CellTypeN <- c(
  "12" = "N1",
  "16" = "N1",
  
  "2"  = "N2",
  "7"  = "N2",
  "18" = "N2",
  
  "0"  = "N3",
  "6"  = "N3",
  "11" = "N3",
  "19" = "N3",
  
  "13" = "N4",
  "20" = "N4",
  
  "1"  = "N5",
  "3"  = "N6",
  "4"  = "N7",
  "8"  = "N8",
  "9"  = "N9",
  
  # Merge into N10
  "14" = "N10",
  "17" = "N10",
  "10" = "N10",
  
  # Renumber subsequent clusters
  "5"  = "N11",
  "15" = "N12"
)

# Create CellTypeN from RNA_snn_res.0.5
neu$CellTypeN <- unname(
  cluster_to_CellTypeN[
    as.character(neu$RNA_snn_res.0.5)
  ]
)

# Set factor level order
neu$CellTypeN <- factor(
  neu$CellTypeN,
  levels = paste0("N", 1:12)
)

# Check cluster-to-cell-type assignments
table(
  original_cluster = neu$RNA_snn_res.0.5,
  CellTypeN = neu$CellTypeN,
  useNA = "ifany"
)

# Check cell counts in the new clusters
table(neu$CellTypeN)

# Set as the active identity
Idents(neu) <- "CellTypeN"

p_cluster <- DimPlot(
  neu,
  reduction = "umap",
  group.by = "CellTypeN",
  label = TRUE,
  label.size = 5,
  repel = TRUE,
  pt.size = 0.5,
  shuffle = TRUE,
  seed = random_seed,
  raster = FALSE
) +
  theme_classic()

out_dir = "/public/home/wutong/Project/PJ6_Brain_Clab/Re/step4_Neu_recluster_SP5C_region"
ggsave(
  filename = file.path(
    out_dir,
    "SP5C_Neu_recluster_cluster_UMAP.pdf"
  ),
  plot = p_cluster,
  width = 6,
  height = 5
)

saveRDS(neu,"/public/home/wutong/Project/PJ6_Brain_Clab/Re/step4_Neu_recluster_SP5C_region/Exc_Nue.rds")

#-----
# DEG
# ============================================================
# Find marker DEGs for each cluster
DefaultAssay(neu) <- "RNA"
Idents(neu) <- neu$CellTypeN

deg_all <- FindAllMarkers(
  object = neu,
  assay = "RNA",
  only.pos = TRUE,
  test.use = "wilcox",
  logfc.threshold = 0.25,
  return.thresh = 1,
  verbose = FALSE
)

saveRDS(deg_all,"/public/home/wutong/Project/PJ6_Brain_Clab/Re/step4_Neu_recluster_SP5C_region/DEG_celltypeN.rds")
write.csv(deg_all,"/public/home/wutong/Project/PJ6_Brain_Clab/Re/step4_Neu_recluster_SP5C_region/DEG_celltypeN.csv",row.names = TRUE)
# NC paper：
# Wilcoxon test
# Bonferroni correction
# avg_log2FC > 0.25
# p_val_adj < 0.05


# Filtering criteria from the NC paper
deg_sig_unique <- deg_all %>%
  tibble::as_tibble() %>%
  dplyr::filter(
    avg_log2FC > 0.25,
    p_val_adj < 0.05
  ) %>%
  dplyr::mutate(
    pct_difference = pct.1 - pct.2,
    rank_score =
      -log10(pmax(p_val, 1e-300)) *
      avg_log2FC
  ) %>%
  
  # Assign each gene only to the cluster in which it is most specific
  dplyr::group_by(gene) %>%
  dplyr::arrange(
    dplyr::desc(pct_difference),
    dplyr::desc(rank_score),
    dplyr::desc(avg_log2FC),
    .by_group = TRUE
  ) %>%
  dplyr::slice_head(n = 1) %>%
  dplyr::ungroup() %>%
  
  # Sort by cluster and rank_score
  dplyr::arrange(
    cluster,
    dplyr::desc(rank_score)
  )
write.csv(deg_sig_unique,"/public/home/wutong/Project/PJ6_Brain_Clab/Re/step4_Neu_recluster_SP5C_region/DEG_celltypeN_filter.csv",row.names = TRUE)
saveRDS(deg_sig_unique,"/public/home/wutong/Project/PJ6_Brain_Clab/Re/step4_Neu_recluster_SP5C_region/DEG_celltypeN_filter.rds")
#Select the top 10 genes for SP5C subdivision
deg_sig_unique %>%
  group_by(cluster) %>%
  slice_max(rank_score, n = 10, with_ties = FALSE) %>%
  ungroup() %>%
  write.csv("/public/home/wutong/Project/PJ6_Brain_Clab/Re/step4_Neu_recluster_SP5C_region/DEG_celltypeN_Top10_markers_by_LogFC_pv.csv",row.names = TRUE)

#-----
library(Seurat)
library(ggplot2)
library(patchwork)

dir.create(outdir, recursive = TRUE, showWarnings = FALSE)

fp <- function(gene, ord) {
  FeaturePlot(
    neu,
    features = gene,
    reduction = "umap",
    order = TRUE,
    min.cutoff = "q05",
    max.cutoff = "q95",
    pt.size = 2,
    cols = c("#440154", "#21908C", "#FDE725"),
    raster = TRUE
  ) 
}

p_Penk     <- fp("Penk", FALSE)
p_Nmbr     <- fp("Nmbr", TRUE)
p_tdTomato <- fp("tdTomato", TRUE)

p_all <- p_Penk | p_Nmbr | p_tdTomato
ggsave(
  file.path(outdir, "Penk_Nmbr_tdTomato_FeaturePlot_v2.pdf"),
  plot = p_all,
  width = 12,
  height = 3.2
)
