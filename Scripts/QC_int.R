#!/usr/bin/env Rscript

suppressPackageStartupMessages({
  library(Seurat)
  library(scDblFinder)
  library(SingleCellExperiment)
  library(ggplot2)
  library(patchwork)
})

set.seed(1234)

# ------------------------------- 参数 ---------------------------------------

data_root <- "/public/home/wutong/Project/PJ6_Brain_Clab/Data"

gtf_file <- "/public/home/wutong/Project/PJ6_Brain_Clab/Ref_gtf/genes.gtf"

sample_dirs <- c(
  SP5C_1 = file.path(data_root, "SP5C_1", "filtered_feature_bc_matrix"),
  SP5C_2 = file.path(data_root, "SP5C_2", "filtered_feature_bc_matrix"),
  SP5C_3 = file.path(data_root, "SP5C_3", "filtered_feature_bc_matrix"),
  SP5C_4 = file.path(data_root, "SP5C_4", "filtered_feature_bc_matrix")
)

out_dir <- "/public/home/wutong/Project/PJ6_Brain_Clab/Re/Step1_QC_int"
dir.create(out_dir, recursive = TRUE, showWarnings = FALSE)

# 不在 GTF 中但必须保留的外源标记基因
extra_gene_names_to_keep <- c("tdTomato")

# 细胞质控参数
min_nfeature <- 500
min_ncount <- 2000
max_percent_mito <- 20

# 整合与聚类参数
n_variable_features <- 1000
n_integration_features <- 1000
n_pcs <- 30
cluster_resolution <- 0.1
rasterize_umap <- FALSE


# ------------------------------- 辅助函数 -----------------------------------

log_message <- function(...) {
  message(
    sprintf(
      "[%s] %s",
      format(Sys.time(), "%F %T"),
      paste0(...)
    )
  )
}

merge_object_list <- function(object_list) {
  if (length(object_list) == 1L) {
    return(object_list[[1]])
  }
  
  merge(
    x = object_list[[1]],
    y = object_list[-1]
  )
}

get_gtf_attribute <- function(x, key) {
  pattern <- paste0(
    "(?:^|;[[:space:]]*)",
    key,
    '[[:space:]]+"([^"]+)"'
  )
  
  matched <- regmatches(
    x,
    regexec(pattern, x, perl = TRUE)
  )
  
  vapply(
    matched,
    function(z) {
      if (length(z) > 1L) z[2] else NA_character_
    },
    character(1)
  )
}

remove_gene_id_version <- function(x) {
  sub("\\.[0-9]+$", "", x)
}


# -------------------------- 读取 GTF 蛋白编码基因 ----------------------------

log_message(
  "Seurat version: ",
  as.character(packageVersion("Seurat"))
)

log_message("Reading GTF: ", gtf_file)

gtf <- read.delim(
  gtf_file,
  header = FALSE,
  sep = "\t",
  quote = "",
  comment.char = "#",
  stringsAsFactors = FALSE
)

if (ncol(gtf) < 9L) {
  stop(
    "GTF 文件少于 9 列：",
    gtf_file,
    call. = FALSE
  )
}

colnames(gtf)[1:9] <- c(
  "chr", "source", "feature", "start", "end",
  "score", "strand", "frame", "attribute"
)

gtf <- gtf[
  gtf$feature == "gene",
  ,
  drop = FALSE
]

gtf$gene_id <- get_gtf_attribute(
  gtf$attribute,
  "gene_id"
)

gtf$gene_biotype <- get_gtf_attribute(
  gtf$attribute,
  "gene_biotype"
)

# 兼容使用 gene_type 而不是 gene_biotype 的 GTF
missing_biotype <- is.na(gtf$gene_biotype)

if (any(missing_biotype)) {
  gtf$gene_biotype[missing_biotype] <- get_gtf_attribute(
    gtf$attribute[missing_biotype],
    "gene_type"
  )
}

gtf$gene_id_clean <- remove_gene_id_version(
  gtf$gene_id
)

protein_coding_ids <- unique(
  gtf$gene_id_clean[
    !is.na(gtf$gene_id_clean) &
      gtf$gene_biotype == "protein_coding"
  ]
)

if (length(protein_coding_ids) == 0L) {
  stop(
    "GTF 中没有找到 gene_biotype/gene_type = protein_coding 的基因",
    call. = FALSE
  )
}

log_message(
  "Protein-coding gene IDs in GTF: ",
  length(protein_coding_ids)
)

rm(gtf)
invisible(gc())


# ---------------------------- 识别基因格式 ----------------------------------

detect_qc_patterns <- function(features, sample_id) {
  if (any(grepl("^MT-", features))) {
    return(list(
      style = "human-like (MT-/RPL/RPS)",
      mito = "^MT-",
      ribo = "^RP[SL]",
      hb = "^HB[^P]",
      remove = "^MT-|^RP[SL]|^ERCC"
    ))
  }
  
  if (any(grepl("^mt-", features))) {
    return(list(
      style = "mouse-like (mt-/Rpl/Rps)",
      mito = "^mt-",
      ribo = "^Rp[sl]",
      hb = "^Hb[^p]",
      remove = "^mt-|^Rp[sl]|^ERCC"
    ))
  }
  
  stop(
    sample_id,
    ": 未检测到以 MT- 或 mt- 开头的线粒体基因。",
    call. = FALSE
  )
}


# ---------------------------- 读取单个样本 ----------------------------------

read_one_sample <- function(data_dir, sample_id) {
  required_files <- c(
    "barcodes.tsv.gz",
    "features.tsv.gz",
    "matrix.mtx.gz"
  )
  
  missing_files <- required_files[
    !file.exists(file.path(data_dir, required_files))
  ]
  
  if (length(missing_files) > 0L) {
    stop(
      "Missing file(s) for ", sample_id, ": ",
      paste(missing_files, collapse = ", "),
      call. = FALSE
    )
  }
  
  counts <- Read10X(
    data.dir = data_dir,
    gene.column = 2,
    unique.features = TRUE
  )
  
  if (is.list(counts)) {
    if (!"Gene Expression" %in% names(counts)) {
      stop(
        sample_id,
        ": 没有找到 Gene Expression 矩阵",
        call. = FALSE
      )
    }
    
    counts <- counts[["Gene Expression"]]
  }
  
  feature_info <- read.delim(
    gzfile(file.path(data_dir, "features.tsv.gz")),
    header = FALSE,
    sep = "\t",
    quote = "",
    stringsAsFactors = FALSE
  )
  
  if (ncol(feature_info) < 2L) {
    stop(
      sample_id,
      ": features.tsv.gz 少于两列",
      call. = FALSE
    )
  }
  
  colnames(feature_info)[1:2] <- c(
    "gene_id",
    "gene_name"
  )
  
  if (ncol(feature_info) >= 3L) {
    colnames(feature_info)[3] <- "feature_type"
  }
  
  # 模拟 Read10X(unique.features = TRUE) 的唯一命名
  feature_info$seurat_name <- make.unique(
    feature_info$gene_name
  )
  
  # 如果包含多种模态，只保留 Gene Expression 注释
  if ("feature_type" %in% colnames(feature_info)) {
    gene_feature_info <- feature_info[
      feature_info$feature_type == "Gene Expression",
      ,
      drop = FALSE
    ]
  } else {
    gene_feature_info <- feature_info
  }
  
  # 兼容不同 Seurat 版本的 make.unique 顺序
  if (!identical(
    rownames(counts),
    gene_feature_info$seurat_name
  )) {
    gene_feature_info$seurat_name <- make.unique(
      gene_feature_info$gene_name
    )
  }
  
  if (!identical(
    rownames(counts),
    gene_feature_info$seurat_name
  )) {
    stop(
      sample_id,
      ": Read10X 行名无法与 features.tsv.gz 对应",
      call. = FALSE
    )
  }
  
  gene_feature_info$gene_id_clean <- remove_gene_id_version(
    gene_feature_info$gene_id
  )
  
  # GTF 中注释为 protein_coding 的基因
  keep_by_annotation <- (
    gene_feature_info$gene_id_clean %in%
      protein_coding_ids
  )
  
  # 不在 GTF 中但必须保留的外源基因
  keep_by_whitelist <- (
    gene_feature_info$gene_name %in%
      extra_gene_names_to_keep
  )
  
  missing_extra_genes <- setdiff(
    extra_gene_names_to_keep,
    gene_feature_info$gene_name
  )
  
  if (length(missing_extra_genes) > 0L) {
    stop(
      sample_id,
      ": features.tsv.gz 中没有找到必须保留的外源基因：",
      paste(missing_extra_genes, collapse = ", "),
      call. = FALSE
    )
  }
  
  extra_features <- gene_feature_info$seurat_name[
    keep_by_whitelist
  ]
  
  # 最终候选：蛋白编码基因或外源白名单基因
  downstream_features <- gene_feature_info$seurat_name[
    keep_by_annotation | keep_by_whitelist
  ]
  
  object <- CreateSeuratObject(
    counts = counts,
    project = sample_id,
    min.cells = 0,
    min.features = 0
  )
  
  object <- RenameCells(
    object,
    add.cell.id = sample_id
  )
  
  object$sample <- sample_id
  object$batch <- sample_id
  object$samplebatch <- sample_id
  
  qc_patterns <- detect_qc_patterns(
    rownames(object),
    sample_id
  )
  
  object@misc$qc_patterns <- qc_patterns
  
  object@misc$downstream_features <- intersect(
    rownames(object),
    downstream_features
  )
  
  object@misc$required_extra_features <- intersect(
    rownames(object),
    extra_features
  )
  
  object@misc$feature_annotation <- gene_feature_info[
    gene_feature_info$seurat_name %in%
      object@misc$downstream_features,
    c(
      "gene_id",
      "gene_name",
      "seurat_name"
    ),
    drop = FALSE
  ]
  
  # 必须在删除基因前计算质控指标
  object[["percent_mito"]] <- PercentageFeatureSet(
    object,
    pattern = qc_patterns$mito
  )
  
  object[["percent_ribo"]] <- PercentageFeatureSet(
    object,
    pattern = qc_patterns$ribo
  )
  
  object[["percent_hb"]] <- PercentageFeatureSet(
    object,
    pattern = qc_patterns$hb
  )
  
  log_message(
    sample_id,
    ": total features = ", nrow(object),
    "; protein-coding + whitelist features = ",
    length(object@misc$downstream_features),
    "; retained external marker = ",
    paste(
      object@misc$required_extra_features,
      collapse = ", "
    ),
    "; naming style = ",
    qc_patterns$style
  )
  
  object
}


# ------------------------------- 单样本质控 ---------------------------------

qc_one_sample <- function(object, sample_id) {
  metadata <- object[[]]
  
  pass_qc <- with(
    metadata,
    nFeature_RNA > min_nfeature &
      nCount_RNA > min_ncount &
      percent_mito < max_percent_mito
  )
  
  object <- subset(
    object,
    cells = rownames(metadata)[pass_qc]
  )
  
  if (ncol(object) == 0L) {
    stop(
      sample_id,
      ": 固定阈值过滤后没有剩余细胞",
      call. = FALSE
    )
  }
  
  n_after_threshold_qc <- ncol(object)
  
  median_nfeature_after_qc <- median(
    object$nFeature_RNA
  )
  
  median_ncount_after_qc <- median(
    object$nCount_RNA
  )
  
  # 蛋白编码基因 + tdT
  retained_genes <- intersect(
    rownames(object),
    object@misc$downstream_features
  )
  
  # 完全删除线粒体、核糖体和 ERCC
  retained_genes <- retained_genes[
    !grepl(
      object@misc$qc_patterns$remove,
      retained_genes
    )
  ]
  
  if (length(retained_genes) == 0L) {
    stop(
      sample_id,
      ": 基因过滤后没有剩余基因",
      call. = FALSE
    )
  }
  
  object <- subset(
    object,
    features = retained_genes
  )
  
  # 强制确认 tdT 没有被删除
  missing_extra_features <- setdiff(
    object@misc$required_extra_features,
    rownames(object)
  )
  
  if (length(missing_extra_features) > 0L) {
    stop(
      sample_id,
      ": 外源标记基因被意外删除：",
      paste(missing_extra_features, collapse = ", "),
      call. = FALSE
    )
  }
  
  # 强制确认 MT、RPL/RPS、ERCC 已完全删除
  remaining_unwanted <- grep(
    object@misc$qc_patterns$remove,
    rownames(object),
    value = TRUE
  )
  
  if (length(remaining_unwanted) > 0L) {
    stop(
      sample_id,
      ": 仍存在未删除的 MT/RPL/RPS/ERCC 基因：",
      paste(remaining_unwanted, collapse = ", "),
      call. = FALSE
    )
  }
  
  log_message(
    sample_id,
    ": cells after threshold QC = ",
    ncol(object),
    "; downstream features = ",
    nrow(object),
    "; retained external marker = ",
    paste(
      object@misc$required_extra_features,
      collapse = ", "
    )
  )
  
  # 先标准化，兼容 Seurat v4/v5 并避免缺少 data layer 的提示
  object <- NormalizeData(
    object,
    verbose = FALSE
  )
  
  # 每个样本独立检测双细胞
  set.seed(1234)
  
  sce <- as.SingleCellExperiment(
    object,
    assay = "RNA"
  )
  
  sce <- scDblFinder::scDblFinder(
    sce,
    verbose = FALSE
  )
  
  object$scDblFinder.score <- sce$scDblFinder.score
  
  object$scDblFinder.class <- as.character(
    sce$scDblFinder.class
  )
  
  n_singlet <- sum(
    object$scDblFinder.class == "singlet"
  )
  
  n_doublet <- sum(
    object$scDblFinder.class == "doublet"
  )
  
  object <- subset(
    object,
    cells = colnames(object)[
      object$scDblFinder.class == "singlet"
    ]
  )
  
  object <- FindVariableFeatures(
    object,
    selection.method = "vst",
    nfeatures = n_variable_features,
    verbose = FALSE
  )
  
  qc_summary <- data.frame(
    sample = sample_id,
    n_cells_before_qc = nrow(metadata),
    min_nFeature_exclusive = min_nfeature,
    min_nCount_exclusive = min_ncount,
    max_percent_mito_exclusive = max_percent_mito,
    n_cells_after_threshold_qc = n_after_threshold_qc,
    median_nFeature_after_qc = median_nfeature_after_qc,
    median_nCount_after_qc = median_ncount_after_qc,
    n_singlet = n_singlet,
    n_doublet = n_doublet,
    doublet_rate = n_doublet / (n_singlet + n_doublet),
    n_cells_final = ncol(object),
    n_downstream_features = nrow(object),
    tdT_retained = (
      "tdTomato" %in% rownames(object)
    ),
    row.names = NULL,
    check.names = FALSE
  )
  
  list(
    object = object,
    summary = qc_summary
  )
}


# ----------------------- nFeature/nCount 质控图 ------------------------------

make_qc_plot_data <- function(object, stage) {
  metadata <- object[[]]
  
  data.frame(
    sample = metadata$sample,
    nFeature_RNA = metadata$nFeature_RNA,
    nCount_RNA = metadata$nCount_RNA,
    stage = stage,
    stringsAsFactors = FALSE
  )
}

save_before_after_qc_plots <- function(
    before_object,
    after_object
) {
  plot_data <- rbind(
    make_qc_plot_data(
      before_object,
      "Before QC"
    ),
    make_qc_plot_data(
      after_object,
      "After QC"
    )
  )
  
  plot_data$stage <- factor(
    plot_data$stage,
    levels = c(
      "Before QC",
      "After QC"
    )
  )
  
  violin_data <- rbind(
    data.frame(
      sample = plot_data$sample,
      stage = plot_data$stage,
      metric = "nFeature_RNA",
      value = plot_data$nFeature_RNA
    ),
    data.frame(
      sample = plot_data$sample,
      stage = plot_data$stage,
      metric = "nCount_RNA",
      value = plot_data$nCount_RNA
    )
  )
  
  violin_plot <- ggplot(
    violin_data,
    aes(
      x = sample,
      y = value,
      fill = sample
    )
  ) +
    geom_violin(
      scale = "width",
      trim = TRUE
    ) +
    facet_grid(
      stage ~ metric,
      scales = "free_y"
    ) +
    labs(
      title = "nFeature/nCount before and after QC",
      x = NULL,
      y = NULL
    ) +
    theme_classic() +
    theme(
      legend.position = "none",
      axis.text.x = element_text(
        angle = 45,
        hjust = 1
      ),
      strip.text = element_text(
        face = "bold",
        size = 11
      )
    )
  
  scatter_plot <- ggplot(
    plot_data,
    aes(
      x = nCount_RNA,
      y = nFeature_RNA,
      color = sample
    )
  ) +
    geom_point(
      size = 0.25,
      alpha = 0.5
    ) +
    facet_wrap(
      ~stage,
      ncol = 1,
      scales = "free"
    ) +
    labs(
      title = "nCount versus nFeature before and after QC",
      x = "nCount_RNA",
      y = "nFeature_RNA",
      color = "Sample"
    ) +
    theme_classic()
  
  ggsave(
    file.path(
      out_dir,
      "QC_nFeature_nCount_before_after_violin.pdf"
    ),
    plot = violin_plot,
    width = 12,
    height = 9
  )
  
  ggsave(
    file.path(
      out_dir,
      "QC_nFeature_nCount_before_after_scatter.pdf"
    ),
    plot = scatter_plot,
    width = 10,
    height = 10
  )
}


# ---------------------------- 数据读取与质控 --------------------------------

log_message("Reading four 10x matrices")

raw_list <- Map(
  f = read_one_sample,
  data_dir = unname(sample_dirs),
  sample_id = names(sample_dirs)
)

names(raw_list) <- names(sample_dirs)

before_qc <- merge_object_list(raw_list)

saveRDS(
  before_qc,
  file.path(
    out_dir,
    "SP5C_before_QC.rds"
  )
)

log_message(
  "Applying fixed QC thresholds and doublet removal"
)

qc_results <- Map(
  f = qc_one_sample,
  object = raw_list,
  sample_id = names(raw_list)
)

names(qc_results) <- names(raw_list)

seurat_list <- lapply(
  qc_results,
  `[[`,
  "object"
)

qc_table <- do.call(
  rbind,
  lapply(
    qc_results,
    `[[`,
    "summary"
  )
)

after_qc <- merge_object_list(
  seurat_list
)

save_before_after_qc_plots(
  before_qc,
  after_qc
)

write.table(
  qc_table,
  file = file.path(
    out_dir,
    "SP5C_QC_summary.tsv"
  ),
  sep = "\t",
  quote = FALSE,
  row.names = FALSE
)

saveRDS(
  seurat_list,
  file.path(
    out_dir,
    "SP5C_after_QC_list.rds"
  )
)

saveRDS(
  after_qc,
  file.path(
    out_dir,
    "SP5C_after_QC_merged.rds"
  )
)


# ------------------------------ RPCA 整合 -----------------------------------

log_message("Selecting integration features")

integration_features <- SelectIntegrationFeatures(
  object.list = seurat_list,
  nfeatures = n_integration_features
)

log_message(
  "Scaling each sample and running PCA"
)

seurat_list <- lapply(
  seurat_list,
  function(object) {
    object <- ScaleData(
      object,
      features = integration_features,
      verbose = FALSE
    )
    
    object <- RunPCA(
      object,
      features = integration_features,
      npcs = n_pcs,
      verbose = FALSE
    )
    
    object
  }
)

log_message("Finding RPCA anchors")

integration_anchors <- FindIntegrationAnchors(
  object.list = seurat_list,
  anchor.features = integration_features,
  reduction = "rpca",
  dims = seq_len(n_pcs)
)

log_message("Integrating four samples")

integrated <- IntegrateData(
  anchorset = integration_anchors,
  dims = seq_len(n_pcs)
)

rm(
  integration_anchors,
  seurat_list
)

invisible(gc())


# -------------------------- 降维、聚类与保存 --------------------------------

log_message(
  "Scaling, dimensional reduction and clustering"
)

DefaultAssay(integrated) <- "integrated"

integrated <- ScaleData(
  integrated,
  verbose = FALSE
)

integrated <- RunPCA(
  integrated,
  npcs = n_pcs,
  verbose = FALSE
)

integrated <- FindNeighbors(
  integrated,
  reduction = "pca",
  dims = 1:10,
  k.param = 60,
  verbose = FALSE
)

integrated <- RunUMAP(
  integrated,
  reduction = "pca",
  dims = 1:5,
  n.neighbors = 200,
  min.dist = 2,
  spread = 1,
  metric = "euclidean",
  seed.use = 919,
  verbose = FALSE
)


cluster_plot <- DimPlot(
  integrated,
  reduction = "umap",
  label = TRUE,
  repel = TRUE
)

integrated <- FindClusters(
  integrated,
  resolution = cluster_resolution,
  random.seed = 1234,
  verbose = FALSE
)



saveRDS(
  integrated,
  file.path(
    out_dir,
    "SP5C_integrated_analyzed_V2.rds"
  )
)


# ------------------------------- UMAP ---------------------------------------

samplebatch_plot <- DimPlot(
  integrated,
  reduction = "umap",
  group.by = "samplebatch"
) +
  ggtitle("Sample batch")

sample_plot <- DimPlot(
  integrated,
  reduction = "umap",
  group.by = "sample"
) +
  ggtitle("Sample")

cluster_plot <- DimPlot(
  integrated,
  reduction = "umap",
  label = TRUE,
  repel = TRUE
) +
  ggtitle(
    paste0(
      "Clusters (resolution = ",
      cluster_resolution,
      ")"
    )
  )

ggsave(
  file.path(
    out_dir,
    paste0(
      "SP5C_UMAP_resolution_",
      cluster_resolution,
      ".pdf"
    )
  ),
  plot = samplebatch_plot +
    sample_plot +
    cluster_plot,
  width = 16,
  height = 4
)


# ---------------------------- cluster QC 图 ---------------------------------

qc_features <- c(
  "nFeature_RNA",
  "nCount_RNA",
  "percent_mito",
  "percent_ribo",
  "percent_hb",
  "scDblFinder.score"
)

cluster_qc_plot <- VlnPlot(
  integrated,
  features = qc_features,
  group.by = "seurat_clusters",
  pt.size = 0,
  ncol = 2
)

ggsave(
  file.path(
    out_dir,
    "cluster_QC_violin.pdf"
  ),
  plot = cluster_qc_plot,
  width = 16,
  height = 12
)

capture.output(
  sessionInfo(),
  file = file.path(
    out_dir,
    "sessionInfo.txt"
  )
)

log_message(
  "Finished. Results written to: ",
  out_dir
)