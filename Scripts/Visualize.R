library(Seurat)
library(dplyr)
library(tibble)
library(ggplot2)

# ============================================================
# 1. 路径、读取数据
# ============================================================

object_file <- paste0(
  "/public/home/wutong/Project/PJ6_Brain_Clab/Re/",
  "step4_Neu_recluster_SP5C_region/Exc_Nue.rds"
)

deg_file <- paste0(
  "/public/home/wutong/Project/PJ6_Brain_Clab/Re/",
  "step4_Neu_recluster_SP5C_region/DEG_celltypeN_filter.rds"
)

out_dir <- paste0(
  "/public/home/wutong/Project/PJ6_Brain_Clab/Re/",
  "Step6_state_tdT_Penk_commen_ratio"
)

dir.create(out_dir, recursive = TRUE, showWarnings = FALSE)

neu <- readRDS(object_file)
deg <- tibble::as_tibble(readRDS(deg_file))

DefaultAssay(neu) <- "RNA"
cluster_levels <- paste0("N", 1:12)

# ============================================================
# 2. 明确区分原始counts、RDS已有data及重新标准化的data
# ============================================================

expr_counts <- GetAssayData(
  neu,
  assay = "RNA",
  slot = "counts"
)

expr_saved_data <- GetAssayData(
  neu,
  assay = "RNA",
  slot = "data"
)

if (!"tdTomato" %in% rownames(expr_counts)) {
  stop("RNA assay的counts中没有找到tdTomato")
}

# 从counts重新进行标准化，避免依赖RDS中data槽的历史处理方式
expr_normalized <- GetAssayData(
  neu_normalized,
  assay = "RNA",
  slot = "data"
)

stopifnot(
  identical(colnames(expr_counts), colnames(expr_normalized))
)

# tdTomato只作阳性/阴性判断：原始UMI至少1个
tdT_positive <- as.numeric(
  expr_counts["tdTomato", ]
) > 0

total_cells <- ncol(expr_counts)
tdT_positive_n <- sum(tdT_positive)
global_tdT_percent <- 100 * tdT_positive_n / total_cells

cat("总细胞数：", total_cells, "\n")
cat("tdT阳性细胞数：", tdT_positive_n, "\n")
cat(
  "总体tdT阳性率：",
  round(global_tdT_percent, 2),
  "%\n"
)

# 直接检查Penk、Abcc9在三种矩阵中>2的细胞数量
check_genes <- intersect(
  c("Penk", "Abcc9"),
  rownames(expr_counts)
)

expression_check <- dplyr::bind_rows(
  lapply(check_genes, function(g) {
    saved <- as.numeric(expr_saved_data[g, ])
    counts <- as.numeric(expr_counts[g, ])
    normalized <- as.numeric(expr_normalized[g, ])
    
    data.frame(
      gene = g,
      
      saved_data_gt2_n =
        sum(saved > 2),
      
      saved_data_gt2_tdT_n =
        sum(saved > 2 & tdT_positive),
      
      raw_counts_gt2_n =
        sum(counts > 2),
      
      raw_counts_gt2_tdT_n =
        sum(counts > 2 & tdT_positive),
      
      LogNormalize_gt2_n =
        sum(normalized > 2),
      
      LogNormalize_gt2_tdT_n =
        sum(normalized > 2 & tdT_positive),
      
      saved_data_equals_counts =
        isTRUE(all.equal(saved, counts))
    )
  })
)

# ============================================================
# 3. 准备DEG排名
# ============================================================

required_cols <- c(
  "cluster", "gene", "p_val", "p_val_adj",
  "avg_log2FC", "pct.1", "pct.2"
)

missing_cols <- setdiff(required_cols, colnames(deg))

if (length(missing_cols) > 0) {
  stop(
    "DEG缺少列：",
    paste(missing_cols, collapse = ", ")
  )
}

deg <- deg %>%
  dplyr::mutate(
    cluster = as.character(cluster),
    gene = as.character(gene)
  )

# 若RDS中已有rank_score，沿用原有排名
if (!"rank_score" %in% colnames(deg)) {
  deg <- deg %>%
    dplyr::mutate(
      rank_score =
        -log10(pmax(p_val, 1e-300)) *
        avg_log2FC
    )
}

# ============================================================
# 4. 每个基因分配给特异性最高的cluster
#    再提取每群rank_score最高的Top1基因
# ============================================================

deg_sig_unique <- deg %>%
  dplyr::filter(
    cluster %in% cluster_levels,
    gene %in% rownames(expr_normalized),
    !is.na(rank_score),
    !is.na(pct.1),
    !is.na(pct.2),
    avg_log2FC > 0.25,
    p_val_adj < 0.05,
    pct.1 > 0.6
  ) %>%
  dplyr::mutate(
    pct_difference = pct.1 - pct.2
  ) %>%
  dplyr::group_by(gene) %>%
  dplyr::arrange(
    dplyr::desc(pct_difference),
    dplyr::desc(rank_score),
    dplyr::desc(avg_log2FC),
    .by_group = TRUE
  ) %>%
  dplyr::slice_head(n = 1) %>%
  dplyr::ungroup()

stopifnot(!anyDuplicated(deg_sig_unique$gene))

top1_marker <- deg_sig_unique %>%
  dplyr::group_by(cluster) %>%
  dplyr::arrange(
    dplyr::desc(rank_score),
    .by_group = TRUE
  ) %>%
  dplyr::slice_head(n = 1) %>%
  dplyr::ungroup() %>%
  dplyr::mutate(
    cluster = factor(cluster, levels = cluster_levels)
  ) %>%
  dplyr::arrange(cluster) %>%
  dplyr::select(
    cluster, gene, rank_score,
    pct.1, pct.2, pct_difference
  )

if (nrow(top1_marker) != 12) {
  stop(
    "只提取到", nrow(top1_marker),
    "个cluster的Top1基因"
  )
}

cat("\n12个cluster的Top1基因：\n")
print(
  tibble::as_tibble(top1_marker),
  n = 12,
  width = 200
)

write.csv(
  top1_marker,
  file.path(out_dir, "Top1_marker_by_cluster.csv"),
  row.names = FALSE
)

# ============================================================
# 5. 标准化表达量>2的tdT共标统计
#    在全部细胞中计算，不限制cluster
# ============================================================

expression_threshold <- 2

top1_threshold2_result <- dplyr::bind_rows(
  lapply(seq_len(nrow(top1_marker)), function(i) {
    
    marker_gene <- top1_marker$gene[i]
    marker_cluster <- as.character(top1_marker$cluster[i])
    
    # 使用刚刚明确执行LogNormalize后得到的表达值
    gene_expression <- as.numeric(
      expr_normalized[marker_gene, ]
    )
    
    DEG_high <- gene_expression > expression_threshold
    double_positive <- tdT_positive & DEG_high
    
    DEG_high_n <- sum(DEG_high)
    double_positive_n <- sum(double_positive)
    
    a <- double_positive_n
    b <- DEG_high_n - double_positive_n
    c <- tdT_positive_n - double_positive_n
    d <- total_cells - a - b - c
    
    fisher_result <- fisher.test(
      matrix(
        c(a, b, c, d),
        nrow = 2,
        byrow = TRUE
      ),
      alternative = "greater"
    )
    
    percent_in_DEG <- if (DEG_high_n > 0) {
      100 * double_positive_n / DEG_high_n
    } else {
      NA_real_
    }
    
    data.frame(
      threshold = expression_threshold,
      cluster = marker_cluster,
      gene = marker_gene,
      total_cells = total_cells,
      tdT_positive_n = tdT_positive_n,
      
      DEG_high_n =
        DEG_high_n,
      
      tdT_DEG_double_positive_n =
        double_positive_n,
      
      # 双阳数 / DEG高表达细胞数
      double_positive_in_DEG_percent =
        percent_in_DEG,
      
      # 双阳数 / 所有tdT阳性细胞数
      double_positive_in_tdT_percent =
        100 * double_positive_n / tdT_positive_n,
      
      global_tdT_percent =
        global_tdT_percent,
      
      fold_over_global =
        percent_in_DEG / global_tdT_percent,
      
      odds_ratio =
        unname(fisher_result$estimate),
      
      p_val =
        fisher_result$p.value
    )
  })
) %>%
  dplyr::mutate(
    p_val_adj = p.adjust(p_val, method = "BH"),
    
    fold_rank = dplyr::min_rank(
      dplyr::desc(fold_over_global)
    ),
    
    cluster = factor(
      cluster,
      levels = cluster_levels
    )
  ) %>%
  dplyr::arrange(fold_rank, cluster)

cat("\nLogNormalize表达量>2的统计结果：\n")
print(
  tibble::as_tibble(top1_threshold2_result),
  n = nrow(top1_threshold2_result),
  width = 200
)

write.csv(
  top1_threshold2_result,
  file.path(
    out_dir,
    "Top1_LogNormalize_gt2_tdT_statistics.csv"
  ),
  row.names = FALSE
)

# ============================================================
# 6. UMAP配色
# ============================================================

all_cluster_cols <- c(
  "N1"  = "#cd1d1d",
  "N2"  = "#17BECF",
  "N3"  = "#8DA0CB",
  "N4"  = "#E78AC3",
  "N5"  = "#A6D854",
  "N6"  = "#FFD92F",
  "N7"  = "#1F78B4",
  "N8"  = "#66C2A5",
  "N9"  = "#6A3D9A",
  "N10" = "#FF7F00",
  "N11" = "#DCDCDC",
  "N12" = "#DCDCDC"
)

# ============================================================
# 7. Fold rank棒棒糖图
# ============================================================

plot_df <- top1_threshold2_result %>%
  dplyr::filter(is.finite(fold_over_global)) %>%
  dplyr::arrange(
    dplyr::desc(fold_over_global),
    cluster
  ) %>%
  dplyr::mutate(
    cluster = as.character(cluster),
    marker = paste(cluster, gene),
    marker = factor(marker, levels = rev(marker)),
    fold_label = sprintf("%.2f×", fold_over_global),
    label_x = pmax(fold_over_global, 0.06)
  )

p_fold_rank <- ggplot(
  plot_df,
  aes(y = marker)
) +
  geom_vline(
    xintercept = 1,
    linetype = "dashed",
    colour = "grey50"
  ) +
  geom_segment(
    aes(
      x = 1,
      xend = fold_over_global,
      yend = marker,
      colour = cluster
    ),
    size = 1
  ) +
  geom_point(
    aes(
      x = fold_over_global,
      fill = cluster
    ),
    shape = 21,
    colour = "black",
    stroke = 0.35,
    size = 3.5
  ) +
  geom_text(
    aes(
      x = label_x,
      label = fold_label
    ),
    vjust = -0.9,
    size = 3.2
  ) +
  scale_colour_manual(values = all_cluster_cols) +
  scale_fill_manual(values = all_cluster_cols) +
  scale_x_continuous(
    limits = c(
      0,
      max(plot_df$fold_over_global) + 0.15
    ),
    expand = c(0, 0)
  ) +
  scale_y_discrete(
    expand = expansion(add = c(0.5, 0.7))
  ) +
  labs(
    x = "Fold enrichment relative to overall tdT+ rate",
    y = NULL,
    title = "tdT colabeling of Top1 DEGs",
    subtitle = "LogNormalize RNA expression > 2",
    caption = "Dashed line = 1 (overall tdT+ rate)"
  ) +
  theme_classic(base_size = 11) +
  theme(
    legend.position = "none",
    axis.text = element_text(colour = "black"),
    plot.title = element_text(face = "bold", hjust = 0.5),
    plot.subtitle = element_text(hjust = 0.5),
    plot.caption = element_text(hjust = 0, size = 8),
    panel.border = element_rect(
      colour = "black",
      fill = NA,
      size = 0.4
    )
  )

p_fold_rank

ggsave(
  file.path(
    out_dir,
    "Top1_LogNormalize_gt2_fold_rank_lollipop.pdf"
  ),
  plot = p_fold_rank,
  width = 6,
  height = 6.5,
  useDingbats = FALSE
)
