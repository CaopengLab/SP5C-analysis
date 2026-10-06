# ==================== Nmbr/Penk positivity and tdT status in neurons ====================

library(Seurat)
library(ggplot2)
library(patchwork)


# ==================== 0. Input and output ====================

input_rds <- "/public/home/wutong/Project/PJ6_Brain_Clab/Re/Step2_annotation/SP5C_integrated_manual_annotation.rds"

out_dir <- "/public/home/wutong/Project/PJ6_Brain_Clab/Re/Step2_annotation/step3_Neu_Nmbr_Penk_tdT"

dir.create(out_dir, recursive = TRUE, showWarnings = FALSE)


# ==================== 1. Load data and extract neurons ====================

integrated <- readRDS(input_rds)

neu <- subset(
  integrated,
  subset = celltype_manual %in% c(
    "Neu"
  )
)

if (ncol(neu) == 0) {
  stop("没有提取到Exc Neu或Inh Neu细胞。")
}

DefaultAssay(neu) <- "RNA"

if (!"umap" %in% Reductions(neu)) {
  stop("Neu对象中没有umap降维结果。")
}


# ==================== 2. Extract RNA counts and check genes ====================
# For Seurat v4.3; define expression positivity as raw RNA counts > 0

rna_counts <- GetAssayData(
  object = neu,
  assay = "RNA",
  slot = "counts"
)

tdt_gene <- "tdTomato"

required_genes <- c("Nmbr", "Penk")
missing_genes <- setdiff(required_genes, rownames(rna_counts))

if (length(missing_genes) > 0) {
  stop(
    "RNA assay中没有找到以下基因：",
    paste(missing_genes, collapse = ", ")
  )
}

message("使用的tdTomato基因名称：", tdt_gene)


# ==================== 3. Define tdT-, Nmbr-, and Penk-positive cells ====================

cell_names <- colnames(neu)
rna_counts <- rna_counts[, cell_names, drop = FALSE]

tdt_positive <- as.numeric(
  rna_counts[tdt_gene, ]
) > 0

nmbr_positive <- as.numeric(
  rna_counts["Nmbr", ]
) > 0

penk_positive <- as.numeric(
  rna_counts["Penk", ]
) > 0


# ==================== 4. Calculate proportions for four groups ====================

n_all_neurons <- length(cell_names)

n_nmbr_positive <- sum(nmbr_positive)
n_nmbr_negative <- n_all_neurons - n_nmbr_positive
n_nmbr_tdt_positive <- sum(nmbr_positive & tdt_positive)
n_nmbr_tdt_negative <- sum(nmbr_positive & !tdt_positive)

n_penk_positive <- sum(penk_positive)
n_penk_negative <- n_all_neurons - n_penk_positive
n_penk_tdt_positive <- sum(penk_positive & tdt_positive)
n_penk_tdt_negative <- sum(penk_positive & !tdt_positive)

percent_nmbr_in_neu <-
  100 * n_nmbr_positive / n_all_neurons

percent_tdt_in_nmbr <- if (n_nmbr_positive > 0) {
  100 * n_nmbr_tdt_positive / n_nmbr_positive
} else {
  NA_real_
}

percent_penk_in_neu <-
  100 * n_penk_positive / n_all_neurons

percent_tdt_in_penk <- if (n_penk_positive > 0) {
  100 * n_penk_tdt_positive / n_penk_positive
} else {
  NA_real_
}


# ==================== 5. Generate pie-chart data and a summary table ====================

make_fraction_data <- function(status, number) {
  total <- sum(number)
  
  percentage <- if (total > 0) {
    100 * number / total
  } else {
    rep(0, length(number))
  }
  
  df <- data.frame(
    status = factor(status, levels = status),
    number = number,
    total = total,
    percentage = percentage
  )
  
  df$label <- sprintf(
    "%s\n%.1f%%\n(n=%d)",
    df$status,
    df$percentage,
    df$number
  )
  
  df
}

nmbr_neu_fraction <- make_fraction_data(
  status = c("Nmbr-", "Nmbr+"),
  number = c(n_nmbr_negative, n_nmbr_positive)
)

nmbr_tdt_fraction <- make_fraction_data(
  status = c("tdT-", "tdT+"),
  number = c(n_nmbr_tdt_negative, n_nmbr_tdt_positive)
)

penk_neu_fraction <- make_fraction_data(
  status = c("Penk-", "Penk+"),
  number = c(n_penk_negative, n_penk_positive)
)

penk_tdt_fraction <- make_fraction_data(
  status = c("tdT-", "tdT+"),
  number = c(n_penk_tdt_negative, n_penk_tdt_positive)
)

add_analysis_name <- function(df, analysis_name) {
  data.frame(
    analysis = analysis_name,
    status = as.character(df$status),
    number = df$number,
    total = df$total,
    percentage = df$percentage
  )
}

result_summary <- rbind(
  add_analysis_name(
    nmbr_neu_fraction,
    "Nmbr status among all Neu"
  ),
  add_analysis_name(
    nmbr_tdt_fraction,
    "tdT status among Nmbr+ Neu"
  ),
  add_analysis_name(
    penk_neu_fraction,
    "Penk status among all Neu"
  ),
  add_analysis_name(
    penk_tdt_fraction,
    "tdT status among Penk+ Neu"
  )
)

print(result_summary)

write.table(
  result_summary,
  file = file.path(
    out_dir,
    "SP5C_Neu_Nmbr_Penk_tdT_summary.tsv"
  ),
  sep = "\t",
  quote = FALSE,
  row.names = FALSE
)


# ==================== 6. Prepare UMAP data ====================

neu_umap <- Embeddings(neu, reduction = "umap")

if (!all(cell_names %in% rownames(neu_umap))) {
  stop("UMAP坐标与Neu细胞名称不能完全匹配。")
}

neu_umap <- neu_umap[cell_names, , drop = FALSE]

if (ncol(neu_umap) < 2) {
  stop("umap降维结果少于2个维度。")
}

umap_df <- data.frame(
  UMAP_1 = neu_umap[, 1],
  UMAP_2 = neu_umap[, 2],
  row.names = cell_names
)

umap_df$nmbr_neu_status <- factor(
  ifelse(nmbr_positive, "Nmbr+", "Nmbr-"),
  levels = c("Nmbr-", "Nmbr+")
)

umap_df$nmbr_tdt_status <- factor(
  ifelse(
    !nmbr_positive,
    "Non-Nmbr+ neuron",
    ifelse(tdt_positive, "Nmbr+ tdT+", "Nmbr+ tdT-")
  ),
  levels = c(
    "Non-Nmbr+ neuron",
    "Nmbr+ tdT-",
    "Nmbr+ tdT+"
  )
)

umap_df$penk_neu_status <- factor(
  ifelse(penk_positive, "Penk+", "Penk-"),
  levels = c("Penk-", "Penk+")
)

umap_df$penk_tdt_status <- factor(
  ifelse(
    !penk_positive,
    "Non-Penk+ neuron",
    ifelse(tdt_positive, "Penk+ tdT+", "Penk+ tdT-")
  ),
  levels = c(
    "Non-Penk+ neuron",
    "Penk+ tdT-",
    "Penk+ tdT+"
  )
)


# Use the same square coordinate limits for all UMAP plots
x_range <- range(umap_df$UMAP_1, na.rm = TRUE)
y_range <- range(umap_df$UMAP_2, na.rm = TRUE)
x_mid <- mean(x_range)
y_mid <- mean(y_range)

square_span <- max(
  diff(x_range),
  diff(y_range)
) * 1.05

if (!is.finite(square_span) || square_span <= 0) {
  square_span <- 1
}

square_xlim <- c(
  x_mid - square_span / 2,
  x_mid + square_span / 2
)

square_ylim <- c(
  y_mid - square_span / 2,
  y_mid + square_span / 2
)


# ==================== 7. UMAP plotting function ====================

umap_theme <- theme(
  axis.line = element_blank(),
  axis.text = element_blank(),
  axis.ticks = element_blank(),
  axis.title = element_blank(),
  panel.background = element_blank(),
  panel.grid = element_blank(),
  panel.border = element_rect(
    colour = "black",
    fill = NA,
    size = 0.8
  ),
  aspect.ratio = 1,
  legend.title = element_blank(),
  legend.text = element_text(
    size = 10,
    colour = "black"
  ),
  plot.title = element_text(
    hjust = 0.5,
    size = 13,
    face = "bold"
  )
)

make_umap <- function(
    df,
    status_column,
    status_levels,
    colours,
    plot_title
) {
  plot_df <- df
  
  plot_df$plot_status <- factor(
    plot_df[[status_column]],
    levels = status_levels
  )
  
  # Plot the gray background first and highlighted cells afterward so positive cells remain visible
  plot_df <- plot_df[
    order(as.numeric(plot_df$plot_status)),
    ,
    drop = FALSE
  ]
  
  ggplot(
    plot_df,
    aes(
      x = UMAP_1,
      y = UMAP_2,
      colour = plot_status
    )
  ) +
    geom_point(
      size = 0.6,
      shape = 16
    ) +
    scale_colour_manual(
      values = colours,
      breaks = status_levels,
      drop = FALSE
    ) +
    ggtitle(plot_title) +
    coord_fixed(
      ratio = 1,
      xlim = square_xlim,
      ylim = square_ylim,
      expand = FALSE
    ) +
    umap_theme +
    guides(
      colour = guide_legend(
        override.aes = list(size = 4)
      )
    )
}


# ==================== 8. Generate four UMAP plots ====================

p_nmbr_neu_umap <- make_umap(
  df = umap_df,
  status_column = "nmbr_neu_status",
  status_levels = c("Nmbr-", "Nmbr+"),
  colours = c(
    "Nmbr-" = "#E5E5E5",
    "Nmbr+" = "#084077"
  ),
  plot_title = sprintf(
    "Nmbr+ among Neu: %.1f%% (%d/%d)",
    percent_nmbr_in_neu,
    n_nmbr_positive,
    n_all_neurons
  )
)

p_nmbr_tdt_umap <- make_umap(
  df = umap_df,
  status_column = "nmbr_tdt_status",
  status_levels = c(
    "Non-Nmbr+ neuron",
    "Nmbr+ tdT-",
    "Nmbr+ tdT+"
  ),
  colours = c(
    "Non-Nmbr+ neuron" = "#E5E5E5",
    "Nmbr+ tdT-" = "#f4d8ae",
    "Nmbr+ tdT+" = "#c31f1d"
  ),
  plot_title = sprintf(
    "tdT+ among Nmbr+ Neu: %.1f%% (%d/%d)",
    percent_tdt_in_nmbr,
    n_nmbr_tdt_positive,
    n_nmbr_positive
  )
)

p_penk_neu_umap <- make_umap(
  df = umap_df,
  status_column = "penk_neu_status",
  status_levels = c("Penk-", "Penk+"),
  colours = c(
    "Penk-" = "#E5E5E5",
    "Penk+" = "#084077"
  ),
  plot_title = sprintf(
    "Penk+ among Neu: %.1f%% (%d/%d)",
    percent_penk_in_neu,
    n_penk_positive,
    n_all_neurons
  )
)

p_penk_tdt_umap <- make_umap(
  df = umap_df,
  status_column = "penk_tdt_status",
  status_levels = c(
    "Non-Penk+ neuron",
    "Penk+ tdT-",
    "Penk+ tdT+"
  ),
  colours = c(
    "Non-Penk+ neuron" = "#E5E5E5",
    "Penk+ tdT-" = "#f4d8ae",
    "Penk+ tdT+" = "#c31f1d"
  ),
  plot_title = sprintf(
    "tdT+ among Penk+ Neu: %.1f%% (%d/%d)",
    percent_tdt_in_penk,
    n_penk_tdt_positive,
    n_penk_positive
  )
)


# ==================== 9. Pie-chart function ====================

make_pie <- function(
    df,
    fill_colours,
    text_colours,
    plot_title
) {
  plot_df <- df[df$number > 0, , drop = FALSE]
  
  if (nrow(plot_df) == 0) {
    return(
      ggplot() +
        annotate(
          "text",
          x = 0,
          y = 0,
          label = "No cells",
          size = 5
        ) +
        labs(title = plot_title) +
        theme_void() +
        theme(
          plot.title = element_text(
            hjust = 0.5,
            face = "bold",
            size = 12
          )
        )
    )
  }
  
  ggplot(
    plot_df,
    aes(
      x = "",
      y = percentage,
      fill = status
    )
  ) +
    geom_col(
      width = 1,
      colour = "black",
      size = 0.4
    ) +
    coord_polar(theta = "y") +
    geom_text(
      aes(
        label = label,
        colour = status
      ),
      position = position_stack(vjust = 0.5),
      size = 3.5,
      lineheight = 1.05
    ) +
    scale_fill_manual(
      values = fill_colours,
      guide = "none",
      drop = FALSE
    ) +
    scale_colour_manual(
      values = text_colours,
      guide = "none",
      drop = FALSE
    ) +
    labs(title = plot_title) +
    theme_void(base_size = 11) +
    theme(
      plot.title = element_text(
        hjust = 0.5,
        face = "bold",
        size = 11
      ),
      plot.margin = margin(5, 5, 5, 5),
      aspect.ratio = 1
    )
}


# ==================== 10. Generate four pie charts ====================

gene_fill_nmbr <- c(
  "Nmbr-" = "#f4d8ae",
  "Nmbr+" = "#084077"
)

gene_text_nmbr <- c(
  "Nmbr-" = "black",
  "Nmbr+" = "white"
)

gene_fill_penk <- c(
  "Penk-" = "#f4d8ae",
  "Penk+" = "#084077"
)

gene_text_penk <- c(
  "Penk-" = "black",
  "Penk+" = "white"
)

tdt_fill <- c(
  "tdT-" = "#D9D9D9",
  "tdT+" = "#c31f1d"
)

tdt_text <- c(
  "tdT-" = "black",
  "tdT+" = "white"
)

p_nmbr_neu_pie <- make_pie(
  df = nmbr_neu_fraction,
  fill_colours = gene_fill_nmbr,
  text_colours = gene_text_nmbr,
  plot_title = "Nmbr status among Neu"
)

p_nmbr_tdt_pie <- make_pie(
  df = nmbr_tdt_fraction,
  fill_colours = tdt_fill,
  text_colours = tdt_text,
  plot_title = "tdT status among Nmbr+ Neu"
)

p_penk_neu_pie <- make_pie(
  df = penk_neu_fraction,
  fill_colours = gene_fill_penk,
  text_colours = gene_text_penk,
  plot_title = "Penk status among Neu"
)

p_penk_tdt_pie <- make_pie(
  df = penk_tdt_fraction,
  fill_colours = tdt_fill,
  text_colours = tdt_text,
  plot_title = "tdT status among Penk+ Neu"
)


# ==================== 11. Combine plots ====================

p_umap_all <- (
  p_nmbr_neu_umap |
    p_nmbr_tdt_umap
) / (
  p_penk_neu_umap |
    p_penk_tdt_umap
) +
  plot_annotation(tag_levels = "A")

p_pie_all <- (
  p_nmbr_neu_pie |
    p_nmbr_tdt_pie
) / (
  p_penk_neu_pie |
    p_penk_tdt_pie
) +
  plot_annotation(tag_levels = "A")


# ==================== 12. Save four UMAP plots ====================

ggsave(
  file.path(out_dir, "SP5C_Neu_Nmbr_status_UMAP.pdf"),
  plot = p_nmbr_neu_umap,
  width = 6,
  height = 5
)

ggsave(
  file.path(out_dir, "SP5C_Nmbr_positive_tdT_status_UMAP.pdf"),
  plot = p_nmbr_tdt_umap,
  width = 6,
  height = 5
)

ggsave(
  file.path(out_dir, "SP5C_Neu_Penk_status_UMAP.pdf"),
  plot = p_penk_neu_umap,
  width = 6,
  height = 5
)

ggsave(
  file.path(out_dir, "SP5C_Penk_positive_tdT_status_UMAP.pdf"),
  plot = p_penk_tdt_umap,
  width = 6,
  height = 5
)

ggsave(
  file.path(out_dir, "SP5C_Neu_Nmbr_Penk_tdT_4UMAP.pdf"),
  plot = p_umap_all,
  width = 12,
  height = 10
)


# ==================== 13. Save four pie charts ====================

ggsave(
  file.path(out_dir, "SP5C_Neu_Nmbr_status_pie.pdf"),
  plot = p_nmbr_neu_pie,
  width = 3,
  height = 3
)

ggsave(
  file.path(out_dir, "SP5C_Nmbr_positive_tdT_status_pie.pdf"),
  plot = p_nmbr_tdt_pie,
  width = 3,
  height = 3
)

ggsave(
  file.path(out_dir, "SP5C_Neu_Penk_status_pie.pdf"),
  plot = p_penk_neu_pie,
  width = 3,
  height = 3
)

ggsave(
  file.path(out_dir, "SP5C_Penk_positive_tdT_status_pie.pdf"),
  plot = p_penk_tdt_pie,
  width = 3,
  height = 3
)

ggsave(
  file.path(out_dir, "SP5C_Neu_Nmbr_Penk_tdT_4pies.pdf"),
  plot = p_pie_all,
  width = 6,
  height = 6
)


# Display the combined plot in an interactive R session
p_umap_all
p_pie_all
