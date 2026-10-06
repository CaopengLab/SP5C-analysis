
#----
# Marker feature plot
library(Seurat)
library(ggplot2)
out_dir <- "/public/home/wutong/Project/PJ6_Brain_Clab/Re/Step2_annotation"
DefaultAssay(integrated) <- "RNA"

markers <- c(
  "Syt1", "Syt4",           # NEU
  "Csf1r", "C1qc",          # MG
  "Ccdc153", "Foxj1",       # EpC
  "Mog", "Mobp",            # OD
  "Slc7a10", "Agt",         # AS
  "Gpr17", "Vcan",          # OPC
  "Cldn5", "Flt1",           # EnC
  "Slc17a6", "Slc32a1",     #Exc Inh
  "Kcnj8","Rgs5",           # Pericyte
  "Col1a1","Dcn"            # Vascular and leptomeningeal cells (VLMCs), Fibroblasts

)

# Check which markers are present
markers_present <- intersect(
  markers,
  rownames(integrated[["RNA"]])
)

setdiff(markers, markers_present)

# Plot all markers
p_marker <- FeaturePlot(
  integrated,
  features = markers_present,
  reduction = "umap",
  ncol = 4,
  order = TRUE,
  min.cutoff = "q05",
  max.cutoff = "q95",
  cols = c("lightgrey", "#d11f21"),
  raster = FALSE
)


ggsave(
  "/public/home/wutong/Project/PJ6_Brain_Clab/Re/Step2_annotation/All_celltype_markers_FeaturePlot.pdf",
  plot = p_marker,
  width = 19,
  height = 19
)

#-------opc anno
DefaultAssay(integrated) <- "RNA"

markers <- c(
  # OPC lineage
  "Pdgfra", "Cspg4", "Olig2", "Sox10",
  "Ptprz1", "Vcan", "Bcan",
  
  # Proliferating OPCs
  "Mki67", "Top2a", "Cdk1", "Ccnb1", "Pcna",
  
  # COP
  "Gpr17", "Bmp4", "Nkx2-2", "Fyn",
  
  # Newly formed oligodendrocytes
  "Enpp6", "Bcas1", "Tcf7l2", "Itpr2", "Myrf",
  
  # Pericytes
  "Pdgfrb", "Rgs5", "Kcnj8", "Abcc9",
  "Des", "Notch3",
  
  # Fibroblasts
  "Col1a1", "Col3a1", "Dcn", "Lum"
)

missing_markers <- setdiff(
  markers,
  rownames(integrated[["RNA"]])
)

if (length(missing_markers) > 0L) {
  message(
    "对象中不存在的基因：",
    paste(missing_markers, collapse = ", ")
  )
}

markers <- intersect(
  markers,
  rownames(integrated[["RNA"]])
)

p_marker <- FeaturePlot(
  integrated,
  features = markers,
  reduction = "umap",
  ncol = 5,
  order = TRUE,
  min.cutoff = "q05",
  max.cutoff = "q95",
  cols = c(
    "lightgrey",
    "red"
  ),
  raster = TRUE,
  combine = TRUE
) +
  plot_annotation(
    title = "OPC lineage and vascular-cell markers"
  )

print(p_marker)

ggsave(
  filename = paste0(
    "/public/home/wutong/Project/PJ6_Brain_Clab/",
    "Re/Step2_annotation/",
    "All_celltype_markers_FeaturePlot_OPC.pdf"
  ),
  plot = p_marker,
  width = 19,
  height = ceiling(length(markers) / 5) * 4,
  units = "in",
  limitsize = FALSE
)

#---- anno
#----
# Celltype annotation
library(Seurat)
library(ggplot2)

out_dir <- paste0(
  "/public/home/wutong/Project/PJ6_Brain_Clab/",
  "Re/Step2_annotation"
)
dir.create(out_dir, recursive = TRUE, showWarnings = FALSE)


cluster_to_celltype <- c(
  "0" = "OD", # "Oligodendrocyte",
  "1" = "Exc Neu", # "Neuron",
  "2" = "Inh Neu", # "Neuron",
  "3" = "AS", 
  "4" = "EnC", #"Endothelial cell",
  "5" = "MG", # "Microglia",
  "6" = "OPC", # "Oligodendrocyte Precursor Cell"
  "7" = "Inh Neu", # "Neuron",
  "8" = "PC", # "Pericytes",
  "9" = "FB" # "Fibroblasts" # Vascular and leptomeningeal cells (VLMCs), Fibroblasts
)

# Add manual annotations based on the original seurat_clusters
integrated$celltype_manual <- unname(
  cluster_to_celltype[
    as.character(integrated$seurat_clusters)
  ]
)

# Set cell-type order
integrated$celltype_manual <- factor(
  integrated$celltype_manual,
  levels = c(
    "Exc Neu",
    "Inh Neu",
    "AS",
    "OD",
    "OPC",
    "MG",
    "EnC",
    "PC",
    "FB"
  )
)

# Check annotation results
table(
  integrated$seurat_clusters,
  integrated$celltype_manual,
  useNA = "ifany"
)

# Set cell type as the active identity
Idents(integrated) <- "celltype_manual"

# Check cell counts for each cell type
print(table(integrated$celltype_manual))

# Save the annotated object
saveRDS(
  integrated,
  file.path(out_dir, "SP5C_integrated_manual_annotation.rds")
)

#----
# Cell-type colors inspired by the NC paper
library(Seurat)
library(ggplot2)

# Keep the original theme
umap_theme <- theme(
  axis.line = element_blank(),
  axis.text.x = element_blank(),
  axis.text.y = element_blank(),
  axis.ticks = element_blank(),
  axis.title.x = element_blank(),
  axis.title.y = element_blank(),
  panel.background = element_blank(),
  panel.border = element_blank(),
  panel.grid.major = element_blank(),
  panel.grid.minor = element_blank()
)

# Use the selected colors
my_cols <- c(
  "Exc Neu" = "#C31F1D",  # red
  "Inh Neu" = "#a4cdea",
  "AS"  = "#fecddd",  # dark purple
  "OD"  = "#483D8B",  # light blue
  "OPC" = "#FADE6B",  # light yellow
  "MG"  = "#E67E22",  # orange
  "EnC" = "#1E8440",  # green
  "PC"  = "#583722",  # teal
  "FB"  = "#4280b3"   # light pink
)

# Get UMAP coordinate ranges
umap_xy <- Embeddings(integrated, reduction = "umap")

x_range <- range(umap_xy[, 1])
y_range <- range(umap_xy[, 2])

x_width  <- diff(x_range)
y_height <- diff(y_range)

# Set the origin and length of the axes in the lower-left corner
x0 <- x_range[1] + 0.06 * x_width
y0 <- y_range[1] + 0.07 * y_height

x1 <- x0 + 0.13 * x_width
y1 <- y0 + 0.13 * y_height

# Plot UMAP
p_umap <- DimPlot(
  integrated,
  reduction = "umap",
  group.by = "celltype_manual",
  cols = my_cols,
  pt.size = 0.2,
  label = FALSE,
  label.size = 6,
  repel = TRUE,
  raster = FALSE,
  shuffle = TRUE,
  seed = 1234
) +
  labs(colour = NULL) +
  umap_theme +
  coord_fixed() +
  theme(
    legend.title = element_blank(),
    legend.text = element_text(
      size = 12,
      colour = "black"
    ),
    legend.key = element_blank()
  )

p_umap

ggsave(
  file.path(out_dir, "SP5C_celltype_UMAP_custom_colors_2.pdf"),
  plot = p_umap,
  width = 6,
  height = 5
)

#

#----
# Bubble Plot
library(Seurat)
library(ggplot2)

out_dir <- "/public/home/wutong/Project/PJ6_Brain_Clab/Re/Step2_annotation"
dir.create(
  out_dir,
  recursive = TRUE,
  showWarnings = FALSE
)

# Plot expression using the RNA assay
DefaultAssay(integrated) <- "RNA"

# Cell-type order: the first type appears at the bottom of the bubble plot
celltype_levels <- c(
  "Exc Neu",
  "Inh Neu",
  "AS",
  "OD",
  "OPC",
  "MG",
  "EnC",
  "PC",
  "FB"
)

integrated$celltype_manual <- factor(
  integrated$celltype_manual,
  levels = celltype_levels
)

# Marker order
marker_order <- c(
  # Neu
  "Syt1", "Syt4", 
  # Exc InH Neu
  "Slc17a6", "Slc32a1",
  # AS
  "Slc7a10", "Agt",
  
  # OD
  "Mog", "Mobp",
  
  # OPC
  "Gpr17", "Vcan",
  
  # MG
  "Csf1r", "C1qc",
  
  # EnC
  "Cldn5", "Flt1",
  
  # PC
  "Kcnj8", "Rgs5",
  
  # FB
  "Col1a1", "Dcn"
)



# Keep only markers present in the RNA assay
markers_present <- intersect(
  marker_order,
  rownames(integrated[["RNA"]])
)

# Show missing markers
missing_markers <- setdiff(
  marker_order,
  markers_present
)

print(missing_markers)

# Use Seurat to calculate mean expression and the proportion of positive cells
dot_data <- DotPlot(
  integrated,
  features = markers_present,
  assay = "RNA",
  group.by = "celltype_manual",
  scale = TRUE,
  col.min = 0,
  col.max = 2
)$data

# Set the marker order on the x-axis
dot_data$features.plot <- factor(
  dot_data$features.plot,
  levels = markers_present
)

# Set the cell-type order on the y-axis
dot_data$id <- factor(
  as.character(dot_data$id),
  levels = celltype_levels
)

# Create a bubble plot inspired by the NC paper
p_dot <- ggplot(
  dot_data,
  aes(
    x = features.plot,
    y = id
  )
) +
  
  geom_point(
    aes(
      size = pct.exp,
      colour = avg.exp.scaled
    )
  ) +
  
  # Expression color; use the selected red for the highest expression
  scale_colour_gradient(
    low = "#FCE8E8",
    high = "#c31f1d",
    limits = c(0, 2),
    breaks = c(0, 1, 2),
    oob = scales::squish,
    name = "Expression"
  ) +
  
  # Proportion of positive cells
  scale_size_area(
    max_size = 10,
    limits = c(0, 75),
    breaks = c(0, 25, 50, 75),
    labels = c("0", "25", "50", "75"),
    oob = scales::squish,
    name = "Percentage"
  ) +
  
  labs(
    title = "Cell classes",
    x = NULL,
    y = NULL
  ) +
  
  theme_classic(base_size = 14) +
  
  theme(
    # Title
    plot.title = element_text(
      hjust = 0.5,
      size = 17,
      face = "plain"
    ),
    
    # Full black border
    panel.border = element_rect(
      colour = "black",
      fill = NA,
      size = 0.8
    ),
    
    # Remove default axis lines
    axis.line = element_blank(),
    
    # Marker names on the x-axis
    axis.text.x = element_text(
      angle = 65,
      hjust = 1,
      vjust = 1,
      size = 12,
      face = "italic",
      colour = "black"
    ),
    
    # Cell types on the y-axis
    axis.text.y = element_text(
      size = 14,
      colour = "black"
    ),
    
    # Keep only x-axis ticks
    axis.ticks.x = element_line(
      colour = "black",
      size = 0.6
    ),
    
    axis.ticks.y = element_blank(),
    
    # Legend
    legend.position = "right",
    legend.box = "vertical",
    legend.spacing.y = grid::unit(0.5, "cm"),
    
    legend.background = element_rect(
      fill = "white",
      colour = "black",
      size = 0.6
    ),
    
    legend.key = element_blank(),
    legend.text = element_text(
      size = 10,
      colour = "black"
    )
  ) +
  
  guides(
    colour = guide_colourbar(
      order = 1,
      title.position = "left",
      title.theme = element_text(
        angle = 90,
        vjust = 0.5,
        size = 12
      ),
      barheight = grid::unit(1.7, "cm"),
      barwidth = grid::unit(0.35, "cm")
    ),
    
    size = guide_legend(
      order = 2,
      title.position = "left",
      title.theme = element_text(
        angle = 90,
        vjust = 0.5,
        size = 12
      ),
      override.aes = list(
        colour = "black"
      )
    )
  )

# Display the plot
p_dot

# Save as PDF
ggsave(
  filename = file.path(
    out_dir,
    "SP5C_celltype_NC_style_BubblePlot.pdf"
  ),
  plot = p_dot,
  width = 8.5,
  height = 5.5
)
