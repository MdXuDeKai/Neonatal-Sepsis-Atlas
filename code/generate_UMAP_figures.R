# ============================================================================
# UMAP 图表生成脚本（单独保存，高分辨率）
# 目标：按照 generate_figures.R 的绘图风格，单独生成所有 UMAP 图
# 格式要求：Times New Roman字体，TIFF和PNG格式，300 DPI
# 保存路径：E:\GBA465\败血症单细胞\figure\figure1\umap图
# ============================================================================

cat("========================================\n")
cat("UMAP 图表生成脚本（单独保存）\n")
cat("按照 generate_figures.R 的绘图风格\n")
cat("========================================\n\n")

# 加载必要的库
required_packages <- c("Seurat", "ggplot2", "dplyr", "RColorBrewer", "viridis")

for (pkg in required_packages) {
  if (!requireNamespace(pkg, quietly = TRUE)) {
    stop(paste("错误：缺少必要的包:", pkg, "\n请先安装:", pkg))
  }
  library(pkg, character.only = TRUE)
}

# 设置工作目录
setwd("E:/GBA465/败血症单细胞")

# 设置图形参数（Times New Roman字体，高分辨率）
# Windows系统使用Times New Roman
if (.Platform$OS.type == "windows") {
  tryCatch({
    windowsFonts(Times = windowsFont("Times New Roman"))
    base_family <- "Times"
  }, error = function(e) {
    base_family <- "serif"  # 如果Times New Roman不可用，使用serif作为备选
    warning("Times New Roman字体不可用，使用serif字体")
  })
} else {
  # Mac/Linux系统
  base_family <- "serif"  # 使用serif作为Times New Roman的近似
}

# 设置默认主题
theme_set(theme_classic(base_family = base_family, base_size = 12))
theme_update(
  text = element_text(family = base_family, size = 12),
  axis.text = element_text(family = base_family, size = 10),
  axis.title = element_text(family = base_family, size = 12),
  legend.text = element_text(family = base_family, size = 10),
  legend.title = element_text(family = base_family, size = 12),
  plot.title = element_text(family = base_family, size = 14, face = "bold")
)

# 创建输出目录
output_dir <- "figure/figure1/umap图"
if (!dir.exists(output_dir)) {
  dir.create(output_dir, recursive = TRUE)
  cat("已创建输出目录:", output_dir, "\n")
}

# 读取Seurat对象
cat("\n正在读取Seurat对象...\n")
if (!file.exists("pbmc_scRNA_merged_seurat.rds")) {
  stop("错误：找不到 pbmc_scRNA_merged_seurat.rds 文件\n",
       "请先运行 Step1_Process_With_Demultiplexing.R 和 apply_hierarchical_annotation.R")
}

pbmc_seurat <- readRDS("pbmc_scRNA_merged_seurat.rds")
cat("成功读取PBMC对象\n")
cat("总细胞数:", ncol(pbmc_seurat), "\n")

# 检查是否有 celltype_minor（来自 apply_hierarchical_annotation.R）
has_minor <- "celltype_minor" %in% colnames(pbmc_seurat@meta.data)
if (has_minor) {
  cat("检测到 celltype_minor，将生成精细亚群图\n")
} else {
  cat("警告：未检测到 celltype_minor，将跳过精细亚群相关图\n")
}

# 确保数据已合并（Seurat v5）
DefaultAssay(pbmc_seurat) <- "RNA"
tryCatch({
  pbmc_seurat <- JoinLayers(pbmc_seurat)
  cat("成功合并所有数据层\n")
}, error = function(e) {
  cat("JoinLayers不适用或已合并，继续...\n")
})

# ============================================================================
# 定义颜色方案（与 generate_figures.R 保持一致）
# ============================================================================
cat("\n=== 设置颜色方案 ===\n")

# celltype_major 颜色
celltype_colors <- c(
  "T_CD4"    = "#1F77B4",
  "T_CD8_NK" = "#FF7F0E",
  "B"        = "#2CA02C",
  "Monocyte" = "#D62728",
  "HSPC"     = "#9467BD",
  "Plasma"   = "#8C564B",
  "Platelet" = "#E377C2",
  "RBC"      = "#7F7F7F"
)

# 如果存在 celltype_minor，添加精细亚群颜色
if (has_minor) {
  minor_colors <- c(
    "T_CD4"            = "#1F77B4",
    "T_CD8_NK"        = "#FF7F0E",
    "B"                = "#2CA02C",
    "Classical_Mono"   = "#D62728",
    "Inflammatory_Mono" = "#FF6B6B",
    "HSPC"             = "#9467BD",
    "Plasma"           = "#8C564B",
    "Platelet"         = "#E377C2",
    "RBC"              = "#7F7F7F"
  )
}

# ============================================================================
# 第一部分：来自 Step1_Process_With_Demultiplexing.R 的 UMAP 图
# ============================================================================
cat("\n========================================\n")
cat("第一部分：基础 UMAP 图（来自 Step1）\n")
cat("========================================\n")

# 1. UMAP by Patient ID
cat("\n1. 生成 UMAP by Patient ID...\n")
p1 <- DimPlot(
  pbmc_seurat,
  reduction = "umap",
  group.by = "patient_id",
  label = TRUE,
  label.size = 4,
  pt.size = 0.3,
  repel = TRUE
) +
  labs(title = "UMAP by Patient ID") +
  theme(
    plot.title = element_text(hjust = 0.5, size = 14, face = "bold", family = base_family),
    legend.position = "right"
  )

ggsave(file.path(output_dir, "UMAP_by_Patient_ID.tiff"), p1,
       width = 8, height = 6, dpi = 300, device = "tiff", compression = "lzw")
ggsave(file.path(output_dir, "UMAP_by_Patient_ID.png"), p1,
       width = 8, height = 6, dpi = 300, device = "png")
cat("  已保存: UMAP_by_Patient_ID.tiff/png\n")

# 2. UMAP by Condition
cat("\n2. 生成 UMAP by Condition...\n")
p2 <- DimPlot(
  pbmc_seurat,
  reduction = "umap",
  group.by = "condition",
  label = TRUE,
  label.size = 4,
  pt.size = 0.3,
  repel = TRUE,
  cols = c("Pre-sepsis" = "#4DAF4A", "Sepsis" = "#E41A1C")
) +
  labs(title = "UMAP by Condition") +
  theme(
    plot.title = element_text(hjust = 0.5, size = 14, face = "bold", family = base_family),
    legend.position = "right"
  )

ggsave(file.path(output_dir, "UMAP_by_Condition.tiff"), p2,
       width = 8, height = 6, dpi = 300, device = "tiff", compression = "lzw")
ggsave(file.path(output_dir, "UMAP_by_Condition.png"), p2,
       width = 8, height = 6, dpi = 300, device = "png")
cat("  已保存: UMAP_by_Condition.tiff/png\n")

# 3. UMAP by Cell Type (celltype_major)
cat("\n3. 生成 UMAP by Cell Type (Major)...\n")
p3 <- DimPlot(
  pbmc_seurat,
  reduction = "umap",
  group.by = "celltype_major",
  label = TRUE,
  label.size = 4,
  pt.size = 0.3,
  repel = TRUE
) +
  scale_color_manual(values = celltype_colors) +
  labs(title = "UMAP by Cell Type (Major)") +
  theme(
    plot.title = element_text(hjust = 0.5, size = 14, face = "bold", family = base_family),
    legend.position = "right"
  )

ggsave(file.path(output_dir, "UMAP_by_CellType_Major.tiff"), p3,
       width = 8, height = 6, dpi = 300, device = "tiff", compression = "lzw")
ggsave(file.path(output_dir, "UMAP_by_CellType_Major.png"), p3,
       width = 8, height = 6, dpi = 300, device = "png")
cat("  已保存: UMAP_by_CellType_Major.tiff/png\n")

# 4. UMAP by Cluster
cat("\n4. 生成 UMAP by Cluster...\n")
p4 <- DimPlot(
  pbmc_seurat,
  reduction = "umap",
  group.by = "seurat_clusters",
  label = TRUE,
  label.size = 4,
  pt.size = 0.3,
  repel = TRUE
) +
  labs(title = "UMAP by Cluster") +
  theme(
    plot.title = element_text(hjust = 0.5, size = 14, face = "bold", family = base_family),
    legend.position = "right"
  )

ggsave(file.path(output_dir, "UMAP_by_Cluster.tiff"), p4,
       width = 8, height = 6, dpi = 300, device = "tiff", compression = "lzw")
ggsave(file.path(output_dir, "UMAP_by_Cluster.png"), p4,
       width = 8, height = 6, dpi = 300, device = "png")
cat("  已保存: UMAP_by_Cluster.tiff/png\n")

# 5. UMAP by Sample
cat("\n5. 生成 UMAP by Sample...\n")
p5 <- DimPlot(
  pbmc_seurat,
  reduction = "umap",
  group.by = "sample",
  label = TRUE,
  label.size = 4,
  pt.size = 0.3,
  repel = TRUE
) +
  labs(title = "UMAP by Sample") +
  theme(
    plot.title = element_text(hjust = 0.5, size = 14, face = "bold", family = base_family),
    legend.position = "right"
  )

ggsave(file.path(output_dir, "UMAP_by_Sample.tiff"), p5,
       width = 8, height = 6, dpi = 300, device = "tiff", compression = "lzw")
ggsave(file.path(output_dir, "UMAP_by_Sample.png"), p5,
       width = 8, height = 6, dpi = 300, device = "png")
cat("  已保存: UMAP_by_Sample.tiff/png\n")

# ============================================================================
# 第二部分：来自 apply_hierarchical_annotation.R 的 UMAP 图
# ============================================================================
if (has_minor) {
  cat("\n========================================\n")
  cat("第二部分：精细亚群 UMAP 图（来自 Hierarchical Annotation）\n")
  cat("========================================\n")
  
  # 6. UMAP by Cell Type Minor
  cat("\n6. 生成 UMAP by Cell Type (Minor)...\n")
  p6 <- DimPlot(
    pbmc_seurat,
    reduction = "umap",
    group.by = "celltype_minor",
    label = TRUE,
    label.size = 4,
    pt.size = 0.3,
    repel = TRUE
  ) +
    scale_color_manual(values = minor_colors) +
    labs(title = "UMAP by Cell Type (Minor)") +
    theme(
      plot.title = element_text(hjust = 0.5, size = 14, face = "bold", family = base_family),
      legend.position = "right"
    )
  
  ggsave(file.path(output_dir, "UMAP_by_CellType_Minor.tiff"), p6,
         width = 8, height = 6, dpi = 300, device = "tiff", compression = "lzw")
  ggsave(file.path(output_dir, "UMAP_by_CellType_Minor.png"), p6,
         width = 8, height = 6, dpi = 300, device = "png")
  cat("  已保存: UMAP_by_CellType_Minor.tiff/png\n")
  
  # 7. Monocyte Subtypes (Classical vs Inflammatory)
  cat("\n7. 生成 Monocyte Subtypes UMAP...\n")
  mono_subset <- pbmc_seurat[, pbmc_seurat$celltype_major == "Monocyte"]
  
  if (ncol(mono_subset) > 0) {
    p7 <- DimPlot(
      mono_subset,
      reduction = "umap",
      group.by = "celltype_minor",
      label = TRUE,
      label.size = 4,
      pt.size = 0.5,
      repel = TRUE,
      cols = c("Classical_Mono" = "#D62728", "Inflammatory_Mono" = "#FF6B6B")
    ) +
      labs(title = "Monocyte Subtypes (Classical vs Inflammatory)") +
      theme(
        plot.title = element_text(hjust = 0.5, size = 14, face = "bold", family = base_family),
        legend.position = "right"
      )
    
    ggsave(file.path(output_dir, "UMAP_Monocyte_Subtypes.tiff"), p7,
           width = 8, height = 6, dpi = 300, device = "tiff", compression = "lzw")
    ggsave(file.path(output_dir, "UMAP_Monocyte_Subtypes.png"), p7,
           width = 8, height = 6, dpi = 300, device = "png")
    cat("  已保存: UMAP_Monocyte_Subtypes.tiff/png\n")
  } else {
    cat("  警告：未找到 Monocyte 细胞，跳过此图\n")
  }
  
  # 8. Inflammatory_Mono Cells (Highlighted)
  cat("\n8. 生成 Inflammatory_Mono Highlighted UMAP...\n")
  inf_cells <- WhichCells(pbmc_seurat, expression = celltype_minor == "Inflammatory_Mono")
  
  if (length(inf_cells) > 0) {
    p8 <- DimPlot(
      pbmc_seurat,
      reduction = "umap",
      cells.highlight = inf_cells,
      cols.highlight = "red",
      cols = "lightgray",
      pt.size = 0.3,
      sizes.highlight = 0.5
    ) +
      labs(title = "Inflammatory_Mono Cells (Highlighted)") +
      theme(
        plot.title = element_text(hjust = 0.5, size = 14, face = "bold", family = base_family),
        legend.position = "right"
      )
    
    ggsave(file.path(output_dir, "UMAP_Inflammatory_Mono_Highlighted.tiff"), p8,
           width = 8, height = 6, dpi = 300, device = "tiff", compression = "lzw")
    ggsave(file.path(output_dir, "UMAP_Inflammatory_Mono_Highlighted.png"), p8,
           width = 8, height = 6, dpi = 300, device = "png")
    cat("  已保存: UMAP_Inflammatory_Mono_Highlighted.tiff/png\n")
  } else {
    cat("  警告：未找到 Inflammatory_Mono 细胞，跳过此图\n")
  }
  
} else {
  cat("\n跳过精细亚群相关图（未检测到 celltype_minor）\n")
}

# ============================================================================
# 总结报告
# ============================================================================
cat("\n========================================\n")
cat("UMAP 图表生成完成！\n")
cat("========================================\n\n")

cat("输出目录:", output_dir, "\n")
cat("生成的文件:\n")

# 列出所有生成的文件
base_plots <- c(
  "UMAP_by_Patient_ID",
  "UMAP_by_Condition",
  "UMAP_by_CellType_Major",
  "UMAP_by_Cluster",
  "UMAP_by_Sample"
)

if (has_minor) {
  minor_plots <- c(
    "UMAP_by_CellType_Minor",
    "UMAP_Monocyte_Subtypes",
    "UMAP_Inflammatory_Mono_Highlighted"
  )
  all_plots <- c(base_plots, minor_plots)
} else {
  all_plots <- base_plots
}

for (plot_name in all_plots) {
  cat("  -", plot_name, ".tiff\n")
  cat("  -", plot_name, ".png\n")
}

cat("\n所有图片已按照 generate_figures.R 的绘图风格生成：\n")
cat("  ✓ Times New Roman 字体\n")
cat("  ✓ theme_classic 主题\n")
cat("  ✓ 300 DPI 高分辨率\n")
cat("  ✓ TIFF 和 PNG 双格式\n")
cat("  ✓ 每张图单独保存，清晰可读\n\n")

cat("脚本执行完成！\n")

