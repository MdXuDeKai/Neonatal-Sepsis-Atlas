# 早产儿败血症单细胞数据 - 图表生成脚本
# 目标：生成符合论文投稿要求的QC、细胞类型注释、比例图
# 格式要求：Times New Roman字体，TIFF和PNG格式，300 DPI

# 加载必要的库
required_packages <- c("Seurat", "ggplot2", "dplyr", "patchwork", "RColorBrewer", "viridis", "scales")

for (pkg in required_packages) {
  if (!requireNamespace(pkg, quietly = TRUE)) {
    cat("正在安装缺失的包:", pkg, "\n")
    install.packages(pkg, dependencies = TRUE)
  }
  library(pkg, character.only = TRUE)
}

# 安装presto以加速FindAllMarkers（可选但强烈推荐）
if (!requireNamespace("presto", quietly = TRUE)) {
  cat("\n检测到presto未安装。presto可以大幅加速FindAllMarkers。\n")
  cat("是否安装presto？(这将需要devtools包)\n")
  cat("正在尝试安装presto...\n")
  tryCatch({
    if (!requireNamespace("devtools", quietly = TRUE)) {
      install.packages("devtools")
    }
    devtools::install_github("immunogenomics/presto")
    cat("presto安装成功！FindAllMarkers将自动使用更快的实现。\n")
  }, error = function(e) {
    warning("presto安装失败，将使用默认的FindAllMarkers实现（可能较慢）\n")
  })
} else {
  cat("presto已安装，FindAllMarkers将使用加速版本。\n")
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

# 读取Seurat对象
cat("正在读取Seurat对象...\n")
pbmc_seurat <- readRDS("pbmc_scRNA_merged_seurat.rds")
tcell_seurat <- readRDS("tcell_EL_scRNA_merged_seurat.rds")

cat("PBMC对象：", ncol(pbmc_seurat), "个细胞\n")
cat("T细胞对象：", ncol(tcell_seurat), "个细胞\n")

# ============================================================================
# 样本名映射函数：将内部样本名映射到 GEO 原始样本名
# ============================================================================
map_to_geo_sample <- function(sample_name, original_sample = NULL) {
  # 如果提供了 original_sample，优先使用
  if (!is.null(original_sample) && !is.na(original_sample)) {
    orig <- as.character(original_sample)
    if (orig == "L163_Pre") return("L163-Pre")
    if (orig == "L163_Sepsis") return("L163-Sepsis")
    if (orig == "L131_S1_presep") return("L131-Pre")
    if (orig == "L131_S2_sep") return("L131-Sepsis")
    if (orig == "EL1") return("EL1")
    if (orig == "EL2") return("EL2")
  }
  
  # 直接匹配样本名
  sample_char <- as.character(sample_name)
  
  # HTO demultiplexing 后的样本（A5_Pre, A11_Pre, B6_Pre 等）
  if (grepl("^(A5|A11|B6)_Pre$", sample_char)) {
    return("L163-Pre")
  }
  if (grepl("^(A5|A11|B6)_Sep$", sample_char)) {
    return("L163-Sepsis")
  }
  
  # L131 样本
  if (sample_char == "L131_S1_presep") return("L131-Pre")
  if (sample_char == "L131_S2_sep") return("L131-Sepsis")
  
  # EL 样本
  if (sample_char == "EL1") return("EL1")
  if (sample_char == "EL2") return("EL2")
  
  # 原始 L163 样本（如果未做 demux）
  if (sample_char == "L163_Pre") return("L163-Pre")
  if (sample_char == "L163_Sepsis") return("L163-Sepsis")
  
  # 如果都不匹配，返回原样（但会显示警告）
  warning("未识别的样本名: ", sample_char)
  return(sample_char)
}

# ============================================================================
# 一、QC & 数据概况图（Figure 1A–B）
# ============================================================================

cat("\n=== 生成Figure 1A: QC小提琴图 ===\n")

# 合并 PBMC 和 T 细胞数据用于 QC 图
cat("合并 PBMC 和 T 细胞数据...\n")

# 检查原始样本名分布
cat("\nPBMC 原始样本名分布:\n")
print(table(pbmc_seurat$sample))
if ("original_sample" %in% colnames(pbmc_seurat@meta.data)) {
  cat("PBMC original_sample 分布:\n")
  print(table(pbmc_seurat$original_sample))
}

cat("\nT 细胞原始样本名分布:\n")
print(table(tcell_seurat$sample))
if ("original_sample" %in% colnames(tcell_seurat@meta.data)) {
  cat("T 细胞 original_sample 分布:\n")
  print(table(tcell_seurat$original_sample))
}

# 为 PBMC 数据添加 GEO 样本名
if ("original_sample" %in% colnames(pbmc_seurat@meta.data)) {
  pbmc_geo_samples <- mapply(
    map_to_geo_sample,
    pbmc_seurat$sample,
    pbmc_seurat$original_sample,
    SIMPLIFY = TRUE
  )
} else {
  pbmc_geo_samples <- sapply(pbmc_seurat$sample, map_to_geo_sample, original_sample = NULL, USE.NAMES = FALSE)
}

# 为 T 细胞数据添加 GEO 样本名
if ("original_sample" %in% colnames(tcell_seurat@meta.data)) {
  tcell_geo_samples <- mapply(
    map_to_geo_sample,
    tcell_seurat$sample,
    tcell_seurat$original_sample,
    SIMPLIFY = TRUE
  )
} else {
  tcell_geo_samples <- sapply(tcell_seurat$sample, map_to_geo_sample, original_sample = NULL, USE.NAMES = FALSE)
}

# Figure 1A: QC小提琴图（PBMC + T 细胞）
qc_features <- c("nFeature_RNA", "nCount_RNA", "percent.mt")
qc_labels <- c("Number of Genes", "Number of UMIs", "Mitochondrial %")

# 准备 PBMC 数据
qc_data_list <- list()
for (i in seq_along(qc_features)) {
  df_pbmc <- data.frame(
    Sample = pbmc_geo_samples,
    Value = pbmc_seurat[[qc_features[i]]][, 1],
    Metric = qc_labels[i]
  )
  qc_data_list[[i]] <- df_pbmc
}

# 准备 T 细胞数据
for (i in seq_along(qc_features)) {
  df_tcell <- data.frame(
    Sample = tcell_geo_samples,
    Value = tcell_seurat[[qc_features[i]]][, 1],
    Metric = qc_labels[i]
  )
  qc_data_list[[length(qc_data_list) + 1]] <- df_tcell
}

qc_data <- do.call(rbind, qc_data_list)

# 定义 GEO 样本名的顺序（按 GSE236099 的顺序）
geo_sample_order <- c("L131-Pre", "L131-Sepsis", "L163-Pre", "L163-Sepsis", "EL1", "EL2")
geo_sample_labels <- c("L131-Pre", "L131-Sepsis", "L163-Pre", "L163-Sepsis", "EL1", "EL2")

# 只保留实际存在的样本
actual_samples <- unique(qc_data$Sample)
geo_sample_order <- geo_sample_order[geo_sample_order %in% actual_samples]
geo_sample_labels <- geo_sample_labels[geo_sample_labels %in% actual_samples]

cat("\n映射后的 GEO 样本名分布:\n")
print(table(qc_data$Sample))
cat("\n将使用以下 GEO 样本名（按顺序）:\n")
cat(paste(geo_sample_labels, collapse = ", "), "\n")

# 设置因子顺序
qc_data$Sample <- factor(qc_data$Sample, 
  levels = geo_sample_order,
  labels = geo_sample_labels
)

# 绘制QC小提琴图
p_qc <- ggplot(qc_data, aes(x = Sample, y = Value, fill = Sample)) +
  geom_violin(trim = FALSE, scale = "width") +
  geom_boxplot(width = 0.1, fill = "white", outlier.size = 0.5) +
  facet_wrap(~ Metric, scales = "free_y", ncol = 3) +
  scale_fill_brewer(palette = "Set2") +
  labs(x = "Sample", y = "Value") +
  theme(
    axis.text.x = element_text(angle = 45, hjust = 1),
    strip.text = element_text(size = 11),
    legend.position = "none"
  )

# 保存Figure 1A（调整宽度以容纳更多样本）
ggsave("Figure_1A_QC_violin.tiff", p_qc, 
       width = 14, height = 4, dpi = 300, device = "tiff", compression = "lzw")
ggsave("Figure_1A_QC_violin.png", p_qc, 
       width = 14, height = 4, dpi = 300, device = "png")

cat("Figure 1A已保存\n")

# Figure 1B: 细胞数 & 基因数总结柱状图
cat("\n=== 生成Figure 1B: 数据概况柱状图 ===\n")

# 合并 PBMC 和 T 细胞数据计算统计信息
# 按 GEO 样本名分组统计

# PBMC 数据统计（按 GEO 样本名分组）
pbmc_stats_list <- list()
for (geo_sample in geo_sample_order) {
  # 找到属于这个 GEO 样本的所有细胞
  cells_in_sample <- which(pbmc_geo_samples == geo_sample)
  if (length(cells_in_sample) > 0) {
    pbmc_stats_list[[geo_sample]] <- data.frame(
      Sample = geo_sample,
      Cells = length(cells_in_sample),
      Median_Genes = median(pbmc_seurat$nFeature_RNA[cells_in_sample])
    )
  }
}

# T 细胞数据统计（按 GEO 样本名分组）
tcell_stats_list <- list()
for (geo_sample in geo_sample_order) {
  # 找到属于这个 GEO 样本的所有细胞
  cells_in_sample <- which(tcell_geo_samples == geo_sample)
  if (length(cells_in_sample) > 0) {
    tcell_stats_list[[geo_sample]] <- data.frame(
      Sample = geo_sample,
      Cells = length(cells_in_sample),
      Median_Genes = median(tcell_seurat$nFeature_RNA[cells_in_sample])
    )
  }
}

# 合并所有统计信息（按 GEO 样本名汇总）
sample_stats_list <- list()
for (geo_sample in geo_sample_order) {
  pbmc_cells <- if (geo_sample %in% names(pbmc_stats_list)) pbmc_stats_list[[geo_sample]]$Cells else 0
  tcell_cells <- if (geo_sample %in% names(tcell_stats_list)) tcell_stats_list[[geo_sample]]$Cells else 0
  total_cells <- pbmc_cells + tcell_cells
  
  if (total_cells > 0) {
    # 计算合并后的中位基因数
    pbmc_genes <- if (pbmc_cells > 0) {
      cells_idx <- which(pbmc_geo_samples == geo_sample)
      pbmc_seurat$nFeature_RNA[cells_idx]
    } else numeric(0)
    
    tcell_genes <- if (tcell_cells > 0) {
      cells_idx <- which(tcell_geo_samples == geo_sample)
      tcell_seurat$nFeature_RNA[cells_idx]
    } else numeric(0)
    
    all_genes <- c(pbmc_genes, tcell_genes)
    median_genes <- median(all_genes)
    
    sample_stats_list[[geo_sample]] <- data.frame(
      Sample = geo_sample,
      Cells = total_cells,
      Median_Genes = median_genes
    )
  }
}

sample_stats <- do.call(rbind, sample_stats_list)

# 设置因子顺序
sample_stats$Sample <- factor(sample_stats$Sample,
  levels = geo_sample_order,
  labels = geo_sample_labels
)

# 绘制细胞数柱状图
p_cells <- ggplot(sample_stats, aes(x = Sample, y = Cells, fill = Sample)) +
  geom_bar(stat = "identity", width = 0.6) +
  geom_text(aes(label = Cells), vjust = -0.5, size = 3.5, family = base_family) +
  scale_fill_brewer(palette = "Set2") +
  labs(x = "Sample", y = "Number of Cells", title = "A") +
  theme(legend.position = "none", axis.text.x = element_text(angle = 45, hjust = 1))

# 绘制中位基因数柱状图
p_genes <- ggplot(sample_stats, aes(x = Sample, y = Median_Genes, fill = Sample)) +
  geom_bar(stat = "identity", width = 0.6) +
  geom_text(aes(label = round(Median_Genes, 0)), vjust = -0.5, size = 3.5, family = base_family) +
  scale_fill_brewer(palette = "Set2") +
  labs(x = "Sample", y = "Median Genes per Cell", title = "B") +
  theme(legend.position = "none", axis.text.x = element_text(angle = 45, hjust = 1))

# 合并图表
p_fig1b <- p_cells | p_genes

# 保存Figure 1B（调整宽度以容纳更多样本）
ggsave("Figure_1B_Data_summary.tiff", p_fig1b, 
       width = 12, height = 4, dpi = 300, device = "tiff", compression = "lzw")
ggsave("Figure_1B_Data_summary.png", p_fig1b, 
       width = 12, height = 4, dpi = 300, device = "png")

cat("Figure 1B已保存\n")

# ============================================================================
# 二、PBMC 主 UMAP：细胞类型注释（Figure 1C–E）
# ============================================================================

cat("\n=== 检查细胞类型注释 ===\n")

# 检查是否已经有正确的 celltype_major（没有 Unknown 或 Unknown 很少）
has_valid_annotation <- FALSE
if ("celltype_major" %in% colnames(pbmc_seurat@meta.data)) {
  unknown_count <- sum(pbmc_seurat$celltype_major == "Unknown", na.rm = TRUE)
  unknown_pct <- 100 * unknown_count / ncol(pbmc_seurat)
  cat("当前 celltype_major 中 Unknown 数量:", unknown_count, "(", round(unknown_pct, 2), "%)\n")
  
  # 如果 Unknown 比例小于 5%，认为已经有有效的注释
  if (unknown_pct < 5) {
    has_valid_annotation <- TRUE
    cat("✓ 检测到有效的细胞类型注释，跳过重新注释步骤\n")
    cat("当前 celltype_major 分布:\n")
    print(table(pbmc_seurat$celltype_major))
    
    # 检查是否有 celltype_minor（双层级注释）
    if ("celltype_minor" %in% colnames(pbmc_seurat@meta.data)) {
      cat("\n✓ 检测到 celltype_minor（双层级注释）\n")
      cat("celltype_minor 分布:\n")
      print(table(pbmc_seurat$celltype_minor))
    }
  }
}

# 如果没有有效注释，才进行重新注释
if (!has_valid_annotation) {
  cat("\n=== 开始细胞类型注释 ===\n")
  
  # 定义经典marker基因（与 process_sepsis_scRNA.R 保持一致）
  # 根据最新专家建议，细化T细胞为T_CD4与T_CD8_NK
  marker_genes <- list(
    "T_CD4"    = c("CD3D", "CD3E", "TRAC", "IL7R", "CCR7", "TCF7", "SELL"),
    "T_CD8_NK" = c("CD3D", "CD3E", "TRAC", "NKG7", "GNLY", "GZMB", "PRF1"),
    "B"        = c("MS4A1", "CD79A", "CD74"),
    "Monocyte" = c("LYZ", "LST1", "S100A8", "S100A9", "CTSS"),
    "HSPC"     = c("CD34", "SPINK2", "KIT", "MKI67", "TOP2A", "TYMS", "STMN1"),
    "Plasma"   = c("MZB1", "XBP1", "SDC1", "JCHAIN"),
    "Platelet" = c("PPBP", "PF4", "NRGN"),
    "RBC"      = c("HBB", "HBA1", "HBA2", "SLC4A1")
  )
  
  # 检查哪些marker基因存在于数据中
  available_markers <- lapply(marker_genes, function(x) {
    x[x %in% rownames(pbmc_seurat)]
  })
  
  cat("可用的marker基因:\n")
  print(available_markers)
  
  # 运行FindAllMarkers（如果还没有运行过）
  if (!"markers" %in% names(pbmc_seurat@misc)) {
    cat("\n正在运行FindAllMarkers（这可能需要一些时间）...\n")
    Idents(pbmc_seurat) <- "seurat_clusters"
    # Seurat v5: 如果当前assay是分层（multi-layer），需要先JoinLayers，否则会出现
    # "data layers are not joined. Please run JoinLayers" 的警告且不会返回DE基因
    DefaultAssay(pbmc_seurat) <- "RNA"
    pbmc_seurat <- JoinLayers(pbmc_seurat)
    cluster_markers <- FindAllMarkers(
      pbmc_seurat,
      only.pos = TRUE,
      min.pct = 0.25,
      logfc.threshold = 0.25,
      verbose = FALSE
    )
    pbmc_seurat@misc$markers <- cluster_markers
  } else {
    cluster_markers <- pbmc_seurat@misc$markers
    cat("使用已保存的marker基因结果\n")
  }
  
  # 根据marker基因表达给cluster分配细胞类型
  cat("\n正在分配细胞类型...\n")
  pbmc_seurat$celltype_major <- "Unknown"
  
  # 计算每个cluster的平均marker表达
  Idents(pbmc_seurat) <- "seurat_clusters"
  clusters <- levels(pbmc_seurat$seurat_clusters)
  
  # 为每个cluster计算marker得分
  cluster_scores <- matrix(0, nrow = length(clusters), ncol = length(marker_genes))
  rownames(cluster_scores) <- clusters
  colnames(cluster_scores) <- names(marker_genes)
  
  # 关键修复：先JoinLayers再进行后续分析（Seurat v5）
  DefaultAssay(pbmc_seurat) <- "RNA"
  pbmc_seurat <- JoinLayers(pbmc_seurat)
  
  # 获取所有可用的marker基因
  all_markers <- unique(unlist(available_markers))
  all_markers <- all_markers[all_markers %in% rownames(pbmc_seurat)]
  
  if (length(all_markers) > 0) {
    # 使用AverageExpression计算平均表达（更稳健，避免Dropout影响）
    cat("正在使用AverageExpression计算每个cluster的marker得分...\n")
    avg_expr <- AverageExpression(
      pbmc_seurat,
      features = all_markers,
      assays = "RNA",
      layer = "data",
      group.by = "seurat_clusters"
    )$RNA
    
    # 为每个cluster计算marker得分
    for (i in seq_along(clusters)) {
      clust <- clusters[i]
      if (clust %in% colnames(avg_expr)) {
        for (j in seq_along(marker_genes)) {
          markers <- available_markers[[j]]
          markers <- markers[markers %in% rownames(avg_expr)]
          if (length(markers) > 0) {
            cluster_scores[i, j] <- mean(avg_expr[markers, clust])
          }
        }
      }
    }
    
    # 为每个cluster分配最高得分的细胞类型
    # 专家建议：降低阈值从0.1到0.01，减少Unknown
    for (i in seq_along(clusters)) {
      if (max(cluster_scores[i, ]) > 0.01) {  # 降低阈值，尽量不给 Unknown
        max_type <- names(marker_genes)[which.max(cluster_scores[i, ])]
        pbmc_seurat$celltype_major[pbmc_seurat$seurat_clusters == clusters[i]] <- max_type
      }
    }
  }
  
  # 对于仍然是Unknown的cluster，使用FindAllMarkers进一步判断
  unknown_count <- sum(pbmc_seurat$celltype_major == "Unknown")
  if (unknown_count > 0 && exists("cluster_markers") && nrow(cluster_markers) > 0) {
    cat("仍有", unknown_count, "个细胞未注释，尝试使用FindAllMarkers结果...\n")
    unknown_clusters <- unique(pbmc_seurat$seurat_clusters[pbmc_seurat$celltype_major == "Unknown"])
    
    for (clust in unknown_clusters) {
      clust_markers <- cluster_markers[cluster_markers$cluster == clust & 
                                        cluster_markers$avg_log2FC > 0.5, ]
      if (nrow(clust_markers) > 0) {
        top_genes <- clust_markers$gene[seq_len(min(5, nrow(clust_markers)))]
        # 根据top marker基因推断类型
        for (type_name in names(marker_genes)) {
          if (any(top_genes %in% marker_genes[[type_name]])) {
            pbmc_seurat$celltype_major[pbmc_seurat$seurat_clusters == clust] <- type_name
            matched_gene <- top_genes[top_genes %in% marker_genes[[type_name]]][1]
            cat("  Cluster", clust, "根据marker", matched_gene, "推断为", type_name, "\n")
            break
          }
        }
      }
    }
  }
} # 结束 if (!has_valid_annotation) 块

cat("\n细胞类型分布:\n")
print(table(pbmc_seurat$celltype_major))

# Figure 1C: UMAP按celltype_major上色
cat("\n=== 生成Figure 1C: UMAP细胞类型图 ===\n")

# 设置颜色（使用更专业的配色方案）
# 注意：根据双层级注释策略，Unknown 已归入 Monocyte，不再单独显示
celltype_colors <- c(
  "T_CD4" = "#E31A1C",      # 红色（T细胞CD4）
  "T_CD8_NK" = "#FB9A99",   # 浅红色（T细胞CD8/NK）
  "T_NK" = "#E31A1C",       # 红色（兼容旧命名）
  "B" = "#1F78B4",          # 蓝色
  "Monocyte" = "#33A02C",   # 绿色
  "DC" = "#FF7F00",         # 橙色
  "HSPC" = "#6A3D9A",       # 紫色
  "Plasma" = "#A6CEE3",     # 浅蓝色（新增：浆细胞）
  "Platelet" = "#B15928",   # 棕色
  "RBC" = "#FF6B6B",        # 粉红色（新增：红细胞）
  "Unknown" = "#CCCCCC"     # 灰色（保留以防万一）
)

# 确保所有出现的celltype都有颜色
all_celltypes <- unique(pbmc_seurat$celltype_major)
missing_colors <- setdiff(all_celltypes, names(celltype_colors))
if (length(missing_colors) > 0) {
  additional_colors <- rainbow(length(missing_colors))
  names(additional_colors) <- missing_colors
  celltype_colors <- c(celltype_colors, additional_colors)
}

p_umap_celltype <- DimPlot(pbmc_seurat, reduction = "umap", group.by = "celltype_major",
                           label = TRUE, label.size = 4, pt.size = 0.3, repel = TRUE) +
  scale_color_manual(values = celltype_colors) +
  labs(title = "Cell Type Annotation") +
  theme(
    plot.title = element_text(hjust = 0.5, size = 14, face = "bold", family = base_family),
    legend.position = "right"
  )

# 保存Figure 1C
ggsave("Figure_1C_UMAP_celltype.tiff", p_umap_celltype, 
       width = 8, height = 6, dpi = 300, device = "tiff", compression = "lzw")
ggsave("Figure_1C_UMAP_celltype.png", p_umap_celltype, 
       width = 8, height = 6, dpi = 300, device = "png")

cat("Figure 1C已保存\n")

# Figure 1D: DotPlot - 细胞类型 vs marker基因
cat("\n=== 生成Figure 1D: Marker基因表达热图 ===\n")

# 准备marker基因列表（每个类型选2-3个最特异的，HSPC使用专家推荐的marker）
# 注意：Monocyte marker 包括经典单核和炎症单核的 marker
dotplot_markers <- c(
  "CD3D", "IL7R", "CCR7",  # T_CD4
  "NKG7", "GNLY", "GZMB",  # T_CD8_NK
  "MS4A1", "CD79A",  # B
  "LYZ", "S100A8", "S100A9", "IL1R2", "FPR2",  # Monocyte（包括炎症单核marker）
  "CD34", "SPINK2", "KIT",  # HSPC（专家建议：SPINK2比CD34更敏感）
  "PPBP", "PF4",  # Platelet
  "MZB1", "JCHAIN", "XBP1",  # Plasma
  "HBB", "HBA1", "SLC4A1"  # RBC
)

# 只保留数据中存在的基因
dotplot_markers <- dotplot_markers[dotplot_markers %in% rownames(pbmc_seurat)]

# 创建DotPlot
Idents(pbmc_seurat) <- "celltype_major"
p_dotplot <- DotPlot(pbmc_seurat, features = dotplot_markers, 
                     cols = c("lightgrey", "red"), dot.scale = 8) +
  RotatedAxis() +
  labs(x = "Marker Genes", y = "Cell Type") +
  theme(
    axis.text.x = element_text(angle = 45, hjust = 1, family = base_family),
    axis.text.y = element_text(family = base_family),
    axis.title = element_text(family = base_family)
  )

# 保存Figure 1D
ggsave("Figure_1D_Marker_dotplot.tiff", p_dotplot, 
       width = 10, height = 6, dpi = 300, device = "tiff", compression = "lzw")
ggsave("Figure_1D_Marker_dotplot.png", p_dotplot, 
       width = 10, height = 6, dpi = 300, device = "png")

cat("Figure 1D已保存\n")

# Figure 1E: 堆叠柱状图 - 各cell type在Pre-sepsis vs Sepsis中的百分比
# 注意：排除EL样本（它们是分选的T细胞，会人为拉低HSPC比例）
cat("\n=== 生成Figure 1E: 细胞类型比例堆叠图 ===\n")
cat("注意：排除EL样本（分选的T细胞），只分析PBMC样本\n")

# 排除EL样本（分选的T细胞）
pbmc_for_prop <- pbmc_seurat[, !(pbmc_seurat$patient_id %in% c("EL"))]

cat("用于比例分析的细胞数（排除EL）:", ncol(pbmc_for_prop), "\n")
cat("HSPC细胞数:", sum(pbmc_for_prop$celltype_major == "HSPC", na.rm = TRUE), "\n")

# 计算比例（只使用PBMC样本）
prop_table <- prop.table(table(pbmc_for_prop$celltype_major, pbmc_for_prop$condition), margin = 2)
prop_df <- as.data.frame(prop_table)
colnames(prop_df) <- c("CellType", "Condition", "Proportion")

# 重新命名condition
prop_df$Condition <- factor(prop_df$Condition,
  levels = c("Pre-sepsis", "Sepsis"),
  labels = c("Pre-sepsis", "Sepsis")
)

# 确保HSPC在图中（即使比例很小）
if (!"HSPC" %in% prop_df$CellType) {
  cat("警告：未检测到HSPC细胞，可能需要调整聚类参数或marker阈值\n")
}

# 绘制堆叠柱状图
p_prop <- ggplot(prop_df, aes(x = Condition, y = Proportion, fill = CellType)) +
  geom_bar(stat = "identity", position = "stack", width = 0.6) +
  scale_fill_manual(values = celltype_colors) +
  labs(x = "Condition", y = "Proportion", fill = "Cell Type",
       subtitle = "PBMC samples only (EL samples excluded as sorted T cells)") +
  scale_y_continuous(labels = scales::percent_format()) +
  theme(
    axis.text = element_text(family = base_family),
    axis.title = element_text(family = base_family),
    legend.text = element_text(family = base_family),
    legend.title = element_text(family = base_family),
    plot.subtitle = element_text(size = 9, hjust = 0.5, face = "italic")
  )

# 保存Figure 1E
ggsave("Figure_1E_Celltype_proportion.tiff", p_prop, 
       width = 6, height = 6, dpi = 300, device = "tiff", compression = "lzw")
ggsave("Figure_1E_Celltype_proportion.png", p_prop, 
       width = 6, height = 6, dpi = 300, device = "png")

cat("Figure 1E已保存\n")

# ============================================================================
# 三、T 细胞 UMAP：Supplement
# ============================================================================

cat("\n=== 生成T细胞Supplement图表 ===\n")

# T细胞marker基因
tcell_markers <- c("CD4", "CD8A", "CCR7", "IL7R", "GZMB", "PRF1", "PDCD1", "LAG3")
tcell_markers <- tcell_markers[tcell_markers %in% rownames(tcell_seurat)]

# Figure S1A: T细胞UMAP按condition上色
p_tcell_umap <- DimPlot(tcell_seurat, reduction = "umap", group.by = "condition",
                        label = FALSE, pt.size = 0.5) +
  labs(title = "T Cell UMAP by Condition") +
  theme(plot.title = element_text(hjust = 0.5, size = 14, face = "bold", family = base_family))

ggsave("Figure_S1A_Tcell_UMAP.tiff", p_tcell_umap, 
       width = 8, height = 6, dpi = 300, device = "tiff", compression = "lzw")
ggsave("Figure_S1A_Tcell_UMAP.png", p_tcell_umap, 
       width = 8, height = 6, dpi = 300, device = "png")

# Figure S1B: T细胞marker基因FeaturePlot
if (length(tcell_markers) > 0) {
  p_tcell_features <- FeaturePlot(tcell_seurat, features = tcell_markers, 
                                 ncol = 4, pt.size = 0.3, combine = FALSE)
  p_tcell_features_combined <- wrap_plots(p_tcell_features, ncol = 4)
  
  ggsave("Figure_S1B_Tcell_markers.tiff", p_tcell_features_combined, 
         width = 16, height = 8, dpi = 300, device = "tiff", compression = "lzw")
  ggsave("Figure_S1B_Tcell_markers.png", p_tcell_features_combined, 
         width = 16, height = 8, dpi = 300, device = "png")
  
  cat("Figure S1B已保存\n")
}

# Figure S1C: CD4/CD8亚群比例条形图
cat("\n=== 生成Figure S1C: T细胞亚群比例图 ===\n")

# 检查是否已有tcell_subset注释（来自process_sepsis_scRNA.R）
if (!"tcell_subset" %in% colnames(tcell_seurat@meta.data) || 
    all(tcell_seurat$tcell_subset == "Unknown")) {
  
  cat("正在进行T细胞亚群注释（先聚类后注释策略）...\n")
  
  # 关键修复：先JoinLayers
  DefaultAssay(tcell_seurat) <- "RNA"
  tcell_seurat <- JoinLayers(tcell_seurat)
  
  # 初始化T细胞亚群
  tcell_seurat$tcell_subset <- "Unknown"
  
  # 使用AverageExpression计算每个cluster的CD4和CD8A平均表达
  Idents(tcell_seurat) <- "seurat_clusters"
  tcell_clusters <- levels(Idents(tcell_seurat))
  
  tcell_markers <- c("CD4", "CD8A", "CD8B")
  tcell_markers <- tcell_markers[tcell_markers %in% rownames(tcell_seurat)]
  
  if (length(tcell_markers) > 0) {
    avg_expr_tcell <- AverageExpression(
      tcell_seurat,
      features = tcell_markers,
      assays = "RNA",
      layer = "data",
      group.by = "seurat_clusters"
    )$RNA
    
    # 以Cluster为单位判断CD4/CD8
    for (clust in tcell_clusters) {
      if (clust %in% colnames(avg_expr_tcell)) {
        cd4_avg <- ifelse("CD4" %in% rownames(avg_expr_tcell), avg_expr_tcell["CD4", clust], 0)
        cd8_avg <- ifelse("CD8A" %in% rownames(avg_expr_tcell), avg_expr_tcell["CD8A", clust], 0)
        if ("CD8B" %in% rownames(avg_expr_tcell)) {
          cd8_avg <- max(cd8_avg, avg_expr_tcell["CD8B", clust])
        }
        
        if (cd4_avg > 0.5 && cd4_avg > cd8_avg * 1.5) {
          tcell_seurat$tcell_subset[tcell_seurat$seurat_clusters == clust] <- "CD4+"
        } else if (cd8_avg > 0.5 && cd8_avg > cd4_avg * 1.5) {
          tcell_seurat$tcell_subset[tcell_seurat$seurat_clusters == clust] <- "CD8+"
        } else if (cd4_avg > 0.3 && cd8_avg > 0.3) {
          tcell_seurat$tcell_subset[tcell_seurat$seurat_clusters == clust] <- "Mixed"
        } else {
          tcell_seurat$tcell_subset[tcell_seurat$seurat_clusters == clust] <- "Other"
        }
      }
    }
  }
} else {
  cat("使用已有的T细胞亚群注释\n")
}

cat("T细胞亚群分布:\n")
print(table(tcell_seurat$tcell_subset))

# 计算比例
tcell_prop <- prop.table(table(tcell_seurat$tcell_subset, tcell_seurat$condition), margin = 2)
tcell_prop_df <- as.data.frame(tcell_prop)
colnames(tcell_prop_df) <- c("Subset", "Condition", "Proportion")

# 绘制条形图
p_tcell_prop <- ggplot(tcell_prop_df, aes(x = Condition, y = Proportion, fill = Subset)) +
  geom_bar(stat = "identity", position = "dodge", width = 0.6) +
  scale_fill_brewer(palette = "Set2") +
  labs(x = "Condition", y = "Proportion", fill = "T Cell Subset") +
  scale_y_continuous(labels = scales::percent_format()) +
  theme(
    axis.text = element_text(family = base_family),
    axis.title = element_text(family = base_family),
    legend.text = element_text(family = base_family),
    legend.title = element_text(family = base_family)
  )

ggsave("Figure_S1C_Tcell_subset_proportion.tiff", p_tcell_prop, 
       width = 6, height = 6, dpi = 300, device = "tiff", compression = "lzw")
ggsave("Figure_S1C_Tcell_subset_proportion.png", p_tcell_prop, 
       width = 6, height = 6, dpi = 300, device = "png")

cat("Figure S1C已保存\n")

# ============================================================================
# 四、HSPC 特别分析（专家建议重点关注）
# ============================================================================

cat("\n=== 生成HSPC特别分析图表 ===\n")

# 检查是否有HSPC
if (sum(pbmc_seurat$celltype_major == "HSPC", na.rm = TRUE) > 0) {
  
  # Figure S1D: HSPC UMAP位置和marker表达
  cat("\n=== 生成Figure S1D: HSPC UMAP和Marker表达 ===\n")
  
  # HSPC marker基因FeaturePlot
  hspc_markers <- c("CD34", "SPINK2", "KIT", "MKI67", "TOP2A")
  hspc_markers <- hspc_markers[hspc_markers %in% rownames(pbmc_seurat)]
  
  if (length(hspc_markers) > 0) {
    # 创建HSPC标识
    pbmc_seurat$is_HSPC <- pbmc_seurat$celltype_major == "HSPC"
    
    # UMAP标注HSPC位置
    p_hspc_umap <- DimPlot(pbmc_seurat, reduction = "umap", group.by = "is_HSPC",
                           cols = c("FALSE" = "lightgrey", "TRUE" = "red"),
                           pt.size = 0.3) +
      labs(title = "HSPC Location in UMAP",
           subtitle = "Red = HSPC cells") +
      theme(plot.title = element_text(hjust = 0.5, size = 14, face = "bold", family = base_family))
    
    # HSPC marker基因FeaturePlot
    p_hspc_features <- FeaturePlot(pbmc_seurat, features = hspc_markers, 
                                  ncol = 3, pt.size = 0.3, combine = FALSE)
    p_hspc_features_combined <- wrap_plots(p_hspc_features, ncol = 3)
    
    # 合并图表
    p_hspc_combined <- p_hspc_umap / p_hspc_features_combined
    
    ggsave("Figure_S1D_HSPC_analysis.tiff", p_hspc_combined, 
           width = 12, height = 10, dpi = 300, device = "tiff", compression = "lzw")
    ggsave("Figure_S1D_HSPC_analysis.png", p_hspc_combined, 
           width = 12, height = 10, dpi = 300, device = "png")
    
    cat("Figure S1D已保存\n")
  }
  
  # Figure S1E: HSPC在不同条件中的比例（排除EL样本）
  cat("\n=== 生成Figure S1E: HSPC比例分析 ===\n")
  
  # 只使用PBMC样本（排除EL）
  pbmc_for_hspc <- pbmc_seurat[, !(pbmc_seurat$patient_id %in% c("EL"))]
  
  # 计算HSPC比例
  hspc_prop <- prop.table(table(pbmc_for_hspc$celltype_major == "HSPC", 
                                 pbmc_for_hspc$condition), margin = 2)
  hspc_prop_df <- data.frame(
    Condition = c("Pre-sepsis", "Sepsis"),
    HSPC_Proportion = as.numeric(hspc_prop["TRUE", ])
  )
  
  # 绘制HSPC比例柱状图
  p_hspc_prop <- ggplot(hspc_prop_df, aes(x = Condition, y = HSPC_Proportion, fill = Condition)) +
    geom_bar(stat = "identity", width = 0.6) +
    geom_text(aes(label = paste0(round(HSPC_Proportion * 100, 3), "%")), 
              vjust = -0.5, size = 4, family = base_family) +
    scale_fill_brewer(palette = "Set1") +
    labs(x = "Condition", y = "HSPC Proportion", 
         title = "HSPC Proportion in PBMC",
         subtitle = "PBMC samples only (EL excluded)") +
    scale_y_continuous(labels = scales::percent_format()) +
    theme(
      axis.text = element_text(family = base_family),
      axis.title = element_text(family = base_family),
      plot.title = element_text(hjust = 0.5, size = 14, face = "bold", family = base_family),
      plot.subtitle = element_text(hjust = 0.5, size = 10, face = "italic"),
      legend.position = "none"
    )
  
  ggsave("Figure_S1E_HSPC_proportion.tiff", p_hspc_prop, 
         width = 6, height = 6, dpi = 300, device = "tiff", compression = "lzw")
  ggsave("Figure_S1E_HSPC_proportion.png", p_hspc_prop, 
         width = 6, height = 6, dpi = 300, device = "png")
  
  cat("Figure S1E已保存\n")
  
} else {
  cat("未检测到HSPC细胞，跳过HSPC特别分析\n")
  cat("建议：检查CD34和SPINK2的表达，可能需要调整聚类参数\n")
}

# ============================================================================
# 五、Monocyte 亚群分析（双层级注释：Classical vs Inflammatory）
# ============================================================================

cat("\n=== 生成Monocyte亚群分析图表 ===\n")

# 检查是否有 celltype_minor 列（双层级注释）
if ("celltype_minor" %in% colnames(pbmc_seurat@meta.data)) {
  
  # 检查是否有 Monocyte 亚群
  mono_minor_types <- unique(pbmc_seurat$celltype_minor[pbmc_seurat$celltype_major == "Monocyte"])
  has_inflammatory <- "Inflammatory_Mono" %in% mono_minor_types
  has_classical <- "Classical_Mono" %in% mono_minor_types
  
  if (has_inflammatory || has_classical) {
    
    cat("检测到Monocyte亚群：", paste(mono_minor_types, collapse = ", "), "\n")
    
    # Figure S1F: Monocyte 亚群 UMAP（Classical vs Inflammatory）
    cat("\n=== 生成Figure S1F: Monocyte亚群UMAP ===\n")
    
    # 提取 Monocyte 细胞
    mono_cells <- pbmc_seurat$celltype_major == "Monocyte"
    if (sum(mono_cells) > 0) {
      mono_subset <- pbmc_seurat[, mono_cells]
      
      # UMAP 按 celltype_minor 上色
      p_mono_umap <- DimPlot(
        mono_subset,
        reduction = "umap",
        group.by = "celltype_minor",
        label = TRUE,
        repel = TRUE,
        pt.size = 0.5
      ) +
        scale_color_manual(
          values = c(
            "Classical_Mono" = "#33A02C",      # 绿色
            "Inflammatory_Mono" = "#FF6B6B"    # 粉红色
          )
        ) +
        labs(
          title = "Monocyte Subtypes",
          subtitle = "Classical vs Inflammatory Monocytes"
        ) +
        theme(
          plot.title = element_text(hjust = 0.5, size = 14, face = "bold", family = base_family),
          plot.subtitle = element_text(hjust = 0.5, size = 11, face = "italic")
        )
      
      ggsave("Figure_S1F_Monocyte_subtypes_UMAP.tiff", p_mono_umap, 
             width = 8, height = 6, dpi = 300, device = "tiff", compression = "lzw")
      ggsave("Figure_S1F_Monocyte_subtypes_UMAP.png", p_mono_umap, 
             width = 8, height = 6, dpi = 300, device = "png")
      
      cat("Figure S1F已保存\n")
      
      # Figure S1G: Inflammatory_Mono marker 基因表达
      cat("\n=== 生成Figure S1G: Inflammatory_Mono Marker表达 ===\n")
      
      # Inflammatory_Mono 的 marker 基因（根据诊断结果）
      inflammatory_markers <- c("IL1R2", "FPR2", "S100A8", "S100A9", "RETN", "S100A12")
      inflammatory_markers <- inflammatory_markers[inflammatory_markers %in% rownames(pbmc_seurat)]
      
      if (length(inflammatory_markers) > 0) {
        # FeaturePlot
        p_inflammatory_features <- FeaturePlot(
          mono_subset,
          features = inflammatory_markers,
          ncol = 3,
          pt.size = 0.3,
          combine = FALSE
        )
        p_inflammatory_combined <- wrap_plots(p_inflammatory_features, ncol = 3)
        
        ggsave("Figure_S1G_Inflammatory_Mono_markers.tiff", p_inflammatory_combined, 
               width = 12, height = 8, dpi = 300, device = "tiff", compression = "lzw")
        ggsave("Figure_S1G_Inflammatory_Mono_markers.png", p_inflammatory_combined, 
               width = 12, height = 8, dpi = 300, device = "png")
        
        cat("Figure S1G已保存\n")
      }
      
      # Figure S1H: Monocyte 亚群比例（按 condition）
      cat("\n=== 生成Figure S1H: Monocyte亚群比例 ===\n")
      
      # 排除EL样本
      mono_for_prop <- mono_subset[, !(mono_subset$patient_id %in% c("EL"))]
      
      if (ncol(mono_for_prop) > 0) {
        # 计算比例
        mono_prop <- prop.table(
          table(mono_for_prop$celltype_minor, mono_for_prop$condition),
          margin = 2
        )
        mono_prop_df <- as.data.frame(mono_prop)
        colnames(mono_prop_df) <- c("Subtype", "Condition", "Proportion")
        
        # 绘制比例图
        p_mono_prop <- ggplot(mono_prop_df, aes(x = Condition, y = Proportion, fill = Subtype)) +
          geom_bar(stat = "identity", position = "stack", width = 0.6) +
          scale_fill_manual(
            values = c(
              "Classical_Mono" = "#33A02C",
              "Inflammatory_Mono" = "#FF6B6B"
            )
          ) +
          labs(
            x = "Condition",
            y = "Proportion",
            fill = "Monocyte Subtype",
            title = "Monocyte Subtype Proportion",
            subtitle = "PBMC samples only (EL excluded)"
          ) +
          scale_y_continuous(labels = scales::percent_format()) +
          theme(
            axis.text = element_text(family = base_family),
            axis.title = element_text(family = base_family),
            legend.text = element_text(family = base_family),
            legend.title = element_text(family = base_family),
            plot.title = element_text(hjust = 0.5, size = 14, face = "bold", family = base_family),
            plot.subtitle = element_text(hjust = 0.5, size = 10, face = "italic")
          )
        
        ggsave("Figure_S1H_Monocyte_subtype_proportion.tiff", p_mono_prop, 
               width = 6, height = 6, dpi = 300, device = "tiff", compression = "lzw")
        ggsave("Figure_S1H_Monocyte_subtype_proportion.png", p_mono_prop, 
               width = 6, height = 6, dpi = 300, device = "png")
        
        cat("Figure S1H已保存\n")
      }
    }
  } else {
    cat("未检测到Monocyte亚群（Classical_Mono 或 Inflammatory_Mono），跳过Monocyte亚群分析\n")
    cat("提示：如果已运行 apply_hierarchical_annotation.R，应该会有 celltype_minor 列\n")
  }
} else {
  cat("未检测到 celltype_minor 列，跳过Monocyte亚群分析\n")
  cat("提示：请先运行 apply_hierarchical_annotation.R 创建双层级注释\n")
}

# 保存更新后的Seurat对象（包含celltype_major注释）
cat("\n正在保存更新后的Seurat对象...\n")
saveRDS(pbmc_seurat, file = "pbmc_scRNA_merged_seurat_annotated.rds")
saveRDS(tcell_seurat, file = "tcell_EL_scRNA_merged_seurat_annotated.rds")

cat("\n=== 所有图表生成完成！===\n")
cat("生成的文件:\n")
cat("主图（Figure 1）:\n")
cat("- Figure_1A_QC_violin.tiff/png\n")
cat("- Figure_1B_Data_summary.tiff/png\n")
cat("- Figure_1C_UMAP_celltype.tiff/png（使用 celltype_major，概览图）\n")
cat("- Figure_1D_Marker_dotplot.tiff/png（包含HSPC marker：CD34, SPINK2, KIT）\n")
cat("- Figure_1E_Celltype_proportion.tiff/png（排除EL样本，使用 celltype_major）\n")
cat("\n补充图（Figure S1）:\n")
cat("- Figure_S1A_Tcell_UMAP.tiff/png\n")
cat("- Figure_S1B_Tcell_markers.tiff/png\n")
cat("- Figure_S1C_Tcell_subset_proportion.tiff/png\n")
if (sum(pbmc_seurat$celltype_major == "HSPC", na.rm = TRUE) > 0) {
  cat("- Figure_S1D_HSPC_analysis.tiff/png（HSPC位置和marker表达）\n")
  cat("- Figure_S1E_HSPC_proportion.tiff/png（HSPC比例分析）\n")
}
if ("celltype_minor" %in% colnames(pbmc_seurat@meta.data)) {
  mono_minor_types <- unique(pbmc_seurat$celltype_minor[pbmc_seurat$celltype_major == "Monocyte"])
  if ("Inflammatory_Mono" %in% mono_minor_types || "Classical_Mono" %in% mono_minor_types) {
    cat("- Figure_S1F_Monocyte_subtypes_UMAP.tiff/png（Monocyte亚群：Classical vs Inflammatory）\n")
    cat("- Figure_S1G_Inflammatory_Mono_markers.tiff/png（Inflammatory_Mono marker表达）\n")
    cat("- Figure_S1H_Monocyte_subtype_proportion.tiff/png（Monocyte亚群比例）\n")
  }
}
cat("\n已保存带注释的Seurat对象:\n")
cat("- pbmc_scRNA_merged_seurat_annotated.rds\n")
cat("- tcell_EL_scRNA_merged_seurat_annotated.rds\n")
cat("\n注意：\n")
cat("- Figure 1 使用 celltype_major（宏观大类），适合概览图\n")
cat("- Figure S1F-H 使用 celltype_minor（精细亚群），展示 Monocyte 内部结构\n")
cat("- 如果未看到 Monocyte 亚群图，请先运行 apply_hierarchical_annotation.R\n")

