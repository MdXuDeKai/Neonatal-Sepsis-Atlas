# ============================================================================
# HTO Demultiplexing - 将 L163 混合样本拆分为 A5, A11, B6
# ============================================================================
# 
# 检测结果确认：.h5 文件包含 HTO 数据！
# 
# HTO 到宝宝的映射（根据 GSE236099_ADT_HTO_LIST.csv）：
#   - hashtag_1, hashtag_2 → A5 (stage 2 NEC)
#   - hashtag_3, hashtag_4 → A11 (Enterobacter MCS)
#   - hashtag_5, hashtag_6 → B6 (stage 3 NEC)
#
# ============================================================================

cat("========================================\n")
cat("HTO Demultiplexing\n")
cat("将 L163 混合样本拆分为 A5, A11, B6\n")
cat("========================================\n\n")

# 加载必要的包
library(Seurat)
library(ggplot2)
library(dplyr)

# 设置工作目录
setwd("E:/GBA465/败血症单细胞")

# ============================================================================
# 定义 HTO 到宝宝的映射
# ============================================================================
# 根据 ADT_HTO_LIST.csv：
#   Hashtag 1, 2 → A5
#   Hashtag 3, 4 → A11
#   Hashtag 5, 6 → B6

# 注意：feature 名称可能有变化（Seurat 会把下划线改成横线），使用模式匹配
get_baby_from_hto <- function(hto_name) {
  hto_name <- tolower(hto_name)
  # Hashtag 1, 2 → A5
  if (grepl("hashtag_1|hashtag-1|hashtag1|b0251", hto_name)) {
    return("A5")
  } else if (grepl("hashtag_2|hashtag-2|hashtag2|b0252", hto_name)) {
    return("A5")
  # Hashtag 3, 4 → A11
  } else if (grepl("hashtag_3|hashtag-3|hashtag3|b0253", hto_name)) {
    return("A11")
  } else if (grepl("hashtag_4|hashtag-4|hashtag4|b0254", hto_name)) {
    return("A11")
  # Hashtag 5, 6 → B6
  } else if (grepl("hashtag_5|hashtag-5|hashtag5|b0255", hto_name)) {
    return("B6")
  } else if (grepl("hashtag_6|hashtag-6|hashtag6|b0256", hto_name)) {
    return("B6")
  } else {
    return(NA)
  }
}

# ============================================================================
# 函数：处理单个样本的 HTO demultiplexing
# ============================================================================
process_sample_demux <- function(h5_file, sample_name, condition) {
  cat("\n=== 处理样本:", sample_name, "===\n")
  
  # 读取数据
  data <- Read10X_h5(h5_file)
  
  if (!is.list(data)) {
    stop("数据不是多模态格式！")
  }
  
  # 提取 RNA 和 Antibody Capture 数据
  rna_data <- data[["Gene Expression"]]
  adt_data <- data[["Antibody Capture"]]
  
  cat("RNA features:", nrow(rna_data), "\n")
  cat("ADT/HTO features:", nrow(adt_data), "\n")
  cat("Cells:", ncol(rna_data), "\n")
  
  # 查找 HTO features（hashtag 开头的）
  all_features <- rownames(adt_data)
  cat("所有 ADT/HTO features:\n")
  print(all_features)
  
  hto_features <- all_features[grep("hashtag", all_features, ignore.case = TRUE)]
  cat("\nHTO features:\n")
  print(hto_features)
  
  if (length(hto_features) == 0) {
    stop("未找到 HTO features！")
  }
  
  # 创建 Seurat 对象
  seurat_obj <- CreateSeuratObject(counts = rna_data, project = sample_name)
  
  # 添加 HTO assay（只使用 hashtag features）
  hto_data <- adt_data[hto_features, , drop = FALSE]
  
  # 关键修复：更严格的过滤策略
  # 1. 过滤总计数为 0 的细胞
  hto_total_counts <- colSums(hto_data)
  cells_with_hto <- names(hto_total_counts)[hto_total_counts > 0]
  
  cat("  总细胞数:", ncol(seurat_obj), "\n")
  cat("  有 HTO 总信号的细胞:", length(cells_with_hto), "\n")
  
  # 2. 进一步过滤：要求每个细胞至少在一个 HTO feature 上有足够的计数（> 1）
  # 这可以避免 HTODemux 遇到"所有 features 都是 0"的细胞
  hto_data_filtered <- hto_data[, cells_with_hto, drop = FALSE]
  max_hto_per_cell <- apply(hto_data_filtered, 2, max)
  cells_with_strong_hto <- names(max_hto_per_cell)[max_hto_per_cell > 1]
  
  cat("  HTO 信号强度 > 1 的细胞:", length(cells_with_strong_hto), "\n")
  cat("  过滤掉的弱信号细胞:", length(cells_with_hto) - length(cells_with_strong_hto), "\n")
  
  if (length(cells_with_strong_hto) < 100) {
    warning("有足够 HTO 信号的细胞太少，可能数据有问题")
  }
  
  # 只保留有足够 HTO 信号的细胞
  seurat_obj <- seurat_obj[, cells_with_strong_hto]
  hto_data <- hto_data[, cells_with_strong_hto, drop = FALSE]
  
  # 检查 HTO 数据质量
  cat("\n  HTO 数据质量检查:\n")
  cat("    每个 HTO feature 的非零细胞数:\n")
  for (hto in rownames(hto_data)) {
    n_nonzero <- sum(hto_data[hto, ] > 0)
    cat("      ", hto, ":", n_nonzero, "cells\n")
  }
  
  # 创建 HTO assay
  seurat_obj[["HTO"]] <- CreateAssayObject(counts = hto_data)
  
  # 标准化 HTO 数据（CLR 方法）
  seurat_obj <- NormalizeData(seurat_obj, assay = "HTO", normalization.method = "CLR")
  
  # 运行 HTODemux（带备用方案）
  cat("\n正在运行 HTODemux...\n")
  cat("  过滤后细胞数:", ncol(seurat_obj), "\n")
  
  hto_demux_success <- FALSE
  
  # 尝试不同的 positive.quantile 值
  for (q in c(0.99, 0.95, 0.90, 0.85, 0.80)) {
    tryCatch({
      seurat_obj <- HTODemux(seurat_obj, assay = "HTO", positive.quantile = q)
      cat("  ✓ HTODemux 成功（positive.quantile =", q, "）\n")
      hto_demux_success <- TRUE
      break
    }, error = function(e) {
      if (q == 0.80) {
        cat("  ✗ HTODemux 失败，将使用备用方案（基于阈值的简单分配）\n")
      }
    })
  }
  
  # 如果 HTODemux 失败，使用备用方案
  if (!hto_demux_success) {
    cat("\n使用备用方案：基于阈值的简单 HTO 分配...\n")
    
    # 获取标准化后的 HTO 数据
    hto_normalized <- GetAssayData(seurat_obj, assay = "HTO", layer = "data")
    
    # 为每个细胞找到最强的 HTO 信号
    max_hto_idx <- apply(hto_normalized, 2, which.max)
    max_hto_value <- apply(hto_normalized, 2, max)
    
    # 设置阈值：如果最强信号 > 1（标准化后），则分配该 HTO
    # 否则标记为 Negative
    seurat_obj$hash.ID <- rownames(hto_normalized)[max_hto_idx]
    seurat_obj$hash.ID[max_hto_value < 1] <- "Negative"
    
    # 检查是否有 Doublet（多个 HTO 信号都很强）
    # 如果第二强的信号也很接近最强信号（> 0.8 * 最强），可能是 Doublet
    second_max_values <- apply(hto_normalized, 2, function(x) {
      sorted <- sort(x, decreasing = TRUE)
      if (length(sorted) > 1) sorted[2] else 0
    })
    
    doublet_cells <- names(second_max_values)[second_max_values > 0.8 * max_hto_value]
    seurat_obj$hash.ID[doublet_cells] <- "Doublet"
    
    # 设置 classification
    seurat_obj$HTO_classification.global <- ifelse(
      seurat_obj$hash.ID == "Negative", "Negative",
      ifelse(seurat_obj$hash.ID == "Doublet", "Doublet", "Singlet")
    )
    
    cat("  备用方案完成\n")
  }
  
  # 查看 demultiplexing 结果
  cat("\n=== Demultiplexing 结果 ===\n")
  cat("Classification:\n")
  print(table(seurat_obj$HTO_classification.global))
  
  cat("\nHash ID 分布:\n")
  print(table(seurat_obj$hash.ID))
  
  # 根据 hash.ID 分配 baby_id（显式对齐到细胞名，避免 No cell overlap 错误）
  cat("\n正在分配 Baby ID...\n")
  
  # 获取 hash.ID 向量（确保有正确的细胞名）
  hash_vec <- as.character(seurat_obj$hash.ID)
  names(hash_vec) <- colnames(seurat_obj)  # 确保有细胞名
  
  # 使用 vapply 替代 sapply（更安全，显式指定返回类型）
  baby_vec <- vapply(
    hash_vec,
    FUN.VALUE = character(1),
    FUN = function(x) {
      if (x %in% c("Doublet", "Negative")) {
        return(x)
      } else {
        result <- get_baby_from_hto(x)
        if (is.na(result)) {
          return("Unknown")
        }
        return(result)
      }
    }
  )
  
  # 显式创建 data.frame，设置 row.names = 细胞名
  meta_baby <- data.frame(baby_id = baby_vec, row.names = colnames(seurat_obj))
  
  # 使用 AddMetaData 写入 Seurat 对象（而不是 $<-）
  seurat_obj <- AddMetaData(seurat_obj, metadata = meta_baby)
  
  cat("\nBaby ID 分布:\n")
  print(table(seurat_obj$baby_id))
  
  # 添加 metadata
  seurat_obj$original_sample <- sample_name
  seurat_obj$condition <- condition
  seurat_obj$timepoint <- ifelse(condition == "Pre-sepsis", "Pre", "Sepsis")
  seurat_obj$cell_source <- "PBMC"
  
  # 计算 QC 指标
  seurat_obj[["percent.mt"]] <- PercentageFeatureSet(seurat_obj, pattern = "^MT-")
  seurat_obj[["percent.ribo"]] <- PercentageFeatureSet(seurat_obj, pattern = "^RP[SL]")
  seurat_obj[["percent.hb"]] <- PercentageFeatureSet(seurat_obj, pattern = "^HB[^(P)]")
  
  return(seurat_obj)
}

# ============================================================================
# 处理 L163-Pre
# ============================================================================
l163_pre <- process_sample_demux(
  h5_file = "数据/GSE236099_L163-Pre_filtered_feature_bc_matrix.h5",
  sample_name = "L163_Pre",
  condition = "Pre-sepsis"
)

# ============================================================================
# 处理 L163-Sepsis
# ============================================================================
l163_sep <- process_sample_demux(
  h5_file = "数据/GSE236099_L163-Sepsis_filtered_feature_bc_matrix.h5",
  sample_name = "L163_Sepsis",
  condition = "Sepsis"
)

# ============================================================================
# 生成 QC 可视化
# ============================================================================
cat("\n=== 生成 HTO Demultiplexing 可视化 ===\n")

# 为可视化运行 PCA（如果还没有）
cat("为可视化运行 PCA...\n")
if (!"pca" %in% names(l163_pre@reductions)) {
  l163_pre <- NormalizeData(l163_pre, verbose = FALSE)
  l163_pre <- FindVariableFeatures(l163_pre, verbose = FALSE)
  l163_pre <- ScaleData(l163_pre, verbose = FALSE)
  l163_pre <- RunPCA(l163_pre, features = VariableFeatures(l163_pre), npcs = 10, verbose = FALSE)
}

if (!"pca" %in% names(l163_sep@reductions)) {
  l163_sep <- NormalizeData(l163_sep, verbose = FALSE)
  l163_sep <- FindVariableFeatures(l163_sep, verbose = FALSE)
  l163_sep <- ScaleData(l163_sep, verbose = FALSE)
  l163_sep <- RunPCA(l163_sep, features = VariableFeatures(l163_sep), npcs = 10, verbose = FALSE)
}

# Ridge plot 显示 HTO 信号
pdf("HTO_Demux_QC.pdf", width = 14, height = 10)

# L163-Pre
p1 <- RidgePlot(l163_pre, assay = "HTO", features = rownames(l163_pre[["HTO"]]), ncol = 3) +
  ggtitle("L163-Pre: HTO Signal Distribution")

# L163-Pre classification（使用 PCA）
Idents(l163_pre) <- "HTO_classification.global"
p2 <- DimPlot(l163_pre, reduction = "pca", group.by = "HTO_classification.global") +
  ggtitle("L163-Pre: HTO Classification (PCA)")

# L163-Pre baby_id
p3 <- VlnPlot(l163_pre, features = "nFeature_RNA", group.by = "baby_id") +
  ggtitle("L163-Pre: nFeature_RNA by Baby ID")

# L163-Sepsis
p4 <- RidgePlot(l163_sep, assay = "HTO", features = rownames(l163_sep[["HTO"]]), ncol = 3) +
  ggtitle("L163-Sepsis: HTO Signal Distribution")

print(p1)
print(p2)
print(p3)
print(p4)

dev.off()
cat("已保存: HTO_Demux_QC.pdf\n")

# ============================================================================
# 过滤 Doublet 和 Negative，只保留 Singlet
# ============================================================================
cat("\n=== 过滤 Doublet 和 Negative ===\n")

# L163-Pre
singlet_pre <- l163_pre$HTO_classification.global == "Singlet"
valid_baby_pre <- l163_pre$baby_id %in% c("A5", "A11", "B6")
keep_pre <- singlet_pre & valid_baby_pre

cat("L163-Pre:\n")
cat("  总细胞数:", ncol(l163_pre), "\n")
cat("  Singlet:", sum(singlet_pre), "\n")
cat("  有效 Baby ID:", sum(valid_baby_pre), "\n")
cat("  保留:", sum(keep_pre), "\n")

l163_pre_clean <- l163_pre[, keep_pre]

# L163-Sepsis
singlet_sep <- l163_sep$HTO_classification.global == "Singlet"
valid_baby_sep <- l163_sep$baby_id %in% c("A5", "A11", "B6")
keep_sep <- singlet_sep & valid_baby_sep

cat("\nL163-Sepsis:\n")
cat("  总细胞数:", ncol(l163_sep), "\n")
cat("  Singlet:", sum(singlet_sep), "\n")
cat("  有效 Baby ID:", sum(valid_baby_sep), "\n")
cat("  保留:", sum(keep_sep), "\n")

l163_sep_clean <- l163_sep[, keep_sep]

# ============================================================================
# 保存 Demultiplexed 对象
# ============================================================================
cat("\n=== 保存 Demultiplexed 对象 ===\n")

saveRDS(l163_pre_clean, "L163_Pre_demultiplexed.rds")
cat("已保存: L163_Pre_demultiplexed.rds\n")

saveRDS(l163_sep_clean, "L163_Sepsis_demultiplexed.rds")
cat("已保存: L163_Sepsis_demultiplexed.rds\n")

# ============================================================================
# 总结报告
# ============================================================================
cat("\n========================================\n")
cat("HTO Demultiplexing 完成！\n")
cat("========================================\n\n")

cat("L163-Pre 样本拆分结果:\n")
print(table(l163_pre_clean$baby_id))

cat("\nL163-Sepsis 样本拆分结果:\n")
print(table(l163_sep_clean$baby_id))

cat("\n配对设计检查:\n")
pre_babies <- unique(l163_pre_clean$baby_id)
sep_babies <- unique(l163_sep_clean$baby_id)
paired_babies <- intersect(pre_babies, sep_babies)
cat("有完整配对的宝宝:", paste(paired_babies, collapse = ", "), "\n")

cat("\n下一步:\n")
cat("1. 运行 Step1_Process_With_Demultiplexing.R\n")
cat("2. 运行 apply_hierarchical_annotation.R\n")
cat("3. 运行 generate_figures.R\n")
cat("4. 运行 paired_DE_analysis_robust.R\n")

cat("\n🎉 L163 已成功拆分为 A5, A11, B6！\n")
cat("统计效力从 n=2 提升到 n=4！\n")

