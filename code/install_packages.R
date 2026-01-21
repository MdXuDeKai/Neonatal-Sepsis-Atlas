# 安装处理单细胞数据所需的R包

# 检查并安装BiocManager（如果需要）
if (!requireNamespace("BiocManager", quietly = TRUE)) {
  install.packages("BiocManager")
}

# 安装CRAN包
cran_packages <- c("Seurat", "hdf5r", "dplyr", "ggplot2", "devtools")

for (pkg in cran_packages) {
  if (!requireNamespace(pkg, quietly = TRUE)) {
    cat("正在安装:", pkg, "\n")
    install.packages(pkg, dependencies = TRUE)
  } else {
    cat(pkg, "已安装\n")
  }
}

# 检查Seurat版本
if (requireNamespace("Seurat", quietly = TRUE)) {
  cat("\nSeurat版本:", packageVersion("Seurat"), "\n")
}

cat("\n所有必需的包安装完成！\n")
cat("如果Seurat安装失败，可以尝试：\n")
cat("remotes::install_github('satijalab/seurat', 'seurat5')\n")


