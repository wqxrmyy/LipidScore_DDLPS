## S33 · Fig6 执行前置：R 包与内存环境盘点
suppressWarnings(suppressMessages({
  ip <- rownames(installed.packages())
}))
pkgs <- c("Seurat","SeuratObject","AUCell","GSEABase","GSVA","CellChat",
          "monocle3","survival","survminer","randomForestSRC","glmnet",
          "ComplexHeatmap","harmony","SingleCellExperiment","msigdbr",
          "estimate","patchwork","ggplot2","Matrix","data.table","BiocParallel")
cat("== R 包 ==\n")
for (p in pkgs) {
  if (p %in% ip) cat(sprintf("  OK      %-20s %s\n", p, as.character(packageVersion(p))))
  else          cat(sprintf("  MISSING %-20s\n", p))
}
cat("\n== 版本 ==\n")
cat("  R ", as.character(getRversion()), "\n", sep = "")
cat("  platform ", R.version$platform, "\n", sep = "")
cat("  n_threads ", parallel::detectCores(), "\n", sep = "")
cat("  libpath:\n"); for (l in .libPaths()) cat("    ", l, "\n", sep = "")
