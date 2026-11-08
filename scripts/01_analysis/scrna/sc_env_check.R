#!/usr/bin/env Rscript
# 单细胞分析工具链盘点
pkgs <- c("Seurat", "SeuratObject", "Matrix", "hdf5r", "AUCell", "GSEABase",
          "harmony", "SingleCellExperiment", "scater", "scran", "scDblFinder",
          "SingleR", "celldex", "ComplexHeatmap", "patchwork", "cowplot",
          "ggplot2", "dplyr", "data.table", "survival", "survminer", "ragg")
cat("R:", R.version.string, "\n")
cat("libpaths:\n"); cat(paste0("  ", .libPaths()), sep = "\n"); cat("\n\n")
cat(sprintf("%-24s %s\n", "package", "version"))
cat(strrep("-", 46), "\n")
for (p in pkgs) {
  v <- tryCatch(if (requireNamespace(p, quietly = TRUE)) as.character(packageVersion(p)) else "MISSING",
                error = function(e) "ERR")
  cat(sprintf("%-24s %s\n", p, v))
}
cat("\nCPU cores:", parallel::detectCores(), "\n")
cat("mem.limit MB:", round(as.numeric(memory.limit()) / 1024, 0), "\n")
