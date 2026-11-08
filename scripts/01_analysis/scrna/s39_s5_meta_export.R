## Step 5 前置：从 s4.rds 导出 metadata，供粒度诊断（Python 侧）使用
suppressPackageStartupMessages({ library(Seurat) })
options(future.globals.maxSize = 8000 * 1024^2)
plan("sequential")   # Windows 铁律：并行在 2GB 级对象上会崩且更慢

p <- file.path(Sys.getenv("LPS_WORK_ROOT", "."), "_work", "s39_seurat.s4.rds")
obj <- readRDS(p)

cat("=== 对象 ===\n")
cat("class:", paste(class(obj), collapse=","), "\n")
cat("dim:", paste(dim(obj), collapse=" x "), "\n")

md <- obj@meta.data
cat("meta ncol =", ncol(md), " nrow =", nrow(md), "\n")
cat("meta cols:", paste(colnames(md), collapse=" | "), "\n")

# reductions 是否含 harmony
cat("reductions:", paste(names(obj@reductions), collapse=" | "), "\n")

# 候选粒度列
cand <- c("patient","tissue","site","subtype_code","L1","L2","L3",
          "celltype_L1","celltype_L2","celltype_L3","grade",
          "percent.mt","nFeature_RNA","nCount_RNA")
keep <- intersect(cand, colnames(md))
cat("keep:", paste(keep, collapse=" | "), "\n")

# 对疑似粒度列做唯一值统计
for (k in keep) {
  v <- md[[k]]
  if (is.numeric(v)) {
    cat(sprintf("  %-14s numeric  min=%.3g med=%.3g max=%.3g\n",
                k, min(v,na.rm=TRUE), median(v,na.rm=TRUE), max(v,na.rm=TRUE)))
  } else {
    u <- unique(as.character(v))
    cat(sprintf("  %-14s factor/char  nuniq=%d\n", k, length(u)))
    if (length(u) <= 12) cat("      vals:", paste(sort(u), collapse=" | "), "\n")
  }
}

out <- data.frame(cell = rownames(md), md[, keep, drop = FALSE],
                  stringsAsFactors = FALSE, check.names = FALSE)
d <- file.path(Sys.getenv("LPS_WORK_ROOT", "."), "_work")
write.csv(out, file.path(d, "s39_s5_meta.csv"), row.names = FALSE, fileEncoding = "utf-8")
cat("wrote s39_s5_meta.csv  rows =", nrow(out), " cols =", ncol(out), "\n")
cat("DONE_S5_META\n")
