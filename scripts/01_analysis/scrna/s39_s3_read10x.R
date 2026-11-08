suppressPackageStartupMessages({library(Seurat); library(Matrix)})
options(future.globals.maxSize = 8 * 1024^3)
W <- file.path(Sys.getenv("LPS_WORK_ROOT", "."), "_work")
D <- file.path(W, "gse221493", "raw10x")
cat("== Read10X from", D, "\n")
t0 <- Sys.time()
m <- Read10X(data.dir = D, gene.column = 2)
cat("dim =", dim(m)[1], "genes x", dim(m)[2], "cells\n")
cat("elapsed =", round(as.numeric(difftime(Sys.time(), t0, units = "secs")), 1), "s\n")
cat("class =", class(m)[1], " nnz =", length(m@x), "\n")
cat("total counts =", format(sum(m@x), big.mark = ","), "\n")

rs <- Matrix::rowSums(m)
o <- order(-rs)
cat("\n--- top 25 genes by total UMI ---\n")
for (i in o[1:25]) cat(sprintf("%-16s %12.0f  detect=%.1f%%\n",
                               rownames(m)[i], rs[i], 100 * mean(m[i, ] > 0)))

cat("\n--- sanity panel ---\n")
pan <- c("ACTB","GAPDH","B2M","MALAT1","ADIPOQ","PLIN1","LEP","FABP4","CD36","PPARG","SCD","FASN")
for (g in pan) {
  if (g %in% rownames(m)) {
    v <- m[g, ]
    cat(sprintf("%-8s nnz=%-8d detect=%6.1f%%  mean=%.4f  max=%d\n",
                g, sum(v > 0), 100 * mean(v > 0), mean(v), max(v)))
  } else cat(g, "MISSING\n")
}

cat("\n--- per-cell depth quantiles ---\n")
cs <- Matrix::colSums(m)
print(round(quantile(cs, c(0, .01, .25, .5, .75, .99, 1)), 1))

out <- file.path(W, "gse221493", "raw_matrix.rds")
cat("\nsaving ->", out, "\n")
t1 <- Sys.time()
saveRDS(m, out, compress = "gzip")
cat("saved,", round(file.size(out) / 1048576, 1), "MB, in",
    round(as.numeric(difftime(Sys.time(), t1, units = "secs")), 1), "s\n")
cat("DONE_R_READ\n")
