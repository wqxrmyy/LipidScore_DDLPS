suppressPackageStartupMessages({
  library(Seurat); library(Matrix); library(future)
})
options(future.globals.maxSize = 16 * 1024^3)
plan("multisession", workers = 8)

W <- file.path(Sys.getenv("LPS_WORK_ROOT", "."), "_work")
LOG <- file.path(W, "s39_s3b_report.txt")
sink(LOG, split = TRUE)
say <- function(...) cat(sprintf(...), "\n")

say("== Step 3b  QC 过滤 + 归一化 + 降维 ==")
say("start: %s", format(Sys.time()))

so <- readRDS(file.path(W, "s39_seurat.qc.rds"))
say("loaded: %d genes x %d cells", nrow(so), ncol(so))

# ---------- 1. QC 过滤 ----------
NF_MIN <- 200
MT_MAX <- 20
n0 <- ncol(so)
md0 <- so@meta.data
before <- table(md0$tissue)

keep <- md0$nFeature_RNA >= NF_MIN & md0$percent.mt <= MT_MAX
say("\n-- 过滤 --")
say("criterion: nFeature_RNA >= %d  AND  percent.mt <= %d", NF_MIN, MT_MAX)
say("kept %d / %d  (%.2f%%), removed %d", sum(keep), n0, 100 * mean(keep), sum(!keep))

so <- subset(so, cells = colnames(so)[keep])
after <- table(so$tissue)
ord <- c("nl", "lipoma", "WDLPS", "DDLPS-WD", "DDLPS-sWD", "DDLPS-DD", "DDLPS-DDK")
ord <- c(ord[ord %in% names(before)], setdiff(names(before), ord))
say("\n%-12s %9s %9s %9s", "tissue", "before", "after", "keep%")
for (t in ord) {
  b <- if (t %in% names(before)) before[[t]] else 0
  a <- if (t %in% names(after)) after[[t]] else 0
  say("%-12s %9d %9d %8.1f%%", t, b, a, 100 * a / b)
}
say("TOTAL        %9d %9d %8.1f%%", n0, ncol(so), 100 * ncol(so) / n0)

pairs_t <- c("P1","P4","P7","P10","P14","P15","P16")
say("\n-- 严格配对 7 例的保留情况 --")
say("%-6s %-12s %8s %8s %8s", "pat", "tissue", "before", "after", "keep%")
for (p in pairs_t) {
  for (t in c("DDLPS-WD", "DDLPS-DD")) {
    b <- sum(md0$patient == p & md0$tissue == t)
    a <- sum(so$patient == p & so$tissue == t)
    if (b > 0) say("%-6s %-12s %8d %8d %7.1f%%", p, t, b, a, 100 * a / b)
  }
}

# ---------- 2. 归一化 ----------
say("\n-- NormalizeData (LogNormalize, scale.factor=1e4) --")
t0 <- Sys.time()
so <- NormalizeData(so, normalization.method = "LogNormalize", scale.factor = 1e4, verbose = FALSE)
say("done in %.1fs", as.numeric(difftime(Sys.time(), t0, units = "secs")))

# ---------- 3. HVG ----------
say("\n-- FindVariableFeatures (vst, n=2000) --")
t0 <- Sys.time()
so <- FindVariableFeatures(so, selection.method = "vst", nfeatures = 2000, verbose = FALSE)
say("done in %.1fs", as.numeric(difftime(Sys.time(), t0, units = "secs")))
hv <- VariableFeatures(so)
say("top 30 HVG: %s", paste(head(hv, 30), collapse = ", "))
score6 <- c("SREBF1", "SCD", "HMGCR", "ACLY", "PPARG", "FASN")
say("评分 6 基因是否入选 HVG: %s", paste(sprintf("%s=%s", score6, score6 %in% hv), collapse = ", "))
markers <- c("FABP4", "CD36", "ADIPOQ", "PLIN1", "LEP", "LPL", "PDGFRA", "DCN", "PTPRC", "PECAM1", "MDM2", "CDK4")
say("marker 入选 HVG: %s", paste(sprintf("%s=%s", markers, markers %in% hv), collapse = ", "))

# ---------- 4. Scale + PCA ----------
say("\n-- ScaleData + RunPCA (npcs=50) --")
t0 <- Sys.time()
so <- ScaleData(so, features = hv, verbose = FALSE)
so <- RunPCA(so, features = hv, npcs = 50, verbose = FALSE)
say("done in %.1fs", as.numeric(difftime(Sys.time(), t0, units = "secs")))
ve <- so[["pca"]]@stdev^2
pct <- 100 * ve / sum(ve)
say("PC1-10 方差解释(%%): %s", paste(sprintf("%.1f", pct[1:10]), collapse = ", "))
say("累计 PC1-20: %.1f%%   PC1-30: %.1f%%", sum(pct[1:20]), sum(pct[1:30]))

# ---------- 5. UMAP ----------
say("\n-- RunUMAP (dims=1:30) --")
t0 <- Sys.time()
so <- RunUMAP(so, dims = 1:30, verbose = FALSE)
say("done in %.1fs", as.numeric(difftime(Sys.time(), t0, units = "secs")))

# ---------- 6. 与官方坐标的一致性对比（L3 kNN 纯度）----------
knn_purity <- function(emb, labels, n = 3000, k = 10, seed = 42) {
  set.seed(seed)
  idx <- sample(seq_along(labels), min(n, length(labels)))
  e <- emb[idx, , drop = FALSE]
  lb <- labels[idx]
  d <- as.matrix(dist(e))
  diag(d) <- Inf
  nn <- t(apply(d, 1, function(x) order(x)[1:k]))
  mean(sapply(seq_len(nrow(nn)), function(i) mean(lb[nn[i, ]] == lb[i])))
}
base <- mean(table(so$L3) * (table(so$L3) / ncol(so)))   # 随机基线 = Σp²
say("\n-- 嵌入一致性：L3 kNN 纯度 (抽样 3000, k=10) --")
say("随机基线 (Σp²) = %.4f", base)
p_off <- knn_purity(cbind(so$umap1, so$umap2), so$L3)
say("官方坐标 (col9/col10) = %.4f  (%.1fx)", p_off, p_off / base)
p_own <- knn_purity(Embeddings(so, "umap"), so$L3)
say("自算 UMAP (dims=1:30) = %.4f  (%.1fx)", p_own, p_own / base)

# 按 L1 也看一遍
base1 <- mean(table(so$L1) * (table(so$L1) / ncol(so)))
p_off1 <- knn_purity(cbind(so$umap1, so$umap2), so$L1)
p_own1 <- knn_purity(Embeddings(so, "umap"), so$L1)
say("L1 层：基线 %.4f | 官方 %.4f | 自算 %.4f", base1, p_off1, p_own1)

# ---------- 7. 落盘 ----------
out <- file.path(W, "s39_seurat.s3.rds")
say("\nsaving -> %s", out)
t0 <- Sys.time()
saveRDS(so, out, compress = "gzip")
say("saved: %.1f MB in %.1fs", file.size(out) / 1048576,
    as.numeric(difftime(Sys.time(), t0, units = "secs")))
say("\nDONE_S3B")
sink()
