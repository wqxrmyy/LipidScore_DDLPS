suppressPackageStartupMessages({library(Seurat); library(Matrix); library(future)})
plan("sequential")

W <- file.path(Sys.getenv("LPS_WORK_ROOT", "."), "_work")
LOG <- file.path(W, "s39_s3b2_report.txt")
sink(LOG, split = TRUE)
say <- function(...) cat(sprintf(...), "\n")
tt <- function(t0) as.numeric(difftime(Sys.time(), t0, units = "secs"))

say("== Step 3b-2  UMAP 与嵌入一致性 ==")
say("start: %s", format(Sys.time()))

so <- readRDS(file.path(W, "s39_seurat.s3.rds"))
say("loaded: %d genes x %d cells, HVG=%d", nrow(so), ncol(so), length(VariableFeatures(so)))

knn_purity <- function(emb, labels, n = 3000, k = 10, seed = 42) {
  set.seed(seed)
  idx <- sample(seq_len(nrow(emb)), min(n, nrow(emb)))
  e <- emb[idx, , drop = FALSE]; lb <- labels[idx]
  d <- as.matrix(dist(e)); diag(d) <- Inf
  nn <- t(apply(d, 1, function(x) order(x)[1:k]))
  mean(sapply(seq_len(nrow(nn)), function(i) mean(lb[nn[i, ]] == lb[i])))
}
base_of <- function(x) { p <- table(x) / length(x); sum(p^2) }

say("\n-- RunUMAP (dims = 1:30, n_threads = 8) --")
t0 <- Sys.time()
so <- RunUMAP(so, dims = 1:30, n_threads = 8, verbose = FALSE)
say("done in %.1fs", tt(t0))

say("\n-- ★ 嵌入一致性：L3 / L1 kNN 纯度（抽样 3000, k = 10）--")
for (lv in c("L1", "L2", "L3")) {
  b  <- base_of(so@meta.data[[lv]])
  po <- knn_purity(cbind(so$umap1, so$umap2), so@meta.data[[lv]])
  pw <- knn_purity(Embeddings(so, "umap"), so@meta.data[[lv]])
  say("%-4s 随机基线 %.4f | 官方坐标 %.4f (%.1fx) | 自算 UMAP %.4f (%.1fx)",
      lv, b, po, po / b, pw, pw / b)
}

say("\n-- 自算 UMAP 的组织分布（每组织质心与离散度）--")
emb <- Embeddings(so, "umap")
ord <- c("nl", "lipoma", "WDLPS", "DDLPS-WD", "DDLPS-sWD", "DDLPS-DD", "DDLPS-DDK")
ord <- c(ord[ord %in% unique(so$tissue)], setdiff(sort(unique(so$tissue)), ord))
say("%-12s %8s %10s %10s %10s", "tissue", "n", "cx", "cy", "spread")
for (t in ord) {
  s <- emb[so$tissue == t, , drop = FALSE]
  say("%-12s %8d %10.2f %10.2f %10.2f", t, nrow(s),
      mean(s[, 1]), mean(s[, 2]), mean(apply(s, 1, function(z) sqrt(sum(z^2)))))
}

out <- file.path(W, "s39_seurat.s3.rds")
say("\nsaving -> %s", out)
t0 <- Sys.time()
saveRDS(so, out, compress = "gzip")
say("saved: %.1f MB in %.1fs", file.size(out) / 1048576, tt(t0))
say("\nDONE_S3B2")
sink()
