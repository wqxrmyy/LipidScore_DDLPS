## =============================================================================
##  S33 · Fig6 6D：TAM 亚群中 mTORC2/PPAR-γ/CD36 的表达 + 与高评分肿瘤细胞的邻域富集
##  ---------------------------------------------------------------------------
##  Q7/Q13 口径：本项目**无空间转录组** → 一律称「UMAP 邻域富集」，**全文禁用 spatial**
##  Q5 口径（不降维）：主口径 = L2 Macrophages（TAM1+TAM3）；敏感性 = + 单核（Monocytes）
##  输入：s39_seurat.s6.rds
##  输出：_work/s33_6d_*.csv / s33_6d_report.txt
## =============================================================================
options(stringsAsFactors = FALSE)
suppressPackageStartupMessages({ library(Seurat); library(Matrix) })

W <- file.path(Sys.getenv("LPS_WORK_ROOT", "."), "_work")
LOG <- file.path(W, "s33_6d_report.txt"); if (file.exists(LOG)) file.remove(LOG)
A <- function(...) { s <- paste0(...); cat(s, "\n", sep = ""); cat(s, "\n", sep = "", file = LOG, append = TRUE) }

A("===============================================================================")
A("S33 6D  TAM 亚群表达 + UMAP 邻域富集   ", format(Sys.time(), "%Y-%m-%d %H:%M:%S"))
A("===============================================================================")

obj <- readRDS(file.path(W, "s39_seurat.s6.rds"))
md <- obj[[]]
A("载入: ", paste(dim(obj), collapse = " x "))
A("")

## ---------- (1) 髓系构成（Q5 口径举证）----------
A("-- [1] 髓系构成 --")
mye <- md[md$L1 == "Myeloid cells", ]
tb <- table(mye$L2); tb3 <- table(mye$L3)
for (k in names(sort(tb, decreasing = TRUE))) A(sprintf("  L2 %-16s %7d", k, tb[k]))
A("  ---- L3（Macrophages 内）----")
mac3 <- table(md$L3[md$L2 == "Macrophages"])
for (k in names(sort(mac3, decreasing = TRUE))) A(sprintf("  L3 %-16s %7d", k, mac3[k]))
A("")

## ---------- (2) 巨噬亚群中轴基因表达 ----------
G <- intersect(c("RICTOR", "MTOR", "CD36", "PPARG", "FABP4",
                 "MRC1", "ARG1", "IL10", "TGFB1", "CD163", "SPP1", "APOE", "TREM2", "LYZ"), rownames(obj))
E <- as.matrix(t(GetAssayData(obj, assay = "RNA", layer = "data")[G, , drop = FALSE]))

A("-- [2] 巨噬亚群轴基因表达（主口径 L2 = Macrophages）--")
mac <- which(md$L2 == "Macrophages")
A(sprintf("  巨噬细胞 n = %d（L3: %s）", length(mac),
          paste(sprintf("%s=%d", names(mac3), as.integer(mac3)), collapse = ", ")))
A(sprintf("  %-8s %10s %10s %10s %10s %10s", "gene", "detect%", "mean", "TAM1_det", "TAM3_det", "肿瘤_mean"))
tum <- which(md$L1 == "Tumor cells")
rows <- list()
for (g in G) {
  d_all <- 100 * mean(E[mac, g] > 0)
  m_all <- mean(E[mac, g])
  d1 <- if (any(md$L3[mac] == "TAM1")) 100 * mean(E[mac[md$L3[mac] == "TAM1"], g] > 0) else NA
  d3 <- if (any(md$L3[mac] == "TAM3")) 100 * mean(E[mac[md$L3[mac] == "TAM3"], g] > 0) else NA
  mt <- mean(E[tum, g])
  A(sprintf("  %-8s %10.2f %10.3f %10.2f %10.2f %10.3f", g, d_all, m_all, d1, d3, mt))
  rows[[length(rows) + 1]] <- data.frame(gene = g, det_mac = d_all, mean_mac = m_all,
                                         det_TAM1 = d1, det_TAM3 = d3, mean_tumor = mt)
}
write.csv(do.call(rbind, rows), file.path(W, "s33_6d_mac_genes.csv"), row.names = FALSE)
A("")

## ---------- (3) 轴基因在巨噬细胞的检出率 vs 全细胞背景 ----------
A("-- [3] 轴基因检出率：巨噬 vs 全体（判断是否巨噬特异）--")
A(sprintf("  %-8s %10s %10s %8s", "gene", "mac%", "all%", "倍数"))
for (g in c("RICTOR", "MTOR", "CD36", "PPARG", "MRC1", "ARG1", "IL10")) {
  dm <- 100 * mean(E[mac, g] > 0); da <- 100 * mean(E[, g] > 0)
  A(sprintf("  %-8s %10.2f %10.2f %8.2f", g, dm, da, dm / max(da, 1e-9)))
}
A("")

## ---------- (4) M2 样标志物打分（Q11：不得称 M2 亚群）----------
A("-- [4] M2 样标志物打分（MRC1/ARG1/IL10，Q11：仅作『M2 样标志物打分』）--")
m2g <- intersect(c("MRC1", "ARG1", "IL10"), rownames(obj))
mm <- as.matrix(GetAssayData(obj, assay = "RNA", layer = "data")[m2g, , drop = FALSE])
z <- t(scale(t(mm)))                       # 逐基因跨细胞 z
m2s <- colMeans(z, na.rm = TRUE)
A(sprintf("  打分：巨噬 mean=%.4f 中位=%.4f | 全体 mean=%.4f", mean(m2s[mac]), median(m2s[mac]), mean(m2s)))
A(sprintf("  TAM1 mean=%.4f | TAM3 mean=%.4f",
          mean(m2s[mac[md$L3[mac] == "TAM1"]]), mean(m2s[mac[md$L3[mac] == "TAM3"]])))
wt <- suppressWarnings(wilcox.test(m2s[mac][md$L3[mac] == "TAM1"], m2s[mac][md$L3[mac] == "TAM3"]))
A(sprintf("  TAM1 vs TAM3 Wilcoxon P = %.4g", wt$p.value))
A("")

## ---------- (5) UMAP 邻域富集（禁止 spatial 措辞）----------
A("===============================================================================")
A("-- [5] UMAP 邻域富集：高评分肿瘤细胞的 k 近邻中巨噬细胞比例 vs 期望 --")
A("===============================================================================")
coord <- as.matrix(md[, c("umap1", "umap2")])
qS <- quantile(md$score_aucell[tum], .75)
hiT <- tum[md$score_aucell[tum] >= qS]
loT <- tum[md$score_aucell[tum] <  qS]
A(sprintf("  肿瘤细胞: 高评分 %d, 低评分 %d（四分位切点 %.4f）", length(hiT), length(loT), qS))

knn_prop <- function(idx, cand, k = 30) {
  ## 在候选池 cand 内计算 idx 的 k 近邻，返回近邻中巨噬细胞比例的均值
  cds <- coord[cand, , drop = FALSE]
  org <- coord[idx, , drop = FALSE]
  out <- numeric(nrow(org))
  step <- 2000
  mac_of_cand <- md$L2[cand] == "Macrophages"
  for (s in seq(1, nrow(org), by = step)) {
    e <- min(s + step - 1, nrow(org))
    blk <- org[s:e, , drop = FALSE]
    ## 距离（分块，避免巨矩阵）
    d <- outer(blk[, 1], cds[, 1], "-")^2 + outer(blk[, 2], cds[, 2], "-")^2
    for (i in seq_len(nrow(blk))) {
      o <- order(d[i, ])[seq_len(k + 1)]        # 首个通常是自身
      o <- setdiff(o, s + i - 1)[seq_len(k)]
      out[s + i - 1] <- mean(mac_of_cand[o])
    }
  }
  out
}

## 用「同患者同组织」的样本池做对照（控制患者/组织混杂）
A("  邻域用『同患者同组织』的细胞池（控制患者与组织混杂）")
res <- list()
run <- function(idx, tag) {
  if (length(idx) < 20) { A(sprintf("  [%s] 过小", tag)); return(NULL) }
  ## 逐患者-组织分层
  key <- paste(md$patient, md$tissue, sep = "_")
  props <- pct_exp <- NULL
  allp <- numeric(0); allexp <- numeric(0)
  for (k in unique(key[idx])) {
    ii <- idx[key[idx] == k]
    if (length(ii) < 5) next
    cand <- which(key == k)
    ## 抽样上限，控制计算量（每层最多 400 个 query）
    if (length(ii) > 400) ii <- sample(ii, 400)
    p <- knn_prop(ii, cand, k = 30)
    exp_mac <- mean(md$L2[cand] == "Macrophages")
    allp <- c(allp, p); allexp <- c(allexp, rep(exp_mac, length(p)))
  }
  if (!length(allp)) { A(sprintf("  [%s] 无有效层", tag)); return(NULL) }
  obs <- mean(allp); exp0 <- mean(allexp)
  tt <- suppressWarnings(t.test(allp, mu = exp0))
  A(sprintf("  [%s] n_query=%d  近邻中巨噬比例 观察=%.4f  期望(池内占比)=%.4f  富集倍数=%.2f  t 检验 P=%.3g",
            tag, length(allp), obs, exp0, obs / exp0, tt$p.value))
  list(tag = tag, n = length(allp), observed = obs, expected = exp0, fold = obs / exp0, p = tt$p.value)
}

A("")
r1 <- run(hiT, "high_score_tumor")
r2 <- run(loT, "low_score_tumor")
A("")
if (!is.null(r1) && !is.null(r2)) {
  A(sprintf("  ★ 高 vs 低 评分肿瘤细胞：近邻巨噬比例差 = %+.4f（%.2f× vs %.2f×）",
            r1$observed - r2$observed, r1$fold, r2$fold))
}
write.csv(do.call(rbind, lapply(list(r1, r2), function(z) if (!is.null(z)) as.data.frame(z))),
          file.path(W, "s33_6d_neighborhood.csv"), row.names = FALSE)
A("")
A("★ 措辞纪律：本节为 UMAP 二维嵌入上的邻域富集，**不是空间邻接**；")
A("  稿件中不得出现 spatial / spatially adjacent 等表述（Q7/Q13）。")
A("")
A("DONE_S33_6D")
