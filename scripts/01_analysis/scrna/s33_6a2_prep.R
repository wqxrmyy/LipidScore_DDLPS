## =============================================================================
##  S33 · 6A-补（pseudobulk 去 dropout）+ 6B 前置导出
##  ---------------------------------------------------------------------------
##  动机：6A 细胞级结果暴露张力 —— 单细胞 dropout（RICTOR 检出 20.4%、CD36 13.5%）
##        使单基因共定位无效（rho≈0.04/0.0015），而 bulk 层 r=0.55 稳健。
##        不降维解法 = pseudobulk 聚合（项目 Step 6 已有实现），并修正共定位口径。
##  同时导出：fig6_data/expr_tumor.rds + cell_meta_tumor.csv（供 6B），含 A−PPARG 评分
##  输出：_work/s33_6a2_report.txt / s33_6a2_pseudobulk.csv / s33_6a2_detect.csv
## =============================================================================
options(stringsAsFactors = FALSE)
suppressPackageStartupMessages({ library(Seurat); library(Matrix) })

W <- file.path(Sys.getenv("LPS_WORK_ROOT", "."), "_work")
D <- file.path(W, "p5_figures/fig6_data"); dir.create(D, showWarnings = FALSE)
LOG <- file.path(W, "s33_6a2_report.txt"); if (file.exists(LOG)) file.remove(LOG)
A <- function(...) { s <- paste0(...); cat(s, "\n", sep = ""); cat(s, "\n", sep = "", file = LOG, append = TRUE) }

A("===============================================================================")
A("S33 6A-补  pseudobulk 去 dropout + 6B 前置导出   ", format(Sys.time(), "%Y-%m-%d %H:%M:%S"))
A("===============================================================================")

t0 <- Sys.time()
obj <- readRDS(file.path(W, "s39_seurat.s6.rds"))
A("载入 s6.rds: ", paste(dim(obj), collapse = " x "), "  ",
  round(as.numeric(difftime(Sys.time(), t0, units = "secs")), 1), " s")
md <- obj[[]]

G <- c("RICTOR", "MTOR", "CD36", "PPARG", "FABP4", "ADIPOQ", "PLIN1",
       "SREBF1", "SCD", "HMGCR", "ACLY", "FASN", "MRC1", "ARG1", "IL10")
G <- intersect(G, rownames(obj))
A("基因: ", paste(G, collapse = ", "))
A("")

## ---------- A−PPARG 敏感性评分（Q3）----------
G5 <- intersect(c("SREBF1", "SCD", "HMGCR", "ACLY", "FASN"), rownames(obj))
m5 <- as.matrix(GetAssayData(obj, assay = "RNA", layer = "data")[G5, , drop = FALSE])
## ★ 逐基因（跨细胞）z-score：scale() 默认按列，须转置 → 与 D18 的 mean(z-score(log2CPM)) 同向
z5 <- t(scale(t(m5)))
md$score_noPPARG <- as.numeric(colMeans(z5, na.rm = TRUE))
A("A−PPARG 评分（5 基因 z 均值）: mean=", round(mean(md$score_noPPARG), 4),
  "  sd=", round(sd(md$score_noPPARG), 4))
A("与 AUCell 评分 Spearman: ",
  round(suppressWarnings(cor(md$score_noPPARG, md$score_aucell, method = "spearman")), 4))
rm(z5); gc(verbose = FALSE)

## ---------- 肿瘤细胞子集 ----------
sel <- which(md$L1 == "Tumor cells")
A("肿瘤细胞: ", length(sel))
md$sample <- paste(md$patient, md$tissue, sep = "_")

X <- GetAssayData(obj, assay = "RNA", layer = "data")[G, sel, drop = FALSE]
A("导出矩阵: ", nrow(X), " 基因 x ", ncol(X), " 细胞")

## ---------- pseudobulk：样本级 / 样本xL2 ----------
A("")
A("===============================================================================")
A("-- pseudobulk：样本级 与 样本xL2 级（评分 vs 轴基因）--")
A("===============================================================================")
expr_sel <- as.matrix(t(GetAssayData(obj, assay = "RNA", layer = "data")[G, sel, drop = FALSE]))
mt <- md[sel, ]

pb <- function(keyvec, tag) {
  key <- factor(keyvec)
  n <- as.vector(table(key))
  sm <- rowsum(mt$score_aucell, group = key, reorder = TRUE)[, 1] / n
  sp <- rowsum(mt$score_noPPARG, group = key, reorder = TRUE)[, 1] / n
  ge <- rowsum(expr_sel, group = key, reorder = TRUE)
  df <- data.frame(key = rownames(ge), n_cells = n, score_aucell = sm, score_noPPARG = sp,
                   stringsAsFactors = FALSE)
  for (g in G) df[[g]] <- ge[, g] / n
  # 检出率（该 key 内表达>0 的细胞占比）
  for (g in G) {
    pos <- tapply(expr_sel[, g] > 0, key, mean)
    df[[paste0("det_", g)]] <- as.numeric(pos)
  }
  df
}

pb_s <- pb(mt$sample, "sample")
A(sprintf("[样本级] %d 个样本（患者x组织）", nrow(pb_s)))
pb_l2 <- pb(paste(mt$sample, mt$L2, sep = "|"), "samplexL2")
A(sprintf("[样本xL2] %d 行", nrow(pb_l2)))
A("")

A("-- ★ 样本级 pseudobulk 相关（n = ", nrow(pb_s), " 样本）--")
A(sprintf("  %-9s %10s %11s %10s %11s %10s", "gene", "rho_aucell", "P", "rho_noPPARG", "P", "detect%"))
for (g in G) {
  r1 <- suppressWarnings(cor(pb_s$score_aucell, pb_s[[g]], method = "spearman"))
  p1 <- suppressWarnings(cor.test(pb_s$score_aucell, pb_s[[g]], method = "spearman")$p.value)
  r2 <- suppressWarnings(cor(pb_s$score_noPPARG, pb_s[[g]], method = "spearman"))
  p2 <- suppressWarnings(cor.test(pb_s$score_noPPARG, pb_s[[g]], method = "spearman")$p.value)
  A(sprintf("  %-9s %10.4f %11.2e %10.4f %11.2e %10.1f", g, r1, p1, r2, p2, 100 * mean(pb_s[[paste0("det_", g)]])))
}
A("")
A("-- ★ 样本xL2 pseudobulk 相关（n = ", nrow(pb_l2), " 行）--")
A(sprintf("  %-9s %10s %11s %10s %11s", "gene", "rho_aucell", "P", "rho_noPPARG", "P"))
for (g in G) {
  r1 <- suppressWarnings(cor(pb_l2$score_aucell, pb_l2[[g]], method = "spearman"))
  p1 <- suppressWarnings(cor.test(pb_l2$score_aucell, pb_l2[[g]], method = "spearman")$p.value)
  r2 <- suppressWarnings(cor(pb_l2$score_noPPARG, pb_l2[[g]], method = "spearman"))
  p2 <- suppressWarnings(cor.test(pb_l2$score_noPPARG, pb_l2[[g]], method = "spearman")$p.value)
  A(sprintf("  %-9s %10.4f %11.2e %10.4f %11.2e", g, r1, p1, r2, p2))
}
A("")
A("★ 注：pseudobulk 使每个点是「样本内肿瘤细胞的均值表达」，规避单细胞 dropout；")
A("  与 bulk（TCGA n=50，r=0.55）口径不同层级，解读时须写明。")
A("")

write.csv(pb_s,  file.path(W, "s33_6a2_pseudobulk_sample.csv"), row.names = FALSE)
write.csv(pb_l2, file.path(W, "s33_6a2_pseudobulk_samplexL2.csv"), row.names = FALSE)
A("写出 s33_6a2_pseudobulk_sample.csv / _samplexL2.csv")
A("")

## ---------- 修正版共定位：只在检出细胞中分层 ----------
A("===============================================================================")
A("-- 修正版共定位（上一版因 75 分位=0 失效）：仅在『检出细胞』内做四格 --")
A("===============================================================================")
det_tab <- list()
for (g in c("RICTOR", "CD36")) {
  v <- expr_sel[, g]
  pos <- v > 0
  A(sprintf("---- %s ----", g))
  A(sprintf("  检出细胞 %d / %d (%.2f%%)", sum(pos), length(pos), 100 * mean(pos)))
  if (sum(pos) < 200) { A("  过小，跳过"); next }
  sv <- mt$score_aucell
  q <- quantile(sv[pos], .75)
  hiS <- sv[pos] >= q
  qg <- quantile(v[pos], .75)
  hiG <- v[pos] >= qg
  A(sprintf("  检出细胞内：评分高 %d (%.1f%%) | %s高 %d (%.1f%%)",
            sum(hiS), 100 * mean(hiS), g, sum(hiG), 100 * mean(hiG)))
  A(sprintf("  交集 %d (%.1f%%) 观察/期望 = %.3f", sum(hiS & hiG), 100 * mean(hiS & hiG),
            mean(hiS & hiG) / (mean(hiS) * mean(hiG))))
  ft <- fisher.test(table(hiS, hiG))
  A(sprintf("  Fisher OR = %.3f  P = %.3g", ft$estimate, ft$p.value))
  det_tab[[g]] <- data.frame(gene = g, n_detect = sum(pos), pct_detect = 100 * mean(pos),
                             n_intersect = sum(hiS & hiG), obs_exp = mean(hiS & hiG) / (mean(hiS) * mean(hiG)),
                             OR = as.numeric(ft$estimate), P = ft$p.value)
  # 评分三分位 vs 检出率（剂量关系）—— 用秩切分，避免并列值使 breaks 不唯一
  rk <- rank(sv, ties.method = "first")
  tt <- cut(rk, breaks = c(0, length(sv) / 3, 2 * length(sv) / 3, length(sv)),
            labels = c("low", "mid", "high"), include.lowest = TRUE)
  dr <- tapply(pos, tt, mean)
  A("  评分三分位下的检出率: ", paste(sprintf("%s=%.1f%%", names(dr), 100 * dr), collapse = "  "))
  A("")
}
if (length(det_tab)) write.csv(do.call(rbind, det_tab), file.path(W, "s33_6a2_detect.csv"), row.names = FALSE)

## ---------- 导出 6B 输入 ----------
A("===============================================================================")
A("-- 导出 6B 输入 --")
A("===============================================================================")
saveRDS(X, file.path(D, "expr_tumor.rds"))
cm <- data.frame(barcode = colnames(X), patient = mt$patient, tissue = mt$tissue,
                 sample = mt$sample, L1 = mt$L1, L2 = mt$L2, L3 = mt$L3,
                 score_aucell = mt$score_aucell, score_ams = mt$score_ams,
                 score_noPPARG = mt$score_noPPARG, bucket = mt$bucket,
                 stringsAsFactors = FALSE)
write.csv(cm, file.path(D, "cell_meta_tumor.csv"), row.names = FALSE, fileEncoding = "UTF-8-BOM")
A("expr_tumor.rds: ", nrow(X), " 基因 x ", ncol(X), " 细胞（log1p CP10K，稀疏）")
A("cell_meta_tumor.csv: ", nrow(cm), " 行（含 score_noPPARG）")
A("")

## ---------- 评分分布核对（供图注）----------
A("-- 评分分布核对（供 6A 图注）--")
A(sprintf("  肿瘤细胞内 AUCell: mean=%.4f  中位=%.4f  P90=%.4f", mean(mt$score_aucell),
          median(mt$score_aucell), quantile(mt$score_aucell, .9)))
A(sprintf("  肿瘤细胞内 A−PPARG: mean=%.4f  中位=%.4f", mean(mt$score_noPPARG), median(mt$score_noPPARG)))
ntum <- length(sel); nall <- ncol(obj)
A(sprintf("  肿瘤细胞占全体比例: %d / %d = %.1f%%", ntum, nall, 100 * ntum / nall))
A("")
A("DONE_S33_6A2")
