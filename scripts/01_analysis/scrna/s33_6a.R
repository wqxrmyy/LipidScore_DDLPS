## =============================================================================
##  S33 · Fig6 6A：评分在单细胞的投影 × RICTOR/CD36 共定位（数字先行）
##  ---------------------------------------------------------------------------
##  依据：D26 §一 Q3（统计只用 RICTOR/CD36 + A−PPARG 敏感性）、Q4（肿瘤细胞主口径 + 敏感性）
##  输入：_work/s39_seurat.s6.rds（含 score_ams / score_aucell / L1-L3 / bucket）
##  输出：_work/s33_6a_*.csv / _work/s33_6a_report.txt
##  纪律：Windows 下 plan("sequential")；不启用 future
## =============================================================================
options(stringsAsFactors = FALSE)
suppressPackageStartupMessages({ library(Seurat); library(Matrix) })

W <- file.path(Sys.getenv("LPS_WORK_ROOT", "."), "_work")
LOG <- file.path(W, "s33_6a_report.txt")
if (file.exists(LOG)) file.remove(LOG)
A <- function(...) { s <- paste0(...); cat(s, "\n", sep = ""); cat(s, "\n", sep = "", file = LOG, append = TRUE) }

A("===============================================================================")
A("S33 6A  评分投影 x RICTOR/CD36 共定位   ", format(Sys.time(), "%Y-%m-%d %H:%M:%S"))
A("===============================================================================")

t0 <- Sys.time()
obj <- readRDS(file.path(W, "s39_seurat.s6.rds"))
A("dim = ", paste(dim(obj), collapse = " x "), "  载入耗时 ",
  round(as.numeric(difftime(Sys.time(), t0, units = "secs")), 1), " s")
md <- obj[[]]
A("meta 列: ", paste(colnames(md), collapse = ", "))
A("")

## ---- 基因（ENSG 已核验；此处按 symbol 取，存在性须查）----
G <- c("RICTOR", "MTOR", "CD36", "PPARG", "FABP4", "ADIPOQ", "PLIN1",
       "SREBF1", "SCD", "HMGCR", "ACLY", "FASN")
miss <- setdiff(G, rownames(obj))
A("基因存在性: ", if (!length(miss)) "全部命中 ✓" else paste0("★缺失 ", paste(miss, collapse = ", ")))
G <- intersect(G, rownames(obj))

expr <- as.matrix(t(GetAssayData(obj, assay = "RNA", layer = "data")[G, , drop = FALSE]))
A("表达矩阵（log1p CP10K）: ", nrow(expr), " 细胞 x ", ncol(expr), " 基因")
A("")

## ---- 细胞池定义（Q4 主口径 + 敏感性）----
md$isTumor <- md$L1 == "Tumor cells"
pools <- list(
  tumor_all   = which(md$isTumor),
  tumor_noAd  = which(md$isTumor & md$L2 != "Adipo diff"),
  tumor_prog  = which(md$isTumor & md$L2 %in% c("Invasion", "Hypoxia", "Stemness", "ECM remodeling", "Angiogenesis", "Cycling cells")),
  adipo_diff  = which(md$L2 == "Adipo diff"),
  all_cells   = seq_len(nrow(md))
)
A("-- 细胞池 --")
for (k in names(pools)) A(sprintf("  %-12s %7d", k, length(pools[[k]])))
A("")

## ---- 评分分布：全体 vs 肿瘤池（回应 Fig5「评分住在稀有 Normal adipocytes」）----
A("-- 评分（AUCell / AddModuleScore）分布 --")
A(sprintf("%-12s %8s %10s %10s %10s %10s %10s", "pool", "n", "auc_mean", "auc_med", "auc_p90", "ams_mean", "ams_med"))
for (k in names(pools)) {
  i <- pools[[k]]
  A(sprintf("%-12s %8d %10.4f %10.4f %10.4f %10.4f %10.4f", k, length(i),
            mean(md$score_aucell[i]), median(md$score_aucell[i]), quantile(md$score_aucell[i], .9),
            mean(md$score_ams[i]), median(md$score_ams[i])))
}
A("")
nk <- function(l2) sum(md$L2 == l2)
A("（对照）Normal adipocytes n = ", nk("Normal adipocytes"),
  "  auc_mean = ", round(mean(md$score_aucell[md$L2 == "Normal adipocytes"]), 4),
  "｜Adipo diff n = ", nk("Adipo diff"),
  "  auc_mean = ", round(mean(md$score_aucell[md$L2 == "Adipo diff"]), 4))
A("")

## ---- 6A 核心：评分 vs 轴基因表达（肿瘤池内，Spearman + 检出率）----
A("===============================================================================")
A("-- 6A 核心：肿瘤细胞内 评分 AUCell vs 轴基因表达 --")
A("===============================================================================")
res <- list()
for (k in c("tumor_all", "tumor_noAd", "tumor_prog", "adipo_diff")) {
  i <- pools[[k]]
  if (length(i) < 50) { A(sprintf("[%s] 过小(%d)，跳过", k, length(i))); next }
  A(sprintf("---- pool = %s (n = %d) ----", k, length(i)))
  A(sprintf("  %-9s %9s %10s %11s %10s %10s", "gene", "detect%", "rho_aucell", "P", "rho_ams", "mean_expr"))
  for (g in G) {
    v <- expr[i, g]
    r1 <- suppressWarnings(cor(md$score_aucell[i], v, method = "spearman"))
    p1 <- suppressWarnings(cor.test(md$score_aucell[i], v, method = "spearman")$p.value)
    r2 <- suppressWarnings(cor(md$score_ams[i], v, method = "spearman"))
    res[[length(res) + 1]] <- data.frame(pool = k, gene = g, n = length(i),
                                         detect_pct = 100 * mean(v > 0),
                                         rho_aucell = r1, p_aucell = p1, rho_ams = r2,
                                         mean_expr = mean(v))
    A(sprintf("  %-9s %9.2f %10.4f %11.2e %10.4f %10.3f", g, 100 * mean(v > 0), r1, p1, r2, mean(v)))
  }
  A("")
}
RES <- do.call(rbind, res)
write.csv(RES, file.path(W, "s33_6a_corr.csv"), row.names = FALSE)
A("写出 s33_6a_corr.csv （", nrow(RES), " 行）")
A("")

## ---- 共定位：高评分 ∩ RICTOR高 ∩ CD36高 ----
A("===============================================================================")
A("-- 共定位：肿瘤细胞中 评分高 & RICTOR高 & CD36高 的交集（各取全池上四分位）--")
A("===============================================================================")
for (k in c("tumor_all", "tumor_noAd")) {
  i <- pools[[k]]
  qs <- quantile(md$score_aucell[i], .75); qr <- quantile(expr[i, "RICTOR"], .75); qc <- quantile(expr[i, "CD36"], .75)
  hiS <- md$score_aucell[i] >= qs; hiR <- expr[i, "RICTOR"] >= qr; hiC <- expr[i, "CD36"] >= qc
  A(sprintf("---- pool = %s (n = %d) ----", k, length(i)))
  A(sprintf("  评分高: %6d (%.1f%%) | RICTOR高: %6d (%.1f%%) | CD36高: %6d (%.1f%%)",
            sum(hiS), 100 * mean(hiS), sum(hiR), 100 * mean(hiR), sum(hiC), 100 * mean(hiC)))
  A(sprintf("  评分高 ∩ RICTOR高: %6d (%.1f%%)   观察/期望 = %.2f", sum(hiS & hiR), 100 * mean(hiS & hiR),
            (mean(hiS & hiR)) / (mean(hiS) * mean(hiR))))
  A(sprintf("  评分高 ∩ CD36高: %6d (%.1f%%)   观察/期望 = %.2f", sum(hiS & hiC), 100 * mean(hiS & hiC),
            (mean(hiS & hiC)) / (mean(hiS) * mean(hiC))))
  A(sprintf("  三者同时: %6d (%.1f%%)   观察/期望 = %.2f", sum(hiS & hiR & hiC), 100 * mean(hiS & hiR & hiC),
            (mean(hiS & hiR & hiC)) / (mean(hiS) * mean(hiR) * mean(hiC))))
  # Fisher 检验 评分高 vs RICTOR高 的富集
  ft <- fisher.test(table(hiS, hiR))
  A(sprintf("  Fisher(评分高 x RICTOR高): OR = %.3f  P = %.3g", ft$estimate, ft$p.value))
  ft2 <- fisher.test(table(hiS, hiC))
  A(sprintf("  Fisher(评分高 x CD36高):   OR = %.3f  P = %.3g", ft2$estimate, ft2$p.value))
  A("")
}

## ---- A−PPARG 敏感性：剔除 PPARG 重算评分，复算与 RICTOR/CD36 的相关 ----
A("===============================================================================")
A("-- Q3 双轨：评分 A−PPARG（剔除 PPARG 后 5 基因 z 均值）复算 --")
A("===============================================================================")
G5 <- c("SREBF1", "SCD", "HMGCR", "ACLY", "FASN")
G5 <- intersect(G5, rownames(obj))
m5 <- as.matrix(GetAssayData(obj, assay = "RNA", layer = "data")[G5, , drop = FALSE])
z <- t(scale(t(m5)))   # ★ 逐基因（跨细胞）z-score，与 D18 的 mean(z-score(log2CPM)) 同向
score_noPP <- colMeans(z, na.rm = TRUE)
i <- pools$tumor_all
for (g in c("RICTOR", "CD36", "PPARG")) {
  if (!g %in% colnames(expr)) next
  r <- suppressWarnings(cor(score_noPP[i], expr[i, g], method = "spearman"))
  p <- suppressWarnings(cor.test(score_noPP[i], expr[i, g], method = "spearman")$p.value)
  r0 <- suppressWarnings(cor(md$score_aucell[i], expr[i, g], method = "spearman"))
  A(sprintf("  tumor_all: A-PPARG vs %-7s rho = %+.4f (P=%.2e)   |   A(AUCell) vs %-7s rho = %+.4f",
            g, r, p, g, r0))
}
A("")
A("注：AUCell 与 z 均值口径不同，此处比的是**方向与显著性是否一致**，不是数值可比。")
A("")

## ---- 导出细胞级子表（供出图 join）----
out <- data.frame(
  barcode = rownames(md), patient = md$patient, tissue = md$tissue, L1 = md$L1, L2 = md$L2,
  bucket = md$bucket, score_aucell = md$score_aucell, score_ams = md$score_ams,
  RICTOR = expr[, "RICTOR"], CD36 = expr[, "CD36"],
  PPARG = if ("PPARG" %in% colnames(expr)) expr[, "PPARG"] else NA_real_,
  MTOR = if ("MTOR" %in% colnames(expr)) expr[, "MTOR"] else NA_real_,
  stringsAsFactors = FALSE
)
write.csv(out, gzfile(file.path(W, "s33_6a_cells.csv.gz")), row.names = FALSE)
A("写出 s33_6a_cells.csv.gz （", nrow(out), " 细胞）")
A("")
A("DONE_S33_6A")
