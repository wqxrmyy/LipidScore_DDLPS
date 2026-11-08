## =============================================================================
##  S33 · 合并补跑：6A-补（CD36 共定位 + 导出）+ 6D-补（样本级共变）
##  ---------------------------------------------------------------------------
##  原因：① 6a2 在 CD36 段因并列值崩溃，6B 输入未导出；
##        ② 6D 的 UMAP 邻域富集观察值 = 0（肿瘤细胞在嵌入空间与巨噬不相邻）
##           → 按「不降维」改为**样本级共变**（更有生物学意义的互锁检验）
##  输入：s39_seurat.s6.rds
##  输出：fig6_data/expr_tumor.rds、cell_meta_tumor.csv；
##        _work/s33_6a3_report.txt、s33_6d2_sample_coupling.csv
## =============================================================================
options(stringsAsFactors = FALSE)
suppressPackageStartupMessages({ library(Seurat); library(Matrix) })

W <- file.path(Sys.getenv("LPS_WORK_ROOT", "."), "_work")
D <- file.path(W, "p5_figures/fig6_data"); dir.create(D, showWarnings = FALSE)
LOG <- file.path(W, "s33_6a3_report.txt"); if (file.exists(LOG)) file.remove(LOG)
A <- function(...) { s <- paste0(...); cat(s, "\n", sep = ""); cat(s, "\n", sep = "", file = LOG, append = TRUE) }

A("===============================================================================")
A("S33 合并补跑  6A-补 + 6D-补   ", format(Sys.time(), "%Y-%m-%d %H:%M:%S"))
A("===============================================================================")
obj <- readRDS(file.path(W, "s39_seurat.s6.rds"))
md <- obj[[]]
A("载入: ", paste(dim(obj), collapse = " x "))
A("")

G <- intersect(c("RICTOR", "MTOR", "CD36", "PPARG", "FABP4",
                 "SREBF1", "SCD", "HMGCR", "ACLY", "FASN",
                 "MRC1", "ARG1", "IL10", "CD163", "SPP1", "APOE", "TREM2", "LYZ"), rownames(obj))
X <- GetAssayData(obj, assay = "RNA", layer = "data")[G, , drop = FALSE]

## ---------- z 评分（逐基因跨细胞）----------
G5 <- intersect(c("SREBF1", "SCD", "HMGCR", "ACLY", "FASN"), rownames(obj))
m5 <- as.matrix(GetAssayData(obj, assay = "RNA", layer = "data")[G5, , drop = FALSE])
z5 <- t(scale(t(m5)))
md$score_noPPARG <- as.numeric(colMeans(z5, na.rm = TRUE))
A("A−PPARG: mean=", round(mean(md$score_noPPARG), 4), " sd=", round(sd(md$score_noPPARG), 4),
  " 与 AUCell rho=", round(suppressWarnings(cor(md$score_noPPARG, md$score_aucell, method = "spearman")), 4))
rm(m5, z5); gc(verbose = FALSE)

## ---------- M2 样打分（MRC1/ARG1/IL10；Q11）----------
m2g <- intersect(c("MRC1", "ARG1", "IL10"), rownames(obj))
mm <- as.matrix(GetAssayData(obj, assay = "RNA", layer = "data")[m2g, , drop = FALSE])
zm <- t(scale(t(mm)))
md$m2like <- as.numeric(colMeans(zm, na.rm = TRUE))
rm(mm, zm); gc(verbose = FALSE)

## =========================================================== (A) CD36 共定位
A("===============================================================================")
A("-- [A] CD36 共定位（检出细胞内四格；承接 6a2 未完成项）--")
A("===============================================================================")
tum <- which(md$L1 == "Tumor cells")
E <- as.matrix(t(X[, tum, drop = FALSE]))
for (g in c("RICTOR", "CD36")) {
  v <- E[, g]; pos <- v > 0
  A(sprintf("---- %s ----", g))
  A(sprintf("  肿瘤细胞内检出 %d / %d (%.2f%%)", sum(pos), length(pos), 100 * mean(pos)))
  if (sum(pos) < 200) { A("  过小，跳过"); next }
  sv <- md$score_aucell[tum][pos]
  q <- quantile(sv, .75); qg <- quantile(v[pos], .75)
  hiS <- sv >= q; hiG <- v[pos] >= qg
  A(sprintf("  检出细胞内 评分高 %d (%.1f%%) | %s高 %d (%.1f%%) | 交集 %d (%.1f%%) 观察/期望=%.3f",
            sum(hiS), 100 * mean(hiS), g, sum(hiG), 100 * mean(hiG),
            sum(hiS & hiG), 100 * mean(hiS & hiG), mean(hiS & hiG) / (mean(hiS) * mean(hiG))))
  ft <- fisher.test(table(hiS, hiG))
  A(sprintf("  Fisher OR = %.3f  P = %.3g", ft$estimate, ft$p.value))
  # 剂量关系（秩三分位；分层基准与 pos 同为「全部肿瘤细胞」）
  sv_all <- md$score_aucell[tum]
  rk <- rank(sv_all, ties.method = "first")
  tt <- cut(rk, breaks = c(0, length(rk) / 3, 2 * length(rk) / 3, length(rk)),
            labels = c("low", "mid", "high"), include.lowest = TRUE)
  dr <- tapply(pos, tt, mean)
  A("  评分三分位下的检出率（分母=全部肿瘤细胞）: ",
    paste(sprintf("%s=%.1f%%", names(dr), 100 * dr), collapse = "  "))
  A("")
}

## =========================================================== (B) 导出 6B 输入
A("===============================================================================")
A("-- [B] 导出 6B 输入 --")
A("===============================================================================")
Xt <- X[, tum, drop = FALSE]
saveRDS(Xt, file.path(D, "expr_tumor.rds"))
cm <- data.frame(barcode = colnames(Xt), patient = md$patient[tum], tissue = md$tissue[tum],
                 sample = paste(md$patient[tum], md$tissue[tum], sep = "_"),
                 L1 = md$L1[tum], L2 = md$L2[tum], L3 = md$L3[tum],
                 score_aucell = md$score_aucell[tum], score_ams = md$score_ams[tum],
                 score_noPPARG = md$score_noPPARG[tum], m2like = md$m2like[tum],
                 bucket = md$bucket[tum], stringsAsFactors = FALSE)
write.csv(cm, file.path(D, "cell_meta_tumor.csv"), row.names = FALSE, fileEncoding = "UTF-8-BOM")
A("expr_tumor.rds: ", nrow(Xt), " 基因 x ", ncol(Xt), " 细胞（log1p CP10K）")
A("cell_meta_tumor.csv: ", nrow(cm), " 行（含 score_noPPARG / m2like）")
A("")

## =========================================================== (C) 6D-补：样本级共变
A("===============================================================================")
A("-- [C] 6D-补：样本级共变（肿瘤细胞评分 vs 巨噬细胞状态）--")
A("===============================================================================")
A("  说明：6D 的 UMAP 邻域富集观察值=0（肿瘤细胞与巨噬在嵌入空间不相邻）→ 几何邻接不可用；")
A("        改测『同一取样内 肿瘤细胞评分均值』与『巨噬细胞 M2 样打分/CD36 表达均值』的共变。")
A("        这是**取样级表达关联**，不是空间邻接（措辞纪律：全文禁用 spatial）。")
A("")
key <- paste(md$patient, md$tissue, sep = "_")
smp <- sort(unique(key))
rows <- list()
for (k in smp) {
  ii <- which(key == k)
  ti <- ii[md$L1[ii] == "Tumor cells"]
  mi <- ii[md$L2[ii] == "Macrophages"]
  if (length(ti) < 20 || length(mi) < 5) next
  rows[[length(rows) + 1]] <- data.frame(
    sample = k, n_tumor = length(ti), n_mac = length(mi),
    tumor_score = mean(md$score_aucell[ti]),
    tumor_score_noPP = mean(md$score_noPPARG[ti]),
    mac_m2like = mean(md$m2like[mi]),
    mac_CD36 = mean(as.matrix(X["CD36", mi])),
    mac_RICTOR = mean(as.matrix(X["RICTOR", mi])),
    mac_PPARG = mean(as.matrix(X["PPARG", mi])),
    stringsAsFactors = FALSE)
}
CP <- do.call(rbind, rows)
A("  可用取样数 = ", nrow(CP), "（要求 肿瘤>=20 细胞且 巨噬>=5 细胞）")
A("")
A(sprintf("  %-14s %7s %7s %11s %11s %11s", "sample", "n_tum", "n_mac", "tum_score", "mac_M2like", "mac_CD36"))
for (i in seq_len(nrow(CP)))
  A(sprintf("  %-14s %7d %7d %11.4f %11.4f %11.4f", CP$sample[i], CP$n_tumor[i], CP$n_mac[i],
            CP$tumor_score[i], CP$mac_m2like[i], CP$mac_CD36[i]))
A("")
pairs_c <- list(
  c("tumor_score", "mac_m2like"), c("tumor_score", "mac_CD36"),
  c("tumor_score", "mac_RICTOR"), c("tumor_score", "mac_PPARG"),
  c("tumor_score_noPP", "mac_m2like"), c("tumor_score_noPP", "mac_CD36"))
A("  -- 取样级 Spearman 相关（n = ", nrow(CP), "）--")
A(sprintf("  %-32s %10s %11s", "pair", "rho", "P"))
cpl <- list()
for (p in pairs_c) {
  r <- suppressWarnings(cor(CP[[p[1]]], CP[[p[2]]], method = "spearman"))
  pv <- suppressWarnings(cor.test(CP[[p[1]]], CP[[p[2]]], method = "spearman")$p.value)
  A(sprintf("  %-32s %10.4f %11.3g", paste(p, collapse = " vs "), r, pv))
  cpl[[length(cpl) + 1]] <- data.frame(x = p[1], y = p[2], rho = r, p = pv, n = nrow(CP))
}
write.csv(CP, file.path(W, "s33_6d2_sample_coupling.csv"), row.names = FALSE)
write.csv(do.call(rbind, cpl), file.path(W, "s33_6d2_coupling_corr.csv"), row.names = FALSE)
A("")
A("写出 s33_6d2_sample_coupling.csv / s33_6d2_coupling_corr.csv")
A("")
A("DONE_S33_6A3")
