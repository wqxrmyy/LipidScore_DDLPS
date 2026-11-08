## =============================================================================
##  S33 · Fig6 6B（重写版）：高/低评分细胞通路活性 + 6D-补（样本级共变）+ 导出
##  ---------------------------------------------------------------------------
##  修正前版三处问题：① 导出仅 15 基因（GSVA 需通路基因并集）；② write.csv 不支持 BOM；
##                    ③ 未做样本级共变（6D 邻域富集观察值=0 → 几何邻接不可用）
##  内存策略：GSVA 前把矩阵子集到「通路基因并集」（避免 49632x36592 dense 爆炸）
##  口径：Q4 四口径；Q3 双轨（含 A−PPARG）；通路 = MSigDB Hallmark + KEGG 主题
##  输出：fig6_data/cell_meta_tumor.csv、fig6_data/panelB_pathway_*.csv、
##        _work/s33_6b_report.txt、s33_6d2_sample_coupling.csv
## =============================================================================
options(stringsAsFactors = FALSE)
suppressPackageStartupMessages({ library(Seurat); library(Matrix); library(GSVA); library(msigdbr); library(data.table) })

W <- file.path(Sys.getenv("LPS_WORK_ROOT", "."), "_work")
D <- file.path(W, "p5_figures/fig6_data"); dir.create(D, showWarnings = FALSE)
LOG <- file.path(W, "s33_6b_report.txt"); if (file.exists(LOG)) file.remove(LOG)
A <- function(...) { s <- paste0(...); cat(s, "\n", sep = ""); cat(s, "\n", sep = "", file = LOG, append = TRUE) }

A("===============================================================================")
A("S33 6B+6D补  通路活性 + 样本级共变   ", format(Sys.time(), "%Y-%m-%d %H:%M:%S"))
A("===============================================================================")
A("GSVA ", as.character(packageVersion("GSVA")), " | msigdbr ", as.character(packageVersion("msigdbr")))
A("")

## ---------- 1. 加载与子集 ----------
t0 <- Sys.time()
obj <- readRDS(file.path(W, "s39_seurat.s6.rds"))
md <- obj[[]]
tum <- which(md$L1 == "Tumor cells")
A("载入: ", paste(dim(obj), collapse = " x "), " | 肿瘤细胞 ", length(tum), "  (",
  round(as.numeric(difftime(Sys.time(), t0, units = "secs")), 1), " s)")
A("")

## 评分（A−PPARG，逐基因跨细胞 z）
G5 <- intersect(c("SREBF1", "SCD", "HMGCR", "ACLY", "FASN"), rownames(obj))
m5 <- as.matrix(GetAssayData(obj, assay = "RNA", layer = "data")[G5, , drop = FALSE])
md$score_noPPARG <- as.numeric(colMeans(t(scale(t(m5))), na.rm = TRUE))
rm(m5); gc(verbose = FALSE)
## M2 样打分
m2g <- intersect(c("MRC1", "ARG1", "IL10"), rownames(obj))
mm <- as.matrix(GetAssayData(obj, assay = "RNA", layer = "data")[m2g, , drop = FALSE])
md$m2like <- as.numeric(colMeans(t(scale(t(mm))), na.rm = TRUE))
rm(mm); gc(verbose = FALSE)
A("A−PPARG: sd=", round(sd(md$score_noPPARG), 4), " 与 AUCell rho=",
  round(suppressWarnings(cor(md$score_noPPARG, md$score_aucell, method = "spearman")), 4))
A("M2 样打分: 巨噬 mean=", round(mean(md$m2like[md$L2 == "Macrophages"]), 4),
  " 全体 mean=", round(mean(md$m2like), 4))
A("")

## 导出 cell_meta_tumor.csv（带 BOM）
cm <- data.frame(barcode = colnames(obj)[tum], patient = md$patient[tum], tissue = md$tissue[tum],
                 sample = paste(md$patient[tum], md$tissue[tum], sep = "_"),
                 L1 = md$L1[tum], L2 = md$L2[tum], L3 = md$L3[tum],
                 score_aucell = md$score_aucell[tum], score_ams = md$score_ams[tum],
                 score_noPPARG = md$score_noPPARG[tum], m2like = md$m2like[tum],
                 bucket = md$bucket[tum])
fwrite(cm, file.path(D, "cell_meta_tumor.csv"), bom = TRUE)
A("写出 cell_meta_tumor.csv（", nrow(cm), " 行，UTF-8-BOM）")
A("")

## ---------- 2. 6D-补：样本级共变 ----------
A("===============================================================================")
A("-- [6D-补] 样本级共变：肿瘤细胞评分 vs 巨噬细胞状态（几何邻接不可用后的升维替代）--")
A("===============================================================================")
key <- paste(md$patient, md$tissue, sep = "_")
Xg <- GetAssayData(obj, assay = "RNA", layer = "data")
rows <- list()
for (k in sort(unique(key))) {
  ii <- which(key == k)
  ti <- ii[md$L1[ii] == "Tumor cells"]; mi <- ii[md$L2[ii] == "Macrophages"]
  if (length(ti) < 20 || length(mi) < 5) next
  rows[[length(rows) + 1]] <- data.frame(
    sample = k, patient = md$patient[ii][1], tissue = md$tissue[ii][1],
    n_tumor = length(ti), n_mac = length(mi),
    tumor_score = mean(md$score_aucell[ti]),
    tumor_score_noPP = mean(md$score_noPPARG[ti]),
    mac_m2like = mean(md$m2like[mi]),
    mac_CD36 = mean(as.matrix(Xg["CD36", mi])),
    mac_RICTOR = mean(as.matrix(Xg["RICTOR", mi])),
    mac_PPARG = mean(as.matrix(Xg["PPARG", mi])),
    stringsAsFactors = FALSE)
}
CP <- do.call(rbind, rows)
A("  可用取样 = ", nrow(CP), "（肿瘤>=20 且 巨噬>=5 细胞）")
fwrite(CP, file.path(W, "s33_6d2_sample_coupling.csv"), bom = TRUE)
A("")
A(sprintf("  %-14s %6s %6s %11s %11s %10s %10s", "sample", "n_tum", "n_mac", "tum_score", "mac_M2like", "mac_CD36", "mac_RICTOR"))
for (i in seq_len(nrow(CP)))
  A(sprintf("  %-14s %6d %6d %11.4f %11.4f %10.4f %10.4f", CP$sample[i], CP$n_tumor[i], CP$n_mac[i],
            CP$tumor_score[i], CP$mac_m2like[i], CP$mac_CD36[i], CP$mac_RICTOR[i]))
A("")
A(sprintf("  %-34s %10s %12s", "pair (取样级 Spearman)", "rho", "P"))
cpl <- list()
for (p in list(c("tumor_score", "mac_m2like"), c("tumor_score", "mac_CD36"),
               c("tumor_score", "mac_RICTOR"), c("tumor_score", "mac_PPARG"),
               c("tumor_score_noPP", "mac_m2like"), c("tumor_score_noPP", "mac_CD36"))) {
  r <- suppressWarnings(cor(CP[[p[1]]], CP[[p[2]]], method = "spearman"))
  pv <- suppressWarnings(cor.test(CP[[p[1]]], CP[[p[2]]], method = "spearman")$p.value)
  A(sprintf("  %-34s %10.4f %12.3g", paste(p, collapse = " vs "), r, pv))
  cpl[[length(cpl) + 1]] <- data.frame(x = p[1], y = p[2], rho = r, p = pv, n = nrow(CP))
}
fwrite(do.call(rbind, cpl), file.path(W, "s33_6d2_coupling_corr.csv"), bom = TRUE)
A("")

## ---------- 3. 通路集 ----------
A("-- 通路集（Hallmark 50 + KEGG 主题）--")
hset <- split(msigdbr(species = "Homo sapiens", collection = "H")$gene_symbol,
              msigdbr(species = "Homo sapiens", collection = "H")$gs_name)
kg <- msigdbr(species = "Homo sapiens", collection = "C2", subcollection = "CP:KEGG_LEGACY")
kset <- split(kg$gene_symbol, kg$gs_name)
kset <- kset[intersect(c("KEGG_MTOR_SIGNALING_PATHWAY", "KEGG_PPAR_SIGNALING_PATHWAY",
                         "KEGG_FATTY_ACID_METABOLISM", "KEGG_OXIDATIVE_PHOSPHORYLATION",
                         "KEGG_ADIPOCYTOKINE_SIGNALING_PATHWAY"), names(kset))]
PWS <- c(hset, kset)
A("  合计 ", length(PWS), " 条（Hallmark ", length(hset), " + KEGG ", length(kset), "）")
## 内存策略：子集到通路基因并集
pw_genes <- unique(unlist(PWS))
pw_genes <- intersect(pw_genes, rownames(obj))
A("  通路基因并集（在库内）= ", length(pw_genes), " 个 → 只对该子集做 GSVA（避免 dense 爆炸）")
Xsub <- Xg[pw_genes, tum, drop = FALSE]
A("  子集矩阵: ", nrow(Xsub), " 基因 x ", ncol(Xsub), " 细胞  nnz=", format(length(Xsub@x), big.mark = ","))
A("")

## ---------- 4. GSVA 四口径 ----------
A("===============================================================================")
A("-- GSVA（kcdf=Gaussian, maxDiff=TRUE）；四口径 + A−PPARG 敏感性 --")
A("===============================================================================")
mdt <- md[tum, ]
run_pool <- function(sel, tag, sc) {
  if (sum(sel) < 100) { A(sprintf("[%s] 细胞过少(%d)", tag, sum(sel))); return(NULL) }
  A(sprintf("---- %s (n = %d) ----", tag, sum(sel)))
  t0 <- Sys.time()
  es <- tryCatch(gsva(gsvaParam(Xsub[, sel, drop = FALSE], PWS, kcdf = "Gaussian", maxDiff = TRUE),
                      verbose = FALSE),
                 error = function(e) { A("  GSVA 错误: ", conditionMessage(e)); NULL })
  if (is.null(es)) return(NULL)
  A("  耗时 ", round(as.numeric(difftime(Sys.time(), t0, units = "secs")), 1), " s")
  sv <- sc[sel]; hi <- sv > median(sv)
  A("  分组: high=%d low=%d", sum(hi), sum(!hi))
  rr <- t(apply(es, 1, function(v) {
    wt <- suppressWarnings(wilcox.test(v[hi], v[!hi]))
    c(delta = mean(v[hi]) - mean(v[!hi]), p = wt$p.value, mean_high = mean(v[hi]), mean_low = mean(v[!hi]))
  }))
  rr <- as.data.frame(rr); rr$pathway <- rownames(rr); rr$pool <- tag
  rr$n_high <- sum(hi); rr$n_low <- sum(!hi); rr$p_adj_BH <- p.adjust(rr$p, "BH")
  rr$FDR_pass <- rr$p_adj_BH < 0.05
  rr <- rr[order(rr$p), ]
  A("  ★ BH(FDR<0.05): ", sum(rr$FDR_pass), " / ", nrow(rr))
  for (i in seq_len(min(6, nrow(rr))))
    A(sprintf("     %-46s %+8.4f  P=%.2e  BH=%.2e", rr$pathway[i], rr$delta[i], rr$p[i], rr$p_adj_BH[i]))
  A("")
  rr
}
prog <- c("Invasion", "Hypoxia", "Stemness", "ECM remodeling", "Angiogenesis", "Cycling cells")
out <- list()
out[["tumor_all"]]      <- run_pool(rep(TRUE, nrow(mdt)), "tumor_all", mdt$score_aucell)
out[["tumor_noAd"]]     <- run_pool(mdt$L2 != "Adipo diff", "tumor_noAd", mdt$score_aucell)
out[["tumor_prog"]]     <- run_pool(mdt$L2 %in% prog, "tumor_prog", mdt$score_aucell)
out[["adipo_diff"]]     <- run_pool(mdt$L2 == "Adipo diff", "adipo_diff", mdt$score_aucell)
A("-- Q3 敏感性：A−PPARG 评分重算 tumor_all --")
out[["tumor_all_noPPARG"]] <- run_pool(rep(TRUE, nrow(mdt)), "tumor_all_A-PPARG", mdt$score_noPPARG)

ALL <- do.call(rbind, out)
if (!is.null(ALL)) {
  fwrite(ALL, file.path(D, "panelB_pathway_activity.csv"), bom = TRUE)
  A("写出 panelB_pathway_activity.csv（", nrow(ALL), " 行）")
  A("")
  A("-- 主题通路摘要（delta，* = BH<0.05）--")
  pls <- unique(ALL$pool)
  thm <- grep("MTOR|PPAR|FATTY_ACID|OXIDATIVE|ADIPOCYTO", ALL$pathway, value = TRUE) |> unique()
  A(sprintf("%-46s %s", "pathway", paste(sprintf("%-24s", pls), collapse = "")))
  for (p in thm) {
    cc <- sapply(pls, function(pl) {
      r <- ALL[ALL$pool == pl & ALL$pathway == p, ]
      if (!nrow(r)) return("--")
      sprintf("%+.3f(P%.1e)%s", r$delta[1], r$p[1], ifelse(r$FDR_pass[1], "*", ""))
    })
    A(sprintf("%-46s %s", p, paste(sprintf("%-24s", cc), collapse = "")))
  }
}
A("")
A("DONE_S33_6B")
