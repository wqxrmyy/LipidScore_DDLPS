## =============================================================================
##  S34 · 6B v5：GSVA（pseudobulk）——★ 修正 v4 的「空转对照」缺陷
##  ---------------------------------------------------------------------------
##  v4 缺陷（自证，非他证）：
##    表达矩阵 present <- intersect(kept, rownames(obj)) 只保留「剔除六基因后」
##    的 4650 基因 → full 基因集里那 6 个基因在矩阵中根本不存在 →
##    g_full 与 g_cln 在数学上恒等 → v4 报告的
##        "Spearman=1.000 | 平均|rho| 0.234 -> 0.234 | 平均降幅 0.000"
##    是**空转对照（no-op）**，不能作为「循环量级=0」的证据。
##  v5 修正：
##    ① 矩阵保留完整并集 present <- intersect(lost_union, rownames(obj))
##       （4656，含评分六基因）→ full / clean 真正可比；
##    ② 新增 E 部分：单列「含评分六基因的基因集」的 full -> clean 变化
##       = 循环落点；并报 |delta| 最大/中位；
##    ③ clean 集读数与 v4 对照（矩阵多 6 个不在 clean 集内的基因，
##       GSVA 逐集只取集内基因 → clean 读数应逐位复现 v4，可作交叉验证）。
##  输出：_work/s34_6b_report.txt / fig6_data/pathway_vs_score.csv
## =============================================================================
options(stringsAsFactors = FALSE)
suppressPackageStartupMessages({ library(Seurat); library(Matrix); library(GSVA); library(msigdbr) })

W <- file.path(Sys.getenv("LPS_WORK_ROOT", "."), "_work")
FIG <- file.path(W, "p5_figures/fig6_data"); dir.create(FIG, showWarnings = FALSE)
LOG <- file.path(W, "s34_6b_report.txt"); if (file.exists(LOG)) file.remove(LOG)
A <- function(...) { s <- paste0(...); cat(s, "\n", sep = ""); cat(s, "\n", sep = "", file = LOG, append = TRUE) }

SCORE6 <- c("SREBF1", "SCD", "HMGCR", "ACLY", "PPARG", "FASN")

A("===============================================================================")
A("S34 6B v5  pathway activity (GSVA pseudobulk) — v4 no-op ctrl FIXED   ", format(Sys.time(), "%Y-%m-%d %H:%M:%S"))
A("===============================================================================")

obj <- readRDS(file.path(W, "s39_seurat.s6.rds")); md <- obj[[]]
A("loaded: ", paste(dim(obj), collapse = " x "))

## ---- A-PPARG（逐基因跨细胞 z）---------------------------------------------
G5 <- intersect(setdiff(SCORE6, "PPARG"), rownames(obj))
M5 <- as.matrix(GetAssayData(obj, assay = "RNA", layer = "data")[G5, , drop = FALSE])
md$score_noPPARG <- as.numeric(Matrix::colMeans(t(scale(t(M5)))))
A(sprintf("score_noPPARG: genes=%d  sd=%.4f  rho vs AUCell=%.4f",
          length(G5), sd(md$score_noPPARG), cor(md$score_noPPARG, md$score_aucell, method = "spearman")))
A("")

## ---- 通路集：full / clean（同一批通路，逐集剔六基因）-----------------------
h <- msigdbr(species = "Homo sapiens", collection = "H")
k <- msigdbr(species = "Homo sapiens", collection = "C2", subcollection = "CP:KEGG_LEGACY")
key <- grep("MTOR|PPAR|FATTY_ACID|OXIDATIVE|ADIPOCYTO|LIPID|INSULIN|GLYCOLYSIS", unique(k$gs_name), value = TRUE)
kk <- k[k$gs_name %in% key, ]
gs <- split(c(h$gene_symbol, kk$gene_symbol), c(h$gs_name, kk$gs_name))
gs <- lapply(gs, unique)
gs_clean_all <- lapply(gs, function(z) setdiff(z, SCORE6))
keepset <- sapply(gs_clean_all, length) >= 5
gs <- gs[keepset]; gs_clean <- gs_clean_all[keepset]
has6 <- sapply(gs, function(z) any(SCORE6 %in% z))
A(sprintf("gene sets: %d（其中含评分六基因者 %d 集）", length(gs), sum(has6)))
for (nm in names(gs)[has6]) A(sprintf("   [%s] 含 %s", nm, paste(intersect(gs[[nm]], SCORE6), collapse = ",")))

lost_union <- unique(unlist(gs)); kept_union <- unique(unlist(gs_clean))
present <- intersect(lost_union, rownames(obj))          ## ★ 4656：保留六基因
A(sprintf("matrix genes: full-union=%d  clean-union=%d  仅full有=%s  -> used=%d",
          length(lost_union), length(kept_union),
          paste(setdiff(lost_union, kept_union), collapse = ","), length(present)))
A("")

## ---- 池定义 ---------------------------------------------------------------
isT <- md$L1 == "Tumor cells"
pools <- list(
  P1_tumor_all  = which(isT),
  P2_tumor_prog = which(isT & md$L2 %in% c("Invasion", "ECM remodeling", "Stemness", "Hypoxia", "Angiogenesis")),
  P3_adipo_diff = which(md$L2 == "Adipo diff"),
  P4_normal_adi = which(md$L1 == "Normal adipocytes")
)
smp <- md$sample
expr <- GetAssayData(obj, assay = "RNA", layer = "data")[present, , drop = FALSE]

tkeys <- unique(smp[pools$P1_tumor_all])
tum_score <- sapply(tkeys, function(z) mean(md$score_aucell[pools$P1_tumor_all][smp[pools$P1_tumor_all] == z]))
tum_score_np <- sapply(tkeys, function(z) mean(md$score_noPPARG[pools$P1_tumor_all][smp[pools$P1_tumor_all] == z]))
names(tum_score) <- names(tum_score_np) <- tkeys

res_list <- list()
for (pn in names(pools)) {
  ii <- pools[[pn]]; keys <- unique(smp[ii])
  M <- sapply(keys, function(z) Matrix::rowMeans(expr[, ii[smp[ii] == z], drop = FALSE]))
  nc <- sapply(keys, function(z) sum(smp[ii] == z))
  keep <- nc >= 5
  A(sprintf("---- %s : cells=%d pseudobulk=%d (>=5 cells: %d) ----", pn, length(ii), length(keys), sum(keep)))
  g_full <- gsva(gsvaParam(M, gs),       verbose = FALSE)   ## ★ 真 full：矩阵含六基因
  g_cln  <- gsva(gsvaParam(M, gs_clean), verbose = FALSE)
  own_score <- sapply(keys, function(z) mean(md$score_aucell[ii][smp[ii] == z]))
  for (gi in rownames(g_cln)) {
    v <- as.numeric(g_cln[gi, ]); vf <- as.numeric(g_full[gi, ])
    ok <- keep & is.finite(v) & is.finite(vf)
    if (sum(ok) < 8) next
    ts <- tum_score[keys[ok]]; tn <- tum_score_np[keys[ok]]
    c1 <- suppressWarnings(cor.test(ts, v[ok],  method = "spearman"))
    c2 <- suppressWarnings(cor.test(tn, v[ok],  method = "spearman"))
    c3 <- suppressWarnings(cor.test(own_score[ok], v[ok], method = "spearman"))
    c4 <- suppressWarnings(cor.test(ts, vf[ok], method = "spearman"))
    res_list[[length(res_list) + 1]] <- data.frame(
      pool = pn, pathway = gi, n = sum(ok),
      rho_tum = unname(c1$estimate), P_tum = c1$p.value,
      rho_tum_noPP = unname(c2$estimate), P_tum_noPP = c2$p.value,
      rho_own = unname(c3$estimate), P_own = c3$p.value,
      rho_tum_FULL = unname(c4$estimate), P_tum_FULL = c4$p.value)
  }
}
tab <- do.call(rbind, res_list)
tab$q_tum <- p.adjust(tab$P_tum, "BH")
tab$q_tum_noPP <- p.adjust(tab$P_tum_noPP, "BH")

A("")
A("===============================================================================")
A("-- A. 循环量级：同一通路的 rho(full, 含六基因) vs rho(clean, 剔六基因) --")
A("===============================================================================")
cc <- cor(tab$rho_tum, tab$rho_tum_FULL, method = "spearman")
A(sprintf("   Spearman(rho_clean, rho_full) = %.4f | 平均|rho| full=%.4f -> clean=%.4f | 平均降幅 %.4f | 最大单集降幅 %.4f",
          cc, mean(abs(tab$rho_tum_FULL)), mean(abs(tab$rho_tum)),
          mean(abs(tab$rho_tum_FULL)) - mean(abs(tab$rho_tum)),
          max(abs(tab$rho_tum_FULL - tab$rho_tum))))

A("")
A("===============================================================================")
A("-- B. clean 通路活性 vs 肿瘤细胞评分（跨细胞室/去循环）--")
A("===============================================================================")
o <- order(-abs(tab$rho_tum))
A(sprintf("%-14s %-46s %3s %9s %10s %10s %9s", "pool", "pathway", "n", "rho_tum", "P", "q(BH)", "rho_FULL"))
for (i in o[seq_len(min(30, nrow(tab)))])
  A(sprintf("%-14s %-46s %3d %+9.4f %10.3g %10.3g %+9.4f",
            tab$pool[i], tab$pathway[i], tab$n[i], tab$rho_tum[i], tab$P_tum[i], tab$q_tum[i], tab$rho_tum_FULL[i]))

A("")
sig <- tab[tab$q_tum < 0.05, ]
A(sprintf("★ FDR<0.05 的通路：%d / %d", nrow(sig), nrow(tab)))
if (nrow(sig)) for (i in seq_len(min(25, nrow(sig))))
  A(sprintf("  * %-14s %-46s rho=%+.3f P=%.3g q=%.3g (full rho=%+.3f)",
            sig$pool[i], sig$pathway[i], sig$rho_tum[i], sig$P_tum[i], sig$q_tum[i], sig$rho_tum_FULL[i]))

A("")
A("-- C. mTOR / PPAR 相关通路专项 --")
sel <- grep("MTOR|PPAR|FATTY_ACID|OXIDATIVE|ADIPOCYTO", tab$pathway)
ss <- tab[sel, ]; ss <- ss[order(ss$P_tum), ]
A(sprintf("%-14s %-46s %3s %9s %10s %10s %9s", "pool", "pathway", "n", "rho_tum", "P", "q(BH)", "rho_FULL"))
for (i in seq_len(min(24, nrow(ss))))
  A(sprintf("%-14s %-46s %3d %+9.4f %10.3g %10.3g %+9.4f",
            ss$pool[i], ss$pathway[i], ss$n[i], ss$rho_tum[i], ss$P_tum[i], ss$q_tum[i], ss$rho_tum_FULL[i]))

A("")
A("-- D. A-PPARG 敏感性（x 轴换成不含 PPARG 的评分，clean 集）--")
o2 <- order(-abs(tab$rho_tum_noPP))
A(sprintf("%-14s %-46s %3s %11s %10s", "pool", "pathway", "n", "rho_noPP", "P"))
for (i in o2[seq_len(min(15, nrow(tab)))])
  A(sprintf("%-14s %-46s %3d %+11.4f %10.3g",
            tab$pool[i], tab$pathway[i], tab$n[i], tab$rho_tum_noPP[i], tab$P_tum_noPP[i]))

A("")
A("-- D2. PPAR/脂肪酸三通路：含 PPARG 评分 vs 去 PPARG 评分 --")
p3 <- c("KEGG_PPAR_SIGNALING_PATHWAY", "KEGG_FATTY_ACID_METABOLISM",
        "KEGG_BIOSYNTHESIS_OF_UNSATURATED_FATTY_ACIDS", "HALLMARK_FATTY_ACID_METABOLISM")
dd <- tab[tab$pathway %in% p3, ]
A(sprintf("%-14s %-46s %3s %11s %11s %9s", "pool", "pathway", "n", "rho(A)", "rho(A-PPARG)", "delta"))
for (i in seq_len(nrow(dd)))
  A(sprintf("%-14s %-46s %3d %+11.4f %+11.4f %+9.4f",
            dd$pool[i], dd$pathway[i], dd$n[i], dd$rho_tum[i], dd$rho_tum_noPP[i],
            dd$rho_tum_noPP[i] - dd$rho_tum[i]))

A("")
A("-- E. ★ 循环落点：含评分六基因的基因集  full -> clean --")
h6 <- tab[tab$pathway %in% names(gs)[has6], ]
h6$d_rho <- h6$rho_tum - h6$rho_tum_FULL
h6 <- h6[order(-abs(h6$d_rho)), ]
A(sprintf("%-14s %-46s %3s %9s %9s %9s", "pool", "pathway", "n", "rho_full", "rho_clean", "delta"))
for (i in seq_len(min(40, nrow(h6))))
  A(sprintf("%-14s %-46s %3d %+9.4f %+9.4f %+9.4f",
            h6$pool[i], h6$pathway[i], h6$n[i], h6$rho_tum_FULL[i], h6$rho_tum[i], h6$d_rho[i]))
A(sprintf("   |delta| 最大 = %.4f   中位 = %.4f   （含六基因集 %d 行 / 全表 %d 行）",
          max(abs(h6$d_rho)), median(abs(h6$d_rho)), nrow(h6), nrow(tab)))

write.csv(tab, file.path(FIG, "pathway_vs_score.csv"), row.names = FALSE)
A("")
A("wrote: pathway_vs_score.csv")
A("DONE")
