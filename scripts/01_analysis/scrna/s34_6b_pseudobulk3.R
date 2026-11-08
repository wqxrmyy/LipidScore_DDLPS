## =============================================================================
##  S34 · 6B v4：GSVA（pseudobulk）——★ 剔除评分六基因以消除循环 + 双评分对照 + BH
##  ---------------------------------------------------------------------------
##  v3 的问题：① 评分六基因（SREBF1/SCD/HMGCR/ACLY/PPARG/FASN）本就属于 PPAR/脂肪酸/脂肪生成
##              等通路基因集 → 与「通路活性 vs 评分」构成**定义性循环**；
##            ② P3/P4 的相关用**同池内**评分自比，进一步放大循环；
##            ③ 未做多重检验校正。
##  本版：① gs_full 与 **gs_clean（逐集剔除评分六基因）** 并列报告，量化循环量级；
##        ② 每个池同时对照「同池评分」与「肿瘤细胞（P1）评分」；
##        ③ BH-FDR 校正。
##  输出：_work/s34_6b_report.txt / fig6_data/pathway_*.csv
## =============================================================================
options(stringsAsFactors = FALSE)
suppressPackageStartupMessages({ library(Seurat); library(Matrix); library(GSVA); library(msigdbr) })

W <- file.path(Sys.getenv("LPS_WORK_ROOT", "."), "_work")
FIG <- file.path(W, "p5_figures/fig6_data"); dir.create(FIG, showWarnings = FALSE)
LOG <- file.path(W, "s34_6b_report.txt"); if (file.exists(LOG)) file.remove(LOG)
A <- function(...) { s <- paste0(...); cat(s, "\n", sep = ""); cat(s, "\n", sep = "", file = LOG, append = TRUE) }

SCORE6 <- c("SREBF1", "SCD", "HMGCR", "ACLY", "PPARG", "FASN")

A("===============================================================================")
A("S34 6B v4  pathway activity (GSVA pseudobulk, overlap-stripped)   ", format(Sys.time(), "%Y-%m-%d %H:%M:%S"))
A("===============================================================================")

obj <- readRDS(file.path(W, "s39_seurat.s6.rds")); md <- obj[[]]
A("loaded: ", paste(dim(obj), collapse = " x "))

## A-PPARG（逐基因跨细胞 z）
G5 <- intersect(setdiff(SCORE6, "PPARG"), rownames(obj))
M5 <- as.matrix(GetAssayData(obj, assay = "RNA", layer = "data")[G5, , drop = FALSE])
md$score_noPPARG <- as.numeric(Matrix::colMeans(t(scale(t(M5)))))
A(sprintf("score_noPPARG: sd=%.4f  rho vs AUCell=%.4f", sd(md$score_noPPARG),
          cor(md$score_noPPARG, md$score_aucell, method = "spearman")))
A("")

## 通路集：full 与 clean（剔除评分六基因）
h <- msigdbr(species = "Homo sapiens", collection = "H")
k <- msigdbr(species = "Homo sapiens", collection = "C2", subcollection = "CP:KEGG_LEGACY")
key <- grep("MTOR|PPAR|FATTY_ACID|OXIDATIVE|ADIPOCYTO|LIPID|INSULIN|GLYCOLYSIS", unique(k$gs_name), value = TRUE)
kk <- k[k$gs_name %in% key, ]
gs <- split(c(h$gene_symbol, kk$gene_symbol), c(h$gs_name, kk$gs_name))
gs <- lapply(gs, unique)
gs_clean <- lapply(gs, function(z) setdiff(z, SCORE6))
gs_clean <- gs_clean[sapply(gs_clean, length) >= 5]
lost <- unique(unlist(gs)) ; kept <- unique(unlist(gs_clean))
A(sprintf("gene sets: full=%d  clean=%d（<5 基因者剔除 %d 集）| 并集 %d -> %d 基因",
          length(gs), length(gs_clean), length(gs) - length(gs_clean), length(lost), length(kept)))
A("")

isT <- md$L1 == "Tumor cells"
pools <- list(
  P1_tumor_all  = which(isT),
  P2_tumor_prog = which(isT & md$L2 %in% c("Invasion", "ECM remodeling", "Stemness", "Hypoxia", "Angiogenesis")),
  P3_adipo_diff = which(md$L2 == "Adipo diff"),
  P4_normal_adi = which(md$L1 == "Normal adipocytes")
)
smp <- md$sample
present <- intersect(kept, rownames(obj))
expr <- GetAssayData(obj, assay = "RNA", layer = "data")[present, , drop = FALSE]

## 肿瘤细胞（P1）评分：所有池共用的「肿瘤端评分」
tkeys <- unique(smp[pools$P1_tumor_all])
tum_score <- sapply(tkeys, function(kk2) mean(md$score_aucell[pools$P1_tumor_all][smp[pools$P1_tumor_all] == kk2]))
tum_score_np <- sapply(tkeys, function(kk2) mean(md$score_noPPARG[pools$P1_tumor_all][smp[pools$P1_tumor_all] == kk2]))
names(tum_score) <- names(tum_score_np) <- tkeys

res_list <- list()
for (pn in names(pools)) {
  ii <- pools[[pn]]; keys <- unique(smp[ii])
  M <- sapply(keys, function(kk2) Matrix::rowMeans(expr[, ii[smp[ii] == kk2], drop = FALSE]))
  nc <- sapply(keys, function(kk2) sum(smp[ii] == kk2))
  keep <- nc >= 5
  A(sprintf("---- %s : cells=%d pseudobulk=%d (>=5 cells: %d) ----", pn, length(ii), length(keys), sum(keep)))
  g_full <- gsva(gsvaParam(M, gs[names(gs) %in% names(gs_clean)]), verbose = FALSE)
  g_cln  <- gsva(gsvaParam(M, gs_clean), verbose = FALSE)
  own_score <- sapply(keys, function(kk2) mean(md$score_aucell[ii][smp[ii] == kk2]))
  own_np    <- sapply(keys, function(kk2) mean(md$score_noPPARG[ii][smp[ii] == kk2]))
  for (gi in rownames(g_cln)) {
    v  <- as.numeric(g_cln[gi, ]); vf <- as.numeric(g_full[gi, ])
    ok <- keep & is.finite(v) & is.finite(vf)
    if (sum(ok) < 8) next
    ts <- tum_score[keys[ok]]; tn <- tum_score_np[keys[ok]]
    c1 <- suppressWarnings(cor.test(ts, v[ok], method = "spearman"))
    c2 <- suppressWarnings(cor.test(tn, v[ok], method = "spearman"))
    c3 <- suppressWarnings(cor.test(own_score[ok], v[ok], method = "spearman"))
    c4 <- suppressWarnings(cor.test(ts, vf[ok], method = "spearman"))
    res_list[[length(res_list) + 1]] <- data.frame(
      pool = pn, pathway = gi, n = sum(ok),
      rho_tum       = unname(c1$estimate), P_tum = c1$p.value,
      rho_tum_noPP  = unname(c2$estimate), P_tum_noPP = c2$p.value,
      rho_own       = unname(c3$estimate), P_own = c3$p.value,
      rho_tum_FULL  = unname(c4$estimate), P_tum_FULL = c4$p.value)
  }
}
tab <- do.call(rbind, res_list)
tab$q_tum <- p.adjust(tab$P_tum, "BH")
tab$q_tum_noPP <- p.adjust(tab$P_tum_noPP, "BH")

## 循环量级：同一通路 full vs clean 的 rho 差
cc <- cor(tab$rho_tum, tab$rho_tum_FULL, method = "spearman")
A("")
A("===============================================================================")
A("-- A. 循环量级量化：剔除评分六基因前(full) vs 后(clean) 的同一通路 rho --")
A("===============================================================================")
A(sprintf("   rho(clean) 与 rho(full) 的 Spearman = %.3f  |  平均 |rho| full=%.3f -> clean=%.3f  |  平均降幅 %.3f",
          cc, mean(abs(tab$rho_tum_FULL)), mean(abs(tab$rho_tum)),
          mean(abs(tab$rho_tum_FULL)) - mean(abs(tab$rho_tum))))

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
A("-- D. A-PPARG 敏感性（rho_tum_noPP，clean 集）--")
o2 <- order(-abs(tab$rho_tum_noPP))
for (i in o2[seq_len(min(15, nrow(tab)))])
  A(sprintf("%-14s %-46s %3d %+9.4f %10.3g %9s",
            tab$pool[i], tab$pathway[i], tab$n[i], tab$rho_tum_noPP[i], tab$P_tum_noPP[i], ""))

write.csv(tab, file.path(FIG, "pathway_vs_score.csv"), row.names = FALSE)
A("")
A("wrote: pathway_vs_score.csv")
A("DONE")
