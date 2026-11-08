## =============================================================================
##  S34 · 6B v2：GSVA 通路活性（pseudobulk 口径）
##  ---------------------------------------------------------------------------
##  起因：s33_6b.R 在 GSVA 阶段内存耗尽被杀（4418 x 36592 密化 >10 GB，无报错退出）
##  解法（不降维）：改为 **样本级 pseudobulk** 后做 GSVA
##        · 统计上更正确：36,592 细胞当独立样本是伪重复（同一患者的细胞不独立）
##        · 方法未换：仍是 GSVA（kcdf=Gaussian, maxDiff=TRUE），仅输入由「细胞」改「患者」
##        · 四口径沿用 6a 的细胞池定义；另报 A−PPARG 敏感性
##  输出：_work/s34_6b_report.txt / fig6_data/pathway_pseudobulk.csv
## =============================================================================
options(stringsAsFactors = FALSE)
suppressPackageStartupMessages({ library(Seurat); library(Matrix); library(GSVA); library(msigdbr) })

W <- file.path(Sys.getenv("LPS_WORK_ROOT", "."), "_work")
FIG <- file.path(W, "p5_figures/fig6_data"); dir.create(FIG, showWarnings = FALSE)
LOG <- file.path(W, "s34_6b_report.txt"); if (file.exists(LOG)) file.remove(LOG)
A <- function(...) { s <- paste0(...); cat(s, "\n", sep = ""); cat(s, "\n", sep = "", file = LOG, append = TRUE) }

A("===============================================================================")
A("S34 6B v2  pathway activity (GSVA, pseudobulk)   ", format(Sys.time(), "%Y-%m-%d %H:%M:%S"))
A("===============================================================================")

obj <- readRDS(file.path(W, "s39_seurat.s6.rds"))
md  <- obj[[]]
A("loaded: ", paste(dim(obj), collapse = " x "))
A("meta cols: ", paste(colnames(md), collapse = ", "))
A("")

## ---- 通路集：Hallmark 50 + KEGG 主题 ----
h <- msigdbr(species = "Homo sapiens", collection = "H")
k <- msigdbr(species = "Homo sapiens", collection = "C2", subcollection = "CP:KEGG_LEGACY")
key <- grep("MTOR|PPAR|FATTY_ACID|OXIDATIVE|ADIPOCYTO|LIPID|INSULIN|GLYCOLYSIS", unique(k$gs_name), value = TRUE)
kk <- k[k$gs_name %in% key, ]
gs <- split(c(h$gene_symbol, kk$gene_symbol), c(h$gs_name, kk$gs_name))
gs <- lapply(gs, unique)
A("pathway sets: Hallmark ", length(unique(h$gs_name)), " + KEGG(theme) ", length(key), " = ", length(gs))
geneset <- unique(unlist(gs))
A("union genes in sets: ", length(geneset))
present <- intersect(geneset, rownames(obj))
A("present in object: ", length(present))
A("")

expr <- GetAssayData(obj, assay = "RNA", layer = "data")[present, , drop = FALSE]

## ---- 细胞池（沿用 6a 定义）----
isT  <- md$L1 == "Tumor cells"
pools <- list(
  tumor_all  = which(isT),
  tumor_noAd = which(isT & md$L2 != "Adipo diff"),
  tumor_prog = which(isT & md$L2 %in% c("Invasion", "Hypoxia", "Stemness", "ECM remodeling",
                                        "Angiogenesis", "Cycling cells")),
  tumor_nonorm = which(isT & md$L2 != "Normal adipocytes")
)
smp <- md$sample
if (is.null(smp)) smp <- md$orig.ident
A("sample column used: ", if (is.null(md$sample)) "orig.ident" else "sample",
  " | levels = ", length(unique(smp)))
A("")

pb_list <- list()
for (pn in names(pools)) {
  ii <- pools[[pn]]
  keys <- unique(smp[ii])
  M <- sapply(keys, function(kk2) Matrix::rowMeans(expr[, ii[smp[ii] == kk2], drop = FALSE]))
  pb_list[[pn]] <- M
  A(sprintf("  pool %-13s cells=%6d  pseudobulk samples=%3d", pn, length(ii), ncol(M)))
}
A("")

## ---- 每 pseudobulk 样本的肿瘤评分（表征该样本肿瘤细胞的评分水位）----
scol <- if ("score_aucell" %in% colnames(md)) "score_aucell" else "score_ams"
scol2 <- "score_noPPARG"
A("score column: ", scol, " | sensitivity: ", scol2)
A("")

out <- list(); meta <- list()
for (pn in names(pb_list)) {
  M <- pb_list[[pn]]
  A("")
  A(sprintf("---- GSVA: %s (genes=%d, samples=%d) ----", pn, nrow(M), ncol(M)))
  t0 <- Sys.time()
  gsva_par <- gsvaParam(M, gs)
  res <- gsva(gsva_par, verbose = FALSE)
  A(sprintf("  done in %.1f s ; dim = %s", as.numeric(difftime(Sys.time(), t0, units = "secs")),
            paste(dim(res), collapse = " x ")))
  out[[pn]] <- res
  keys <- colnames(M)
  meta[[pn]] <- data.frame(
    pseudobulk = paste(pn, keys, sep = "|"),
    pool = pn, sample = keys,
    n_cells = sapply(keys, function(kk2) sum(smp[pools[[pn]]] == kk2)),
    score  = sapply(keys, function(kk2) mean(md[[scol]][pools[[pn]]][smp[pools[[pn]]] == kk2])),
    score_noPP = sapply(keys, function(kk2) mean(md[[scol2]][pools[[pn]]][smp[pools[[pn]]] == kk2])),
    row.names = NULL)
}

## ---- 通路活性 vs 肿瘤评分（样本级相关）----
A("")
A("===============================================================================")
A("-- pathway activity vs tumor-cell score (sample level, per pool) --")
A("===============================================================================")
allm <- do.call(rbind, meta)
tab <- data.frame()
for (pn in names(out)) {
  R <- out[[pn]]; m <- allm[allm$pool == pn, ]
  idx <- match(m$sample, colnames(R))
  ok <- !is.na(idx)
  for (gi in rownames(R)) {
    v <- as.numeric(R[gi, idx[ok]]); s <- m$score[ok]
    if (length(v) < 6 || sd(v) == 0) next
    ct <- suppressWarnings(cor.test(s, v, method = "spearman"))
    sp <- suppressWarnings(cor.test(m$score_noPP[ok], v, method = "spearman"))
    tab <- rbind(tab, data.frame(pool = pn, pathway = gi, n = length(v),
                                 rho = unname(ct$estimate), P = ct$p.value,
                                 rho_noPPARG = unname(sp$estimate), P_noPPARG = sp$p.value))
  }
}
KEY <- grep("MTOR|PPAR|FATTY_ACID|OXIDATIVE|ADIPOCYTO|LIPID|GLYCOLYSIS|INSULIN", tab$pathway)
tab2 <- tab[KEY, ]; tab2 <- tab2[order(-abs(tab2$rho)), ]
A(sprintf("%-13s %-46s %4s %8s %10s %10s %10s", "pool", "pathway", "n", "rho", "P", "rho_noPP", "P_noPP"))
for (i in seq_len(min(40, nrow(tab2))))
  A(sprintf("%-13s %-46s %4d %+8.4f %10.3g %+10.4f %10.3g",
            tab2$pool[i], tab2$pathway[i], tab2$n[i], tab2$rho[i], tab2$P[i],
            tab2$rho_noPPARG[i], tab2$P_noPPARG[i]))

w <- do.call(rbind, lapply(names(out), function(pn)
  data.frame(pool = pn, pathway = rownames(out[[pn]]), out[[pn]], check.names = FALSE)))
write.csv(w, file.path(FIG, "pathway_pseudobulk.csv"), row.names = FALSE)
write.csv(allm, file.path(FIG, "pathway_pseudobulk_meta.csv"), row.names = FALSE)
write.csv(tab, file.path(FIG, "pathway_vs_score.csv"), row.names = FALSE)
A("")
A("wrote: pathway_pseudobulk.csv / pathway_pseudobulk_meta.csv / pathway_vs_score.csv")
A("DONE")
