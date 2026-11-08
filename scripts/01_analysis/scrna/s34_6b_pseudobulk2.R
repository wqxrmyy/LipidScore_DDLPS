## =============================================================================
##  S34 · 6B v3：GSVA 通路活性（pseudobulk 口径，修正池定义 + 自算 A-PPARG）
##  ---------------------------------------------------------------------------
##  v2 问题：① score_noPPARG 列在对象中不存在（6a2 仅内存计算）→ 此处按 D18 口径
##           re-compute：逐基因跨细胞 z（t(scale(t(mat)))，非 scale() 的按列！）；
##           ② 池名未命中：Tumor cells 的 L2 实为 {Adipo diff, Invasion, ECM remodeling,
##           Stemness, Hypoxia, Angiogenesis}，合计 12881+23711=36592（逐类核对通过）
##  输出：_work/s34_6b_report.txt / fig6_data/pathway_*.csv
## =============================================================================
options(stringsAsFactors = FALSE)
suppressPackageStartupMessages({ library(Seurat); library(Matrix); library(GSVA); library(msigdbr) })

W <- file.path(Sys.getenv("LPS_WORK_ROOT", "."), "_work")
FIG <- file.path(W, "p5_figures/fig6_data"); dir.create(FIG, showWarnings = FALSE)
LOG <- file.path(W, "s34_6b_report.txt"); if (file.exists(LOG)) file.remove(LOG)
A <- function(...) { s <- paste0(...); cat(s, "\n", sep = ""); cat(s, "\n", sep = "", file = LOG, append = TRUE) }

A("===============================================================================")
A("S34 6B v3  pathway activity (GSVA pseudobulk)   ", format(Sys.time(), "%Y-%m-%d %H:%M:%S"))
A("===============================================================================")

obj <- readRDS(file.path(W, "s39_seurat.s6.rds"))
md  <- obj[[]]
A("loaded: ", paste(dim(obj), collapse = " x "))
A("")

## ---- L2 逐类核对（不得凭印象）----
tb <- table(md$L1, md$L2)
A("-- Tumor cells 的 L2 类 --")
tumL2 <- names(which(colSums(tb[, , drop = FALSE]) > 0))
tcL2 <- colnames(tb)[tb["Tumor cells", ] > 0]
for (k in tcL2) A(sprintf("  %-18s %7d", k, tb["Tumor cells", k]))
A(sprintf("  合计 = %d  (L1 Tumor cells n = %d)  %s", sum(tb["Tumor cells", tcL2]),
          sum(md$L1 == "Tumor cells"),
          if (sum(tb["Tumor cells", tcL2]) == sum(md$L1 == "Tumor cells")) "OK" else "MISMATCH"))
A("")

## ---- 自算 A-PPARG（D18 口径：逐基因跨细胞 z）----
G5 <- intersect(c("SREBF1", "SCD", "HMGCR", "ACLY", "FASN"), rownames(obj))
A("A-PPARG genes: ", paste(G5, collapse = ", "), " (n=", length(G5), ")")
M5 <- as.matrix(GetAssayData(obj, assay = "RNA", layer = "data")[G5, , drop = FALSE])
Z5 <- t(scale(t(M5)))                     # 逐基因跨细胞标准化（关键：非按列）
stopifnot(all(is.finite(colMeans(Z5))))
md$score_noPPARG <- as.numeric(Matrix::colMeans(Z5))
A(sprintf("score_noPPARG: mean=%.4f sd=%.4f | vs score_aucell rho=%.4f | vs ScoreA1 rho=%.4f",
          mean(md$score_noPPARG), sd(md$score_noPPARG),
          cor(md$score_noPPARG, md$score_aucell, method = "spearman"),
          cor(md$score_noPPARG, md$ScoreA1, method = "spearman")))
A("")

## ---- 通路集 ----
h <- msigdbr(species = "Homo sapiens", collection = "H")
k <- msigdbr(species = "Homo sapiens", collection = "C2", subcollection = "CP:KEGG_LEGACY")
key <- grep("MTOR|PPAR|FATTY_ACID|OXIDATIVE|ADIPOCYTO|LIPID|INSULIN|GLYCOLYSIS", unique(k$gs_name), value = TRUE)
kk <- k[k$gs_name %in% key, ]
gs <- split(c(h$gene_symbol, kk$gene_symbol), c(h$gs_name, kk$gs_name))
gs <- lapply(gs, unique)
A("pathway sets: Hallmark ", length(unique(h$gs_name)), " + KEGG(theme) ", length(key), " = ", length(gs))

## ---- 四个池（均为「Tumor cells 内部」或明确标注）----
isT <- md$L1 == "Tumor cells"
pools <- list(
  P1_tumor_all  = which(isT),
  P2_tumor_prog = which(isT & md$L2 %in% c("Invasion", "ECM remodeling", "Stemness", "Hypoxia", "Angiogenesis")),
  P3_adipo_diff = which(md$L2 == "Adipo diff"),
  P4_normal_adi = which(md$L1 == "Normal adipocytes")
)
present <- intersect(unique(unlist(gs)), rownames(obj))
A("union genes in sets: ", length(unique(unlist(gs))), " | present: ", length(present))
expr <- GetAssayData(obj, assay = "RNA", layer = "data")[present, , drop = FALSE]
smp <- md$sample
A("")
for (pn in names(pools)) A(sprintf("  %-14s cells=%6d", pn, length(pools[[pn]])))
A("")

out <- list(); meta <- list()
for (pn in names(pools)) {
  ii <- pools[[pn]]
  keys <- unique(smp[ii])
  M <- sapply(keys, function(kk2) Matrix::rowMeans(expr[, ii[smp[ii] == kk2], drop = FALSE]))
  nc <- sapply(keys, function(kk2) sum(smp[ii] == kk2))
  A(sprintf("---- GSVA %s : genes=%d pseudobulk=%d (cells/pb: min=%d med=%.0f max=%d) ----",
            pn, nrow(M), ncol(M), min(nc), median(nc), max(nc)))
  r <- gsva(gsvaParam(M, gs), verbose = FALSE)
  out[[pn]] <- r
  meta[[pn]] <- data.frame(pool = pn, sample = keys, n_cells = as.integer(nc),
                           score = sapply(keys, function(kk2) mean(md$score_aucell[ii][smp[ii] == kk2])),
                           score_noPPARG = sapply(keys, function(kk2) mean(md$score_noPPARG[ii][smp[ii] == kk2])),
                           row.names = NULL)
  A(sprintf("     dim=%s  cells>=5 的 pseudobulk 数=%d", paste(dim(r), collapse = " x "), sum(nc >= 5)))
}

## ---- 通路 vs 评分 ----
A("")
A("===============================================================================")
A("-- pathway activity vs tumor-cell score (Spearman, sample level) --")
A("===============================================================================")
tab <- data.frame()
for (pn in names(out)) {
  R <- out[[pn]]; m <- meta[[pn]]
  keep <- m$n_cells >= 5
  if (sum(keep) < 6) { A(sprintf("  [%s] 有效 pseudobulk 仅 %d 个 → 从略", pn, sum(keep))); next }
  mm <- m[keep, ]; idx <- match(mm$sample, colnames(R))
  for (gi in rownames(R)) {
    v <- as.numeric(R[gi, idx])
    if (sd(v) == 0) next
    c1 <- suppressWarnings(cor.test(mm$score, v, method = "spearman"))
    c2 <- suppressWarnings(cor.test(mm$score_noPPARG, v, method = "spearman"))
    tab <- rbind(tab, data.frame(pool = pn, pathway = gi, n = length(v),
                                 rho = unname(c1$estimate), P = c1$p.value,
                                 rho_noPPARG = unname(c2$estimate), P_noPPARG = c2$p.value))
  }
}
KEY <- grep("MTOR|PPAR|FATTY_ACID|OXIDATIVE|ADIPOCYTO|LIPID|GLYCOLYSIS|INSULIN", tab$pathway)
t2 <- tab[KEY, ]; t2 <- t2[order(-abs(t2$rho)), ]
A(sprintf("%-14s %-44s %3s %8s %10s %10s %10s", "pool", "pathway", "n", "rho", "P", "rho_noPP", "P_noPP"))
for (i in seq_len(min(45, nrow(t2))))
  A(sprintf("%-14s %-44s %3d %+8.4f %10.3g %+10.4f %10.3g",
            t2$pool[i], t2$pathway[i], t2$n[i], t2$rho[i], t2$P[i], t2$rho_noPPARG[i], t2$P_noPPARG[i]))
A("")
A(sprintf("显著项（P<0.05）合计 %d / %d", sum(tab$P < 0.05), nrow(tab)))
if (sum(tab$P < 0.05)) {
  sig <- tab[tab$P < 0.05, ]; sig <- sig[order(sig$P), ]
  for (i in seq_len(min(20, nrow(sig))))
    A(sprintf("  * %-14s %-44s rho=%+.3f P=%.3g", sig$pool[i], sig$pathway[i], sig$rho[i], sig$P[i]))
}

write.csv(do.call(rbind, lapply(names(out), function(pn)
  data.frame(pool = pn, pathway = rownames(out[[pn]]), out[[pn]], check.names = FALSE))),
  file.path(FIG, "pathway_pseudobulk.csv"), row.names = FALSE)
write.csv(do.call(rbind, meta), file.path(FIG, "pathway_pseudobulk_meta.csv"), row.names = FALSE)
write.csv(tab, file.path(FIG, "pathway_vs_score.csv"), row.names = FALSE)
A("")
A("wrote: pathway_pseudobulk.csv / pathway_pseudobulk_meta.csv / pathway_vs_score.csv")
A("DONE")
