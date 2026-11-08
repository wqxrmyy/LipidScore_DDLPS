## =============================================================================
##  S34 · 6E：slingshot 拟时序 —— 肿瘤细胞内的「脂质身份」轨迹
##  ---------------------------------------------------------------------------
##  为什么不是 monocle3：monocle3 未安装（s34_6c_probe.txt）；slingshot 2.20 在位，
##  且 D26 §三把 slingshot 列为备选首选（Bioc、轻）。口径差异须显式声明：
##  slingshot = 簇水平 MST + 曲线拟合；monocle3 = 单细胞 DDRTree。二者都非「速率」模型。
##  ★ 轨迹假设的争议（必须并列，不得择一）：
##    V1（主）  : 肿瘤细胞全 6 个 L2 类，start = Adipo diff（脂肪样分化表型为高分化端）
##    V2（敏感性）: 排除 Adipo diff（5 类），start = Hypoxia
##    理由：Adipo diff 源注释 L1 = Tumor cells，但在非肿瘤组织同样高丰度
##          （nl 15.48% / lipoma 16.12% / WDLPS 16.65%，仅 DD 塌缩至 3.20%，见 LIM-P3-015）
##          → 它是否属肿瘤谱系本身存疑，故 V1/V2 并列。
##  降维：主用已有 harmony 校正后的 UMAP（与 Fig5 同空间，可对照）；
##        另算 1000-HVG 的 20 维 PCA 作敏感性（避免 2 维 UMAP 距离失真）。
##  输出：s34_6e_report.txt / fig6_data/slingshot_*.csv
## =============================================================================
options(stringsAsFactors = FALSE)
suppressPackageStartupMessages({
  library(Seurat); library(Matrix); library(slingshot)
  library(SingleCellExperiment); library(irlba)
})
future_ok <- requireNamespace("future", quietly = TRUE)
if (future_ok) future::plan("sequential")   ## ★ Windows 下 multisession 会崩

W <- file.path(Sys.getenv("LPS_WORK_ROOT", "."), "_work")
FIG <- file.path(W, "p5_figures/fig6_data"); dir.create(FIG, showWarnings = FALSE)
LOG <- file.path(W, "s34_6e_report.txt"); if (file.exists(LOG)) file.remove(LOG)
A <- function(...) { s <- paste0(...); cat(s, "\n", sep = ""); cat(s, "\n", sep = "", file = LOG, append = TRUE) }
P <- function(x, ...) { o <- capture.output(print(x, ...)); for (l in o) A(l) }

set.seed(20260924)
A("===============================================================================")
A("S34 6E slingshot pseudotime (lipid-identity trajectory, tumor compartment)   ", format(Sys.time(), "%Y-%m-%d %H:%M:%S"))
A("===============================================================================")

obj <- readRDS(file.path(W, "s39_seurat.s6.rds")); md <- obj[[]]
A("loaded: ", paste(dim(obj), collapse = " x "))
A("colData cols: ", paste(colnames(md), collapse = ", "))

TUM <- which(md$L1 == "Tumor cells")
A(sprintf("tumor cells = %d", length(TUM)))
A("-- tumor L2 counts --");  P(table(md$L2[TUM]))
A("-- tumor tissue counts --"); P(table(md$tissue[TUM]))
A("-- tumor L2 x tissue --"); P(table(md$L2[TUM], md$tissue[TUM]))

cells <- colnames(obj)[TUM]
coldf <- data.frame(
  clus   = factor(md$L2[TUM]),
  tissue = factor(md$tissue[TUM]),
  sample = factor(md$sample[TUM]),
  score  = as.numeric(md$score_aucell[TUM]),
  score1 = as.numeric(md$ScoreA1[TUM]),
  umap1  = as.numeric(md$umap1[TUM]),
  umap2  = as.numeric(md$umap2[TUM]),
  row.names = cells
)
A(sprintf("score_aucell in tumor: mean=%.4f sd=%.4f range=[%.3f, %.3f]",
          mean(coldf$score), sd(coldf$score), min(coldf$score), max(coldf$score)))

## ---- 敏感性用的 20 维 PCA（1000 HVG，irlba 直接算，避免稠密化全矩阵）-----
A("")
A("-- 0. PCA (1000 HVG, 20 PCs) for sensitivity --")
hv <- tryCatch({
  sub <- subset(obj, cells = cells)
  sub <- FindVariableFeatures(sub, nfeatures = 1000, verbose = FALSE)
  VariableFeatures(sub)
}, error = function(e) { A(paste("   HVG failed:", conditionMessage(e))); NULL })
pc <- NULL
if (!is.null(hv)) {
  M <- as.matrix(GetAssayData(obj, assay = "RNA", layer = "data")[hv, cells, drop = FALSE])
  A(sprintf("   HVG matrix: %d x %d (%.0f MB)", nrow(M), ncol(M), as.numeric(object.size(M)) / 1024^2))
  pc <- irlba::prcomp_irlba(t(M), n = 20, center = TRUE, scale. = FALSE)
  A(sprintf("   PCA done: %d x %d ; var explained (PC1-5) = %s", nrow(pc$x), ncol(pc$x),
            paste(sprintf("%.3f", (pc$sdev^2 / sum(pc$sdev^2))[1:5]), collapse = ", ")))
  rm(M); gc(verbose = FALSE)
}
if (exists("sub")) rm(sub)
gc(verbose = FALSE)

## ---- 通用：跑一个 slingshot 版本 -------------------------------------------
run_ss <- function(tag, keep_clus, start_clus, rd_mat, rd_name) {
  A("")
  A(paste0("=============== ", tag, " ==============="))
  A(sprintf("clusters kept = %s | start = %s | reducedDim = %s", paste(keep_clus, collapse = "/"),
            start_clus, rd_name))
  ii <- which(as.character(coldf$clus) %in% keep_clus)
  A(sprintf("cells = %d", length(ii)))
  P(table(coldf$clus[ii]))
  if (length(ii) < 200) { A("   !! too few cells, skip"); return(NULL) }
  cd <- coldf[ii, , drop = FALSE]
  rd <- rd_mat[ii, , drop = FALSE]
  sce <- SingleCellExperiment(assays = list(dummy = matrix(0, 1, nrow(cd))))
  reducedDim(sce, rd_name) <- rd
  colData(sce)$clus <- droplevels(cd$clus)
  sce <- slingshot(sce, clusterLabels = "clus", reducedDim = rd_name, start.clus = start_clus)
  pt <- slingPseudotime(sce)
  A(sprintf("lineages = %d ; pseudotime dims = %s", ncol(pt), paste(dim(pt), collapse = " x ")))
  A("-- lineage cluster order --")
  for (l in seq_len(ncol(pt))) {
    ok <- is.finite(pt[, l])
    if (sum(ok) < 10) next
    mm <- tapply(pt[ok, l], cd$clus[ok], median)
    A(sprintf("   L%d : %s", l, paste(sprintf("%s(%.2f)", names(mm), mm), collapse = " -> ")))
  }
  A("-- pseudotime range per lineage --")
  P(apply(pt, 2, function(z) round(range(z, na.rm = TRUE), 3)))
  ## ★ 沿轨迹的评分趋势（核心读数）
  A("-- ★ score_aucell vs pseudotime (Spearman) --")
  for (l in seq_len(ncol(pt))) {
    ok <- is.finite(pt[, l])
    if (sum(ok) < 30) next
    ct <- suppressWarnings(cor.test(pt[ok, l], cd$score[ok], method = "spearman"))
    ## 首/末五分位均值
    q <- quantile(pt[ok, l], c(0.2, 0.8))
    m1 <- mean(cd$score[ok][pt[ok, l] <= q[1]]); m2 <- mean(cd$score[ok][pt[ok, l] >= q[2]])
    A(sprintf("   L%d : rho=%+.4f  P=%.3g  n=%d  | bottom20%%=%.4f  top20%%=%.4f  Δ=%+.4f",
              l, unname(ct$estimate), ct$p.value, sum(ok), m1, m2, m2 - m1))
  }
  ## 组织分布（WD/sWD/DD）
  A("-- tissue distribution along lineage (median pseudotime) --")
  for (l in seq_len(ncol(pt))) {
    ok <- is.finite(pt[, l]); if (sum(ok) < 30) next
    mm <- tapply(pt[ok, l], droplevels(cd$tissue[ok]), median)
    nn <- table(cd$tissue[ok])
    A(sprintf("   L%d : %s", l, paste(sprintf("%s(med=%.2f,n=%d)", names(mm), mm, nn[names(mm)]), collapse = " | ")))
  }
  ## 沿轨迹的六基因表达趋势
  G <- intersect(c("SREBF1", "SCD", "HMGCR", "ACLY", "PPARG", "FASN", "CD36", "FABP4",
                   "ADIPOQ", "PLIN1", "LPL", "MKI67", "TOP2A"), rownames(obj))
  E <- as.matrix(GetAssayData(obj, assay = "RNA", layer = "data")[G, cells[ii], drop = FALSE])
  A("-- gene expression vs pseudotime (Spearman rho, lineage with strongest |rho|) --")
  for (g in G) {
    best <- NULL
    for (l in seq_len(ncol(pt))) {
      ok <- is.finite(pt[, l]); if (sum(ok) < 30) next
      r <- suppressWarnings(cor(pt[ok, l], E[g, ok], method = "spearman"))
      if (is.null(best) || abs(r) > abs(best$r)) best <- list(r = r, l = l, p = NA)
    }
    if (!is.null(best)) A(sprintf("   %-8s best L%d rho=%+.4f", g, best$l, best$r))
  }
  list(sce = sce, pt = pt, cd = cd, tag = tag)
}

V1 <- run_ss("V1 全 6 类 start=Adipo diff", 
             c("Adipo diff", "Invasion", "ECM remodeling", "Stemness", "Hypoxia", "Angiogenesis"),
             "Adipo diff", cbind(coldf$umap1, coldf$umap2), "UMAP")
V2 <- run_ss("V2 去 Adipo diff (5 类) start=Hypoxia",
             c("Invasion", "ECM remodeling", "Stemness", "Hypoxia", "Angiogenesis"),
             "Hypoxia", cbind(coldf$umap1, coldf$umap2), "UMAP")
if (!is.null(pc)) {
  run_ss("V3 敏感性：PCA20 空间，全 6 类 start=Adipo diff",
         c("Adipo diff", "Invasion", "ECM remodeling", "Stemness", "Hypoxia", "Angiogenesis"),
         "Adipo diff", pc$x, "PCA")
}

## ---- 导出 pseudotime 供 P5 出图 -------------------------------------------
if (!is.null(V1)) {
  pt1 <- V1$pt
  out <- data.frame(cell = rownames(V1$cd), clus = V1$cd$clus, tissue = V1$cd$tissue,
                    score = V1$cd$score, umap1 = V1$cd$umap1, umap2 = V1$cd$umap2)
  for (l in seq_len(ncol(pt1))) out[[paste0("pt_L", l)]] <- pt1[, l]
  write.csv(out, file.path(FIG, "slingshot_pseudotime_V1.csv"), row.names = FALSE)
  A("")
  A("wrote: slingshot_pseudotime_V1.csv")
}
saveRDS(list(V1 = V1, V2 = V2), file.path(W, "s34_6e_slingshot.rds"))
A("saved: s34_6e_slingshot.rds")
A("")
A("DONE")
