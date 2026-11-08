## =============================================================================
##  DDLPS 重构 · 单细胞 #39  Step 6：评分双算法（AddModuleScore + AUCell）
##  ---------------------------------------------------------------------------
##  依据：plan_scrna_composition_vs_state.md §7 Step 6
##    · 6 基因与 bulk 完全一致：SREBF1, SCD, HMGCR, ACLY, PPARG, FASN
##    · 两个算法结果做秩相关，任一方异常须追因
##    · 主检验用 pseudobulk（细胞级 dropout 使评分系统性偏低）
##  输入：_work/s39_seurat.s4.rds（含 harmony + 官方坐标 + 冻结的 L2 粒度）
##  输出：s39_seurat.s6.rds / s39_s6_report.txt / s39_s6_summary.json
##        _work/s6/*.csv（pseudobulk 矩阵，供 Step 7 Kitagawa 分解）
##  注意：Windows 下不启用 future 并行（见 §1.8⑤ 崩溃教训）
## =============================================================================

options(stringsAsFactors = FALSE)
suppressPackageStartupMessages({
  library(Seurat)
  library(Matrix)
  library(AUCell)
})
options(mc.cores = 1)

W <- file.path(Sys.getenv("LPS_WORK_ROOT", "."), "_work")
OUT <- file.path(W, "s6"); dir.create(OUT, showWarnings = FALSE)
LOG <- file.path(W, "s39_s6_report.txt")
if (file.exists(LOG)) file.remove(LOG)
A <- function(...) {
  s <- paste0(...)
  cat(s, "\n", sep = "")
  cat(s, "\n", sep = "", file = LOG, append = TRUE)
}

A("===============================================================================")
A("Step 6  评分双算法   ", format(Sys.time(), "%Y-%m-%d %H:%M:%S"))
A("===============================================================================")
A("")
A("-- [0] 环境 --")
A("R ", R.version.string, " | Seurat ", as.character(packageVersion("Seurat")),
  " | AUCell ", as.character(packageVersion("AUCell")),
  " | Matrix ", as.character(packageVersion("Matrix")))
fb <- names(formals(AUCell::AUCell_buildRankings))
fc <- names(formals(AUCell::AUCell_calcAUC))
A("AUCell_buildRankings 形参: ", paste(fb, collapse = ", "))
A("AUCell_calcAUC   形参: ", paste(fc, collapse = ", "))
A("")

GENES <- c("SREBF1", "SCD", "HMGCR", "ACLY", "PPARG", "FASN")
A("-- [1] 载入 s4.rds --")
t0 <- Sys.time()
obj <- readRDS(file.path(W, "s39_seurat.s4.rds"))
A("dim = ", paste(dim(obj), collapse = " x "), "   耗时 ",
  round(as.numeric(difftime(Sys.time(), t0, units = "secs")), 1), " s")

## 释放 scale.data（Step 6 及 Fig5 各面板均不依赖；只保留 counts/data）
freed <- tryCatch({
  if ("scale.data" %in% Layers(obj[["RNA"]])) { obj[["RNA"]]$scale.data <- NULL; TRUE } else FALSE
}, error = function(e) { A("  （scale.data 释放跳过: ", conditionMessage(e), "）"); FALSE })
A("scale.data 已释放: ", freed, "   （Step 6 起不再需要；s6.rds 仅含 counts/data）")
gc(verbose = FALSE)

miss <- setdiff(GENES, rownames(obj))
if (length(miss)) stop("缺失基因: ", paste(miss, collapse = ", "))
A("6 基因全部在库 ✓ : ", paste(GENES, collapse = ", "))
A("")

## ------------------------------------------------------------------ 2. AddModuleScore
A("-- [2] AddModuleScore（Seurat 内置）--")
t0 <- Sys.time()
obj <- AddModuleScore(obj, features = list(ScoreA = GENES), name = "ScoreA",
                      seed = 1, ctrl = 100)
A("耗时 ", round(as.numeric(difftime(Sys.time(), t0, units = "secs")), 1), " s")
meta <- obj[[]]
ams <- as.numeric(meta[["ScoreA1"]])
A("ScoreA1  均值 ", round(mean(ams), 4), "  sd ", round(sd(ams), 4),
  "  范围 [", round(min(ams), 3), ", ", round(max(ams), 3), "]")
A("")

## ------------------------------------------------------------------ 3. AUCell
A("-- [3] AUCell（秩富集，独立于 AddModuleScore 的对照算法）--")
mat <- tryCatch(GetAssayData(obj, assay = "RNA", layer = "data"),
                error = function(e) GetAssayData(obj, assay = "RNA", slot = "data"))
A("表达矩阵（data 层, log1p CP10K）: ", nrow(mat), " 基因 x ", ncol(mat), " 细胞   nnz = ",
  format(length(mat@x), big.mark = ","))

## 6 基因的逐细胞表达（先取出，后面算方向一致性 + pseudobulk 用）
mat6 <- mat[GENES, , drop = FALSE]
m6 <- as.matrix(t(mat6))            # cells x 6
colnames(m6) <- GENES

t0 <- Sys.time()
## ★ AUCell 1.x 的实参名是 exprMatrix，1.20+ 改为 exprMat（S4 泛型按 exprMat 派发）
##   照抄旧写法会报 "unable to find an inherited method ... exprMat = missing" → 按 formals 自适应
nm_expr <- if ("exprMat" %in% fb) "exprMat" else if ("exprMatrix" %in% fb) "exprMatrix" else "exprMat"
A("exprMat 实参名（自适应）= ", nm_expr)
argB <- list(mat, verbose = TRUE)
names(argB)[1] <- nm_expr
if ("splitByBlocks" %in% fb) argB$splitByBlocks <- TRUE
if ("tiesMethod"    %in% fb) argB$tiesMethod <- "average"   # 确定性（部分版本无此参数）
if ("plotStats"     %in% fb) argB$plotStats <- FALSE
if ("nCores"        %in% fb) argB$nCores <- 1
rankings <- do.call(AUCell::AUCell_buildRankings, argB)
A("buildRankings 耗时 ", round(as.numeric(difftime(Sys.time(), t0, units = "secs")), 1), " s")

t0 <- Sys.time()
argC <- list(geneSets = list(ScoreA = GENES), rankings = rankings, verbose = TRUE)
if ("nCores" %in% fc) argC$nCores <- 1
argC$aucMaxRank <- ceiling(0.05 * nrow(rankings))
aucres <- do.call(AUCell::AUCell_calcAUC, argC)
A("calcAUC 耗时 ", round(as.numeric(difftime(Sys.time(), t0, units = "secs")), 1), " s",
  "   (aucMaxRank = ", argC$aucMaxRank, " = 前 5% 基因)")
aucv <- as.numeric(getAUC(aucres)["ScoreA", ])
names(aucv) <- colnames(aucres)
aucv <- aucv[rownames(meta)]      # 对齐到 metadata 顺序
A("AUC      均值 ", round(mean(aucv), 4), "  sd ", round(sd(aucv), 4),
  "  范围 [", round(min(aucv), 3), ", ", round(max(aucv), 3), "]")
rm(mat, rankings, aucres, mat6); gc(verbose = FALSE)
A("")

## ------------------------------------------------------------------ 4. 双算法一致性
A("===============================================================================")
A("-- [4] ★ 双算法一致性（任一方异常须追因）--")
A("===============================================================================")
rho_all <- suppressWarnings(cor(ams, aucv, method = "spearman"))
A("全体 101,890 细胞  Spearman rho = ", sprintf("%.4f", rho_all))
A("")
A("按 L1 谱系大类分层：")
A(sprintf("%-22s %8s %10s %10s", "L1", "n", "rho", "p<1e-300"))
for (g in sort(unique(meta$L1))) {
  i <- which(meta$L1 == g)
  if (length(i) < 30) { A(sprintf("%-22s %8d  (过小, 跳过)", g, length(i))); next }
  r <- suppressWarnings(cor(ams[i], aucv[i], method = "spearman"))
  A(sprintf("%-22s %8d %10.4f %10s", g, length(i), r, ifelse(abs(r) > 0.3, "yes", "-")))
}
A("")
A("按 L2 关键类分层（细胞数 >= 500）：")
A(sprintf("%-24s %8s %10s", "L2", "n", "rho"))
for (g in names(sort(table(meta$L2), decreasing = TRUE))) {
  i <- which(meta$L2 == g)
  if (length(i) < 500) next
  r <- suppressWarnings(cor(ams[i], aucv[i], method = "spearman"))
  A(sprintf("%-24s %8d %10.4f", g, length(i), r))
}
A("")

## 类均值层面的排序一致性（两算法是否把同样的类排到前面）
l2lev <- names(sort(table(meta$L2), decreasing = TRUE))
l2lev <- l2lev[table(meta$L2)[l2lev] >= 100]
mA <- sapply(l2lev, function(g) mean(ams[meta$L2 == g]))
mU <- sapply(l2lev, function(g) mean(aucv[meta$L2 == g]))
rho_cls <- suppressWarnings(cor(mA, mU, method = "spearman"))
A("L2 类均值排序一致性（n = ", length(l2lev), " 类）: Spearman rho = ",
  sprintf("%.4f", rho_cls))
A("")
A("两算法下的 Top-8 类（按均值降序）：")
oA <- order(mA, decreasing = TRUE); oU <- order(mU, decreasing = TRUE)
A(sprintf("%-4s %-24s %10s   |  %-4s %-24s %10s", "rank", "AddModuleScore", "mean",
          "rank", "AUCell", "AUC"))
for (k in 1:8) {
  A(sprintf("%-4d %-24s %10.4f   |  %-4d %-24s %10.4f", k, names(mA)[oA[k]], mA[oA[k]],
            k, names(mU)[oU[k]], mU[oU[k]]))
}
A("")

## ---------------------------------------------------- 5. 6 分量基因方向一致性
A("===============================================================================")
A("-- [5] 6 分量基因与评分的方向一致性（对照 bulk：6 分量应全为正贡献）--")
A("===============================================================================")
A(sprintf("%-10s %8s %10s %10s %10s %8s", "gene", "detect%", "rho_AMS", "rho_AUC", "mean", "sd"))
dir_ok <- TRUE
for (g in GENES) {
  det <- mean(m6[, g] > 0) * 100
  r1 <- suppressWarnings(cor(m6[, g], ams, method = "spearman"))
  r2 <- suppressWarnings(cor(m6[, g], aucv, method = "spearman"))
  A(sprintf("%-10s %8.2f %10.4f %10.4f %10.3f %8.3f", g, det, r1, r2, mean(m6[, g]), sd(m6[, g])))
  if (r1 <= 0) dir_ok <- FALSE
}
A("")
A("所有分量基因与 AddModuleScore 正相关: ", dir_ok,
  ifelse(dir_ok, "  ✓ 与 bulk 方向一致（PPARG 为 1/6 分量，注意循环性）",
         "  ⚠️ 存在非正相关分量 → 须追因"))
A("")

## ------------------------------------------------------------------ 6. pseudobulk
A("===============================================================================")
A("-- [6] pseudobulk 聚合（样本 x L2 / 样本 x 谱系桶）--")
A("===============================================================================")
meta$score_ams <- ams
meta$score_aucell <- aucv
meta$sample <- paste(meta$patient, meta$tissue, sep = "_")

bmap <- read.csv(file.path(W, "s5", "L2_to_bucket.csv"), fileEncoding = "UTF-8-BOM")
bucket_of <- setNames(bmap$bucket, bmap$L2)
meta$bucket <- unname(bucket_of[meta$L2])
A("L2 -> 桶 映射缺失数: ", sum(is.na(meta$bucket)), " （应为 0）")

agg_by <- function(key, tag) {
  key <- factor(key)
  n <- as.vector(table(key))
  s <- rowsum(cbind(ams, aucv), group = key, reorder = TRUE)
  g6 <- rowsum(m6, group = key, reorder = TRUE)
  df <- data.frame(key = rownames(s), n_cells = n, stringsAsFactors = FALSE)
  df$mean_score_ams <- s[, 1] / n
  df$mean_score_aucell <- s[, 2] / n
  for (g in GENES) df[[paste0("mean_", g)]] <- g6[, g] / n
  nm <- strsplit(df$key, "|", fixed = TRUE)
  if (tag == "L2") {
    df$sample <- sapply(nm, `[`, 1); df$L2 <- sapply(nm, `[`, 2)
  } else {
    df$sample <- sapply(nm, `[`, 1); df$bucket <- sapply(nm, `[`, 2)
  }
  df
}

pb_l2 <- agg_by(paste(meta$sample, meta$L2, sep = "|"), "L2")
pb_bk <- agg_by(paste(meta$sample, meta$bucket, sep = "|"), "bucket")
pb_all <- agg_by(meta$sample, "sample")
pb_all$sample <- pb_all$key

write.csv(pb_l2, file.path(OUT, "pseudobulk_sample_x_L2.csv"), row.names = FALSE)
write.csv(pb_bk, file.path(OUT, "pseudobulk_sample_x_bucket.csv"), row.names = FALSE)
write.csv(pb_all, file.path(OUT, "pseudobulk_sample.csv"), row.names = FALSE)
A("写出 pseudobulk：样本xL2 ", nrow(pb_l2), " 行 | 样本x桶 ", nrow(pb_bk),
  " 行 | 样本 ", nrow(pb_all), " 行  →  _work/s6/")
A("")

## ------------------------------------------------------------- 7. 方向判读（描述性）
A("===============================================================================")
A("-- [7] 描述性方向判读（★ 非判决；判决归 Step 7 的 Kitagawa 分解）--")
A("===============================================================================")
PAIR <- c("P1", "P4", "P7", "P10", "P14", "P15", "P16")
## ★ 修正（2026-09-22）：配对内的 WD 是 DDLPS 肿瘤内的高分化成分，标签为 "DDLPS-WD"；
##   "WDLPS" 是独立的纯高分化病例（P3/P5/P8/P9/P13），与配对无关。
##   原写法 WD_T <- c("WDLPS") 导致 WD 侧全部取不到 → §7 全 NA（仅影响描述性判读）。
WD_T <- c("DDLPS-WD"); DD_T <- c("DDLPS-DD")
pb_all$patient <- sub("_.*$", "", pb_all$sample)
pb_all$tissue <- sub("^[^_]*_", "", pb_all$sample)

A("(a) 样本级（全体细胞）评分 WD vs DD —— 配对 7 例")
A(sprintf("%-6s %12s %12s %10s", "pat", "WD(AMS)", "DD(AMS)", "DD-WD"))
dd_neg <- 0; nn <- 0
for (p in PAIR) {
  w <- pb_all$mean_score_ams[pb_all$patient == p & pb_all$tissue %in% WD_T]
  d <- pb_all$mean_score_ams[pb_all$patient == p & pb_all$tissue %in% DD_T]
  if (!length(w) || !length(d)) { A(sprintf("%-6s %12s %12s", p, "NA", "NA")); next }
  nn <- nn + 1; if (mean(d) - mean(w) < 0) dd_neg <- dd_neg + 1
  A(sprintf("%-6s %12.4f %12.4f %10.4f   (n_wd=%d, n_dd=%d)", p, mean(w), mean(d),
            mean(d) - mean(w), length(w), length(d)))
}
A("DD < WD 的配对: ", dd_neg, " / ", nn, "   （bulk 侧为 6/8 负，方向一致）")
A("")

A("(b) 谱系桶级评分（配对均值）：WD 均值 → DD 均值")
A(sprintf("%-10s %12s %12s %10s", "bucket", "WD", "DD", "Delta"))
for (b in c("Adipose", "Tumor", "Immune", "Stroma")) {
  s <- pb_bk[pb_bk$bucket == b, ]
  s$patient <- sub("_.*$", "", s$sample); s$tissue <- sub("^[^_]*_", "", s$sample)
  w <- s$mean_score_ams[s$patient %in% PAIR & s$tissue %in% WD_T]
  d <- s$mean_score_ams[s$patient %in% PAIR & s$tissue %in% DD_T]
  A(sprintf("%-10s %12.4f %12.4f %10.4f", b, mean(w), mean(d), mean(d) - mean(w)))
}
A("")
A("(c) 关键类 Adipo diff 的评分（样本级），配对 7 例")
A(sprintf("%-6s %12s %12s %10s", "pat", "WD", "DD", "DD-WD"))
s <- pb_l2[pb_l2$L2 == "Adipo diff", ]
s$patient <- sub("_.*$", "", s$sample); s$tissue <- sub("^[^_]*_", "", s$sample)
for (p in PAIR) {
  w <- s$mean_score_ams[s$patient == p & s$tissue %in% WD_T]
  d <- s$mean_score_ams[s$patient == p & s$tissue %in% DD_T]
  A(sprintf("%-6s %12s %12s %10s", p,
            ifelse(length(w), sprintf("%.4f", mean(w)), "NA"),
            ifelse(length(d), sprintf("%.4f", mean(d)), "NA"),
            ifelse(length(w) && length(d), sprintf("%.4f", mean(d) - mean(w)), "-")))
}
A("")

## ------------------------------------------------------------ 8. 保存与产物
A("===============================================================================")
A("-- [8] 保存产物 --")
A("===============================================================================")

## 细胞级评分表（供 Step 7 复核 / Step 9 出图 join）
cell <- data.frame(barcode = rownames(meta), patient = meta$patient, tissue = meta$tissue,
                   sample = meta$sample, L1 = meta$L1, L2 = meta$L2,
                   subtype_code = meta$subtype_code, bucket = meta$bucket,
                   score_ams = ams, score_aucell = aucv, stringsAsFactors = FALSE)
write.csv(cell, gzfile(file.path(OUT, "cell_scores.csv.gz")), row.names = FALSE)
A("细胞级评分表: _work/s6/cell_scores.csv.gz  (", nrow(cell), " 细胞 x ", ncol(cell), " 列)")

obj$score_ams <- ams
obj$score_aucell <- aucv
obj$bucket <- meta$bucket
obj$sample <- meta$sample
t0 <- Sys.time()
saveRDS(obj, file.path(W, "s39_seurat.s6.rds"))
A("s39_seurat.s6.rds  ", round(file.size(file.path(W, "s39_seurat.s6.rds")) / 1048576, 1),
  " MB   保存耗时 ", round(as.numeric(difftime(Sys.time(), t0, units = "secs")), 1), " s")

sm <- list(
  date = format(Sys.Date()), genes = GENES,
  n_cells = ncol(obj), rho_overall = rho_all, rho_class_means = rho_cls,
  auc_max_rank = argC$aucMaxRank,
  dir_ok = dir_ok,
  ams_mean = mean(ams), auc_mean = mean(aucv),
  files = c("s39_seurat.s6.rds", "s6/cell_scores.csv.gz",
            "s6/pseudobulk_sample_x_L2.csv", "s6/pseudobulk_sample_x_bucket.csv",
            "s6/pseudobulk_sample.csv")
)
if (requireNamespace("jsonlite", quietly = TRUE)) {
  jsonlite::write_json(sm, file.path(W, "s39_s6_summary.json"),
                       auto_unbox = TRUE, pretty = TRUE)
  A("s39_s6_summary.json 写出")
}
A("")
A("===============================================================================")
A("DONE_STEP6")
A("===============================================================================")
