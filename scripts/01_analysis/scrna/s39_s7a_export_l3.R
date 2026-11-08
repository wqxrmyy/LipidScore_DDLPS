## =============================================================================
##  DDLPS 重构 · 单细胞 #39  Step 7a：导出「样本 × L3」pseudobulk
##  ---------------------------------------------------------------------------
##  目的：Step 7 的粒度敏感性需要 L1 / L2 / L3 三个粒度。
##        Step 6 落盘的 cell_scores.csv.gz 只含 L1 / L2 / subtype_code（**无 L3**），
##        故从 s6.rds 的 metadata 直接聚合 L3 —— 不给主分析留缺口。
##  输入：_work/s39_seurat.s6.rds
##  输出：_work/s6/pseudobulk_sample_x_L3.csv  /  s39_s7a_report.txt
## =============================================================================

options(stringsAsFactors = FALSE)
suppressPackageStartupMessages({ library(Seurat) })
options(mc.cores = 1)

W   <- file.path(Sys.getenv("LPS_WORK_ROOT", "."), "_work")
OUT <- file.path(W, "s6")
LOG <- file.path(W, "s39_s7a_report.txt")
if (file.exists(LOG)) file.remove(LOG)
A <- function(...) {
  s <- paste0(...); cat(s, "\n", sep = "")
  cat(s, "\n", sep = "", file = LOG, append = TRUE)
}

A("===============================================================================")
A("Step 7a  导出 L3 pseudobulk   ", format(Sys.time(), "%Y-%m-%d %H:%M:%S"))
A("===============================================================================")

obj <- readRDS(file.path(W, "s39_seurat.s6.rds"))
A("dim = ", paste(dim(obj), collapse = " x "))
meta <- obj[[]]
A("meta 列: ", paste(colnames(meta), collapse = ", "))

need <- c("patient", "tissue", "L1", "L2", "L3", "score_ams", "score_aucell")
miss <- setdiff(need, colnames(meta))
if (length(miss)) stop("metadata 缺列: ", paste(miss, collapse = ", "))
A("必需列齐全 ✓")

meta$sample <- paste(meta$patient, meta$tissue, sep = "_")
A("样本数 = ", length(unique(meta$sample)))

## 断言：配对 7 例齐全（沿用 Step 6c 的防护 —— 分组标签必须实测）
PAIR <- c("P1", "P4", "P7", "P10", "P14", "P15", "P16")
tb <- table(meta$patient, meta$tissue)
A("")
A("-- 配对完整性断言 --")
ok <- 0L
for (p in PAIR) {
  has_wd <- ("DDLPS-WD" %in% colnames(tb)) && tb[p, "DDLPS-WD"] > 0
  has_dd <- ("DDLPS-DD" %in% colnames(tb)) && tb[p, "DDLPS-DD"] > 0
  ok <- ok + (has_wd && has_dd)
  A(sprintf("  %-4s WD=%-5s DD=%-5s  %s", p,
            ifelse(has_wd, "ok", "MISS"), ifelse(has_dd, "ok", "MISS"),
            ifelse(has_wd && has_dd, "齐全", "★缺失")))
}
A("完整配对: ", ok, " / ", length(PAIR))
if (ok != length(PAIR)) stop("配对不完整，终止")

agg <- function(col) {
  kf <- factor(paste(meta$sample, meta[[col]], sep = "|"))
  n  <- as.vector(table(kf))
  s  <- rowsum(cbind(ams = meta$score_ams, auc = meta$score_aucell),
               group = kf, reorder = TRUE)
  nm <- strsplit(rownames(s), "|", fixed = TRUE)
  df <- data.frame(
    key = rownames(s),
    n_cells = n,
    mean_score_ams = s[, "ams"] / n,
    mean_score_aucell = s[, "auc"] / n,
    sample = sapply(nm, `[`, 1),
    cls = sapply(nm, `[`, 2),
    stringsAsFactors = FALSE
  )
  names(df)[ncol(df)] <- col
  df
}

d3 <- agg("L3")
write.csv(d3, file.path(OUT, "pseudobulk_sample_x_L3.csv"), row.names = FALSE)
A("")
A("样本xL3: ", nrow(d3), " 行 | L3 类数 = ", length(unique(d3$L3)),
  "  ->  s6/pseudobulk_sample_x_L3.csv")

## 对照：L1 / bucket 的行数应与既往一致（L1 需另聚合确认）
d1 <- agg("L1")
A("对照 样本xL1: ", nrow(d1), " 行 | L1 类数 = ", length(unique(d1$L1)))
A("对照 样本xL2（既有文件）应由 Python 端复核")
A("")
A("DONE_STEP7A")
