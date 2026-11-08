# S22 · D18 后：用 A 评分（6 基因 z 分均值）重算 TCGA-SARC 的 Harrell C
# 对照：同时输出 C 评分（02_feature_scores.csv）下的 Harrell C，确认 A 与 C 的差异
suppressWarnings(suppressMessages({
  library(survival)
}))
setwd(file.path(Sys.getenv("LPS_DATA_ROOT", "."), "1"))
set.seed(20260917)

sink(file.path(Sys.getenv("LPS_WORK_ROOT", "."), "_work", "s22_harrell_A_out.txt"))

# ---- 读 A 评分（TCGA_Score.csv）----
a <- read.csv("data/raw/TCGA_Score.csv", stringsAsFactors = FALSE)
rownames(a) <- a$Sample
cat("A 评分样本数:", nrow(a), "\n")

# ---- 读 C 评分（02_feature_scores.csv）作对照 ----
cdf <- read.csv("data/results/02_feature_scores.csv", row.names = 1, check.names = FALSE)
cat("C 评分样本数:", nrow(cdf), "\n")

# ---- 读临床 ----
clin <- read.csv("data/results/03_clinical_outcomes.csv", stringsAsFactors = FALSE)
rownames(clin) <- clin$Sample_ID
tc <- clin[clin$Cohort == "TCGA", , drop = FALSE]
tc <- tc[!is.na(tc$OS_time) & !is.na(tc$OS_status), , drop = FALSE]
cat("TCGA 临床样本数:", nrow(tc), " 事件:", sum(tc$OS_status), "\n")

# ---- 对齐 A 评分 ----
common_a <- intersect(rownames(tc), rownames(a))
cat("A 与临床对齐样本数:", length(common_a), "\n")
tc$score_A <- a[rownames(tc), "Score"]
# 对齐 C 评分
tc$score_C <- cdf[rownames(tc), "Lipid_Metabolism_Score"]
tc$t_yr <- tc$OS_time / 365

# 去掉评分 NA 的样本
tcA <- tc[!is.na(tc$score_A), , drop = FALSE]
tcC <- tc[!is.na(tc$score_C), , drop = FALSE]
cat("A 基准可用 n =", nrow(tcA), " 事件 =", sum(tcA$OS_status), "\n")
cat("C 基准可用 n =", nrow(tcC), " 事件 =", sum(tcC$OS_status), "\n")

harrell <- function(t, e, s, label) {
  fit <- survConcordance(Surv(t, e) ~ s)
  cval <- fit$concordance
  se <- tryCatch({
    cf <- concordance(Surv(t, e) ~ s)
    sqrt(cf$var)
  }, error = function(e) NA_real_)
  cat(sprintf("  %-28s C = %.4f  (SE %.4f, 95%%CI %.3f-%.3f)\n",
              label, cval, se, cval - 1.96 * se, cval + 1.96 * se))
  invisible(c(cval = cval, se = se))
}

cat("\n==================== Harrell C（D18 后，A 基准）====================\n")
rA <- harrell(tcA$t_yr, tcA$OS_status, tcA$score_A, "TCGA-SARC OS (A 评分)")
cat("\n==================== 对照：C 评分 ====================\n")
rC <- harrell(tcC$t_yr, tcC$OS_status, tcC$score_C, "TCGA-SARC OS (C 评分)")

# 风险方向确认（reverse 约定见统计铁律 1）
cat("\n==================== 方向确认 ====================\n")
# 高分 vs 低分：C-index > 0.5 表示高分风险更高（若 score 越大 OS 越短）
m <- median(tcA$score_A)
hi <- tcA$score_A > m
km_hi <- survfit(Surv(t_yr, OS_status) ~ 1, data = tcA[hi, ])
km_lo <- survfit(Surv(t_yr, OS_status) ~ 1, data = tcA[!hi, ])
cat("A 评分中位数分组的 KM（高分 vs 低分），中位 OS：",
    ifelse(is.null(summary(km_hi)$table["median"]), "NA", round(summary(km_hi)$table["median"], 2)),
    " vs ",
    ifelse(is.null(summary(km_lo)$table["median"]), "NA", round(summary(km_lo)$table["median"], 2)),
    " 年\n")

# 落盘
out <- data.frame(
  score_basis = c("A", "C"),
  n = c(nrow(tcA), nrow(tcC)),
  events = c(sum(tcA$OS_status), sum(tcC$OS_status)),
  harrell_C = c(rA["cval"], rC["cval"]),
  se = c(rA["se"], rC["se"])
)
write.csv(out, file.path(Sys.getenv("LPS_WORK_ROOT", "."), "_work", "s22_harrell_A.csv"),
          row.names = FALSE, fileEncoding = "UTF-8")
cat("\n已写出 s22_harrell_A.csv\n")
sink()
cat("DONE\n")
