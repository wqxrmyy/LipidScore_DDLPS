# ============================================================
# S4 · 判别度指标口径核查
#   1) 真 Harrell C-index（survival::concordance）
#   2) 复现手稿的朴素固定时点 AUC（Figure6 的做法，无 IPCW）
#   3) 正确的时间依赖 AUC（timeROC，IPCW）
#   队列：TCGA-SARC (OS, n=50) / GSE30929 all (DFS, n=140) / GSE30929 DDLPS (DFS, n=40)
# ============================================================
suppressWarnings(suppressMessages({
  library(survival)
  library(pROC)
}))
setwd(file.path(Sys.getenv("LPS_DATA_ROOT", "."), "1"))
set.seed(20260917)

cat("==================== 环境 ====================\n")
cat("R:", R.version.string, "\n")
for (p in c("survival", "pROC", "timeROC", "survivalROC", "riskRegression")) {
  cat(sprintf("  %-16s %s\n", p, requireNamespace(p, quietly = TRUE)))
}
options(width = 200)

# ---------- 1. TCGA ----------
scores <- read.csv("data/results/02_feature_scores.csv", row.names = 1)
clinical <- read.csv("data/results/03_clinical_outcomes.csv", row.names = 1)
tc <- clinical[clinical$Cohort == "TCGA", , drop = FALSE]
tc <- tc[!is.na(tc$OS_time) & !is.na(tc$OS_status), , drop = FALSE]
tc$score <- scores[rownames(tc), "Lipid_Metabolism_Score"]
tc$t_yr <- tc$OS_time / 365
cat("\n==================== TCGA-SARC (OS) ====================\n")
cat("n =", nrow(tc), " 事件 =", sum(tc$OS_status), " 随访(年) 中位 =",
    round(median(tc$t_yr), 2), " 最大 =", round(max(tc$t_yr), 2), "\n")

harrell <- function(t, e, s, label) {
  fit <- survConcordance(Surv(t, e) ~ s)
  cval <- fit$concordance
  # 标准误：用 concordance() 的 variance（仅 survival>=3.0 可用）
  se <- tryCatch({
    cf <- concordance(Surv(t, e) ~ s)
    sqrt(cf$var)
  }, error = function(e) NA_real_)
  cat(sprintf("  [真 Harrell C] %-22s C = %.4f  (SE %.4f, 95%%CI %.3f-%.3f)\n",
              label, cval, se, cval - 1.96 * se, cval + 1.96 * se))
  invisible(cval)
}

c_tcga <- harrell(tc$t_yr, tc$OS_status, tc$score, "TCGA-SARC OS")

# ---------- 2. 复现手稿的朴素时点 AUC ----------
naive_auc <- function(time, status, score, yr) {
  st <- ifelse(time <= yr & status == 1, 1, ifelse(time > yr, 0, NA))
  ok <- !is.na(st)
  if (sum(ok) > 10 && length(unique(st[ok])) == 2) {
    r <- roc(st[ok], score[ok], quiet = TRUE)
    return(list(auc = as.numeric(auc(r)), n = sum(ok), npos = sum(st[ok] == 1)))
  }
  list(auc = NA, n = sum(ok), npos = NA)
}
cat("\n  --- 复现手稿 Figure6A 的朴素时点 AUC（无 IPCW）---\n")
naive_res <- list()
for (yr in c(1, 3, 5)) {
  r <- naive_auc(tc$t_yr, tc$OS_status, tc$score, yr)
  naive_res[[as.character(yr)]] <- r
  cat(sprintf("    %d-year: AUC = %.3f   (纳入 %d 例, 事件 %s)\n", yr, r$auc, r$n,
              ifelse(is.na(r$npos), "NA", r$npos)))
}

# ---------- 3. 正确的时间依赖 AUC ----------
cat("\n  --- 时间依赖 AUC（IPCW, timeROC）---\n")
td_res <- NULL
if (requireNamespace("timeROC", quietly = TRUE)) {
  tr <- timeROC::timeROC(T = tc$t_yr, delta = tc$OS_status,
                         marker = tc$score, cause = 1,
                         times = c(1, 3, 5), iid = TRUE)
  cat(sprintf("    AUC(t=1) = %.4f  (SE %.4f)\n", tr$AUC[1], tr$inference$vect_sd_1[1]))
  cat(sprintf("    AUC(t=3) = %.4f  (SE %.4f)\n", tr$AUC[2], tr$inference$vect_sd_1[2]))
  cat(sprintf("    AUC(t=5) = %.4f  (SE %.4f)\n", tr$AUC[3], tr$inference$vect_sd_1[3]))
  td_res <- tr$AUC
} else {
  cat("    timeROC 未安装，跳过\n")
}

# ---------- 4. GSE30929 ----------
cat("\n==================== GSE30929 (DFS) ====================\n")
g <- readRDS(file.path(Sys.getenv("LPS_DATA_ROOT", "."), "2", "raw", "bulk", "GSE30929_clinical_all.rds"))
cat("总样本 =", nrow(g), "  亚型:", paste(names(table(g$subtype)), table(g$subtype), sep = "=", collapse = " "), "\n")
cat("随访(月) 中位 =", round(median(g$OS_time), 1), " 最大 =", round(max(g$OS_time), 1),
    " 事件 =", sum(g$OS_status), "\n")

c_g_all <- harrell(g$OS_time, g$OS_status, g$LipidScore, "GSE30929 全队列 DFS")
gd <- g[g$subtype == "dedifferentiated", , drop = FALSE]
cat("  DDLPS 子集 n =", nrow(gd), " 事件 =", sum(gd$OS_status), "\n")
c_g_dd <- harrell(gd$OS_time, gd$OS_status, gd$LipidScore, "GSE30929 DDLPS DFS")

cat("\n==================== 汇总对照 ====================\n")
cat(sprintf("%-30s %-12s %-12s %-12s\n", "队列", "手稿值", "真HarrellC", "朴素1yrAUC"))
cat(sprintf("%-30s %-12s %-12.4f %-12s\n", "TCGA-SARC (OS, n=50)", "0.720",
            c_tcga, sprintf("%.3f", naive_res[["1"]]$auc)))
cat(sprintf("%-30s %-12s %-12.4f %-12s\n", "GSE30929 all (DFS, n=140)", "0.629", c_g_all, "-"))
cat(sprintf("%-30s %-12s %-12.4f %-12s\n", "GSE30929 DDLPS (DFS, n=40)", "0.580", c_g_dd, "-"))
if (!is.null(td_res)) {
  cat(sprintf("\n手稿的 1/3/5 年 AUC (0.72 / 0.68 / 0.65)\n  朴重复现: %.3f / %.3f / %.3f\n  正确 IPCW: %.4f / %.4f / %.4f\n",
              naive_res[["1"]]$auc, naive_res[["3"]]$auc, naive_res[["5"]]$auc,
              td_res[1], td_res[2], td_res[3]))
}

# ---------- 5. 落盘 ----------
out <- data.frame(
  cohort = c("TCGA-SARC_OS", "GSE30929_DFS", "GSE30929_DDLPS_DFS"),
  n = c(nrow(tc), nrow(g), nrow(gd)),
  events = c(sum(tc$OS_status), sum(g$OS_status), sum(gd$OS_status)),
  manuscript_value = c(0.720, 0.629, 0.580),
  harrell_C = c(c_tcga, c_g_all, c_g_dd)
)
write.csv(out, file.path(Sys.getenv("LPS_WORK_ROOT", "."), "_work", "s4_cindex.csv"),
          row.names = FALSE, fileEncoding = "UTF-8")
td <- data.frame(horizon_yr = c(1, 3, 5),
                 naive_auc = c(naive_res[["1"]]$auc, naive_res[["3"]]$auc, naive_res[["5"]]$auc),
                 ipcw_auc = if (is.null(td_res)) NA else as.numeric(td_res))
write.csv(td, file.path(Sys.getenv("LPS_WORK_ROOT", "."), "_work", "s4_auc_tcga.csv"),
          row.names = FALSE, fileEncoding = "UTF-8")
cat("\n已写出 s4_cindex.csv / s4_auc_tcga.csv\nDONE\n")
