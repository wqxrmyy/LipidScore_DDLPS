# s150_gse30929_verify.R
# 目的：对 GSE30929 用同一份本机文件，实测四件事，判定现稿 Table S6 的数值
#      与 S4 文档的数值哪一套可复现。
#   (1) 全队列 median-split HR (High vs Low) + Harrell C
#   (2) DDLPS 亚组 (subtype == dedifferentiated) HR + C
#   (3) 评分区分 WD / DD 的 AUC
#   (4) High 组的亚型构成
# 产物：_work/s148_out/s150_gse30929_verify.csv（幂等覆盖）
# 口径：C-index 一律“风险约定”(reverse=TRUE, 高分=风险高)

suppressPackageStartupMessages({
  library(survival)
  library(pROC)
})

IN  <- file.path(Sys.getenv("LPS_DATA_ROOT", "."), "2", "raw", "bulk", "GSE30929_Score_all.csv")
OUT <- file.path(Sys.getenv("LPS_WORK_ROOT", "."), "_work", "s148_out", "s150_gse30929_verify.csv")

d <- read.csv(IN, stringsAsFactors = FALSE)
cat("=== 输入 ===\n")
cat("file :", IN, "\n")
cat("rows :", nrow(d), "\n")
cat("cols :", paste(names(d), collapse = ", "), "\n")
cat("subtype 取值:\n"); print(table(d$subtype, useNA = "ifany"))
cat("ScoreGroup 取值:\n"); print(table(d$ScoreGroup, useNA = "ifany"))
cat("train_test 取值:\n"); print(table(d$train_test, useNA = "ifany"))
cat("\n")

# 统一列名
d$time   <- as.numeric(d$OS_time)
d$status <- as.numeric(d$OS_status)
d$score  <- as.numeric(d$LipidScore)

fitset <- function(x, label) {
  x <- x[is.finite(x$time) & is.finite(x$status) & is.finite(x$score), ]
  ev <- sum(x$status == 1)
  grp <- ifelse(x$score >= median(x$score), "High", "Low")
  grp <- factor(grp, levels = c("Low", "High"))
  cx  <- coxph(Surv(time, status) ~ grp, data = x)
  sm  <- summary(cx)
  hr  <- unname(sm$conf.int["grpHigh", "exp(coef)"])
  lo  <- unname(sm$conf.int["grpHigh", "lower .95"])
  hi  <- unname(sm$conf.int["grpHigh", "upper .95"])
  p   <- unname(sm$coefficients["grpHigh", "Pr(>|z|)"])
  cc  <- concordance(Surv(time, status) ~ score, data = x, reverse = TRUE)
  cat(sprintf("%-16s n=%3d ev=%2d | HR(High vs Low)=%.4f (%.3f-%.3f) P=%.4f | C=%.4f (%.4f-%.4f)\n",
              label, nrow(x), ev, hr, lo, hi, p,
              cc$concordance, cc$concordance - 1.96 * sqrt(cc$var),
              cc$concordance + 1.96 * sqrt(cc$var)))
  data.frame(layer = label, n = nrow(x), events = ev,
             HR_high_vs_low = round(hr, 4), lo = round(lo, 3), hi = round(hi, 3),
             P = signif(p, 4), C = round(cc$concordance, 4),
             C_lo = round(cc$concordance - 1.96 * sqrt(cc$var), 3),
             C_hi = round(cc$concordance + 1.96 * sqrt(cc$var), 3),
             stringsAsFactors = FALSE)
}

cat("=== (1)(2) 各层 HR / C ===\n")
res <- list()
res[["all"]] <- fitset(d, "GSE30929 all")
if ("dedifferentiated" %in% d$subtype)
  res[["dd"]] <- fitset(d[d$subtype == "dedifferentiated", ], "GSE30929 DDLPS")
if (any(grepl("^well", d$subtype, ignore.case = TRUE)))
  res[["wd"]] <- fitset(d[grepl("^well", d$subtype, ignore.case = TRUE), ], "GSE30929 WD")
R <- do.call(rbind, res)

cat("\n=== (3) 评分区分 WD vs DD 的 AUC ===\n")
sub <- d[d$subtype %in% c("dedifferentiated", "well-differentiated"), ]
sub$isDD <- as.integer(sub$subtype == "dedifferentiated")
if (length(unique(sub$isDD)) == 2) {
  r_dd <- pROC::roc(sub$isDD, sub$score, quiet = TRUE)
  cat(sprintf("score -> DD vs WD      AUC = %.4f (n=%d)\n", as.numeric(pROC::auc(r_dd)), nrow(sub)))
  r_wd <- pROC::roc(1 - sub$isDD, sub$score, quiet = TRUE)
  cat(sprintf("score -> WD vs rest    AUC = %.4f\n", as.numeric(pROC::auc(r_wd))))
  cat(sprintf("(1-AUC_dd) 等价 WD 侧 = %.4f\n", 1 - as.numeric(pROC::auc(r_dd))))
}

cat("\n=== (4) High/Low 组的亚型构成 ===\n")
d$grp <- ifelse(d$score >= median(d$score), "High", "Low")
print(table(d$grp, d$subtype))

cat("\n=== (5) High/Low 组事件与中位随访 ===\n")
for (g in c("High", "Low")) {
  z <- d[d$grp == g, ]
  cat(sprintf("%s: n=%d events=%d\n", g, nrow(z), sum(z$status == 1)))
}

dir.create(dirname(OUT), recursive = TRUE, showWarnings = FALSE)
if (file.exists(OUT)) file.remove(OUT)
write.csv(R, OUT, row.names = FALSE)
cat("\n落盘:", OUT, "\n")
