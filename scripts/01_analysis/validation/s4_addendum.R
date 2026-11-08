# S4 附：符号约定校验 + GSE30929 方向核查
suppressWarnings(suppressMessages({ library(survival) }))
setwd(file.path(Sys.getenv("LPS_DATA_ROOT", "."), "1"))
options(width = 200)

# 用 concordancefit 的 reverse 参数直接暴露两种约定
h2 <- function(t, e, x) {
  a <- survival::concordancefit(Surv(t, e), x, reverse = FALSE)$concordance
  b <- survival::concordancefit(Surv(t, e), x, reverse = TRUE)$concordance
  c(forward = a, reverse = b)
}

cat("========== A. 符号约定基准测试 ==========\n")
set.seed(1)
n <- 400
risk <- rnorm(n)
Tsurv <- rexp(n, rate = exp(0.8 * risk))          # risk 越大 -> 生存越短
E <- rbinom(n, 1, 0.8)
Tobs <- ifelse(E == 1, Tsurv, Tsurv * 1.5)
r <- h2(Tobs, E, risk)
cat(sprintf("  构造: risk 越大 -> 生存越短 (HR = exp(0.8) = %.2f)\n", exp(0.8)))
cat(sprintf("  reverse=FALSE  C = %.3f\n", r[1]))
cat(sprintf("  reverse=TRUE   C = %.3f\n", r[2]))
if (r[1] > r[2]) {
  cat("  => reverse=FALSE 把高分当【高风险】: C>0.5 表示 高分->生存短 方向一致\n")
  CONVENTION <- "risk"
} else {
  cat("  => reverse=FALSE 把高分当【长生存】\n")
  CONVENTION <- "survival"
}
cat("  本脚本采用约定:", CONVENTION, "\n")

cat("\n========== B. TCGA-SARC ==========\n")
scores <- read.csv("data/results/02_feature_scores.csv", row.names = 1)
clinical <- read.csv("data/results/03_clinical_outcomes.csv", row.names = 1)
tc <- clinical[clinical$Cohort == "TCGA", , drop = FALSE]
tc <- tc[!is.na(tc$OS_time) & !is.na(tc$OS_status), , drop = FALSE]
tc$score <- scores[rownames(tc), "Lipid_Metabolism_Score"]
tc$s_z <- as.numeric(scale(tc$score))
tc$t_yr <- tc$OS_time / 365
co <- coxph(Surv(t_yr, OS_status) ~ s_z, data = tc)
s <- summary(co)
cat(sprintf("  n=%d 事件=%d\n", nrow(tc), sum(tc$OS_status)))
cat(sprintf("  Cox 每 +1SD: HR = %.3f (%.3f-%.3f)  P = %.4f\n",
            exp(coef(co)), exp(confint(co))[1], exp(confint(co))[2], s$coef[5]))
rt <- h2(tc$t_yr, tc$OS_status, tc$s_z)
cat(sprintf("  Harrell C: reverse=FALSE %.4f | reverse=TRUE %.4f\n", rt[1], rt[2]))
c_tcga <- if (CONVENTION == "risk") rt[1] else rt[2]
cat(sprintf("  => 按风险约定, TCGA Harrell C = %.4f\n", c_tcga))

cat("\n========== C. GSE30929 方向核查（关键） ==========\n")
g <- readRDS(file.path(Sys.getenv("LPS_DATA_ROOT", "."), "2", "raw", "bulk", "GSE30929_clinical_all.rds"))
g$s_z <- as.numeric(scale(g$LipidScore))
co2 <- coxph(Surv(OS_time, OS_status) ~ s_z, data = g)
s2 <- summary(co2)
cat(sprintf("  全队列 n=%d 事件=%d 随访(月)中位 %.1f\n", nrow(g), sum(g$OS_status), median(g$OS_time)))
cat(sprintf("  Cox 每 +1SD: HR = %.3f (%.3f-%.3f)  P = %.4f\n",
            exp(coef(co2)), exp(confint(co2))[1], exp(confint(co2))[2], s2$coef[5]))
rg <- h2(g$OS_time, g$OS_status, g$s_z)
cat(sprintf("  Harrell C: reverse=FALSE %.4f | reverse=TRUE %.4f\n", rg[1], rg[2]))
c_g <- if (CONVENTION == "risk") rg[1] else rg[2]
cat(sprintf("  => 按风险约定, GSE30929 全队列 Harrell C = %.4f  %s\n", c_g,
            if (c_g > 0.5) "(方向与叙事一致)" else "(方向与叙事【相反】：高分预后更好)"))

g$grp <- factor(ifelse(g$LipidScore > median(g$LipidScore), "High", "Low"), levels = c("Low", "High"))
sd_ <- survdiff(Surv(OS_time, OS_status) ~ grp, data = g)
cat(sprintf("  中位分组 KM log-rank P = %.4f\n", 1 - pchisq(sd_$chisq, 1)))
cat("  各组中位/总体生存:\n")
print(survfit(Surv(OS_time, OS_status) ~ grp, data = g))
gd <- g[g$subtype == "dedifferentiated", , drop = FALSE]
gd$s_z <- as.numeric(scale(gd$LipidScore))
co3 <- coxph(Surv(OS_time, OS_status) ~ s_z, data = gd)
cat(sprintf("\n  DDLPS 子集 n=%d 事件=%d: Cox 每+1SD HR = %.3f  P = %.4f\n",
            nrow(gd), sum(gd$OS_status), exp(coef(co3)), summary(co3)$coef[5]))
rd <- h2(gd$OS_time, gd$OS_status, gd$s_z)
c_dd <- if (CONVENTION == "risk") rd[1] else rd[2]
cat(sprintf("  Harrell C = %.4f\n", c_dd))

cat("\n========== D. 手稿值溯源（0.629 / 0.580 / 0.72） ==========\n")
naive <- function(time, status, sc, yr) {
  st <- ifelse(time <= yr & status == 1, 1, ifelse(time > yr, 0, NA))
  ok <- !is.na(st)
  if (sum(ok) > 10 && length(unique(st[ok])) == 2) {
    r <- pROC::roc(st[ok], sc[ok], quiet = TRUE)
    return(sprintf("%.3f (n=%d, 阳性=%d)", as.numeric(pROC::auc(r)), sum(ok), sum(st[ok] == 1)))
  }
  "-"
}
cat("  [GSE30929 全队列] 朴素时点 AUC:\n")
for (yr in c(2, 3, 5, 10)) cat(sprintf("     %2d yr: %s\n", yr, naive(g$OS_time, g$OS_status, g$LipidScore, yr)))
cat("  [GSE30929 DDLPS] 朴素时点 AUC:\n")
for (yr in c(2, 3, 5, 10)) cat(sprintf("     %2d yr: %s\n", yr, naive(gd$OS_time, gd$OS_status, gd$LipidScore, yr)))
if (requireNamespace("timeROC", quietly = TRUE)) {
  tr <- timeROC::timeROC(T = g$OS_time, delta = g$OS_status, marker = g$LipidScore,
                         cause = 1, times = c(2, 3, 5, 10), iid = FALSE)
  cat(sprintf("  [GSE30929 全队列] IPCW AUC (2/3/5/10yr) = %s\n", paste(sprintf("%.3f", tr$AUC), collapse = " / ")))
  tr2 <- timeROC::timeROC(T = gd$OS_time, delta = gd$OS_status, marker = gd$LipidScore,
                          cause = 1, times = c(2, 3, 5, 10), iid = FALSE)
  cat(sprintf("  [GSE30929 DDLPS]  IPCW AUC (2/3/5/10yr) = %s\n", paste(sprintf("%.3f", tr2$AUC), collapse = " / ")))
}

cat("\n========== E. 汇总 ==========\n")
res <- data.frame(
  cohort = c("TCGA-SARC (OS)", "GSE30929 all (DFS)", "GSE30929 DDLPS (DFS)"),
  n = c(nrow(tc), nrow(g), nrow(gd)),
  events = c(sum(tc$OS_status), sum(g$OS_status), sum(gd$OS_status)),
  manuscript = c(0.720, 0.629, 0.580),
  harrell_C = c(c_tcga, c_g, c_dd),
  direction_ok = c(c_tcga > 0.5, c_g > 0.5, c_dd > 0.5)
)
print(res)
write.csv(res, file.path(Sys.getenv("LPS_WORK_ROOT", "."), "_work", "s4_cindex.csv"),
          row.names = FALSE, fileEncoding = "UTF-8")
cat("\nDONE\n")
