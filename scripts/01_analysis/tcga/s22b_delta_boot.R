# S22b · A vs C 的 Harrell C 差异：配对 bootstrap（同一批 50 例重抽，保持两评分相关性）
# 目的：回答「A 比 C 高 0.072，这个差是真的吗」；为 Fig4D 提供 ΔC 与其 95% CI
suppressWarnings(suppressMessages(library(survival)))
setwd(file.path(Sys.getenv("LPS_DATA_ROOT", "."), "1"))
set.seed(20260918)

a   <- read.csv("data/raw/TCGA_Score.csv", stringsAsFactors = FALSE)
rownames(a) <- a$Sample
cdf <- read.csv("data/results/02_feature_scores.csv", row.names = 1, check.names = FALSE)
clin <- read.csv("data/results/03_clinical_outcomes.csv", stringsAsFactors = FALSE)
rownames(clin) <- clin$Sample_ID

tc <- clin[clin$Cohort == "TCGA", , drop = FALSE]
tc <- tc[!is.na(tc$OS_time) & !is.na(tc$OS_status), , drop = FALSE]
tc$score_A <- a[rownames(tc), "Score"]
tc$score_C <- cdf[rownames(tc), "Lipid_Metabolism_Score"]
tc$t_yr <- tc$OS_time / 365
d <- tc[!is.na(tc$score_A) & !is.na(tc$score_C), , drop = FALSE]

cat("配对可用 n =", nrow(d), " 事件 =", sum(d$OS_status), "\n\n")

cz <- function(t, e, s) {
  if (length(unique(s)) < 2) return(NA_real_)
  survConcordance(Surv(t, e) ~ s)$concordance
}

cA <- cz(d$t_yr, d$OS_status, d$score_A)
cC <- cz(d$t_yr, d$OS_status, d$score_C)
cat(sprintf("原始：A = %.4f   C = %.4f   Δ(A-C) = %+.4f\n\n", cA, cC, cA - cC))

# ---- 配对 bootstrap：同一次重抽同时算 A 与 C ----
B <- 5000
n <- nrow(d)
dA <- dC <- numeric(B)
for (b in 1:B) {
  idx <- sample.int(n, n, replace = TRUE)
  db <- d[idx, , drop = FALSE]
  dA[b] <- cz(db$t_yr, db$OS_status, db$score_A)
  dC[b] <- cz(db$t_yr, db$OS_status, db$score_C)
}
ok <- is.finite(dA) & is.finite(dC)
dD <- (dA - dC)[ok]
cat("有效 bootstrap 次数 =", length(dD), "/", B, "\n")

ci <- quantile(dD, c(0.025, 0.975), na.rm = TRUE)
# 双侧 P：Δ 分布跨越 0 的比例（按偏离方向取双尾）
p_lo <- mean(dD <= 0); p_hi <- mean(dD >= 0)
pval <- 2 * min(p_lo, p_hi)
pval <- min(pval, 1)

cat(sprintf("ΔC = %+.4f   95%%CI %+.4f ~ %+.4f   双侧 P = %.3f\n", cA - cC, ci[1], ci[2], pval))

# 各自 C 的 bootstrap CI（与解析 SE 对拍，验证 bootstrap 口径正常）
ciA <- quantile(dA[ok], c(0.025, 0.975)); ciC <- quantile(dC[ok], c(0.025, 0.975))
cat(sprintf("A bootstrap CI %.3f-%.3f  (解析 0.562-0.817)\n", ciA[1], ciA[2]))
cat(sprintf("C bootstrap CI %.3f-%.3f  (解析 0.508-0.727)\n", ciC[1], ciC[2]))

# 两评分本身的相关性（解释为何要配对）
cat(sprintf("\nscore A 与 score C 的 Pearson r = %.4f\n", cor(d$score_A, d$score_C)))

out <- data.frame(
  item = c("C_A", "C_C", "delta", "delta_lo", "delta_hi", "p_boot",
           "C_A_lo", "C_A_hi", "C_C_lo", "C_C_hi", "B", "n", "events", "r_AC"),
  value = c(cA, cC, cA - cC, ci[1], ci[2], pval,
            ciA[1], ciA[2], ciC[1], ciC[2], length(dD), nrow(d), sum(d$OS_status),
            cor(d$score_A, d$score_C))
)
write.csv(out,
          file.path(Sys.getenv("LPS_WORK_ROOT", "."), "_work", "s22b_delta_boot.csv"),
          row.names = FALSE, fileEncoding = "UTF-8")
cat("\n已写出 s22b_delta_boot.csv\n")
