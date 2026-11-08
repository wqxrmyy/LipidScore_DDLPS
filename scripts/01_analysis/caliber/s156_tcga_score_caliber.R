# s156 · TCGA-SARC C 值口径定案 —— 0.690 vs 0.6175 到底是"分歧"还是"两套评分"？
# 背景：S4 文档记 TCGA OS Harrell C = 0.6175（CI 0.508-0.727）；
#       现稿/ s149 复算记 0.690（CI 0.562-0.817）。
#       线索：SI §"Rationale for score selection" 里 C 评分 = 0.618, CI 0.508-0.727 —— 与 0.6175 逐位吻合。
# 目的：把三套评分（A / B / C）在同一个 OS 端点、同一份临床表上各算 C 与 HR，一次性定案。
# 口径：与 S4 / s149 严格一致 —— concordance(Surv~score, reverse=TRUE)，风险约定（高分=风险高）。
# 幂等：全部覆盖写。
# 日期：2026-09-25

suppressMessages(library(survival))

base   <- file.path(Sys.getenv("LPS_DATA_ROOT", "."), "1", "data")
outdir <- file.path(Sys.getenv("LPS_WORK_ROOT", "."), "_work", "s148_out")
dir.create(outdir, showWarnings = FALSE, recursive = TRUE)

clin <- read.csv(file.path(base, "raw/TCGA_clinical_real.csv"),
                 fileEncoding = "UTF-8-BOM", check.names = FALSE)
clin$key <- substr(clin$barcode, 1, 12)

read_score <- function(path, col, label, aggfun = mean) {
  d <- read.csv(path, fileEncoding = "UTF-8-BOM", check.names = FALSE)
  names(d)[1] <- "Sample"
  names(d)[which(names(d) == col)] <- "Score"
  d$key <- substr(d$Sample, 1, 12)
  a <- aggregate(Score ~ key, data = d, FUN = aggfun)
  a$score_label <- label
  a
}

S <- list(
  A = read_score(file.path(base, "raw/TCGA_Score.csv"), "Score", "A (6-gene mean z; main text)"),
  B = read_score(file.path(base, "results/02_feature_scores_37genes.csv"),
                 "Lipid_Metabolism_Score", "B (37-gene principal component)"),
  C = read_score(file.path(base, "results/02_feature_scores.csv"),
                 "Lipid_Metabolism_Score", "C (alternative Lipid_Metabolism_Score)")
)

cat("== 评分文件与匹配情况 ==\n")
for (k in names(S)) {
  cat(sprintf("  %s : 聚合后 key 数 = %d\n", k, nrow(S[[k]])))
}

# ---- 在 OS 上逐套计算 ------------------------------------------------------
res <- list()
for (k in names(S)) {
  m   <- merge(clin, S[[k]][, c("key", "Score")], by = "key", all = FALSE)
  ok  <- !is.na(m$OS_time) & !is.na(m$OS_status) & m$OS_time > 0
  mm  <- m[ok, ]
  n   <- nrow(mm); ev <- sum(mm$OS_status == 1)
  Sv  <- Surv(mm$OS_time, mm$OS_status)
  s   <- mm$Score
  sz  <- (s - mean(s)) / sd(s)

  cx  <- coxph(Sv ~ sz)
  sm  <- summary(cx)
  cc  <- concordance(Sv ~ s, reverse = TRUE)   # 风险约定
  Cse <- sqrt(cc$var)

  res[[k]] <- data.frame(
    score   = S[[k]]$score_label[1],
    n = n, events = ev,
    HR_perSD = round(sm$conf.int[1, 1], 3),
    CI_low   = round(sm$conf.int[1, 3], 3),
    CI_high  = round(sm$conf.int[1, 4], 3),
    P_cox    = signif(sm$coef[1, 5], 3),
    C_index  = round(cc$concordance, 4),
    C_CI_low = round(cc$concordance - 1.96 * Cse, 4),
    C_CI_high= round(cc$concordance + 1.96 * Cse, 4),
    direction= ifelse(cc$concordance > 0.5, "higher score = worse OS", "REVERSED"),
    stringsAsFactors = FALSE
  )
}
R <- do.call(rbind, res)
rownames(R) <- NULL

cat("\n== TCGA-SARC, endpoint = OS, 风险约定 ==\n")
print(R, row.names = FALSE)

write.csv(R, file.path(outdir, "s156_tcga_score_caliber.csv"),
          row.names = FALSE, fileEncoding = "UTF-8")

# ---- 与两份文档记载的数值对照 ---------------------------------------------
cat("\n== 对照 ==\n")
cat(sprintf("  S4 文档记载        : C = 0.6175 (CI 0.508-0.727), HR/SD = 1.450 (1.006-2.092), P = 0.047\n"))
cat(sprintf("  现稿 / s149        : C = 0.690  (CI 0.562-0.817), HR/SD = 1.98  (1.30-3.00),   P = 0.0013\n"))
for (k in names(S)) {
  r <- res[[k]]
  cat(sprintf("  实测 %s : C = %.4f (CI %.4f-%.4f), HR/SD = %.3f (%.3f-%.3f), P = %s\n",
              k, r$C_index, r$C_CI_low, r$C_CI_high, r$HR_perSD, r$CI_low, r$CI_high, format(r$P_cox)))
}

# ---- A 与 C 之间的 bootstrap 配对比较（复现 SI 的 0.072 / P=0.33 说法）------
mAc <- merge(merge(clin, S$A[, c("key","Score")], by="key"),
                   S$C[, c("key","Score")], by="key", suffixes=c("_A","_C"))
ok  <- !is.na(mAc$OS_time) & !is.na(mAc$OS_status) & mAc$OS_time > 0
mAc <- mAc[ok, ]
cat(sprintf("\n== A 与 C 可比样本数 = %d（事件 %d）==\n", nrow(mAc), sum(mAc$OS_status == 1)))
set.seed(20260925)
B <- 5000
Sv <- Surv(mAc$OS_time, mAc$OS_status)
ca <- concordance(Sv ~ mAc$Score_A, reverse = TRUE)$concordance
cc2<- concordance(Sv ~ mAc$Score_C, reverse = TRUE)$concordance
d  <- numeric(B)
for (i in seq_len(B)) {
  idx <- sample(nrow(mAc), replace = TRUE)
  Svi <- Surv(mAc$OS_time[idx], mAc$OS_status[idx])
  if (sum(mAc$OS_status[idx] == 1) < 3) { d[i] <- NA; next }
  d[i] <- concordance(Svi ~ mAc$Score_A[idx], reverse = TRUE)$concordance -
          concordance(Svi ~ mAc$Score_C[idx], reverse = TRUE)$concordance
}
d <- d[!is.na(d)]
cat(sprintf("  C_A = %.4f | C_C = %.4f | 差 = %.4f | 95%% CI %.4f ~ %.4f | P = %.2f\n",
            ca, cc2, ca - cc2, quantile(d, .025), quantile(d, .975),
            mean(abs(d - mean(d)) >= abs(ca - cc2))))

cat("\n[OK] s156 done\n")
