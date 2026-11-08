# s151_gse30929_subgroup_A.R
# 目的：在现稿采用的 A 基准（6 基因 z 分均值，与 s144 完全同源）上，
#       重算 GSE30929 的亚组读数 —— 因为 s150/CSV 那套（r=0.776 的旧向量）
#       不能与现稿混用。
# 产出：_work/s148_out/s151_gse30929_subgroup_A.csv（幂等覆盖）
# 口径：C-index 风险约定（reverse = TRUE，高分 = 风险高）

suppressPackageStartupMessages(library(survival))
B   <- file.path(Sys.getenv("LPS_DATA_ROOT", "."), "2", "raw", "bulk")
OUT <- file.path(Sys.getenv("LPS_WORK_ROOT", "."), "_work", "s148_out", "s151_gse30929_subgroup_A.csv")

e  <- readRDS(file.path(B, "GSE30929_6genes_expr.rds"))   # 6 x 140
cl <- readRDS(file.path(B, "GSE30929_clinical_all.rds"))  # 140 x 8

z <- t(scale(t(e)))
A <- colMeans(z, na.rm = TRUE)
cl$A_score <- as.numeric(A[cl$sample])

cat("=== 基准 ===\n")
cat(sprintf("A(重建) vs rds LipidScore: r = %.4f  [与 s144 一致即视为同基准]\n",
            cor(cl$A_score, cl$LipidScore, use = "complete.obs")))

fitset <- function(x, label) {
  x <- x[is.finite(x$OS_time) & is.finite(x$OS_status) & is.finite(x$A_score), ]
  ev  <- sum(x$OS_status == 1)
  grp <- factor(ifelse(x$A_score > median(x$A_score), "High", "Low"), levels = c("Low", "High"))
  cx  <- coxph(Surv(OS_time, OS_status) ~ grp, data = x)
  s   <- summary(cx)
  hr  <- unname(s$conf.int["grpHigh", "exp(coef)"])
  lo  <- unname(s$conf.int["grpHigh", "lower .95"])
  hi  <- unname(s$conf.int["grpHigh", "upper .95"])
  p   <- unname(s$coefficients["grpHigh", "Pr(>|z|)"])
  cc  <- concordance(Surv(OS_time, OS_status) ~ A_score, data = x, reverse = TRUE)
  se  <- sqrt(cc$var)
  cat(sprintf("%-18s n=%3d ev=%2d | HR=%.4f (%.3f-%.3f) P=%.4f | C=%.4f (%.3f-%.3f)\n",
              label, nrow(x), ev, hr, lo, hi, p, cc$concordance,
              cc$concordance - 1.96 * se, cc$concordance + 1.96 * se))
  data.frame(layer = label, n = nrow(x), events = ev,
             HR_high_vs_low = round(hr, 4), lo = round(lo, 3), hi = round(hi, 3),
             P = signif(p, 4), C = round(cc$concordance, 4),
             C_lo = round(cc$concordance - 1.96 * se, 3),
             C_hi = round(cc$concordance + 1.96 * se, 3),
             stringsAsFactors = FALSE)
}

cat("\n=== A 基准下的分层读数 ===\n")
res <- list()
res[["all"]] <- fitset(cl, "A: all (n=140)")
res[["dd"]]  <- fitset(cl[cl$subtype == "dedifferentiated", ], "A: DDLPS (n=40)")
res[["wd"]]  <- fitset(cl[grepl("well", cl$subtype, ignore.case = TRUE), ], "A: WD (n=52)")
R <- do.call(rbind, res)

cat("\n=== 亚组 x 分组 ===\n")
cl$grp <- factor(ifelse(cl$A_score > median(cl$A_score), "High", "Low"), levels = c("Low", "High"))
print(table(cl$subtype, cl$grp))

cat("\n=== 对比：同一亚组在旧 CSV 向量上（仅供诊断，不入稿） ===\n")
old <- read.csv(file.path(B, "GSE30929_Score_all.csv"), stringsAsFactors = FALSE)
olddd <- old[old$subtype == "dedifferentiated", ]
cxo <- coxph(Surv(as.numeric(OS_time), as.numeric(OS_status)) ~
               factor(ifelse(LipidScore > median(LipidScore), "High", "Low"),
                      levels = c("Low", "High")), data = olddd)
cco <- concordance(Surv(as.numeric(olddd$OS_time), as.numeric(olddd$OS_status)) ~ olddd$LipidScore,
                   reverse = TRUE)
cat(sprintf("旧CSV 向量 DDLPS: HR=%.4f  C=%.4f   (s150 复现值 HR=1.9071 C=0.5527)\n",
            unname(summary(cxo)$conf.int["factor(ifelse(LipidScore > median(LipidScore), \"High\", \"Low\"), levels = c(\"Low\", \"High\"))High", "exp(coef)"]),
            cco$concordance))

dir.create(dirname(OUT), recursive = TRUE, showWarnings = FALSE)
if (file.exists(OUT)) file.remove(OUT)
write.csv(R, OUT, row.names = FALSE)
cat("\n落盘:", OUT, "\n")
