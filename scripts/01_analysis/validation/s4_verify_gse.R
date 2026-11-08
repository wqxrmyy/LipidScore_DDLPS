suppressWarnings(suppressMessages({ library(survival); library(pROC) }))
out_dir <- file.path(Sys.getenv("LPS_WORK_ROOT", "."), "_work")

# ============ 1. 多因素模型变体：手稿 2.85 (1.47-5.52) P=0.002 能否复现 ============
g <- readRDS(file.path(Sys.getenv("LPS_DATA_ROOT", "."), "2", "raw", "bulk", "GSE30929_clinical_all.rds"))
g$grp <- factor(ifelse(g$LipidScore > median(g$LipidScore), "High", "Low"), levels = c("Low", "High"))
g$s_z <- as.numeric(scale(g$LipidScore))
g$grade <- ifelse(g$subtype %in% c("dedifferentiated", "pleomorphic"), "High-grade", "Low-grade")

cat("===== 1. 多因素模型变体（手稿值 2.85 / 1.47-5.52 / P=0.002）=====\n")
cat("  参照: High vs Low; 若手稿是 Low vs High 则取倒数\n")
models <- list(
  "grp + subtype"          = Surv(OS_time, OS_status) ~ grp + subtype,
  "grp + grade"            = Surv(OS_time, OS_status) ~ grp + grade,
  "grp + subtype + grade"  = Surv(OS_time, OS_status) ~ grp + subtype + grade,
  "score_z + subtype"      = Surv(OS_time, OS_status) ~ s_z + subtype,
  "grp + subtype + tissue" = NULL
)
for (nm in c("grp + subtype", "grp + grade", "grp + subtype + grade", "score_z + subtype")) {
  f <- models[[nm]]
  co <- tryCatch(coxph(f, data = g), error = function(e) NULL)
  if (is.null(co)) { cat(sprintf("  %-24s 拟合失败\n", nm)); next }
  s <- summary(co)
  rn <- grep("^grp|^s_z", rownames(s$coef))
  hr <- exp(s$coef[rn[1], 1]); lo <- exp(s$conf.int[rn[1], 3]); hi <- exp(s$conf.int[rn[1], 4]); p <- s$coef[rn[1], 5]
  cat(sprintf("  %-24s HR=%.4f (%.3f-%.3f) P=%.5f | 倒数 %.3f (%.2f-%.2f)\n",
              nm, hr, lo, hi, p, 1 / hr, 1 / hi, 1 / lo))
}

cat("\n===== 2. 评分是否为组织学类型的代理变量 =====\n")
g$isDD <- as.integer(g$subtype == "dedifferentiated")
g$isWD <- as.integer(g$subtype == "well-differentiated")
auc_dd <- as.numeric(auc(roc(g$isDD, g$LipidScore, quiet = TRUE)))
auc_wd <- as.numeric(auc(roc(g$isWD, g$LipidScore, quiet = TRUE)))
cat(sprintf("  AUC(评分 判别 DDLPS 与否)   = %.3f\n", auc_dd))
cat(sprintf("  AUC(评分 判别 WD 与否)      = %.3f\n", auc_wd))
cat("  => 0.5 表示与组织学无关；偏离越远表示评分越像组织学代理\n")
cat("\n  TCGA 侧对照（全为 DDLPS，无亚型混杂）:\n")
sc <- read.csv(file.path(Sys.getenv("LPS_DATA_ROOT", "."), "1", "data", "results", "02_feature_scores.csv"), row.names = 1)
cl <- read.csv(file.path(Sys.getenv("LPS_DATA_ROOT", "."), "1", "data", "results", "03_clinical_outcomes.csv"), row.names = 1)
tc <- cl[cl$Cohort == "TCGA", ]
cat("    TCGA 亚型分布:", paste(names(table(tc$Subtype)), table(tc$Subtype), sep = "=", collapse = " "), "\n")

cat("\n===== 3. 各亚组 HR（原方向 High vs Low）与手稿值对照 =====\n")
sub <- function(d, lab) {
  d$grp <- factor(ifelse(d$LipidScore > median(d$LipidScore), "High", "Low"), levels = c("Low", "High"))
  if (length(unique(d$grp)) < 2 || sum(d$OS_status) < 3) { cat(sprintf("  %-22s 不可估\n", lab)); return(invisible()) }
  co <- coxph(Surv(OS_time, OS_status) ~ grp, data = d)
  s <- summary(co)
  hr <- exp(coef(co))[1]; lo <- exp(s$conf.int[1, 3]); hi <- exp(s$conf.int[1, 4])
  # 每 +1SD 连续
  d$s_z <- as.numeric(scale(d$LipidScore))
  co2 <- coxph(Surv(OS_time, OS_status) ~ s_z, data = d)
  cat(sprintf("  %-22s n=%3d ev=%2d | 分组 HR=%.3f (%.2f-%.2f) P=%.4f | 每SD HR=%.3f P=%.4f\n",
              lab, nrow(d), sum(d$OS_status), hr, lo, hi, s$coef[1, 5],
              exp(coef(co2))[1], summary(co2)$coef[1, 5]))
}
sub(g[g$subtype == "pleomorphic", ], "pleomorphic")
sub(g[g$subtype == "dedifferentiated", ], "DDLPS")
sub(g[g$subtype == "well-differentiated", ], "well-differentiated")
cat("  手稿 3.3: pleomorphic HR=0.48 P=0.09 ; DDLPS HR=0.52 P=0.13\n")
cat("  手稿 3.4 三分位: 见 Fig 7D（图中 High 曲线在上方）\n")

cat("\n===== 4. 关键结论汇总 =====\n")
cat("  · 全队列 HR: 图 7A 显示 High 组生存更好(保护) -> 真实 HR(High vs Low)=0.3218\n")
cat("    手稿文本写 3.11 = 1/0.3218 的精确倒数，P/CI 完全一致 -> 方向被颠倒\n")
cat("  · 校正亚型后本机 HR=0.591 P=0.177（不显著），手稿称 HR=2.85 P=0.002 -> 不可复现\n")
cat("  · 评分在 GSE30929 中高度富集于 WD（High 组 46/52 为 WD）-> 组织学代理\n")
cat("DONE\n")
