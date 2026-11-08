## =============================================================================
##  S33 · Fig6 6F：评分 A x CD36 四层 KM + 多因素（Q6' = 4 层，不建新模型）
##  ---------------------------------------------------------------------------
##  依据：D26 §一 Q6'；LIM-P3-001（n=50/事件 23，不建 RSF+LASSO 新模型）
##  输入：TCGA_Score.csv / TCGA_clinical_real.csv / TCGA_LPS_real_log2cpm.rds
##  基因：ENSG 定位（已核验）RICTOR=ENSG00000133030 CD36=ENSG00000135218
##  输出：_work/s33_6f_*.csv / s33_6f_report.txt / fig6_data/km_*.pdf|png
##  注：★ CD36 单变量与 OS 无关（P=0.456，见 D26 §4.5）→ CD36 维度是否有区分须实测
## =============================================================================
options(stringsAsFactors = FALSE)
suppressPackageStartupMessages({ library(survival); library(survminer) })

Root <- file.path(Sys.getenv("LPS_DATA_ROOT", "."), "1")
W <- file.path(Sys.getenv("LPS_WORK_ROOT", "."), "_work")
FIG <- file.path(W, "p5_figures/fig6_data"); dir.create(FIG, showWarnings = FALSE)
LOG <- file.path(W, "s33_6f_report.txt"); if (file.exists(LOG)) file.remove(LOG)
A <- function(...) { s <- paste0(...); cat(s, "\n", sep = ""); cat(s, "\n", sep = "", file = LOG, append = TRUE) }

A("===============================================================================")
A("S33 6F  评分 A x CD36 四层 KM（Q6'=4 层）  ", format(Sys.time(), "%Y-%m-%d %H:%M:%S"))
A("===============================================================================")

sc <- read.csv(file.path(Root, "data/raw/TCGA_Score.csv"))
cl <- read.csv(file.path(Root, "data/raw/TCGA_clinical_real.csv"), check.names = FALSE)
x  <- readRDS(file.path(Root, "data/raw/TCGA_LPS_real_log2cpm.rds"))
rn <- sub("\\..*$", "", rownames(x))
ENS <- c(RICTOR = "ENSG00000133030", CD36 = "ENSG00000135218", MTOR = "ENSG00000198793")
stopifnot(all(ENS %in% rn))
E <- as.data.frame(t(x[match(ENS, rn), , drop = FALSE])); colnames(E) <- names(ENS)
E$Sample <- rownames(E)

df <- merge(sc, E, by = "Sample")
df <- merge(df, cl[, c("barcode", "OS_time", "OS_status", "paper_purity", "paper_CN.cluster")],
            by.x = "Sample", by.y = "barcode")
df <- df[!is.na(df$OS_time) & df$OS_time > 0 & !is.na(df$OS_status), ]
A("n = ", nrow(df), "   事件 = ", sum(df$OS_status == 1))
A("")

## ---------------------------------------------------------------- 分层口径
A("-- 分组口径（中位切分；RICTOR/MTOR 检查并列导致的失衡）--")
med <- function(v) median(v)
for (v in c("Score", "CD36", "RICTOR", "MTOR")) {
  n_gt <- sum(df[[v]] > med(df[[v]])); n_le <- sum(df[[v]] <= med(df[[v]]))
  A(sprintf("  %-8s median=%.4f  > median: %2d   <= median: %2d %s", v, med(df[[v]]), n_gt, n_le,
            if (abs(n_gt - n_le) > 2) "  ★ 并列导致失衡" else ""))
}
A("")

df$S <- ifelse(df$Score  > med(df$Score),  "Score_high", "Score_low")
df$C <- ifelse(df$CD36   > med(df$CD36),   "CD36_high",  "CD36_low")
df$R <- ifelse(df$RICTOR > med(df$RICTOR), "RICTOR_high","RICTOR_low")
df$S <- factor(df$S, levels = c("Score_low", "Score_high"))
df$C <- factor(df$C, levels = c("CD36_low", "CD36_high"))
df$R <- factor(df$R, levels = c("RICTOR_low", "RICTOR_high"))

## ---------------------------------------------------------------- KM 工具
km_tab <- function(g, lab, tag) {
  A(sprintf("---- %s ----", lab))
  f <- as.formula(paste("Surv(OS_time, OS_status) ~", g))
  fit <- survfit(f, data = df)
  st <- summary(fit)$table          # 行名形如 "S=Score_low"；列含 records/events/median
  if (is.null(dim(st))) st <- t(as.matrix(st))
  med_os <- as.numeric(st[, "median"])
  tab <- data.frame(stratum = sub("^[^=]*=", "", rownames(st)),
                    n = as.integer(st[, "records"]),
                    events = as.integer(st[, "events"]),
                    median_OS = round(med_os, 1), row.names = NULL)
  print(tab, row.names = FALSE); A(paste(capture.output(print(tab, row.names = FALSE)), collapse = "\n"))
  lr <- survdiff(f, data = df)
  pv <- 1 - pchisq(lr$chisq, length(lr$n) - 1)
  A(sprintf("  log-rank chisq = %.3f  df = %d  P = %.4g", lr$chisq, length(lr$n) - 1, pv))
  A("")
  list(tab = tab, p = pv, chisq = lr$chisq)
}

A("===============================================================================")
A("-- [1] 基线：评分 A 2 层（现稿已用）--")
A("===============================================================================")
r2 <- km_tab("S", "Score A (2 strata)", "s2")

A("===============================================================================")
A("-- [2] ★ 主分析：评分 A x CD36 四层（Q6' 裁定）--")
A("===============================================================================")
df$SC <- factor(paste(df$S, df$C, sep = " | "),
                levels = c("Score_low | CD36_low", "Score_low | CD36_high",
                           "Score_high | CD36_low", "Score_high | CD36_high"))
r4 <- km_tab("SC", "Score A x CD36 (4 strata)", "s4")

A("===============================================================================")
A("-- [3] 对照（因 CD36 单变量不显著）：评分 A x RICTOR 四层 --")
A("===============================================================================")
df$SR <- factor(paste(df$S, df$R, sep = " | "),
                levels = c("Score_low | RICTOR_low", "Score_low | RICTOR_high",
                           "Score_high | RICTOR_low", "Score_high | RICTOR_high"))
r4b <- km_tab("SR", "Score A x RICTOR (4 strata)", "s4b")

A("===============================================================================")
A("-- [4] 多因素 Cox（连续变量；QR 检查）--")
A("===============================================================================")
z <- function(v) as.numeric(scale(v))
for (vars in list(c("Score"), c("Score", "CD36"), c("Score", "RICTOR"), c("Score", "CD36", "RICTOR"))) {
  d <- df
  for (v in vars) d[[paste0("z_", v)]] <- z(d[[v]])
  fs <- paste("Surv(OS_time, OS_status) ~", paste(paste0("z_", vars), collapse = " + "))
  fit <- coxph(as.formula(fs), data = d)
  A(sprintf("  %s", fs))
  s <- summary(fit)
  for (i in seq_len(nrow(s$conf.int))) {
    A(sprintf("    %-12s HR=%.3f  95%%CI [%.3f, %.3f]  P=%.4g",
              rownames(s$conf.int)[i], s$conf.int[i, 1], s$conf.int[i, 3], s$conf.int[i, 4], s$coefficients[i, 5]))
  }
  A(sprintf("    (n=%d, events=%d, C-index=%.4f, LR test P=%.4g)",
            fit$n, fit$nevent, unname(s$concordance[1]), s$logtest["pvalue"]))
  A("")
}

A("===============================================================================")
A("-- [5] ★ 四层分层的分层效力（每格 n / 事件）与最小格判定 --")
A("===============================================================================")
for (g in c("SC", "SR")) {
  t <- table(df[[g]], df$OS_status)
  nn <- tapply(df$OS_status, df[[g]], length); ev <- tapply(df$OS_status, df[[g]], sum)
  A(sprintf("  [%s] 最小 n = %d, 最小事件 = %d  → %s", g, min(nn), min(ev),
            if (min(ev) < 5) "★ 效力偏紧，须在图中标注每格 n/事件" else "可用"))
}
A("")

## ---------------------------------------------------------------- 出图
A("-- [6] 出图（KM；先出数字版，风格统一于 Fig5 主题）--")
mk <- function(g, lab, fn, cols) {
  f <- as.formula(paste("Surv(OS_time, OS_status) ~", g))
  fit <- survfit(f, data = df)
  ## ★ survminer 会回溯 fit$call$formula；函数内以符号 f 建 fit 会丢引用 → 显式回填
  fit$call$formula <- f
  fit$call$data <- quote(df)
  nl <- length(fit$strata)
  if (length(cols) < nl) cols <- rep(cols, length.out = nl)
  p <- ggsurvplot(fit, data = df, pval = TRUE, risk.table = TRUE, conf.int = FALSE,
                  palette = cols[seq_len(nl)],
                  legend.title = "", legend.labs = sub(".*=", "", names(fit$strata)),
                  xlab = "Time (days)", ylab = "Overall survival", censor.size = 2,
                  ggtheme = theme_classic(base_size = 7), tables.height = 0.28,
                  fontsize = 2.8)
  pdf(file.path(FIG, paste0(fn, ".pdf")), width = 2.6, height = 2.7)
  print(p, newpage = FALSE); dev.off()
  png(file.path(FIG, paste0(fn, ".png")), width = 2.6, height = 2.7, units = "in", res = 300)
  print(p, newpage = FALSE); dev.off()
  A("  ", fn, ".pdf/.png 写出")
}
mk("S",  "Score A (2)", "km_score2", c("#2C6E9B", "#B00000"))
mk("SC", "Score x CD36 (4)", "km_score_cd36_4", c("#2C6E9B", "#7FB3D5", "#E08A7A", "#B00000"))
mk("SR", "Score x RICTOR (4)", "km_score_rictor_4", c("#2C6E9B", "#7FB3D5", "#E08A7A", "#B00000"))
A("")

## ---------------------------------------------------------------- 导出
write.csv(r2$tab, file.path(W, "s33_6f_km_2strata.csv"), row.names = FALSE)
write.csv(r4$tab, file.path(W, "s33_6f_km_score_cd36_4.csv"), row.names = FALSE)
write.csv(r4b$tab, file.path(W, "s33_6f_km_score_rictor_4.csv"), row.names = FALSE)
write.csv(df[, c("Sample", "Score", "CD36", "RICTOR", "MTOR", "S", "C", "R", "SC", "SR",
                 "OS_time", "OS_status", "paper_purity", "paper_CN.cluster")],
          file.path(W, "s33_6f_data.csv"), row.names = FALSE)
A("写出 s33_6f_km_*.csv / s33_6f_data.csv")
A("")
A(sprintf("★ 汇总：2层 P=%.3g | 4层(S x CD36) P=%.3g | 4层(S x RICTOR) P=%.3g",
          r2$p, r4$p, r4b$p))
A("DONE_S33_6F")
