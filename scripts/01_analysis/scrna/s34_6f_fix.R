## =============================================================================
##  S34 · 6F 重跑（RICTOR ID 修正）：ENSG00000164327(chr5)  vs  旧 ID(MPRIP,chr17)
##  ---------------------------------------------------------------------------
##  背景：s33_6f.R 使用 RICTOR=ENSG00000133030，该 ID 实为 MPRIP（Ensembl 反查证实）
##        -> D26 的 6F「评分 x RICTOR」实为「评分 x MPRIP」，须重算
##  口径：Q6'=4 层；不建新模型（LIM-P3-001）；中位切分；C-index 报点估计
##  输出：_work/s34_6f_fix_report.txt；fig6_data/km_*_FIX_*.png|pdf
##  注：ASCII-only（Windows R 非 ASCII 曾致脚本失败）
## =============================================================================
options(stringsAsFactors = FALSE)
suppressPackageStartupMessages({ library(survival); library(survminer) })

Root <- file.path(Sys.getenv("LPS_DATA_ROOT", "."), "1")
W    <- file.path(Sys.getenv("LPS_WORK_ROOT", "."), "_work")
FIG  <- file.path(W, "p5_figures/fig6_data"); dir.create(FIG, showWarnings = FALSE)
LOG  <- file.path(W, "s34_6f_fix_report.txt"); if (file.exists(LOG)) file.remove(LOG)
A <- function(...) { s <- paste0(...); cat(s, "\n", sep = ""); cat(s, "\n", sep = "", file = LOG, append = TRUE) }

A("===============================================================================")
A("S34 6F FIX  score x RICTOR(corrected) 4-strata KM   ", format(Sys.time(), "%Y-%m-%d %H:%M:%S"))
A("===============================================================================")

sc <- read.csv(file.path(Root, "data/raw/TCGA_Score.csv"))
cl <- read.csv(file.path(Root, "data/raw/TCGA_clinical_real.csv"), check.names = FALSE)
x  <- readRDS(file.path(Root, "data/raw/TCGA_LPS_real_log2cpm.rds"))
rn <- sub("\\..*$", "", rownames(x))

ENS <- c(RICTOR = "ENSG00000164327",   # corrected (chr5)
         MPRIP  = "ENSG00000133030",   # what the old ID actually was (chr17)
         CD36   = "ENSG00000135218",
         MTOR   = "ENSG00000198793",
         RPTOR  = "ENSG00000141564")   # mTORC1-committed negative control
miss <- setdiff(ENS, rn)
A("gene availability: ", if (!length(miss)) "all present OK" else paste("MISSING:", paste(miss, collapse = ",")))
stopifnot(!length(miss))
E <- as.data.frame(t(x[match(ENS, rn), , drop = FALSE])); colnames(E) <- names(ENS)
E$Sample <- rownames(E)

df <- merge(sc, E, by = "Sample")
df <- merge(df, cl[, c("barcode", "OS_time", "OS_status")], by.x = "Sample", by.y = "barcode")
df <- df[!is.na(df$OS_time) & df$OS_time > 0 & !is.na(df$OS_status), ]
A("n = ", nrow(df), " | events = ", sum(df$OS_status == 1))
A("")

med <- function(v) median(v)
A("-- median split balance --")
for (v in c("Score", "RICTOR", "MPRIP", "CD36", "MTOR", "RPTOR")) {
  a <- sum(df[[v]] > med(df[[v]])); b <- sum(df[[v]] <= med(df[[v]]))
  A(sprintf("  %-7s med=%.4f  >med:%2d  <=med:%2d%s", v, med(df[[v]]), a, b,
            if (abs(a - b) > 2) "   <-- tie imbalance" else ""))
}
A("")

df$S <- factor(ifelse(df$Score  > med(df$Score),  "Score_high", "Score_low"), c("Score_low", "Score_high"))
df$C <- factor(ifelse(df$CD36   > med(df$CD36),   "CD36_high",  "CD36_low"),  c("CD36_low", "CD36_high"))
df$R <- factor(ifelse(df$RICTOR > med(df$RICTOR), "RICTOR_high","RICTOR_low"),c("RICTOR_low","RICTOR_high"))
df$M <- factor(ifelse(df$MPRIP  > med(df$MPRIP),  "MPRIP_high", "MPRIP_low"), c("MPRIP_low","MPRIP_high"))
df$SR <- factor(paste(df$S, df$R, sep = " x "))
df$SC <- factor(paste(df$S, df$C, sep = " x "))
df$SM <- factor(paste(df$S, df$M, sep = " x "))

## ---------------------------------------------------------------- KM
km_tab <- function(g, lab) {
  A(sprintf("---- %s ----", lab))
  f <- as.formula(paste("Surv(OS_time, OS_status) ~", g))
  fit <- survfit(f, data = df)
  ev <- tapply(df$OS_status, df[[g]], sum); nn <- tapply(df$OS_status, df[[g]], length)
  med_t <- summary(fit)$table
  o <- data.frame(stratum = names(nn), n = as.integer(nn), events = as.integer(ev),
                  median_OS = round(as.numeric(med_t[, "median"]), 1), row.names = NULL)
  print(o, row.names = FALSE)
  for (i in seq_len(nrow(o)))
    A(sprintf("  %-24s n=%2d  events=%2d  medianOS=%s", o$stratum[i], o$n[i], o$events[i],
              ifelse(is.na(o$median_OS[i]), "NA", format(o$median_OS[i]))))
  sd <- survdiff(f, data = df)
  p <- 1 - pchisq(sd$chisq, length(sd$n) - 1)
  A(sprintf("  log-rank chisq=%.3f  df=%d  P=%.5f", sd$chisq, length(sd$n) - 1, p))
  A("")
  invisible(list(tab = o, p = p))
}

s2 <- km_tab("S", "score A (2 strata)")
s2r <- km_tab("R", "RICTOR corrected (2 strata)")
s2c <- km_tab("C", "CD36 (2 strata)")
s2m <- km_tab("M", "MPRIP = old-ID (2 strata)")
r4 <- km_tab("SR", "score A x RICTOR(corrected) (4 strata)")
c4 <- km_tab("SC", "score A x CD36 (4 strata)")
m4 <- km_tab("SM", "score A x MPRIP(old ID) (4 strata)")

## ---------------------------------------------------------------- uni/multi Cox
A("-- univariate Cox (per 1 SD of log2CPM) --")
A(sprintf("  %-10s %-9s %-22s %-9s %s", "var", "HR", "95%CI", "P", "C-index"))
uv <- list()
for (v in c("Score", "RICTOR", "MPRIP", "CD36", "MTOR", "RPTOR")) {
  d <- df; d$z <- as.numeric(scale(d[[v]]))
  f <- coxph(Surv(OS_time, OS_status) ~ z, data = d)
  s <- summary(f)
  uv[[v]] <- c(HR = s$conf.int[1, "exp(coef)"], lo = s$conf.int[1, "lower .95"],
               hi = s$conf.int[1, "upper .95"], P = s$coef[1, "Pr(>|z|)"], C = s$concordance[1])
  A(sprintf("  %-10s %-9.3f %-22s %-9.4g %.4f", v, uv[[v]]["HR"],
            sprintf("[%.3f, %.3f]", uv[[v]]["lo"], uv[[v]]["hi"]), uv[[v]]["P"], uv[[v]]["C"]))
}
A("")
A("-- multivariable Cox (z-scaled) --")
mv <- function(vars, lab) {
  d <- df
  for (v in vars) d[[paste0("z_", v)]] <- as.numeric(scale(d[[v]]))
  f <- as.formula(paste("Surv(OS_time, OS_status) ~", paste(paste0("z_", vars), collapse = " + ")))
  s <- summary(coxph(f, data = d))
  A(sprintf("  %-26s C-index=%.4f", lab, s$concordance[1]))
  for (i in seq_len(nrow(s$coef)))
    A(sprintf("      %-14s HR=%.3f [%.3f, %.3f]  P=%.4g", rownames(s$coef)[i],
              s$conf.int[i, 1], s$conf.int[i, 3], s$conf.int[i, 4], s$coef[i, 5]))
  s$concordance[1]
}
c_score <- mv("Score", "Score")
c_sc    <- mv(c("Score", "CD36"), "Score + CD36")
c_sr    <- mv(c("Score", "RICTOR"), "Score + RICTOR(corrected)")
c_smr   <- mv(c("Score", "MPRIP"), "Score + MPRIP(old ID)")
c_all   <- mv(c("Score", "CD36", "RICTOR"), "Score + CD36 + RICTOR(corrected)")
A("")
A(sprintf("C-index: Score=%.4f | Score+CD36=%.4f | Score+RICTOR(corr)=%.4f | Score+MPRIP(old)=%.4f | Score+CD36+RICTOR(corr)=%.4f",
          c_score, c_sc, c_sr, c_smr, c_all))

## ---------------------------------------------------------------- figures
mk <- function(g, lab, fn, cols) {
  f <- as.formula(paste("Surv(OS_time, OS_status) ~", g))
  fit <- survfit(f, data = df)
  fit$call$formula <- f
  fit$call$data <- quote(df)
  p <- ggsurvplot(fit, data = df, pval = TRUE, risk.table = TRUE, conf.int = FALSE,
                  palette = cols, legend.title = "", legend.labs = sub(".*=", "", names(fit$strata)),
                  xlab = "Time (days)", ylab = "Overall survival", censor.size = 2,
                  ggtheme = theme_classic(base_size = 7), tables.height = 0.30, fontsize = 2.8)
  pdf(file.path(FIG, paste0(fn, ".pdf")), width = 2.7, height = 2.8); print(p, newpage = FALSE); dev.off()
  png(file.path(FIG, paste0(fn, ".png")), width = 2.7, height = 2.8, units = "in", res = 300)
  print(p, newpage = FALSE); dev.off()
  A(sprintf("  wrote %s.pdf/.png", fn))
}
A("")
A("-- figures --")
mk("SR", "score x RICTOR corrected", "km_score_rictorFIX_4", c("#4C72B0", "#DD8452", "#55A868", "#C44E52"))
mk("SC", "score x CD36", "km_score_cd36_4", c("#4C72B0", "#DD8452", "#55A868", "#C44E52"))
mk("R", "RICTOR corrected 2 strata", "km_rictorFIX_2", c("#4C72B0", "#C44E52"))

A("")
A("DONE")
