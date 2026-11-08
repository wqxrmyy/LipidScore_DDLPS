suppressWarnings(suppressMessages({ library(survival); library(pROC) }))
OUT <- file.path(Sys.getenv("LPS_WORK_ROOT", "."), "_work")

hr_ci <- function(co, term) {
  s <- summary(co)
  i <- match(term, rownames(s$coef))
  if (is.na(i)) stop(paste5("term not found:", term, paste(rownames(s$coef), collapse = "|")))
  ci <- exp(confint(co))
  out <- unname(c(exp(coef(co))[i], ci[i, 1], ci[i, 2], s$coef[i, 5]))
  out
}
paste5 <- function(...) paste(...)

g <- readRDS(file.path(Sys.getenv("LPS_DATA_ROOT", "."), "2", "raw", "bulk", "GSE30929_clinical_all.rds"))
g$grp <- factor(ifelse(g$LipidScore > median(g$LipidScore), "High", "Low"), levels = c("Low", "High"))
g$s_z <- as.numeric(scale(g$LipidScore))
g$grade <- ifelse(g$subtype %in% c("dedifferentiated", "pleomorphic"), "High-grade", "Low-grade")

rows <- list()
add <- function(scope, model, r, note = "") {
  rows[[length(rows) + 1]] <<- data.frame(
    scope = I(scope), model = I(model),
    HR = round(r["hr"], 4), CI_low = round(r["lo"], 3), CI_high = round(r["hi"], 3),
    P = signif(r["p"], 4),
    HR_reciprocal = round(1 / r[1], 4),
    CI_reciprocal = sprintf("%.2f-%.2f", 1 / r[3], 1 / r[2]),
    note = I(note), stringsAsFactors = FALSE)
}

cat("=========== 1) GSE30929 各层 HR（原方向 High vs Low）===========\n")
subsets <- list(
  "全队列 (Fig 7A)" = g,
  "DDLPS (Fig 7B)" = g[g$subtype == "dedifferentiated", ],
  "pleomorphic" = g[g$subtype == "pleomorphic", ],
  "well-differentiated" = g[g$subtype == "well-differentiated", ],
  "high-grade" = g[g$grade == "High-grade", ],
  "low-grade" = g[g$grade == "Low-grade", ]
)
for (nm in names(subsets)) {
  d <- subsets[[nm]]
  d$grp <- factor(ifelse(d$LipidScore > median(d$LipidScore), "High", "Low"), levels = c("Low", "High"))
  co <- coxph(Surv(OS_time, OS_status) ~ grp, data = d)
  r <- hr_ci(co, "grpHigh")
  add(nm, "Cox: grp", r, sprintf("n=%d ev=%d", nrow(d), sum(d$OS_status)))
  cat(sprintf("  %-22s HR=%.4f (%.3f-%.3f) P=%.5f | 倒数 %.3f (%s)\n",
              nm, r[1], r[2], r[3], r[4], 1 / r[1], sprintf("%.2f-%.2f", 1 / r[3], 1 / r[2])))
}

cat("\n=========== 2) 多因素模型变体 ===========\n")
for (nm in c("grp + subtype", "grp + grade", "score_z + subtype")) {
  f <- switch(nm,
    "grp + subtype" = Surv(OS_time, OS_status) ~ grp + subtype,
    "grp + grade" = Surv(OS_time, OS_status) ~ grp + grade,
    "score_z + subtype" = Surv(OS_time, OS_status) ~ s_z + subtype)
  co <- coxph(f, data = g)
  term <- if (grepl("^grp", nm)) "grpHigh" else "s_z"
  r <- hr_ci(co, term)
  add("全队列", nm, r, "手稿称 2.85(1.47-5.52) P=0.002")
  cat(sprintf("  %-20s HR=%.4f (%.3f-%.3f) P=%.5f | 倒数 %.3f (%s)\n",
              nm, r[1], r[2], r[3], r[4], 1 / r[1], sprintf("%.2f-%.2f", 1 / r[3], 1 / r[2])))
}

cat("\n=========== 3) 评分=组织学代理 的判别力 ===========\n")
a_dd <- as.numeric(auc(roc(as.integer(g$subtype == "dedifferentiated"), g$LipidScore, quiet = TRUE)))
a_wd <- as.numeric(auc(roc(as.integer(g$subtype == "well-differentiated"), g$LipidScore, quiet = TRUE)))
cat(sprintf("  AUC(score -> DDLPS) = %.3f\n  AUC(score -> WD)   = %.3f\n", a_dd, a_wd))
tab <- table(g$ScoreGroup, g$subtype)
cat("  High 组亚型构成:\n"); print(tab["High", ])
cat("  Low  组亚型构成:\n"); print(tab["Low", ])

res <- do.call(rbind, rows)
res$manuscript <- c("3.11 (1.64-5.88) P=0.00049", "0.52 (0.23-1.22) P=0.13",
                    "0.48 P=0.09", "-", "-", "-",
                    "2.85 (1.47-5.52) P=0.002", "-", "-")
write.csv(res, file.path(OUT, "s4_gse30929_hrs.csv"), row.names = FALSE, fileEncoding = "UTF-8")
cat("\n已写出 s4_gse30929_hrs.csv\nDONE\n")
