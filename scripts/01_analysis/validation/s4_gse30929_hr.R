suppressWarnings(suppressMessages(library(survival)))
g <- readRDS(file.path(Sys.getenv("LPS_DATA_ROOT", "."), "2", "raw", "bulk", "GSE30929_clinical_all.rds"))

show <- function(d, lab) {
  d$grp <- factor(ifelse(d$LipidScore > median(d$LipidScore), "High", "Low"),
                  levels = c("Low", "High"))
  if (length(unique(d$grp)) < 2 || sum(d$OS_status) < 3) {
    cat(sprintf("  %-24s 不可估\n", lab)); return(invisible())
  }
  co <- coxph(Surv(OS_time, OS_status) ~ grp, data = d)
  s <- summary(co)
  hr <- exp(coef(co))[1]; lo <- exp(confint(co))[1, 1]; hi <- exp(confint(co))[1, 2]
  cat(sprintf("  %-24s n=%3d ev=%2d  HR(High vs Low)=%.4f (%.3f-%.3f) P=%.5f | 倒数 HR=%.3f (%.2f-%.2f) P同\n",
              lab, nrow(d), sum(d$OS_status), hr, lo, hi, s$coef[1, 5], 1 / hr, 1 / hi, 1 / lo))
}

cat("===== GSE30929 中位分组 Cox（High vs Low）=====\n")
show(g, "全队列")
show(g[g$subtype == "dedifferentiated", ], "DDLPS")
show(g[g$subtype == "pleomorphic", ], "pleomorphic")
show(g[g$subtype == "well-differentiated", ], "well-differentiated")
show(g[g$subtype == "myxoid", ], "myxoid")
show(g[g$subtype %in% c("myxoid", "myxoid/round cell"), ], "myxoid+roundcell")
g$grade <- ifelse(g$subtype %in% c("dedifferentiated", "pleomorphic"), "High-grade", "Low-grade")
show(g[g$grade == "High-grade", ], "High-grade")
show(g[g$grade == "Low-grade", ], "Low-grade")

cat("\n===== 与手稿值对照 =====\n")
cat("  手稿 3.3 节/Fig 7A 图注: 全队列 HR=3.11 (1.64-5.88) P=0.00049\n")
cat("  手稿 3.3 节:           多因素校正亚型 HR=2.85 (1.47-5.52) P=0.002\n")
g$grp <- factor(ifelse(g$LipidScore > median(g$LipidScore), "High", "Low"),
                levels = c("Low", "High"))
co2 <- coxph(Surv(OS_time, OS_status) ~ grp + subtype, data = g)
s2 <- summary(co2)
cat(sprintf("  本机 多因素(grp+subtype): HR=%.4f (%.3f-%.3f) P=%.5f | 倒数 %.3f (%.2f-%.2f)\n",
            exp(coef(co2))[1], exp(confint(co2))[1, 1], exp(confint(co2))[1, 2], s2$coef[1, 5],
            1 / exp(coef(co2))[1], 1 / exp(confint(co2))[1, 2], 1 / exp(confint(co2))[1, 1]))

cat("\n===== rds 自带 ScoreGroup 与亚型交叉 =====\n")
print(table(g$ScoreGroup, g$subtype))
cat("\n===== 手稿 3.3 各亚组值 vs 本机 =====\n")
cat("  手稿: pleomorphic HR=0.48 P=0.09 ; DDLPS HR=0.52 P=0.13 （被写成 consistent trends）\n")
cat("DONE\n")
