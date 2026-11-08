# ============================================================
# Fig 4A · DD vs WD 评分箱线（官方口径 score_all53）
# 官方值：AUC = 0.2491 (P = 0.0131)；Δ = -0.469；Welch P = 0.025
# 数据：fig4_score.csv（prep_fig4_gse_official.py 产出并断言校验）
# ============================================================
suppressWarnings(suppressMessages(library(ggplot2)))
source(file.path(Sys.getenv("LPS_FIG_ROOT", "."), "utils", "colors.R"))
source(file.path(Sys.getenv("LPS_FIG_ROOT", "."), "utils", "themes.R"))

OUT <- file.path(Sys.getenv("LPS_FIG_ROOT", "."), "out")

d <- read.csv(file.path(Sys.getenv("LPS_FIG_ROOT", "."), "fig4_score.csv"),
              stringsAsFactors = FALSE)
d <- d[d$group %in% c("DD", "WD"), ]
d$group <- factor(d$group, levels = c("WD", "DD"))

# 官方统计值（d11 B5_selfcheck，prep 已断言复算一致）
auc <- 0.249; auc_p <- 0.0131; delta <- -0.469; welch_p <- 0.025

p <- ggplot(d, aes(x = group, y = score, fill = group)) +
  geom_boxplot(outlier.shape = NA, alpha = 0.6, width = 0.5) +
  geom_jitter(width = 0.12, size = 1.3, alpha = 0.8, color = "grey20") +
  scale_fill_manual(values = c(WD = "#2980B9", DD = "#E67E22")) +
  labs(x = "", y = "Lipid metabolism score (A)",
       title = "DD vs WD (GSE221492, paired design)",
       subtitle = sprintf("AUC = %.3f (P = %.4g)\nmean delta = %.3f, Welch P = %.3f",
                          auc, auc_p, delta, welch_p)) +
  theme_publication() +
  theme(legend.position = "none")

export_fig(p, file.path(OUT, "fig4a.pdf"), w = 4.2, h = 4.2)
cat("fig4a done\n")
