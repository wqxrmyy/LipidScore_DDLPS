# ============================================================
# Fig S1A · 评分分布（GSE21122 验证队列）
# 数据：GSE21122_Score.csv（原 Fig S1B，2026-09-28 重构统一为 Fig. S1 A）
# ============================================================
suppressWarnings(suppressMessages(library(ggplot2)))
source(file.path(Sys.getenv("LPS_FIG_ROOT", "."), "utils", "colors.R"))
source(file.path(Sys.getenv("LPS_FIG_ROOT", "."), "utils", "themes.R"))

OUT <- file.path(Sys.getenv("LPS_FIG_ROOT", "."), "out")

load(file.path(Sys.getenv("LPS_FIG_ROOT", "."), "p5_data.RData"))
s <- score_GSE21122$Score

p <- ggplot(data.frame(Score = s), aes(x = Score)) +
  geom_histogram(aes(y = after_stat(density)), bins = 18,
                 fill = "#2980B9", alpha = 0.55, color = "grey20") +
  stat_function(fun = dnorm, args = list(mean = mean(s), sd = sd(s)),
                color = "#E74C3C", size = 1) +
  labs(x = "Lipid metabolism score (A)", y = "Density",
       title = "Score distribution (GSE21122, n = 158)") +
  theme_publication()

export_fig(p, file.path(OUT, "figS1a.pdf"), w = 4, h = 3.5)
cat("figS1a done\n")
