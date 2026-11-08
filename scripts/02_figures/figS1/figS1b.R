# ============================================================
# Fig S1B · 评分分布（GSE159659 验证队列）
# 数据：GSE159659_Score.csv（原 Fig S1C，2026-09-28 重构统一为 Fig. S1 B）
# ============================================================
suppressWarnings(suppressMessages(library(ggplot2)))
source(file.path(Sys.getenv("LPS_FIG_ROOT", "."), "utils", "colors.R"))
source(file.path(Sys.getenv("LPS_FIG_ROOT", "."), "utils", "themes.R"))

OUT <- file.path(Sys.getenv("LPS_FIG_ROOT", "."), "out")

load(file.path(Sys.getenv("LPS_FIG_ROOT", "."), "p5_data.RData"))
s <- score_GSE159659$Score

p <- ggplot(data.frame(Score = s), aes(x = Score)) +
  geom_histogram(aes(y = after_stat(density)), bins = 12,
                 fill = "#2980B9", alpha = 0.55, color = "grey20") +
  stat_function(fun = dnorm, args = list(mean = mean(s), sd = sd(s)),
                color = "#E74C3C", size = 1) +
  labs(x = "Lipid metabolism score (A)", y = "Density",
       title = "Score distribution (GSE159659, n = 45)") +
  theme_publication()

export_fig(p, file.path(OUT, "figS1b.pdf"), w = 4, h = 3.5)
cat("figS1b done\n")
