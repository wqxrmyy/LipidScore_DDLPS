# ============================================================
# Fig 1E · 评分 A 分布（直方图 + 密度）
# ============================================================
suppressWarnings(suppressMessages(library(ggplot2)))
source(file.path(Sys.getenv("LPS_FIG_ROOT", "."), "utils", "colors.R"))
source(file.path(Sys.getenv("LPS_FIG_ROOT", "."), "utils", "themes.R"))

load(file.path(Sys.getenv("LPS_FIG_ROOT", "."), "p5_data.RData"))
OUT <- file.path(Sys.getenv("LPS_FIG_ROOT", "."), "out")

df <- data.frame(Score = score_A$Score)

p <- ggplot(df, aes(x = Score)) +
  geom_histogram(aes(y = ..density..), bins = 15, fill = "#8E44AD", alpha = 0.4, color = "grey30") +
  geom_density(color = "#8E44AD", size = 1) +
  geom_vline(xintercept = median(df$Score), linetype = "dashed", color = "red") +
  labs(x = "Lipid metabolism score (A)", y = "Density", title = "Score distribution (TCGA, n=50)") +
  theme_publication()

export_fig(p, file.path(OUT, "fig1e.pdf"), w = 4, h = 3.5)
cat("fig1e done\n")
