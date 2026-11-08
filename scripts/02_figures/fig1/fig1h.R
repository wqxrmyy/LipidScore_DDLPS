# ============================================================
# Fig 1H · 三队列评分分布对比（箱线）
# TCGA(50) + GSE21122(158) + GSE159659(45) 的 A 评分
# ============================================================
suppressWarnings(suppressMessages(library(ggplot2)))
source(file.path(Sys.getenv("LPS_FIG_ROOT", "."), "utils", "colors.R"))
source(file.path(Sys.getenv("LPS_FIG_ROOT", "."), "utils", "themes.R"))

load(file.path(Sys.getenv("LPS_FIG_ROOT", "."), "p5_data.RData"))
OUT <- file.path(Sys.getenv("LPS_FIG_ROOT", "."), "out")

# 三队列评分合并（Score 列名需统一）
score_TCGA <- data.frame(Cohort = "TCGA", Score = score_A$Score)
score_21122 <- data.frame(Cohort = "GSE21122", Score = score_GSE21122[[2]])
score_159659 <- data.frame(Cohort = "GSE159659", Score = score_GSE159659[[2]])

df <- rbind(score_TCGA, score_21122, score_159659)
df$Cohort <- factor(df$Cohort, levels = c("TCGA", "GSE21122", "GSE159659"))

p <- ggplot(df, aes(x = Cohort, y = Score, fill = Cohort)) +
  geom_boxplot(outlier.shape = NA, alpha = 0.6) +
  geom_jitter(width = 0.15, size = 0.4, alpha = 0.3) +
  scale_fill_manual(values = c(TCGA = "#E74C3C", GSE21122 = "#2980B9", GSE159659 = "#27AE60")) +
  labs(x = "", y = "Lipid metabolism score", title = "Score across cohorts") +
  theme_publication() +
  theme(legend.position = "none")

export_fig(p, file.path(OUT, "fig1h.pdf"), w = 4.5, h = 3.5)
cat("fig1h done\n")
