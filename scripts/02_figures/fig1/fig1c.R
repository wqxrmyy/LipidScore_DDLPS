# ============================================================
# Fig 1C · 6 基因 z 分热图（TCGA n=50）
# 数据：expr6（6 基因 log2CPM，已 z 分）
# ============================================================
suppressWarnings(suppressMessages(library(ggplot2)))
suppressWarnings(suppressMessages(library(reshape2)))
source(file.path(Sys.getenv("LPS_FIG_ROOT", "."), "utils", "colors.R"))
source(file.path(Sys.getenv("LPS_FIG_ROOT", "."), "utils", "themes.R"))

load(file.path(Sys.getenv("LPS_FIG_ROOT", "."), "p5_data.RData"))
OUT <- file.path(Sys.getenv("LPS_FIG_ROOT", "."), "out")

# 每基因 z 分（跨样本）
z <- t(scale(t(expr6)))  # 6 x 50
colnames(z) <- colnames(expr6)
rownames(z) <- rownames(expr6)

# 样本按评分排序
ord <- order(score_A[colnames(z), "Score"])
z <- z[, ord, drop = FALSE]

melt_z <- melt(z)
names(melt_z) <- c("Gene", "Sample", "Z")

p <- ggplot(melt_z, aes(x = Sample, y = Gene, fill = Z)) +
  geom_tile() +
  scale_fill_gradient2(low = "#2980B9", mid = "white", high = "#E74C3C", midpoint = 0) +
  labs(x = "Samples (ordered by score)", y = "Gene", fill = "z-score") +
  theme_publication() +
  theme(axis.text.x = element_blank(), axis.ticks.x = element_blank())

export_fig(p, file.path(OUT, "fig1c.pdf"), w = 8, h = 3)
cat("fig1c done\n")
