# ============================================================
# Fig S2A · 分量基因 DD vs WD 小提琴（8 基因并列，官方 log2CPM 口径）
# 数据：figS3b_gse_data.csv（prep_figS3b_gse.py，raw counts -> log2CPM，
#       CD36 均值与 Table S5 完全一致：DD 7.228 / WD 10.536）
# （原 Fig S3B，2026-09-28 重构统一为 Fig. S2 A）
# ============================================================
suppressWarnings(suppressMessages(library(ggplot2)))
source(file.path(Sys.getenv("LPS_FIG_ROOT", "."), "utils", "colors.R"))
source(file.path(Sys.getenv("LPS_FIG_ROOT", "."), "utils", "themes.R"))

OUT <- file.path(Sys.getenv("LPS_FIG_ROOT", "."), "out")

d <- read.csv(file.path(Sys.getenv("LPS_FIG_ROOT", "."), "figS3b_gse_data.csv"),
              stringsAsFactors = FALSE, check.names = FALSE)
genes8 <- c("SREBF1", "SCD", "HMGCR", "ACLY", "PPARG", "FASN", "CD36", "CD8A")
d <- d[d$gene %in% genes8 & d$group %in% c("DD", "WD"), ]
d$gene <- factor(d$gene, levels = genes8)
d$group <- factor(d$group, levels = c("WD", "DD"))

p <- ggplot(d, aes(x = group, y = expr, fill = group)) +
  geom_violin(alpha = 0.6, trim = TRUE, color = "grey20") +
  geom_boxplot(width = 0.12, outlier.shape = NA, fill = "white", alpha = 0.4) +
  facet_wrap(~ gene, scales = "free_y", nrow = 4) +
  scale_fill_manual(values = c(WD = "#2980B9", DD = "#E67E22"), name = "") +
  labs(x = "", y = "log2CPM", title = "Component-gene expression: DD vs WD") +
  theme_publication() +
  theme(strip.text = element_text(size = 6.5, face = "bold"))

export_fig(p, file.path(OUT, "figS2a.pdf"), w = 8, h = 5)
cat("figS2a done\n")
