# ============================================================
# Fig S1C · 关键证据：6 基因两两相关矩阵热图（TCGA n=50）
# 数据：TCGA_LPS_real_log2cpm.rds -> expr6（6 基因 log2CPM）
# （原 Fig S1D，2026-09-28 重构统一为 Fig. S1 C）
# ============================================================
suppressWarnings(suppressMessages(library(ggplot2)))
suppressWarnings(suppressMessages(library(reshape2)))
source(file.path(Sys.getenv("LPS_FIG_ROOT", "."), "utils", "colors.R"))
source(file.path(Sys.getenv("LPS_FIG_ROOT", "."), "utils", "themes.R"))

OUT <- file.path(Sys.getenv("LPS_FIG_ROOT", "."), "out")

load(file.path(Sys.getenv("LPS_FIG_ROOT", "."), "p5_data.RData"))

cmat <- cor(t(expr6), use = "pairwise.complete.obs")
# 上三角掩码
cmat[upper.tri(cmat, diag = TRUE)] <- NA
m <- melt(cmat, na.rm = TRUE)
names(m) <- c("G1", "G2", "r")

p <- ggplot(m, aes(x = G1, y = G2, fill = r)) +
  geom_tile(color = "white", size = 0.5) +
  geom_text(aes(label = sprintf("%.2f", r)), size = 3, color = "grey20") +
  scale_fill_gradient2(low = "#2980B9", mid = "white", high = "#E74C3C",
                       midpoint = 0, limits = c(-1, 1), name = "r") +
  coord_equal() +
  labs(x = "", y = "", title = "Component-gene inter-correlation (TCGA)") +
  theme_publication() +
  theme(axis.text.x = element_text(angle = 45, hjust = 1))

export_fig(p, file.path(OUT, "figS1c.pdf"), w = 5.5, h = 4)
cat("figS1c done\n")
