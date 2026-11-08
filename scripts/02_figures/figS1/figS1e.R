# ============================================================
# Fig S1E · 评分 vs 各分量基因散点矩阵（6 基因并列，2×3）
# 来源：原 Fig 1F（正文从未引用；1G 信息已在 Figure 3F 与 Table S2）
# 迁移：原 Fig S5A，2026-09-28 重构统一为 Fig. S1 E
# ============================================================
suppressWarnings(suppressMessages(library(ggplot2)))
suppressWarnings(suppressMessages(library(reshape2)))
source(file.path(Sys.getenv("LPS_FIG_ROOT", "."), "utils", "colors.R"))
source(file.path(Sys.getenv("LPS_FIG_ROOT", "."), "utils", "themes.R"))

load(file.path(Sys.getenv("LPS_FIG_ROOT", "."), "p5_data.RData"))
OUT <- file.path(Sys.getenv("LPS_FIG_ROOT", "."), "out")

# 构建长表：每个样本 × 每个基因 的表达 + 评分
genes <- rownames(expr6)
samples <- colnames(expr6)
df <- data.frame()
for (g in genes) {
  df <- rbind(df, data.frame(
    Gene = g,
    Score = score_A[samples, "Score"],
    Expr = as.numeric(expr6[g, samples]),
    r = cor(score_A[samples, "Score"], as.numeric(expr6[g, samples])),
    p = cor.test(score_A[samples, "Score"], as.numeric(expr6[g, samples]))$p.value
  ))
}
df$lab <- sprintf("%s\nr=%.2f, %s", df$Gene, df$r, sapply(df$p, sig_label))

p <- ggplot(df, aes(x = Score, y = Expr)) +
  geom_point(size = 0.8, alpha = 0.6, color = "#8E44AD") +
  geom_smooth(method = "lm", se = TRUE, color = "#E74C3C", size = 0.7) +
  facet_wrap(~ lab, scales = "free_y", ncol = 3) +
  labs(x = "Lipid metabolism score (A)", y = "log2(CPM)") +
  theme_publication()

export_fig(p, file.path(OUT, "figS1e.pdf"), w = 6.5, h = 3.9)
cat("figS1e done\n")
