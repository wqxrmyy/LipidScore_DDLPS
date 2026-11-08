# ============================================================
# Fig 2A · CN 簇 vs 评分 A（箱线 + 散点叠加 + 统计标注）
# 支柱 2 主图：主效应 r=0.5354, P=6.17e-05
# ============================================================
suppressWarnings(suppressMessages(library(ggplot2)))
source(file.path(Sys.getenv("LPS_FIG_ROOT", "."), "utils", "colors.R"))
source(file.path(Sys.getenv("LPS_FIG_ROOT", "."), "utils", "themes.R"))

load(file.path(Sys.getenv("LPS_FIG_ROOT", "."), "p5_data.RData"))
OUT <- file.path(Sys.getenv("LPS_FIG_ROOT", "."), "out")

# 构建数据
d <- data.frame(
  cluster = factor(clin_real[score_A$Sample, "paper_CN.cluster"], levels = c("2", "3", "6")),
  Score = score_A$Score
)
d <- d[!is.na(d$cluster), ]

# Kruskal-Wallis P（与 d12 一致 H=12.02 P=0.00245）
kw <- kruskal.test(Score ~ cluster, data = d)

p <- ggplot(d, aes(x = cluster, y = Score, fill = cluster)) +
  geom_boxplot(outlier.shape = NA, alpha = 0.55, width = 0.55) +
  geom_jitter(width = 0.12, size = 1.1, alpha = 0.75, color = "grey20") +
  scale_fill_manual(values = c("2" = "#2980B9", "3" = "#7F8C8D", "6" = "#E74C3C")) +
  labs(x = "Copy-number cluster", y = "Lipid metabolism score (A)",
       title = "Score by CN cluster",
       # v3：半幅 3.25 in / 7 pt 副标题安全上限约 54 字符，原句 64 字符被设备边缘裁断
       subtitle = sprintf("KW H = %.2f, P = %.3g; r = 0.535, P = 6.2e-05", kw$statistic, kw$p.value)) +
  theme_publication() +
  ## ★★ S39（2026-09-24）：右余量原仅 14 px（实测最右墨迹 = x 轴标题
  ##   "Copy-number cluster" 右端 x=959）；下距 4 -> 6 pt 给轴标签呼吸。
  ##   theme_publication 默认 margin(4,4,4,4) 在半幅 974 px 下左右各仅 17 px。
  theme(legend.position = "none", plot.margin = margin(4, 9, 6, 4))

export_fig(p, file.path(OUT, "fig2a.pdf"), w = 4.5, h = 4)
cat("fig2a done\n")
