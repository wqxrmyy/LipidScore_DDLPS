# ============================================================
# Fig 1A · 研究设计流程图（3 队列 + 三支柱）
# 纯 R 绘制示意图（无外部数据依赖）
# 注意：全英文标签（cairo 设备无法渲染中文，避免乱码）
# ============================================================
suppressWarnings(suppressMessages(library(ggplot2)))
source(file.path(Sys.getenv("LPS_FIG_ROOT", "."), "utils", "colors.R"))
source(file.path(Sys.getenv("LPS_FIG_ROOT", "."), "utils", "themes.R"))

OUT <- file.path(Sys.getenv("LPS_FIG_ROOT", "."), "out")

d <- data.frame(
  x = c(1, 1, 1, 2, 2, 2, 3, 3, 3),
  y = c(3, 2, 1, 3, 2, 1, 3, 2, 1),
  label = c("TCGA-SARC\n(n = 50, DDLPS)",
            "GSE221492\n(n = 53, paired WD/DD)",
            "GSE30929 / 21122 / 159659\n(expression reference)",
            "Lipid metabolism score (A)\nmean z-score of 6 genes",
            "Lipid metabolism score (A)\nmean z-score of 6 genes",
            "Lipid metabolism score (A)\nmean z-score of 6 genes",
            "Pillar 2\nCopy-number cluster\nassociation (Card D)",
            "Pillar 1\nDedifferentiation:\nidentity-gene loss (Card B)",
            "Pillar 3\nCD36 association\n(Card C)"),
  col = c("#E74C3C", "#2980B9", "#7F8C8D",
          "#8E44AD", "#8E44AD", "#8E44AD",
          "#E74C3C", "#2980B9", "#E67E22")
)

p <- ggplot(d, aes(x, y)) +
  geom_rect(aes(xmin = x - 0.42, xmax = x + 0.42, ymin = y - 0.42, ymax = y + 0.42),
            fill = d$col, alpha = 0.88, color = "white", linewidth = 0.4) +
  geom_text(aes(label = label), size = 2.35, lineheight = 0.95,
            color = "white", fontface = "bold") +
  # 箭头：队列 → 评分
  geom_segment(aes(x = 1.44, y = 3, xend = 1.56, yend = 3),
               arrow = arrow(length = unit(0.14, "cm")), linewidth = 0.7) +
  geom_segment(aes(x = 1.44, y = 2, xend = 1.56, yend = 2),
               arrow = arrow(length = unit(0.14, "cm")), linewidth = 0.7) +
  geom_segment(aes(x = 1.44, y = 1, xend = 1.56, yend = 1),
               arrow = arrow(length = unit(0.14, "cm")), linewidth = 0.7) +
  # 箭头：评分 → 支柱
  geom_segment(aes(x = 2.44, y = 3, xend = 2.56, yend = 3),
               arrow = arrow(length = unit(0.14, "cm")), linewidth = 0.7) +
  geom_segment(aes(x = 2.44, y = 2, xend = 2.56, yend = 2),
               arrow = arrow(length = unit(0.14, "cm")), linewidth = 0.7) +
  geom_segment(aes(x = 2.44, y = 1, xend = 2.56, yend = 1),
               arrow = arrow(length = unit(0.14, "cm")), linewidth = 0.7) +
  # 列标题
  annotate("text", x = 1, y = 3.62, label = "Cohorts", fontface = "bold", size = 3, color = "grey25") +
  annotate("text", x = 2, y = 3.62, label = "Score", fontface = "bold", size = 3, color = "grey25") +
  annotate("text", x = 3, y = 3.62, label = "Three pillars", fontface = "bold", size = 3, color = "grey25") +
  xlim(0.4, 3.6) + ylim(0.35, 3.8) +
  theme_void() +
  theme(plot.title = element_text(hjust = 0.5, face = "bold", size = 12))

export_fig(p, file.path(OUT, "fig1a.pdf"), w = 6.2, h = 3.2)
cat("fig1a done\n")
