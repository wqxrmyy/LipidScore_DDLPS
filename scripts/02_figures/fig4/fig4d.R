# ============================================================
# Fig 4D · DD vs 各参照组的判别（AUC 棒棒糖，官方 B5_selfcheck 值）
# 官方：DD_vs_WD 0.2491 | DD_vs_anyWD 0.3129 | DD_vs_lipoma 0.1176 | DD_vs_fat 0.1882
# 数据：fig4_discrimination.csv（prep_fig4_gse_official.py 从 d11 JSON 原值导出）
# ============================================================
suppressWarnings(suppressMessages(library(ggplot2)))
source(file.path(Sys.getenv("LPS_FIG_ROOT", "."), "utils", "colors.R"))
source(file.path(Sys.getenv("LPS_FIG_ROOT", "."), "utils", "themes.R"))

OUT <- file.path(Sys.getenv("LPS_FIG_ROOT", "."), "out")

df <- read.csv(file.path(Sys.getenv("LPS_FIG_ROOT", "."), "fig4_discrimination.csv"),
               stringsAsFactors = FALSE)
# v3：把最长的 y 标签折成两行。原因：ggplot 标题默认对齐「绘图区」左缘，
#     y 轴标签占宽多少，标题就少多少。半幅 3.19 in 下，单行
#     "DD vs any WD (WD+sWD+WDLPS)"（≈1.7 in）会把标题挤到只剩 1.4 in → 被右缘裁断。
#     折行后标签宽降到 ≈0.9 in，标题可用 ≈2.3 in（"DD component discrimination" ≈1.8 in）✓
#     注意：levels 沿用原逻辑（rev(CSV 顺序)），只替换显示文本，不改排序。
lvl <- sub("DD vs any WD \\(WD\\+sWD\\+WDLPS\\)",
           "DD vs any WD\n(WD+sWD+WDLPS)", as.character(df$comparison))
df$comparison <- factor(lvl, levels = rev(lvl))
# 半幅下数值标签宽度受限，压掉过长前缀（x 轴标题已声明 AUC）
df$lab <- sprintf("%.3f, P = %.2g", df$auc, df$p)

p <- ggplot(df, aes(x = auc, y = comparison)) +
  geom_vline(xintercept = 0.5, linetype = "dashed", color = "grey40") +
  geom_segment(aes(x = 0.5, xend = auc, y = comparison, yend = comparison),
               color = "#2980B9", size = 1.1) +
  geom_point(size = 4, color = "#E74C3C") +
  geom_text(aes(label = lab), hjust = -0.06, vjust = -0.8, size = 2.4) +
  coord_cartesian(xlim = c(0, 1)) +
  labs(x = "AUC (low score predicts DD)", y = "", title = "DD component discrimination",
       subtitle = "Lower score = closer to DD") +
  theme_publication() +
  # v4（2026-09-25，S65 续）：只调折行标签的行内间距，不动版面。
  #   加高到 653 px 后类别间距已有 105 px（不再相撞），但折行标签自身两行
  #   仍按默认 lineheight=1.0 -> 两行基线距 29 px、墨高 26 px、净空仅 6 px，
  #   视觉上与类别间 66 px 的大留白反差明显，像"两个不同类别"。
  #   1.15 后基线距 33 px、净空 10 px，与上下留白更均衡。字号仍 7 pt（未动）。
  theme(axis.text.y = element_text(lineheight = 1.15))

export_fig(p, file.path(OUT, "fig4d.pdf"), w = 6.2, h = 3.5)
cat("fig4d done\n")
