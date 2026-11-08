# ============================================================
# Fig 4G · 三队列方向统一示意（WD 高 / DD 低：S8+S13+S17 三线统一）
# 条形：各队列/对比中的标准化方向指标
# ============================================================
suppressWarnings(suppressMessages(library(ggplot2)))
source(file.path(Sys.getenv("LPS_FIG_ROOT", "."), "utils", "colors.R"))
source(file.path(Sys.getenv("LPS_FIG_ROOT", "."), "utils", "themes.R"))

OUT <- file.path(Sys.getenv("LPS_FIG_ROOT", "."), "out")

# 三条独立线索的方向证据（v2.1 narrative 第一节 + d11）
df <- data.frame(
  # v3：整行面板高 516 px（1.72 in），绘图区 ≈1.27 in → 每行 23 pt。
  #     原三行标签（≈24 pt）会互相压字 → 压到两行。
  Evidence = c("GSE221492\nDD vs WD (delta -0.47)",
               "GSE30929\nhigh group (46/52 WD)",
               "Identity genes\nCD36/PPARG/FABP4",
               "TCGA CN cluster\nscore up (aggressive)"),
  Value = c(-0.469, -1, -1, 1),  # 方向：负=DD/去分化侧丢失；正=TCGA 侵袭侧升高
  Type = c("Paired", "Cross", "Genes", "Cross"),
  Dir = c("Lost in DD", "High = WD-like", "Lost in DD", "High = aggressive")
)
df$Evidence <- factor(df$Evidence, levels = df$Evidence)
df$hjust <- ifelse(df$Value < 0, 1.1, -0.1)

p <- ggplot(df, aes(x = Value, y = Evidence, fill = Dir)) +
  geom_vline(xintercept = 0, color = "grey30") +
  geom_bar(stat = "identity", width = 0.55, alpha = 0.85) +
  scale_fill_manual(values = c("Lost in DD" = "#E67E22", "High = WD-like" = "#2980B9",
                               "High = aggressive" = "#E74C3C"), name = "") +
  labs(x = "Direction of association (standardized)", y = "",
       # v3：标题原 53 字符 / 副标题原 67 字符，在整行 1949 px 下虽可容纳，
       #     但为与其余面板保持同一字号口径，按 8 pt 标题（≈3.19 in 上限）压缩
       title = "Three lines converge: identity loss in DD",
       subtitle = "WD-high / DD-low across cohorts") +
  theme_publication() +
  # v4（2026-09-25，S65 续）：三处调整，均不改内容、不改版式。
  #   ① axis.text.y 7 -> 6.4 pt、lineheight 0.95 -> 0.9：
  #      病灶：绘图区仅 246 px / 4 个 y 类别 -> 类别间距 246/4.2 = 58.6 px，
  #      而折行标签两行的完整占高 = em 29 px + 行距 27.7 px = 57 px -> 净空仅 1.6 px，
  #      放大裁图实测四对标签全部紧贴（TCGA CN cluster/score up (aggressive) 等）。
  #      6.4 pt + 0.9 后两行占 26.7 + 23.0 = 49.7 px -> 净空涨到 ~9 px。
  #      ★ 6.4 pt 仍高于期刊下限 6 pt。
  #   ② plot.margin 上下 4 -> 2 pt：把省下的高度让给绘图区（改善 4.77:1 的极端扁比例）。
  #   ③ plot.margin 右 4 -> 12 pt：用户报告「拼板右侧留的空间不够」——实测图例文字
  #      右端距画布仅 34 px（0.11 in），加大右边距给图例呼吸空间。
  theme(axis.text.y = element_text(size = 6.4, lineheight = 0.9),
        plot.margin = margin(2, 12, 2, 2),
        plot.title = element_text(size = 8, face = "bold", margin = margin(b = 1)),
        plot.subtitle = element_text(size = 7, lineheight = 1.0, margin = margin(t = 0)))

export_fig(p, file.path(OUT, "fig4g.pdf"), w = 6.8, h = 4)
cat("fig4g done\n")
