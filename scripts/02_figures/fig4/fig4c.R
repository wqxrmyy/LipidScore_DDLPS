# ============================================================
# Fig 4C · 分量基因 DD vs WD（配对 Δ，z 分量纲，棒棒糖）
# 身份基因在 DD 内显著丢失；MDM2/CDK4 上调为对照
# ============================================================
suppressWarnings(suppressMessages(library(ggplot2)))
source(file.path(Sys.getenv("LPS_FIG_ROOT", "."), "utils", "colors.R"))
source(file.path(Sys.getenv("LPS_FIG_ROOT", "."), "utils", "themes.R"))

load(file.path(Sys.getenv("LPS_FIG_ROOT", "."), "p5_gse221492.RData"))
OUT <- file.path(Sys.getenv("LPS_FIG_ROOT", "."), "out")

d <- gse_meta[gse_meta$Group %in% c("DD", "WD"), ]
tab <- table(d$Patient, d$Group)
paired_pat <- rownames(tab)[tab[, "DD"] == 1 & tab[, "WD"] == 1]

# 每基因跨全 53 样本 z 分（统一量纲）
zall <- t(scale(t(as.matrix(gse_genes))))

genes <- rownames(zall)
rows <- NULL
for (g in genes) {
  ex <- as.numeric(zall[g, ])
  names(ex) <- gse_meta$Sample
  dd <- ex[gse_meta$Sample[gse_meta$Group == "DD" & gse_meta$Patient %in% paired_pat]]
  wd <- ex[gse_meta$Sample[gse_meta$Group == "WD" & gse_meta$Patient %in% paired_pat]]
  tt <- suppressWarnings(t.test(dd, wd, paired = TRUE))
  rows <- rbind(rows, data.frame(Gene = g, delta = mean(dd - wd), p = tt$p.value))
}
rows$Gene <- factor(rows$Gene, levels = rows$Gene[order(rows$delta)])
rows$sig <- sapply(rows$p, sig_label)
rows$col <- ifelse(rows$p < 0.05 & rows$delta < 0, "#E67E22",
                   ifelse(rows$p < 0.05 & rows$delta > 0, "#C0392B", "#7F8C8D"))

p <- ggplot(rows, aes(x = delta, y = Gene)) +
  geom_vline(xintercept = 0, color = "grey40", linetype = "dotted") +
  geom_segment(aes(x = 0, xend = delta, y = Gene, yend = Gene, color = col), size = 1.1) +
  geom_point(aes(color = col), size = 3.5) +
  geom_text(aes(label = sig), hjust = -0.35, size = 3, color = "grey20") +
  # v4（2026-09-25，S65 续）：图例文字折成两行 —— 用「不改一个字」的方式修 x 轴标题被裁。
  #   病灶：成品 PNG 左缘 3 px 内有 16 个墨点（触边），标题左端 "P" 被 cairo 设备裁掉，
  #     成品实际显示为 "aired mean difference (DD - WD), gene z-score"。
  #   根因（gtable 实测列宽）：guide-box-right 图例占 434 px（45% 画布宽，
  #     因为 "Gained in DD (P < 0.05)" 一行就有 ~350 px），绘图区被压到仅 309 px；
  #     而 x 轴标题约 750 px 居中于 309 px -> 两侧各溢出 220 px，
  #     左端 191-220 = -29 px 落到画布外 -> 被裁。
  #   ★ 为什么不把图例移到底部：底部图例占 74.5 px 高，会把 14 个基因的行距
  #     从 31.5 px 压到 25 px < 7 pt 文字的 em 行高 29 px -> 会在标签处制造压字。
  #   ★ 为什么不用 hjust=0 左对齐：会与同图其余面板的居中轴标题风格不一致。
  #   ★ 为什么不用「x 轴标题折行」：那一版实测把绘图区从 440 压到 404 px
  #     （净空 11 -> 6.7 px），是用退化换修复。改图例折行则不占额外宽度、
  #     不占额外高度（图例高 125 px << 绘图区 440 px），是纯收益：
  #     图例宽 -> ~250 px，绘图区 309 -> ~493 px，宽高比 0.70 -> ~1.12:1，
  #     标题 750 px 居中于 493 px -> 溢出收窄到 ~128 px，左端 63 px > 0 ✓。
  #   ★ 文字内容一字未改，只在原标签中间加换行符。
  scale_color_identity(name = NULL,
                       breaks = c("#E67E22", "#C0392B", "#7F8C8D"),
                       labels = c("Lost in DD\n(P < 0.05)", "Gained in DD\n(P < 0.05)", "NS"),
                       guide = "legend") +
  labs(x = "Paired mean difference (DD - WD), gene z-score", y = "",
       title = "Identity genes lost in DD component",
       # v3：原句 55 字符，半幅会被裁断
       subtitle = "Within-patient paired comparison (8 patients)") +
  theme_publication() +
  # v3：14 个基因类别 + 底部图例会把绘图区压到每行 6.9 pt（< 7 pt 轴标签 → 压字）。
  #     图例移到右侧（本图 x 轴 0 右侧留白充足），把 ~45 px 高度让给绘图区；
  #     行高同时按内容密度重排（PANEL_BOX fig4c = .half(0.72) → 700 px，绘图区 ≈9.2 pt/行）。
  theme(legend.position = "right", legend.text = element_text(size = 7),
        legend.key.height = unit(0.32, "cm"),
        axis.text.y = element_text(size = 7))

export_fig(p, file.path(OUT, "fig4c.pdf"), w = 6, h = 4.8)
cat("fig4c done\n")
