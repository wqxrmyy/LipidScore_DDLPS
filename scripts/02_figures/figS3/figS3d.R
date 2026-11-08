# ============================================================
# Fig S3D · 通讯层：TAM<->肿瘤的抑制性通讯图谱（CellChat）
#   （原 Fig 6C，2026-09-28 重构统一为 Fig. S3 D）
#   上栏：单一棒棒糖图 —— 双向各 top 8（共 16 行）。用连续 y 轴 + 手动刻度，
#         因为 CD99_CD99 在两个方向都出现，factor 无法区分同名行。
#         水平虚线分隔两个方向；SPP1 家族标红、THBS2-CD36 标青。
#   ★ 本面板不检验脂质轴：CellChatDB 不含 CD36/PPARG/FABP4 的脂质介导互作
#   数据：figS6_data/figS6c_lr_pairs.csv + figS6c_hub.csv
# ============================================================
suppressWarnings(suppressMessages({library(ggplot2); library(patchwork); library(scales)}))
source(file.path(file.path(Sys.getenv("LPS_FIG_ROOT", "."), "utils"), "themes.R"))
source(file.path(file.path(Sys.getenv("LPS_FIG_ROOT", "."), "utils"), "figS6_common.R"))

lr <- read6("figS6c_lr_pairs.csv")
stopifnot(nrow(lr) == 50)
TOP <- 8
sel <- do.call(rbind, lapply(c("TAM3->tumour", "tumour->TAM3"), function(d) {
  z <- lr[lr$dir == d, ]
  head(z[order(-z$weight), ], TOP)
}))
sel$hl <- "other"
sel$hl[sel$flag == "exclusive"] <- "SPP1 family"
sel$hl[sel$pair == "THBS2_CD36"] <- "THBS2-CD36"
sel$hl <- factor(sel$hl, levels = c("SPP1 family", "THBS2-CD36", "other"))
HL_SHOWN <- intersect(c("SPP1 family", "THBS2-CD36", "other"), unique(as.character(sel$hl)))
stopifnot(length(HL_SHOWN) >= 1)
cat("   S3D 图例项: ", paste(HL_SHOWN, collapse = ", "), "\n", sep = "")
HL_COL <- c("SPP1 family" = "#C0392B", "THBS2-CD36" = "#1B6B7A",
            "other" = "#B3B3B3")
HL_LAB <- c("SPP1 family" = "SPP1 fam.", "THBS2-CD36" = "THBS2-CD36",
            "other" = "other pairs")

sel <- sel[order(match(sel$dir, c("TAM3->tumour", "tumour->TAM3")), -sel$weight), ]
n1 <- sum(sel$dir == "TAM3->tumour")
## ★★ 2026-09-29（第 3 轮，走廊加宽）：原布局把「小标题 + 虚线」同时塞进
##   组间 2.0 单位的走廊，实测（300 dpi，行距 26.6 px/单位）虚线 y=357.8 落在
##   小标题墨迹 y=338.5-358.5 之内 —— 二者贴死，即 PI 说的「堆叠在一起」。
##   修正：走廊 2.0 -> 2.6 单位（组1 上移 0.6，Y_TOP 同步 18.8 -> 19.5），
##   并把顺序改成 行 -> 虚线 -> 小标题 -> 行（小标题挂到自己那组的头上）。
##   实测余量（新行距 ≈25.7 px/单位）：小标题↔下行 8.2px、↔虚线 9.2px、虚线↔上行 9.2px。
sel$ypos <- c(seq(17.6, 17.6 - n1 + 1), seq(8, 8 - TOP + 1))
GAP <- 9.85
Y_TOP <- 19.5; Y_BOT <- 0.55
XMAX <- 0.80
cat(sprintf("S3D up: rows=%d (group1=%d) | y %.1f..%.1f\n",
            nrow(sel), n1, min(sel$ypos), max(sel$ypos)))

pC1 <- ggplot(sel, aes(weight, ypos)) +
  geom_hline(yintercept = GAP, linetype = 2, colour = "grey30", linewidth = 0.55) +
  geom_segment(aes(x = 0, xend = weight, y = ypos, yend = ypos),
               colour = "grey65", linewidth = 0.4) +
  geom_point(aes(fill = hl), shape = 21, size = 2.2, stroke = 0.5, colour = "grey20") +
  geom_text(aes(label = sprintf("%.3f", weight)), hjust = -0.28, size = 2.0,
            colour = "grey20") +
  annotate("text", x = 0.005, y = 18.95, label = "TAM3 -> tumour", hjust = 0,
           size = 2.35, fontface = "bold", colour = "grey15") +
  annotate("text", x = 0.005, y = 9.1, label = "tumour -> TAM3", hjust = 0,
           size = 2.35, fontface = "bold", colour = "grey15") +
  scale_fill_manual(values = HL_COL, name = NULL, drop = FALSE,
                    breaks = HL_SHOWN,
                    labels = HL_LAB[HL_SHOWN]) +
  scale_x_continuous(limits = c(0, XMAX), breaks = c(0, 0.2, 0.4, 0.6),
                     expand = c(0, 0)) +
  scale_y_continuous(breaks = sel$ypos, labels = sel$pair,
                     limits = c(Y_BOT, Y_TOP), expand = c(0, 0)) +
  labs(x = "Communication probability", y = NULL,
       title = "Top 8 LR pairs per direction") +
  theme_publication(base_size = 7.2) +
  theme(axis.text.y = element_text(size = 5.9),
        axis.text.x = element_text(size = 6.0),
        axis.title.x = element_text(size = 6.2),
        panel.grid.major.y = element_blank(),
        plot.title = element_text(size = 6.6),
        plot.title.position = "plot",
        legend.position = "bottom",
        legend.key.height = unit(6, "pt"), legend.key.width = unit(6, "pt"),
        legend.spacing.x = unit(1.5, "pt"),
        legend.text = element_text(size = 5.2),
        legend.margin = margin(0, 0, 0, 0),
        ## 2026-09-29 沿革：底距 2→9pt（PI：两幅子子图堆叠在一起）→ **9→3pt**（PI 复看：
        ##   「D 面板嫌空」）。★ 让间隙的正解不是堆 margin —— margin 是从各自面板里
        ##   扣地，扣得越多两块内容越矮；正确做法是**收 margin + 调 heights 配比**。
        plot.margin = margin(3, 5, 3, 3))

## ---------------------------------------------------------------- 下栏
hub <- read6("figS6c_hub.csv")
stopifnot(nrow(hub) == 2)
hl <- rbind(
  data.frame(subtype = hub$subtype, metric = "outgoing", value = hub$outgoing),
  data.frame(subtype = hub$subtype, metric = "incoming", value = hub$incoming))
hl$metric  <- factor(hl$metric, levels = c("outgoing", "incoming"))
hl$subtype <- factor(hl$subtype, levels = c("TAM1", "TAM3"))
PER <- sprintf("per-cell: %.4f vs %.4f", hub$per_cell_out[1], hub$per_cell_out[2])

pC2 <- ggplot(hl, aes(metric, value, fill = subtype)) +
  geom_col(position = position_dodge(width = 0.68), width = 0.6) +
  geom_text(aes(label = sprintf("%.2f", value)),
            position = position_dodge(width = 0.68), vjust = -0.42,
            size = 2.15, colour = "grey15") +
  scale_fill_manual(values = TAM_COL, name = NULL) +
  scale_y_continuous(limits = c(0, 9.4), expand = c(0, 0)) +
  labs(x = NULL, y = "Total signalling", title = "TAM3 is the hub",
       subtitle = wrap_txt(PER, 30)) +
  theme_publication(base_size = 7.2) +
  theme(axis.text.x = element_text(size = 6.8),
        axis.text.y = element_text(size = 6.2),
        legend.position = "none",
        plot.title = element_text(size = 7.2),
        plot.subtitle = element_text(size = 5.9, colour = "grey30"),
        ## 2026-09-29 沿革：顶距 3→10pt（与上栏合计让出间隙）→ **10→3pt**（PI：嫌空）。
        plot.margin = margin(3, 3, 2, 3))

## heights 配比 2.45:1 → **2.1:1**（2026-09-29，PI：「下面柱状图显得矮了」）。
##   净效果：pC2 槽 325→361 px（+36），再叠加其顶距回收的 7pt≈29px，
##   下栏绘图区合计约 +65px（柱高 ~125→~190px）；pC1 槽 795→759（-36）但底距
##   回收 6pt≈25px，两者基本对冲 -> 上栏行距 25.7→~25.4 px/单位，16 行仍不挤。
## ★ 槽位不能靠 ggsave 探针外推：实测 wrap_elements 在 ggsave 下会溢出
##   （pC1 978 + pC2 284 > 1120），只能以 save6 实渲染后回读为准。
p <- figS6_stack(pC1, pC2, heights = c(2.1, 1))
save6(p, "figS3d")

cat(sprintf("   hub: TAM1 out=%.3f in=%.3f | TAM3 out=%.3f in=%.3f\n",
            hub$outgoing[1], hub$incoming[1], hub$outgoing[2], hub$incoming[2]))
cat(sprintf("   max weight: TAM3->tumour %.4f / tumour->TAM3 %.4f\n",
            max(sel$weight[sel$dir == "TAM3->tumour"]),
            max(sel$weight[sel$dir == "tumour->TAM3"])))
cat(sprintf("   THBS2_CD36 = %.4f\n", lr$weight[lr$pair == "THBS2_CD36"]))
