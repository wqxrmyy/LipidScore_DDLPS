# ============================================================
# Fig 5C · 点图：6 个分量基因 × L2 细胞类
#   ★ 2026-09-22 合规改造（PI 裁定，见 S25 §9）：
#     ① y 轴按「同一轴不混类型」（PI 规范）**分子轴** —— 3 个子轴，按表型口径重组：
#        Adipocytic(2) / Immune & stromal(7) / Tumour programs(4)
#     ② `Adipo diff` 归 **Adipocytic** 侧（PI 方案 X，表型口径）→ 与面板 D 的桶色
#        （Adipocytic 金）自洽，消除 S25 §8.6 的「轴分组 vs 色标」同图矛盾
#     ③ x 轴 6 个基因名由 45° 改**水平**（0°）：腾出纵向空间，避免分轴后标签粘连
#        （实测 45° 时出现 3 处两行粘连，墨高 45–50 px；改 0° 后 13 行全 20–25 px）
#   分组定义唯一真相源：fig5_common.R 的 AXIS_OF_L2 / AXIS_LEV（C/D 共用，防漂移）
#   粒度：L2（Step 5 冻结的主粒度；与 Step 7 分解同粒度）
#   编码：点大小 = 表达细胞比例%；点色 = 平均表达（按基因跨类 z 标准化，
#         即 Seurat DotPlot 默认 scale=TRUE 口径 —— 6 基因绝对量级相差 100 倍，
#         不标准化则色阶被 SCD 支配、其余 5 基因不可辨）
#   数据来源：s39_s9c_L2genes.R（Seurat 归一化 data 层，expm1 均值）
#   ★ 2026-09-23 标签排版 **C 版**（PI 裁定；依据 _work/s29_evidence_pack/label_policy_proposal.md）：
#     ④ x 轴 6 个基因标签墨宽曾与位距等宽（位距 74.0 px < 墨宽 75 px）→ 视觉串粘。
#        处置只用 **L0 结构回收**（图例间距/边距/色条宽）+ **L1 措辞压缩**（图例标题），
#        不旋转、不缩字号、不动高度、不动画布：
#          `mean expr.(z per gene)` -> `z-score(per gene)`；`% expressing` -> `% expr.`
#          色条宽 6 -> 5.5 pt；legend.spacing.x/.margin/.box.margin/.box.spacing 收紧
#        实测位距 74.0 -> **90.0 px**；但**逐对净间隙**复测仅 **4/5 合格**
#        （SREBF1|SCD 15 / SCD|HMGCR 14 / **HMGCR|ACLY 5** / ACLY|PPARG 13 / PPARG|FASN 12 px，
#         判据 d_min = 0.35 em ≈ 9.5 px）——旧记「6 个标签间隙 12–15 px」系 12 px
#        切词阈值造成的**误判**，已作废（详见 _work/s30_evidence_pack/fig5c_label_closeout.md §二）
#        图注为完整句表述（dot size = percentage expressing / colour = mean expression
#        z-standardised per gene），**不受图例标题缩写影响 → 图注字数零变化**
#   ★ 2026-09-23 标签收口 **D3 版**（PI 裁定 = C + 色标标题收为一行；**仅动本面板**）：
#     ⑤ `z-score\n(per gene)` -> `z-score`（标题由两行收为**一行**）→ 右块宽度
#        不再被两行标题绑住，回收宽度全部让给绘图区：位距 90.0 -> **96.0 px**，
#        逐对间隙 20 / 19 / **10** / 17 / 18 px → **5/5 合格**。
#        定义不丢失：图注 (C) 已完整给出「dot colour = per-gene z-standardised
#        mean expression」→ **图注字数零变化**（505 词不动）。
#     ⑥ 未采 D1（类名折行，100.0 px / 间隙 15 px）：须同步重出面板 D 且压紧 y 行距；
#        未采 D2（x 标签两行错位）：标签样式将与 Fig1–4 不一致。
#     ⑦ 已知代价：最小间隙 10.0 px vs 判据 9.5 px，**余量仅 0.5 px**（登记于 s30 §四）。
#   ⚠️ 面板数 / 字母数 / 行数 / 尺寸(974×516 px) / 字号 全部不变
# ============================================================
suppressWarnings(suppressMessages({library(ggplot2)}))
source(file.path(Sys.getenv("LPS_FIG_ROOT", "."), "utils", "colors.R"))
source(file.path(Sys.getenv("LPS_FIG_ROOT", "."), "utils", "themes.R"))
source(file.path(Sys.getenv("LPS_FIG_ROOT", "."), "utils", "fig5_common.R"))

GENES <- c("SREBF1", "SCD", "HMGCR", "ACLY", "PPARG", "FASN")
sel   <- fig5_l2_classes()
tab   <- sel$tab[order(-sel$tab$mean_score), ]
cls   <- tab$L2
stopifnot(all(cls %in% names(AXIS_OF_L2)))
cat(sprintf("纳入 %d 类，覆盖 %.1f%% 细胞；未纳入 %d 类: %s\n",
            length(cls), sel$coverage, length(sel$omitted), paste(sel$omitted, collapse = ", ")))
cat("  分轴（按类数）: ",
    paste(sapply(AXIS_LEV, function(a) sprintf("%s=%d", a, sum(AXIS_OF_L2[cls] == a))), collapse = " | "), "\n")

d <- read.csv(file.path(FIG5_DATA, "L2_genes.csv"), stringsAsFactors = FALSE)
d <- d[d$gene %in% GENES & d$L2 %in% cls, ]
stopifnot(nrow(d) == length(GENES) * length(cls))
d$gene <- factor(d$gene, levels = GENES)
lev    <- rev(unname(L2_SHORT[cls]))              # rev：ggplot y 轴自下而上
d$cls   <- factor(unname(L2_SHORT[d$L2]), levels = lev)
d$axis  <- factor(unname(AXIS_OF_L2[d$L2]), levels = AXIS_LEV)
d$z     <- ave(d$mean_exp, d$gene, FUN = function(x) (x - mean(x)) / sd(x))

pct_at <- c(10, 30, 60)

p <- ggplot(d, aes(x = gene, y = cls)) +
  geom_point(aes(size = pct_exp, fill = z), shape = 21, colour = "grey25", stroke = 0.25) +
  scale_size_continuous(name = "% expr.", range = c(0.5, 3.6), breaks = pct_at,
                        limits = c(0, max(d$pct_exp))) +
  scale_fill_gradientn(name = "z-score",
                       colours = c("#2166AC", "#F7F7F7", "#B2182B"),
                       limits = c(-2.6, 2.6), oob = scales::squish,
                       guide = guide_colourbar(barwidth = unit(5.5, "pt"),
                                               barheight = unit(52, "pt"),
                                               ticks.colour = "black",
                                               frame.colour = "black")) +
  facet_grid(axis ~ ., scales = "free_y", space = "free_y", switch = "y") +
  labs(x = NULL, y = NULL, title = "Marker genes (L2 classes)") +
  theme_publication() +
  theme(panel.grid.major = element_line(colour = "grey93", linewidth = 0.25),
        axis.line = element_blank(),
        axis.ticks = element_blank(),
        axis.text.x = element_text(size = 6.5, angle = 0, hjust = 0.5, vjust = 1),
        axis.text.y = element_text(size = 6.5),
        strip.placement = "outside",
        strip.background = element_rect(fill = "grey94", colour = "grey45", linewidth = 0.3),
        strip.text.y = element_text(size = 6.2, angle = 0, margin = margin(0.6, 1.6, 0.6, 1.6)),
        panel.spacing.y = unit(1.5, "pt"),
        legend.position = "right", legend.box = "vertical",
        legend.key.size = unit(7, "pt"),
        legend.text = element_text(size = 6.2),
        legend.title = element_text(size = 6.2),
        legend.spacing.y = unit(1, "pt"),
        legend.spacing.x = unit(1, "pt"),
        legend.margin = margin(0, 0, 0, 0),
        legend.box.margin = margin(0, 0, 0, 0),
        legend.box.spacing = unit(1, "pt"),
        plot.margin = margin(2, 2, 2, 2))

export_fig(p, file.path(FIG5_OUT, "fig5c.pdf"), w = 3.2, h = 1.7)
cat("fig5c done | 类顺序（自下而上）:", paste(lev, collapse = " < "), "\n")
cat("  子轴自顶向下:", paste(AXIS_LEV, collapse = " > "), "\n")
