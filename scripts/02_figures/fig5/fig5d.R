# ============================================================
# Fig 5D · 评分 A 在各 L2 细胞类中的位置
#   ★ 2026-09-22 合规改造（PI 裁定，见 S25 §9）：
#     ① y 轴按「同一轴不混类型」（PI 规范）**分子轴** —— 与 5C 同三轴、同类序
#        Adipocytic(2) / Immune & stromal(7) / Tumour programs(4)
#     ② `Adipo diff` 归 **Adipocytic** 侧（PI 方案 X）→ 与桶色标（Adipocytic 金）
#        自洽；这是 S25 §8.6 指出的「轴分组 vs 色标」矛盾的修复
#     ③ x 轴三子轴**共享同一评分刻度**（facet 默认 scales=fixed 对 x 生效）
#        → 跨子轴可横向对读评分水平（若放开 x 会使轴间不可比）
#   分组定义唯一真相源：fig5_common.R 的 AXIS_OF_L2 / AXIS_LEV
#   图型：IQR 横条 + 均值点（点估计），非传统四分位箱线
#     ① 与 Step 7 分解口径一致：Kitagawa 的 m 用「类内**均值**」→ 本图也画均值
#        （若画中位数，则本图与分解用的量不是同一个，读者会误对）
#     ② 单细胞评分重尾（Normal adipocytes 的 5–95% 分位跨 −0.09…0.96），
#        画全须会把其余 12 类压扁 → 只画 IQR（剔除极尾）
#     ⚠️ 相对冻结版「箱线/密度」的偏离：中位数→均值、须→仅 IQR
#        已在 LAYOUT_SPEC_v3 §8.5 与方案书 §1.13 登记
#   ⚠️ 面板数 / 字母数 / 行数 / 尺寸(974×516 px) / 字号 全部不变
# ============================================================
suppressWarnings(suppressMessages({library(ggplot2)}))
source(file.path(Sys.getenv("LPS_FIG_ROOT", "."), "utils", "colors.R"))
source(file.path(Sys.getenv("LPS_FIG_ROOT", "."), "utils", "themes.R"))
source(file.path(Sys.getenv("LPS_FIG_ROOT", "."), "utils", "fig5_common.R"))

sel <- fig5_l2_classes()
tab <- sel$tab[order(-sel$tab$mean_score), ]
cls <- tab$L2
stopifnot(all(cls %in% names(AXIS_OF_L2)))

cs <- read.csv(file.path(W5, "s6", "cell_scores.csv.gz"), stringsAsFactors = FALSE)
stopifnot(nrow(cs) == 101890)
d  <- cs[cs$L2 %in% cls, ]
stopifnot(nrow(d) == sum(tab$n))

qs <- do.call(rbind, lapply(cls, function(k) {
  x <- d$score_ams[d$L2 == k]
  data.frame(L2 = k, n = length(x), mean = mean(x), sd = sd(x),
             q25 = quantile(x, .25), med = median(x), q75 = quantile(x, .75))
}))
qs$bucket <- tab$bucket[match(qs$L2, tab$L2)]
qs$cls    <- factor(unname(L2_SHORT[qs$L2]), levels = rev(unname(L2_SHORT[cls])))
qs$axis   <- factor(unname(AXIS_OF_L2[qs$L2]), levels = AXIS_LEV)
qs$lo <- qs$mean - 1.96 * qs$sd / sqrt(qs$n)
qs$hi <- qs$mean + 1.96 * qs$sd / sqrt(qs$n)
print(qs[order(-qs$mean), c("axis", "L2", "n", "mean", "med", "q25", "q75")], row.names = FALSE, digits = 3)
cat(sprintf("均值跨度 [%.3f, %.3f]（轴范围据此定）\n", min(qs$mean), max(qs$mean)))

BLAB <- c(Adipose = "Adipocytic", Tumor = "Tumour", Immune = "Immune", Stroma = "Stroma")
stopifnot(all(qs$bucket %in% names(BUCKET_COL)))
qs$blab <- factor(BLAB[qs$bucket], levels = BLAB)
# 一致性自检：Adipo diff 的桶色标必须与它的子轴同侧（方案 X 的核心断言）
ad <- qs[qs$L2 == "Adipo diff", ]
stopifnot(as.character(ad$axis) == "Adipocytic", as.character(ad$blab) == "Adipocytic")
cat("自检通过：Adipo diff 的子轴 = Adipocytic，桶色 = Adipocytic（方案 X 自洽）\n")

p <- ggplot(qs, aes(y = cls)) +
  geom_vline(xintercept = 0, colour = "grey55", linetype = "dashed", linewidth = 0.35) +
  geom_errorbarh(aes(xmin = q25, xmax = q75, colour = blab), height = 0.42, linewidth = 2.1) +
  geom_errorbarh(aes(xmin = lo, xmax = hi), height = 0.13, colour = "grey20", linewidth = 0.4) +
  geom_point(aes(x = mean), size = 1.15, colour = "grey10") +
  scale_colour_manual(values = setNames(unname(BUCKET_COL[names(BLAB)]), unname(BLAB)), name = NULL) +
  facet_grid(axis ~ ., scales = "free_y", space = "free_y", switch = "y") +
  labs(x = "Score A (bar = IQR; dot = mean, 95% CI)", y = NULL,
       title = "Score across L2 classes") +
  theme_publication() +
  theme(axis.line.y = element_blank(), axis.ticks.y = element_blank(),
        axis.text.y = element_text(size = 6.5),
        axis.text.x = element_text(size = 6.5),
        axis.title.x = element_text(size = 6.2),
        strip.placement = "outside",
        strip.background = element_rect(fill = "grey94", colour = "grey45", linewidth = 0.3),
        strip.text.y = element_text(size = 6.2, angle = 0, margin = margin(0.6, 1.6, 0.6, 1.6)),
        panel.spacing.y = unit(2.5, "pt"),
        legend.position = "right", legend.box = "vertical",
        legend.key.size = unit(7, "pt"), legend.text = element_text(size = 6.2),
        legend.spacing.y = unit(0.5, "pt"),
        plot.margin = margin(2, 2, 2, 2))

export_fig(p, file.path(FIG5_OUT, "fig5d.pdf"), w = 3.2, h = 1.7)
cat("fig5d done | 子轴自顶向下:", paste(AXIS_LEV, collapse = " > "), "\n")
