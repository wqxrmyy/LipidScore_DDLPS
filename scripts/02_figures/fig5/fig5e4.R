# ============================================================
# Fig 5E4 · 组成层关键证据：脂肪系占比的跨组织梯度
#   组序（左->右按分化程度）：DDLPS-DD < DDLPS-WD / DDLPS-sWD / WDLPS / lipoma / nl
#   参照比较组 = 脂肪组织（nl 瘤旁脂肪、lipoma），用于判断 DD 的塌缩幅度
#   ⚠️ nl 为瘤旁脂肪（peritumoral fat）：与肿瘤样本同患者、同后腹膜位点（见 P1 取证 / LIM-P3-016）
#      → 只作比较组，**不作正常上限**（距肿瘤远近未知，可能已受肿瘤微环境影响）
#   编码：空心点 = 各样本；实心点 = DDLPS-DD（与 E1–E3 的「实心 = DD」一致）
#         短竖线 = 组中位数
#   ⚠️ 未纳入本面板的组织：DDLPS-DDK（单样本、注释含义不明）—— 见 prep 报告，
#      属显式披露而非静默剔除
# ============================================================
suppressWarnings(suppressMessages({library(ggplot2)}))
source(file.path(Sys.getenv("LPS_FIG_ROOT", "."), "utils", "colors.R"))
source(file.path(Sys.getenv("LPS_FIG_ROOT", "."), "utils", "themes.R"))
source(file.path(Sys.getenv("LPS_FIG_ROOT", "."), "utils", "fig5_common.R"))

d <- fig5_read("adipose_gradient.csv")
GRP <- c("DDLPS-DD", "DDLPS-WD", "DDLPS-sWD", "WDLPS", "lipoma", "nl")
d$grp <- factor(d$grp, levels = GRP)
d$gy  <- as.numeric(d$grp)          # 全部用数值 y（混用离散/连续会让 ggplot 训练失败）
d <- d[order(d$gy), ]
med <- aggregate(cbind(frac, gy) ~ grp, data = d, FUN = median)
med$isDD <- med$grp == "DDLPS-DD"
d$isDD   <- d$grp == "DDLPS-DD"
dd_med <- med$frac[med$isDD]; wd_med <- med$frac[med$grp == "DDLPS-WD"]
fold   <- wd_med / dd_med
ann    <- sprintf("median DD %.3f\nvs WD %.3f (%.1fx)", dd_med, wd_med, fold)
# 防溢出：1.62 in 面板的绘图区约 1.1 in，5.8 pt 下每行 ≈ 24 字符封顶
stopifnot(all(nchar(strsplit(ann, "\n")[[1]]) <= 22))

p <- ggplot(d, aes(y = gy, x = frac)) +
  geom_segment(data = med, aes(x = frac, xend = frac, y = gy - 0.30, yend = gy + 0.30,
                               colour = isDD),
               linewidth = 0.9, show.legend = FALSE) +
  geom_point(aes(fill = isDD), shape = 21, size = 1.7, stroke = 0.45,
             colour = KEY_COL, position = position_jitter(height = 0.11, seed = 1),
             show.legend = FALSE) +
  scale_fill_manual(values = c(`FALSE` = "white", `TRUE` = KEY_COL)) +
  scale_colour_manual(values = c(`FALSE` = "grey60", `TRUE` = KEY_COL)) +
  scale_x_continuous(limits = c(0, 0.52), breaks = c(0, 0.25, 0.5)) +
  scale_y_continuous(breaks = seq_along(GRP), labels = GRP, limits = c(0.45, 7.95)) +
  labs(x = "Adipocytic lineage fraction", y = NULL,
       title = "Fraction by tissue") +
  annotate("text", x = 0.012, y = 7.9, label = ann, hjust = 0, vjust = 1,
           size = 2.16, colour = "grey20", lineheight = 0.95) +
  theme_publication() +
  theme(axis.line.y = element_blank(), axis.ticks.y = element_blank(),
        axis.text.y = element_text(size = 6.2),
        axis.text.x = element_text(size = 6.2),
        axis.title.x = element_text(size = 6.2),
        plot.margin = margin(2, 2, 2, 2))

export_fig(p, file.path(FIG5_OUT, "fig5e4.pdf"), w = 1.62, h = 1.61)
cat(sprintf("fig5e4 done | DD 中位 %.4f vs WD 中位 %.4f | %.2fx\n", dd_med, wd_med, fold))
