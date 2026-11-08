# ============================================================
# Fig 3D · Harrell C 对比（Score A vs Alternative score C）
# 图型：水平点估计 + 95% CI（forest 式），非柱状
#   理由① 柱状图截断纵轴会夸大 0.07 的差；不截断则两柱等高、差异不可见
#   理由② 点+CI 天然规避「数值标签与统计线/误差棒遮挡」，标签置 CI 右端外侧
# 数值来源：
#   s22_harrell_A_out.txt (L9/L12) -> s22_harrell_A.csv   官方 Harrell C 与解析 SE
#   s22b_delta_boot.csv            -> ΔC 及其配对 bootstrap 95%CI / 双侧 P
# ============================================================
suppressWarnings(suppressMessages(library(ggplot2)))
source(file.path(Sys.getenv("LPS_FIG_ROOT", "."), "utils", "colors.R"))
source(file.path(Sys.getenv("LPS_FIG_ROOT", "."), "utils", "themes.R"))

OUT <- file.path(Sys.getenv("LPS_FIG_ROOT", "."), "out")
WK  <- file.path(Sys.getenv("LPS_WORK_ROOT", "."), "_work")

# ---- 官方 C 与 SE ----
# 注意：两份 CSV 均为 utf-8-sig（带 BOM），必须显式按 UTF-8-BOM 读，
#       否则首列名会变成 "\ufeffscore_basis" 导致取值为空
hc <- read.csv(file.path(WK, "s22_harrell_A.csv"),
               stringsAsFactors = FALSE, fileEncoding = "UTF-8-BOM")
names(hc)[1] <- "score_basis"
hc$harrell_C <- as.numeric(hc$harrell_C)
hc$se        <- as.numeric(hc$se)

# ---- ΔC 的配对 bootstrap（免硬编码：全部从脚本输出读入）----
bt <- read.csv(file.path(WK, "s22b_delta_boot.csv"),
               stringsAsFactors = FALSE, fileEncoding = "UTF-8-BOM")
gv <- function(k) as.numeric(bt$value[bt$item == k])

dC   <- gv("delta"); dlo <- gv("delta_lo"); dhi <- gv("delta_hi"); pbt <- gv("p_boot")
B    <- gv("B");     nn <- gv("n");        ev  <- gv("events")

# 防回归：官方值须与 s22_harrell_A_out.txt 全等；Δ 须与解析值一致；结论须仍为不显著
stopifnot(all(!is.na(hc$harrell_C)),
          abs(hc$harrell_C[hc$score_basis == "A"] - 0.6898) < 1e-6,
          abs(hc$harrell_C[hc$score_basis == "C"] - 0.6175) < 1e-6,
          abs(dC - (hc$harrell_C[hc$score_basis == "A"] -
                    hc$harrell_C[hc$score_basis == "C"])) < 5e-3,
          pbt > 0.05)   # 若差异变为显著，本图注释与本脚本假设失效 -> 必须人工复核

df <- data.frame(
  lab = factor(c("Score A\n(this study)", "Alternative\nscore (C)"),
               levels = c("Alternative\nscore (C)", "Score A\n(this study)")),
  C   = c(hc$harrell_C[hc$score_basis == "A"], hc$harrell_C[hc$score_basis == "C"]),
  se  = c(hc$se[hc$score_basis == "A"],        hc$se[hc$score_basis == "C"])
)
df$lo <- df$C - 1.96 * df$se
df$hi <- df$C + 1.96 * df$se
df$yc <- as.numeric(df$lab)

COL <- setNames(c("#7F8C8D", "#E74C3C"), levels(df$lab))

# 面板内注释：两行，起点置于 0.5 参考线右侧（x = 0.535），避开参考线
ann1 <- sprintf("Difference (A vs C) = %.3f", dC)
ann2 <- sprintf("95%% CI %.3f to %.3f, P = %.2f", dlo, dhi, pbt)

p <- ggplot(df, aes(y = yc, color = lab)) +
  # 随机基线
  geom_vline(xintercept = 0.5, linetype = "dashed", color = "grey45", size = 0.5) +
  # 95% CI 横线 + 两端帽（用 segment 画，规避 geom_errorbarh 的版本兼容问题）
  geom_segment(aes(x = lo, xend = hi, yend = yc), size = 0.75) +
  geom_segment(aes(x = lo, xend = lo, yend = yc - 0.11), size = 0.75) +
  geom_segment(aes(x = hi, xend = hi, yend = yc + 0.11), size = 0.75) +
  geom_point(aes(x = C), size = 4.2) +
  # 点估计标签置 CI 右端外侧，hjust = 0 -> 与误差棒、参考线均无重叠
  geom_text(aes(x = hi + 0.022, label = sprintf("%.3f", C)),
            hjust = 0, size = 3.5, fontface = "bold", show.legend = FALSE) +
  # 差异注释置于面板下部空白区，且整体位于参考线右侧
  # 0929 压矮重排（682->415）：原 y=0.76/0.50 与灰色 CI 条（下帽至 0.89）相压，
  #   下移注释块并把 ylim 下缘扩到 0.26 让出空间
  # 0929 二次修正：0.58/0.38 行距 0.2 单位 < 字高，两行相互挤压。
  #   ylim 下缘扩到 0.12 腾出竖向空间，行距拉开到 0.30 单位（>字高）
  annotate("text", x = 0.535, y = 0.60, label = ann1, hjust = 0, size = 3.0, color = "grey20") +
  annotate("text", x = 0.535, y = 0.30, label = ann2, hjust = 0, size = 3.0, color = "grey20") +
  scale_color_manual(values = COL) +
  scale_x_continuous(breaks = seq(0.5, 0.9, by = 0.1)) +
  scale_y_continuous(breaks = c(1, 2),
                     labels = c("Alternative\nscore (C)", "Score A\n(this study)")) +
  coord_cartesian(xlim = c(0.45, 0.98), ylim = c(0.12, 2.22)) +
  labs(x = "Harrell C (OS, TCGA-SARC; n = 50, 23 events)",
       y = "",
       title = "Discrimination of the score",
       subtitle = "Error bars = 95% CI; dashed line = chance (0.5)") +
  theme_publication() +
  theme(legend.position = "none")

export_fig(p, file.path(OUT, "fig3d.pdf"), w = 4.8, h = 3.5)
# df 行序：row1 = Score A，row2 = Alternative score (C)
cat(sprintf("fig3d done | A=%.4f (%.3f-%.3f)  C=%.4f (%.3f-%.3f)  dC=%+.4f (%.4f-%.4f) P=%.3f | n=%d events=%d B=%.0f\n",
            df$C[1], df$lo[1], df$hi[1], df$C[2], df$lo[2], df$hi[2], dC, dlo, dhi, pbt, nn, ev, B))
