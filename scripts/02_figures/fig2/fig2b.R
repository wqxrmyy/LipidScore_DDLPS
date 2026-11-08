# ============================================================
# Fig 2B · 留一法稳健性（50 次 LOO 的 r 分布直方图）
# 数据：d12_p4_robustness.json -> D_leave_one_out
# ============================================================
suppressWarnings(suppressMessages(library(ggplot2)))
source(file.path(Sys.getenv("LPS_FIG_ROOT", "."), "utils", "colors.R"))
source(file.path(Sys.getenv("LPS_FIG_ROOT", "."), "utils", "themes.R"))

OUT <- file.path(Sys.getenv("LPS_FIG_ROOT", "."), "out")

# 读 d12 结果（含 50 次 LOO 的 r）
d12 <- jsonlite::fromJSON(file.path(Sys.getenv("LPS_WORK_ROOT", "."), "_work", "d12_p4_robustness.json"))

# 重建完整 LOO 分布：d12 JSON 只存 worst5/best3 + 汇总；用主脚本重算（快）
load(file.path(Sys.getenv("LPS_FIG_ROOT", "."), "p5_data.RData"))

numf <- function(x) suppressWarnings(as.numeric(as.character(x)))
cn <- numf(clin_real[score_A$Sample, "paper_CN.cluster"])
ok <- !is.na(cn)
x <- cn[ok]; y <- score_A$Score[ok]

r_loo <- sapply(seq_along(x), function(i) cor(x[-i], y[-i]))
p_loo <- sapply(seq_along(x), function(i) cor.test(x[-i], y[-i])$p.value)

df <- data.frame(r = r_loo, sig = p_loo < 0.05)
n_sig <- sum(df$sig)

p <- ggplot(df, aes(x = r)) +
  geom_histogram(bins = 12, fill = "#2980B9", alpha = 0.75, color = "grey20") +
  geom_vline(xintercept = min(r_loo), linetype = "dashed", color = "#E74C3C", size = 0.8) +
  annotate("text", x = min(r_loo), y = Inf, hjust = -0.05, vjust = 1.2,
           label = sprintf("min r = %.3f\n%d/%d still P<0.05", min(r_loo), n_sig, nrow(df)),
           size = 3.2, color = "#E74C3C") +
  labs(x = "Pearson r (leave-one-out)", y = "Count",
       title = "Leave-one-out robustness (n = 50)") +
  ## ★★ S39：右距 4 -> 9 pt 改善右余量（原 14 px）；上距 4 -> 8 pt
  ##   容纳 annotate(y=Inf, vjust=1.2) 的顶部标注
  theme_publication() +
  theme(plot.margin = margin(8, 9, 6, 4))

export_fig(p, file.path(OUT, "fig2b.pdf"), w = 4.5, h = 3.5)
cat("fig2b done\n")
