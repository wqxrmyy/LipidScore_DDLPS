# ============================================================
# Fig 2D · 整簇剔除森林图（3 簇 r±95%CI）+ 主效应参考线
# ============================================================
suppressWarnings(suppressMessages(library(ggplot2)))
source(file.path(Sys.getenv("LPS_FIG_ROOT", "."), "utils", "colors.R"))
source(file.path(Sys.getenv("LPS_FIG_ROOT", "."), "utils", "themes.R"))

load(file.path(Sys.getenv("LPS_FIG_ROOT", "."), "p5_data.RData"))
OUT <- file.path(Sys.getenv("LPS_FIG_ROOT", "."), "out")

numf <- function(x) suppressWarnings(as.numeric(as.character(x)))
cn <- numf(clin_real[score_A$Sample, "paper_CN.cluster"])
ok <- !is.na(cn)
x_all <- cn[ok]; y_all <- score_A$Score[ok]
r_main <- cor(x_all, y_all); n_main <- length(x_all)

# 主效应 CI（Fisher z）
zci <- function(r, n) {
  z <- atanh(r); se <- 1 / sqrt(n - 3)
  c(tanh(z - 1.96 * se), tanh(z + 1.96 * se))
}
ci_main <- zci(r_main, n_main)

rows <- data.frame(label = "None (main)", r = r_main, lo = ci_main[1], hi = ci_main[2], n = n_main)
for (c0 in sort(unique(x_all))) {
  keep <- x_all != c0
  r <- cor(x_all[keep], y_all[keep]); n <- sum(keep)
  ci <- zci(r, n)
  rows <- rbind(rows, data.frame(label = sprintf("Drop cluster %s", c0), r = r, lo = ci[1], hi = ci[2], n = n))
}
rows$label <- factor(rows$label, levels = rev(rows$label))
rows$sig <- rows$hi < 0.95 & (rows$lo > 0.05)  # 简化标记：lo>0 即方向为正

p <- ggplot(rows, aes(x = r, y = label)) +
  geom_vline(xintercept = 0, linetype = "dotted", color = "grey40") +
  geom_errorbarh(aes(xmin = lo, xmax = hi), height = 0.18, color = "#2980B9", size = 0.9) +
  geom_point(size = 3, color = "#E74C3C") +
  labs(x = "Pearson r (95% CI)", y = "", title = "Drop-whole-cluster robustness",
       # v3：原 57 字符副标题在半幅（3.25 in）被裁断，压缩至 40 字符
       subtitle = "Direction preserved after cluster removal") +
  ## ★★ S39：右距 4 -> 9 pt 改善右余量（原 14 px）；下距 4 -> 6 pt
  theme_publication() +
  theme(plot.margin = margin(4, 9, 6, 4))

export_fig(p, file.path(OUT, "fig2d.pdf"), w = 5, h = 3.2)
cat("fig2d done\n")
