# ============================================================
# Fig 2C · 置换检验零分布（20000 次直方图 + 观测值红线）
# 观测 r=0.5354, P_perm=0.0001
# ============================================================
suppressWarnings(suppressMessages(library(ggplot2)))
source(file.path(Sys.getenv("LPS_FIG_ROOT", "."), "utils", "colors.R"))
source(file.path(Sys.getenv("LPS_FIG_ROOT", "."), "utils", "themes.R"))

load(file.path(Sys.getenv("LPS_FIG_ROOT", "."), "p5_data.RData"))
OUT <- file.path(Sys.getenv("LPS_FIG_ROOT", "."), "out")

set.seed(20260917)
numf <- function(x) suppressWarnings(as.numeric(as.character(x)))
cn <- numf(clin_real[score_A$Sample, "paper_CN.cluster"])
ok <- !is.na(cn)
x <- cn[ok]; y <- score_A$Score[ok]
r_obs <- cor(x, y)

NP <- 20000
r_perm <- replicate(NP, cor(sample(x), y))

df <- data.frame(r = r_perm)
cnt <- sum(abs(r_perm) >= abs(r_obs))
p_perm <- (cnt + 1) / (NP + 1)

p <- ggplot(df, aes(x = r)) +
  geom_histogram(bins = 60, fill = "#7F8C8D", alpha = 0.8, color = "white") +
  geom_vline(xintercept = r_obs, color = "#E74C3C", size = 1.1) +
  annotate("text", x = r_obs, y = Inf, hjust = 1.08, vjust = 1.1,
           label = sprintf("observed r = %.3f\nP(perm) = %.4g", r_obs, p_perm),
           size = 3.2, color = "#E74C3C") +
  coord_cartesian(xlim = c(-0.55, 0.68)) +
  labs(x = "Pearson r under permutation", y = "Count",
       title = sprintf("Permutation test (%d permutations)", NP)) +
  ## ★★ S39：右距 4 -> 9 pt（原右余量 14 px）；上距 4 -> 8 pt
  ##   容纳 annotate(y=Inf, vjust=1.1) 的顶部标注
  theme_publication() +
  theme(plot.margin = margin(8, 9, 6, 4))

export_fig(p, file.path(OUT, "fig2c.pdf"), w = 4.5, h = 3.5)
cat("fig2c done\n")
