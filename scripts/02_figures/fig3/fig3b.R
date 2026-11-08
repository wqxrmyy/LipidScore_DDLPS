# ============================================================
# Fig 3B · CD36 vs 评分散点（GSE221492 外部复现，官方口径）
# 数据：figS4d_gse_data.csv（prep_figS4d_gse.py 复刻 d11 官方管线：
#       raw counts -> log2CPM -> DDLPS n=36 内 z 评分；CD36 log2CPM）
# 官方值：Spearman rho = 0.5086, P = 0.001538（d11 Aprime_replication.CD36）
# 注：P 值直接用官方 JSON 值（R cor.test 的 AS89 近似给 0.0018，与 scipy 略异）
# ============================================================
suppressWarnings(suppressMessages(library(ggplot2)))
source(file.path(Sys.getenv("LPS_FIG_ROOT", "."), "utils", "colors.R"))
source(file.path(Sys.getenv("LPS_FIG_ROOT", "."), "utils", "themes.R"))

OUT <- file.path(Sys.getenv("LPS_FIG_ROOT", "."), "out")

df <- read.csv(file.path(Sys.getenv("LPS_FIG_ROOT", "."), "figS4d_gse_data.csv"),
               stringsAsFactors = FALSE)
rho_off <- 0.5086; p_off <- 0.001538
stopifnot(nrow(df) == 36)

p <- ggplot(df, aes(x = score, y = cd36_log2cpm)) +
  geom_point(size = 2, alpha = 0.8, color = "#E67E22") +
  geom_smooth(method = "lm", se = TRUE, color = "#2980B9", size = 0.8) +
  labs(x = "Lipid metabolism score (A)", y = "CD36 expression (log2 CPM)",
       title = "CD36 association replicated (GSE221492)",
       subtitle = sprintf("Spearman rho = %.3f, P = %.4g (n = %d)", rho_off, p_off, nrow(df))) +
  theme_publication()

export_fig(p, file.path(OUT, "fig3b.pdf"), w = 4.2, h = 3.5)
cat("fig3b done\n")
