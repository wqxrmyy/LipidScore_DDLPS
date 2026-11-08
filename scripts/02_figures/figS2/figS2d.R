# ============================================================
# Fig S2D · 2 标志物（非分量）vs 评分森林图（TCGA, n = 50）
#   由原 figS2a 拆出（2026-09-28 重构：2 标志物 -> Fig. S2 D）
#   数据：d12_p4_robustness.json -> score_vs_genes（CD36/CD8A 键）
# ============================================================
suppressWarnings(suppressMessages(library(ggplot2)))
suppressWarnings(suppressMessages(library(jsonlite)))
source(file.path(Sys.getenv("LPS_FIG_ROOT", "."), "utils", "colors.R"))
source(file.path(Sys.getenv("LPS_FIG_ROOT", "."), "utils", "themes.R"))

OUT <- file.path(Sys.getenv("LPS_FIG_ROOT", "."), "out")

d12 <- fromJSON(file.path(Sys.getenv("LPS_WORK_ROOT", "."), "_work", "d12_p4_robustness.json"))
gv <- d12$score_vs_genes

# 构造森林图数据：每个变量的 r 与 95%CI（Fisher z 变换近似）
build_ci <- function(r, n) {
  z <- atanh(r); se <- 1 / sqrt(n - 3)
  lo <- tanh(z - 1.96 * se); hi <- tanh(z + 1.96 * se)
  c(lo = lo, hi = hi)
}

MARKERS <- c("CD36", "CD8A")

forest_plot <- function(genes, title) {
  rows <- lapply(genes, function(g) {
    r <- gv[[g]]$r; n <- gv[[g]]$n
    ci <- build_ci(r, n)
    data.frame(Var = g, r = r, lo = ci["lo"], hi = ci["hi"], n = n,
               P = gv[[g]]$p, stringsAsFactors = FALSE)
  })
  df <- do.call(rbind, rows)
  df <- df[order(df$r), ]
  df$Var <- factor(df$Var, levels = df$Var)
  df$sig <- ifelse(df$P < 0.05, "P < 0.05", "ns")

  ggplot(df, aes(x = r, y = Var, color = sig)) +
    geom_vline(xintercept = 0, linetype = "dashed", color = "grey50") +
    geom_errorbarh(aes(xmin = lo, xmax = hi), height = 0.25, size = 0.8) +
    geom_point(size = 2.5) +
    scale_color_manual(values = c(`P < 0.05` = "#E74C3C", ns = "#7F8C8D"), name = "") +
    labs(x = "Pearson r (vs score A)", y = "", title = title) +
    theme_publication()
}

export_fig(forest_plot(MARKERS, "Markers, not components (TCGA, n = 50)"),
           file.path(OUT, "figS2d.pdf"))
cat("figS2d done\n")
