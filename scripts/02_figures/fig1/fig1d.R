# ============================================================
# Fig 1D · 6 基因表达小提琴（多基因并列）
# ============================================================
suppressWarnings(suppressMessages(library(ggplot2)))
suppressWarnings(suppressMessages(library(reshape2)))
source(file.path(Sys.getenv("LPS_FIG_ROOT", "."), "utils", "colors.R"))
source(file.path(Sys.getenv("LPS_FIG_ROOT", "."), "utils", "themes.R"))

load(file.path(Sys.getenv("LPS_FIG_ROOT", "."), "p5_data.RData"))
OUT <- file.path(Sys.getenv("LPS_FIG_ROOT", "."), "out")

melt_expr <- melt(expr6)
names(melt_expr) <- c("Gene", "Sample", "log2CPM")

p <- ggplot(melt_expr, aes(x = Gene, y = log2CPM, fill = Gene)) +
  geom_violin(scale = "width", alpha = 0.7) +
  geom_boxplot(width = 0.15, outlier.shape = NA) +
  geom_jitter(width = 0.1, size = 0.6, alpha = 0.5) +
  scale_fill_manual(values = c(SREBF1 = "#7F8C8D", SCD = "#E74C3C", HMGCR = "#2980B9",
                               ACLY = "#E67E22", PPARG = "#27AE60", FASN = "#8E44AD")) +
  labs(x = "", y = "log2(CPM)", title = "Six-gene expression (TCGA, n=50)") +
  theme_publication() +
  theme(legend.position = "none", axis.text.x = element_text(angle = 45, hjust = 1))

export_fig(p, file.path(OUT, "fig1d.pdf"), w = 5, h = 4)
cat("fig1d done\n")
