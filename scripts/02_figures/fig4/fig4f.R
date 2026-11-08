# ============================================================
# Fig 4F · 身份基因 DD vs WD 小提琴（多基因并列；官方 log2CPM）
# 数据：fig4_genes.csv（prep_fig4_gse_official.py，CD36 7.228/10.536 = Table S5）
# ============================================================
suppressWarnings(suppressMessages(library(ggplot2)))
source(file.path(Sys.getenv("LPS_FIG_ROOT", "."), "utils", "colors.R"))
source(file.path(Sys.getenv("LPS_FIG_ROOT", "."), "utils", "themes.R"))

OUT <- file.path(Sys.getenv("LPS_FIG_ROOT", "."), "out")

sig_label <- function(p) {
  if (p < 0.001) "***" else if (p < 0.01) "**" else if (p < 0.05) "*" else "ns"
}

d <- read.csv(file.path(Sys.getenv("LPS_FIG_ROOT", "."), "fig4_genes.csv"),
              stringsAsFactors = FALSE)
show_genes <- c("CD36", "PPARG", "FABP4", "PLIN1", "LPL", "ADIPOQ")
d <- d[d$gene %in% show_genes & d$group %in% c("DD", "WD"), ]
d$group <- factor(d$group, levels = c("WD", "DD"))

pp <- sapply(show_genes, function(g) {
  d0 <- d[d$gene == g, ]
  suppressWarnings(wilcox.test(expr ~ group, data = d0)$p.value)
})
gene_lab <- paste0(show_genes, "\n", sapply(pp, sig_label))
names(gene_lab) <- show_genes
d$GeneLab <- gene_lab[as.character(d$gene)]
d$GeneLab <- factor(d$GeneLab, levels = gene_lab)

p <- ggplot(d, aes(x = group, y = expr, fill = group)) +
  geom_violin(scale = "width", alpha = 0.6) +
  geom_boxplot(width = 0.18, outlier.shape = NA, alpha = 0.9) +
  facet_wrap(~ GeneLab, scales = "free_y", ncol = 3) +
  scale_fill_manual(values = c(WD = "#2980B9", DD = "#E67E22")) +
  labs(x = "", y = "log2 CPM", title = "Identity genes: WD vs DD (GSE221492)") +
  theme_publication() +
  theme(legend.position = "bottom")

export_fig(p, file.path(OUT, "fig4f.pdf"), w = 7, h = 5)
cat("fig4f done\n")
