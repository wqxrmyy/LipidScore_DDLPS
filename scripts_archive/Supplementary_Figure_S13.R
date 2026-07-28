library(GEOquery)

# 然后逐段执行 Figure S13 代码
gse159659 <- getGEO(filename = "raw/GSE159659/GSE159659_series_matrix.txt.gz",
                    GSEMatrix = TRUE, getGPL = FALSE)
pheno_159659 <- pData(gse159659)
expr_6genes_159659 <- readRDS("raw/GSE159659_6genes_expr.rds")

library(ggplot2)
library(ggpubr)
library(tidyr)

plot_data <- data.frame(
  Sample = colnames(expr_6genes_159659),
  Tissue = pheno_159659$`tissue:ch1`[match(colnames(expr_6genes_159659), 
                                           pheno_159659$geo_accession)],
  t(expr_6genes_159659)
)

plot_long <- pivot_longer(plot_data, 
                          cols = all_of(c("FASN", "ACLY", "SCD", "SREBF1", "HMGCR", "PPARG")),
                          names_to = "Gene", values_to = "Expression")

p_box <- ggplot(plot_long, aes(x = Tissue, y = Expression, fill = Tissue)) +
  geom_boxplot(outlier.shape = NA, alpha = 0.7) +
  geom_jitter(width = 0.15, size = 0.6, alpha = 0.5) +
  facet_wrap(~ Gene, scales = "free_y", ncol = 3) +
  stat_compare_means(ref.group = "adipose tissue", label = "p.signif", hide.ns = TRUE) +
  scale_fill_manual(values = c("#4DAF4A", "#E41A1C", "#377EB8")) +
  theme_bw(base_size = 12) +
  theme(
    legend.position = "bottom",
    axis.text.x = element_text(angle = 30, hjust = 1),
    strip.text = element_text(face = "bold")
  ) +
  labs(
    title = "Six-Gene Expression: DDLPS vs WDLPS vs Normal Adipose (GSE159659)",
    x = "",
    y = "Expression (log2)",
    fill = ""
  )

cairo_pdf("results/figures/extended/FigSX_GSE159659_boxplot.pdf", 
          width = 12, height = 8, family = "Arial")
print(p_box)
dev.off()

cat("Figure S13 已保存\n")