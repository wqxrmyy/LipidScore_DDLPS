# ============================================================
# Figure SX: 六基因相关性矩阵（TCGA 发现队列）
# ============================================================
library(corrplot)

tcga_6genes <- readRDS("raw/TCGA_6genes_expr.rds")

cor_matrix <- cor(t(tcga_6genes), method = "spearman")

cairo_pdf("results/figures/extended/FigSX_gene_correlation.pdf", 
          width = 7, height = 6, family = "Arial")
corrplot(cor_matrix, 
         method = "color",
         type = "upper",
         order = "hclust",
         addCoef.col = "black",
         number.cex = 0.9,
         tl.col = "black",
         tl.cex = 1.1,
         col = colorRampPalette(c("#377EB8", "white", "#E41A1C"))(200),
         title = "Six-Gene Spearman Correlation (TCGA-SARC)",
         mar = c(1, 1, 2, 1))
dev.off()

cat("六基因相关性矩阵已保存\n")