# ============================================================
# Stage 3 Multiomics: ssGSEA Immune Signatures + Immune Landscape
# 功能：功能性免疫特征分析、免疫排斥表型验证、免疫相关性热图
# ============================================================

setwd('C:/Users/Tao25/Documents/Liposarcoma/2')
library(GSVA); library(org.Hs.eg.db); library(ggplot2); library(pheatmap)

# ---- 1. ssGSEA 免疫特征分析 ----
tcga_expr <- readRDS('raw/TCGA_LPS_real_log2cpm.rds')

# ENSG → Gene Symbol 转换
ensg_ids <- gsub('\\..*', '', rownames(tcga_expr))
symbol_map <- mapIds(org.Hs.eg.db, keys = ensg_ids, column = 'SYMBOL',
                     keytype = 'ENSEMBL', multiVals = 'first')
expr_df <- as.data.frame(tcga_expr)
expr_df$symbol <- symbol_map[ensg_ids]
expr_df <- expr_df[!is.na(expr_df$symbol), ]
expr_df <- expr_df[order(rowMeans(expr_df[, -ncol(expr_df)]), decreasing = TRUE), ]
expr_df <- expr_df[!duplicated(expr_df$symbol), ]
rownames(expr_df) <- expr_df$symbol; expr_df$symbol <- NULL
tcga_symbol <- as.matrix(expr_df)

# 五特征免疫基因集
immune_genesets <- list(
  T_cell_exhaustion = c('PDCD1', 'CTLA4', 'HAVCR2', 'LAG3', 'TIGIT'),
  Cytotoxic = c('GZMA', 'GZMB', 'PRF1', 'GNLY', 'NKG7'),
  IFN_gamma = c('IFNG', 'STAT1', 'CXCL9', 'CXCL10', 'CXCL11'),
  MDSC = c('CD33', 'CD14', 'IL4I1', 'ARG1'),
  Treg = c('FOXP3', 'IL2RA', 'IKZF2', 'CTLA4')
)
for (sig in names(immune_genesets)) {
  immune_genesets[[sig]] <- intersect(immune_genesets[[sig]], rownames(tcga_symbol))
}
immune_genesets <- immune_genesets[lengths(immune_genesets) >= 2]

# ssGSEA 计算
param <- gsvaParam(tcga_symbol, immune_genesets)
immune_scores <- gsva(param, verbose = FALSE)
immune_scores <- t(immune_scores)

# 合并脂评分
tcga_score <- read.csv('raw/TCGA_Score.csv')
tcga_score$Sample <- gsub('-', '.', tcga_score$Sample)
common <- intersect(rownames(immune_scores), tcga_score$Sample)
immune_df <- data.frame(Sample = common, LipidScore = tcga_score$Score[match(common, tcga_score$Sample)], immune_scores[common, ])
immune_df$Group <- ifelse(immune_df$LipidScore > median(immune_df$LipidScore), 'High', 'Low')

# 相关性 + 组间比较
cat('\nLipidScore vs Immune Signatures:\n')
for (sig in names(immune_genesets)) {
  cor_res <- cor.test(immune_df$LipidScore, immune_df[[sig]], method = 'spearman')
  cat(sprintf('  %s: R=%.3f, P=%.4f\n', sig, cor_res$estimate, cor_res$p.value))
}
write.csv(immune_df, 'results/tables/extended/immune_scores_vs_lipid.csv', row.names = FALSE)

# ---- 2. 免疫相关性热图 ----
cor_matrix <- cor(immune_df[, -c(1,2)], method = 'spearman')
cairo_pdf('results/figures/extended/FigS14_immune_correlation.pdf', width = 8, height = 6)
pheatmap(cor_matrix, color = colorRampPalette(c('#377EB8', 'white', '#E41A1C'))(100),
         display_numbers = TRUE, number_format = '%.3f',
         main = 'LipidScore vs Immune Signatures (TCGA-SARC)', fontsize = 11)
dev.off()

cat('\nStage 3 Multiomics 完成\n')
