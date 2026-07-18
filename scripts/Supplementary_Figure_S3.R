# ============================================
# Supplementary Figure S3: Comprehensive Analysis
# 优化版 - 消除透明度警告
# ============================================

rm(list = ls())
setwd("C:/Users/Tao25/Documents/Liposarcoma/1")

library(ggplot2)
library(pheatmap)
library(reshape2)
library(patchwork)
library(AnnotationDbi)
library(org.Hs.eg.db)
library(ggplotify)
library(survival)

# ========== 颜色定义 ==========
dark_blue <- "#073791"
bright_red <- "#e5051f"
light_green <- "#89c38e"
lake_blue <- "#007bc4"
indigo <- "#6053a0"
warm_orange <- "#ffa700"
gray_mid <- "#999999"

dir.create("figures", recursive = TRUE, showWarnings = FALSE)
dir.create("tables", recursive = TRUE, showWarnings = FALSE)

# ========== 数据加载 ==========
scores <- read.csv("data/results/02_feature_scores.csv", row.names = 1)
clinical <- read.csv("data/results/03_clinical_outcomes.csv", row.names = 1)
log2cpm_raw <- readRDS("data/raw/TCGA_LPS_real_log2cpm.rds")

# 基因名转换
ensembl_ids <- gsub("\\.[0-9]+$", "", rownames(log2cpm_raw))
gene_symbols <- mapIds(org.Hs.eg.db, keys = ensembl_ids, column = "SYMBOL",
                       keytype = "ENSEMBL", multiVals = "first")
keep <- !is.na(gene_symbols)
log2cpm_filtered <- log2cpm_raw[keep, ]
rownames(log2cpm_filtered) <- gene_symbols[keep]

log2cpm <- aggregate(log2cpm_filtered, by = list(rownames(log2cpm_filtered)), FUN = mean)
rownames(log2cpm) <- log2cpm$Group.1
log2cpm$Group.1 <- NULL
log2cpm <- as.matrix(log2cpm)

# 样本过滤
common_samples <- intersect(colnames(log2cpm), rownames(scores))
common_samples <- intersect(common_samples, rownames(clinical))
log2cpm <- log2cpm[, common_samples]
scores <- scores[common_samples, ]
clinical <- clinical[common_samples, ]

valid_idx <- !is.na(clinical$OS_time) & !is.na(clinical$OS_status)
log2cpm <- log2cpm[, valid_idx]
scores <- scores[valid_idx, ]
clinical <- clinical[valid_idx, ]

# 分组
median_score <- median(scores$Lipid_Metabolism_Score)
group <- ifelse(scores$Lipid_Metabolism_Score > median_score, "High", "Low")
group <- factor(group, levels = c("Low", "High"))

# ============================================
# 子图 A: 免疫浸润热图
# ============================================
immune_genesets <- list(
  "CD8+ T" = c("CD8A", "CD8B"),
  "CD4+ T" = c("CD4", "CD40LG"),
  "Treg" = c("FOXP3", "IL2RA"),
  "B" = c("CD19", "MS4A1"),
  "NK" = c("NCAM1", "KLRK1"),
  "Macrophage" = c("CD68", "CD163")
)

immune_scores <- matrix(NA, nrow = length(immune_genesets), ncol = ncol(log2cpm))
rownames(immune_scores) <- names(immune_genesets)
colnames(immune_scores) <- colnames(log2cpm)

for(i in 1:length(immune_genesets)) {
  genes <- immune_genesets[[i]]
  available_genes <- genes[genes %in% rownames(log2cpm)]
  if(length(available_genes) >= 1) {
    immune_scores[i, ] <- colMeans(log2cpm[available_genes, , drop = FALSE], na.rm = TRUE)
  }
}
immune_scores <- immune_scores[!apply(is.na(immune_scores), 1, all), ]

annot_col <- data.frame(Group = group)
rownames(annot_col) <- colnames(immune_scores)

pA <- as.ggplot(function() {
  pheatmap(immune_scores,
           scale = "row",
           color = colorRampPalette(c(dark_blue, "white", bright_red))(100),
           annotation_col = annot_col,
           annotation_colors = list(Group = c(Low = dark_blue, High = bright_red)),
           show_colnames = FALSE,
           fontsize = 7,
           main = "A",
           border_color = NA,
           legend = TRUE)
}) + theme(plot.title = element_text(hjust = 0.5, face = "bold", size = 12, color = dark_blue))

# ============================================
# 子图 B: 免疫浸润箱线图
# ============================================
immune_long <- melt(immune_scores)
colnames(immune_long) <- c("CellType", "Sample", "Score")
immune_long$Group <- rep(group, each = nrow(immune_scores))

# 计算p值
cell_types <- unique(immune_long$CellType)
p_values <- numeric(length(cell_types))
names(p_values) <- cell_types

for(ct in cell_types) {
  sub <- immune_long[immune_long$CellType == ct, ]
  if(length(unique(sub$Group)) == 2) {
    test_result <- tryCatch({
      wilcox.test(Score ~ Group, data = sub)$p.value
    }, error = function(e) NA)
    p_values[ct] <- ifelse(is.null(test_result), NA, test_result)
  } else {
    p_values[ct] <- NA
  }
}

p_labels <- ifelse(is.na(p_values), "",
                   ifelse(p_values < 0.001, "***",
                          ifelse(p_values < 0.01, "**",
                                 ifelse(p_values < 0.05, "*", "ns"))))

p_labels_df <- data.frame(
  CellType = cell_types,
  label = p_labels,
  y = max(immune_long$Score, na.rm = TRUE) * 1.08
)
p_labels_df <- p_labels_df[p_labels_df$label != "ns" & p_labels_df$label != "", ]

pB <- ggplot(immune_long, aes(x = CellType, y = Score, fill = Group)) +
  geom_boxplot(alpha = 0.75, outlier.size = 0.4, position = position_dodge(0.8), width = 0.7) +
  scale_fill_manual(values = c(Low = dark_blue, High = bright_red)) +
  labs(x = "", y = "Mean Expression (log2CPM)", title = "B") +
  theme_bw(base_size = 8) +
  theme(axis.text.x = element_text(angle = 45, hjust = 1, size = 8),
        axis.text.y = element_text(size = 8),
        plot.title = element_text(hjust = 0.5, face = "bold", size = 12, color = dark_blue),
        legend.position = "bottom",
        legend.title = element_text(size = 8),
        legend.text = element_text(size = 7),
        legend.key.size = unit(0.5, "cm"))

if(nrow(p_labels_df) > 0) {
  pB <- pB + geom_text(data = p_labels_df, 
                       aes(x = CellType, y = y, label = label),
                       size = 4, color = "black", inherit.aes = FALSE)
}

# ============================================
# 子图 C: 药物敏感性预测
# ============================================
drug_map <- list(
  "FASN inhibitors" = "FASN",
  "ACLY inhibitors" = "ACLY",
  "SCD inhibitors" = "SCD",
  "CDK4/6 inhibitors" = c("CDK4", "CDK6"),
  "Aurora kinase inhibitors" = c("AURKA", "AURKB"),
  "mTOR inhibitors" = c("MTOR", "RICTOR", "RPTOR"),
  "Immune checkpoint" = c("PDCD1", "CTLA4", "LAG3", "TIGIT")
)

drug_res <- data.frame()
for(drug in names(drug_map)) {
  genes <- drug_map[[drug]]
  avail <- genes[genes %in% rownames(log2cpm)]
  if(length(avail) == 0) next
  if(length(avail) == 1) {
    high <- log2cpm[avail, group == "High"]
    low <- log2cpm[avail, group == "Low"]
  } else {
    high <- colMeans(log2cpm[avail, group == "High", drop = FALSE])
    low <- colMeans(log2cpm[avail, group == "Low", drop = FALSE])
  }
  fc <- mean(high) / mean(low)
  pv <- tryCatch(wilcox.test(high, low)$p.value, error = function(e) 1)
  sens <- ifelse(fc > 1.3 & pv < 0.05, "High-score sensitive",
                 ifelse(fc < 0.77 & pv < 0.05, "Low-score sensitive", "Similar"))
  drug_res <- rbind(drug_res, data.frame(Drug = drug, FoldChange = fc, Sensitivity = sens, P_value = pv))
}

drug_res <- drug_res[order(drug_res$FoldChange, decreasing = TRUE), ]
drug_res$Drug <- factor(drug_res$Drug, levels = rev(drug_res$Drug))
drug_res$p_label <- ifelse(drug_res$Sensitivity != "Similar", "*", "")

pC <- ggplot(drug_res, aes(x = Drug, y = FoldChange, fill = Sensitivity)) +
  geom_bar(stat = "identity", width = 0.7) + 
  coord_flip() +
  scale_fill_manual(values = c("High-score sensitive" = bright_red,
                               "Low-score sensitive" = dark_blue,
                               "Similar" = gray_mid)) +
  geom_hline(yintercept = 1, linetype = "dashed", color = "gray50", linewidth = 0.5) +
  geom_text(aes(label = paste0(round(FoldChange, 2), p_label)), 
            hjust = ifelse(drug_res$FoldChange > 1, -0.15, 1.15), 
            size = 2.8, color = "gray20") +
  labs(x = "", y = "Fold Change (High/Low)", title = "C") +
  theme_bw(base_size = 8) +
  theme(plot.title = element_text(hjust = 0.5, face = "bold", size = 12, color = dark_blue),
        axis.text.y = element_text(size = 7),
        axis.text.x = element_text(size = 8),
        legend.position = "bottom",
        legend.title = element_text(size = 8),
        legend.text = element_text(size = 7),
        legend.key.size = unit(0.4, "cm"))

# ============================================
# 子图 D: 基因相关性网络
# ============================================
network_genes <- c("FASN", "ACLY", "SCD", "SREBF1", "HMGCR", "PPARG")
available_net <- network_genes[network_genes %in% rownames(log2cpm)]

if(length(available_net) >= 2) {
  expr_net <- t(log2cpm[available_net, ])
  cor_mat <- cor(expr_net, method = "spearman")
  cor_melt <- melt(cor_mat)
  colnames(cor_melt) <- c("Gene1", "Gene2", "Correlation")
  
  pD <- ggplot(cor_melt, aes(x = Gene1, y = Gene2, fill = Correlation)) +
    geom_tile(color = "white", linewidth = 0.3) +
    scale_fill_gradient2(low = dark_blue, high = bright_red, mid = "white",
                         midpoint = 0, limit = c(-1, 1), 
                         name = "Spearman\ncorrelation") +
    geom_text(aes(label = round(Correlation, 2)), size = 2.8, color = "gray20") +
    labs(x = "", y = "", title = "D") +
    theme_bw(base_size = 8) +
    theme(axis.text.x = element_text(angle = 45, hjust = 1, size = 8),
          axis.text.y = element_text(size = 8),
          plot.title = element_text(hjust = 0.5, face = "bold", size = 12, color = dark_blue),
          legend.position = "bottom",
          legend.title = element_text(size = 7),
          legend.text = element_text(size = 6),
          legend.key.width = unit(0.6, "cm"),
          panel.grid = element_blank())
} else {
  pD <- ggplot() + 
    annotate("text", x = 0.5, y = 0.5, label = "Insufficient genes for correlation", size = 4) +
    labs(title = "D") +
    theme_bw() +
    theme(plot.title = element_text(hjust = 0.5, face = "bold", size = 12, color = dark_blue))
}

# ============================================
# 组合 Figure S3 (2x2 布局)
# ============================================
figure_S3 <- (pA | pB) / (pC | pD) +
  plot_layout(heights = c(1, 1.1), widths = c(1.1, 1)) +
  plot_annotation(
    title = "Supplementary Figure S3: Multi-Omics Landscape of Lipid Metabolism",
    subtitle = "A-B) Immune infiltration analysis | C) Drug sensitivity prediction | D) Gene correlation network",
    caption = "Note: *p < 0.05, **p < 0.01, ***p < 0.001; ns, not significant. Drug sensitivity based on target gene expression.",
    theme = theme(
      plot.title = element_text(hjust = 0.5, face = "bold", size = 14, color = dark_blue),
      plot.subtitle = element_text(hjust = 0.5, size = 10, color = indigo),
      plot.caption = element_text(hjust = 0, size = 8, color = "gray40")
    )
  )

# ========== 保存多种格式（优化版） ==========

# 1. PDF格式（矢量图，适合论文发表）- 无透明度警告
ggsave("figures/Supplementary_Figure_S3.pdf", 
       figure_S3, width = 12, height = 11, dpi = 300, device = "pdf")
cat("✅ PDF已保存: figures/Supplementary_Figure_S3.pdf\n")

# 2. PNG格式（位图，适合预览）
ggsave("figures/Supplementary_Figure_S3.png", 
       figure_S3, width = 12, height = 11, dpi = 300, device = "png")
cat("✅ PNG已保存: figures/Supplementary_Figure_S3.png\n")

# 3. TIFF格式（高分辨率位图，期刊常用）
ggsave("figures/Supplementary_Figure_S3.tiff", 
       figure_S3, width = 12, height = 11, dpi = 300, 
       device = "tiff", compression = "lzw")
cat("✅ TIFF已保存: figures/Supplementary_Figure_S3.tiff\n")

# 4. 高质量TIFF（600 dpi，用于印刷）
ggsave("figures/Supplementary_Figure_S3_600dpi.tiff", 
       figure_S3, width = 12, height = 11, dpi = 600, 
       device = "tiff", compression = "lzw")
cat("✅ 高质量TIFF (600dpi)已保存: figures/Supplementary_Figure_S3_600dpi.tiff\n")

# 5. EPS格式（矢量图）- 先转为不支持透明度的格式
# 使用cairo_ps设备避免透明度警告
cairo_ps(file = "figures/Supplementary_Figure_S3.eps", 
         width = 12, height = 11, onefile = TRUE)
print(figure_S3)
dev.off()
cat("✅ EPS已保存: figures/Supplementary_Figure_S3.eps\n")

# 6. SVG格式（矢量图，可编辑）
ggsave("figures/Supplementary_Figure_S3.svg", 
       figure_S3, width = 12, height = 11, dpi = 300, device = "svg")
cat("✅ SVG已保存: figures/Supplementary_Figure_S3.svg\n")

cat("\n========================================\n")
cat("✅ 所有格式已保存完成！\n")
cat("📁 输出路径: figures/\n")
cat("📊 布局: 2行x2列 (A|B / C|D)\n")
cat("📝 格式: PDF, PNG, TIFF (300dpi), TIFF (600dpi), EPS, SVG\n")
cat("========================================\n")