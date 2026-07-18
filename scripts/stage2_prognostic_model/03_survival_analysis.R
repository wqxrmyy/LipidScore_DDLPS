# ============================================================
# Figure 2: Lipid Metabolism Landscape
# 出版级规范版本（修复下标出界错误）
# 子图: A, B, C, D, E, F, G, H
# ============================================================

rm(list = ls())
setwd("C:/Users/Tao25/Documents/Liposarcoma/1")

library(ggplot2)
library(patchwork)
library(ggpubr)
library(reshape2)
library(ggridges)
library(grid)
library(AnnotationDbi)
library(org.Hs.eg.db)

# ============================================================
# 1. 颜色系统
# ============================================================

dark_blue <- "#0F4C81"
soft_pink <- "#E69F8C"
light_green <- "#8DAA7D"
warm_orange <- "#E69F00"
bright_red <- "#C93312"
lake_blue <- "#6BAED6"
indigo <- "#6053a0"
gray_light <- "#E8E8E8"

dir.create("figures", recursive = TRUE, showWarnings = FALSE)

# ============================================================
# 2. 数据读取
# ============================================================

scores <- read.csv("data/results/02_feature_scores.csv", row.names = 1)
clinical <- read.csv("data/results/03_clinical_outcomes.csv", row.names = 1)
tcga_raw <- read.csv("data/raw/TCGA_clinical_real.csv", row.names = 1)

common_rows <- intersect(rownames(scores), rownames(clinical))
scores <- scores[common_rows, ]
clinical <- clinical[common_rows, ]
valid_idx <- !is.na(scores$Lipid_Metabolism_Score)
scores_valid <- scores[valid_idx, ]
clinical_valid <- clinical[valid_idx, ]

tcga_idx <- clinical_valid$Cohort == "TCGA"
tcga_data <- clinical_valid[tcga_idx, ]
tcga_scores <- scores_valid[tcga_idx, ]

tcga_data$Gender <- tcga_raw[rownames(tcga_data), "gender"]
tcga_data$Gender <- ifelse(tcga_data$Gender == "male", "Male", "Female")

# ============================================================
# 3. 图 A: 亚型密度分布
# ============================================================

main_subtypes <- c("Dedifferentiated", "DDLPS", "Well_differentiated", "Myxoid", "Pleomorphic")
score_df <- data.frame(Subtype = clinical_valid$Subtype, Score = scores_valid$Lipid_Metabolism_Score)
score_df <- score_df[score_df$Subtype %in% main_subtypes, ]
score_df$Subtype <- factor(score_df$Subtype, levels = rev(main_subtypes))

pA <- ggplot(score_df, aes(x = Score, y = Subtype, fill = Subtype)) +
  geom_density_ridges(alpha = 0.7, scale = 1.2, rel_min_height = 0.01) +
  scale_fill_manual(values = c(soft_pink, light_green, warm_orange, lake_blue, indigo)) +
  labs(x = "Lipid Score", y = "", title = "A") +
  theme_bw(base_size = 10) +
  theme(
    legend.position = "none",
    plot.title = element_text(hjust = 0.5, face = "plain", size = 11, color = dark_blue,
                              family = "serif", margin = margin(b = 2)),
    axis.text.x = element_text(size = 9, color = "#333333", family = "serif"),
    axis.text.y = element_text(size = 9, color = "#333333", family = "serif"),
    axis.title.x = element_text(size = 10, color = "#333333", family = "serif"),
    panel.grid.minor = element_blank(),
    panel.border = element_rect(color = "gray70", linewidth = 0.5)
  )

# ============================================================
# 4. 图 B: 队列密度分布
# ============================================================

cohort_df <- data.frame(
  Cohort = factor(clinical_valid$Cohort, levels = rev(c("TCGA", "GSE21122", "GSE159659"))),
  Score = scores_valid$Lipid_Metabolism_Score
)

pB <- ggplot(cohort_df, aes(x = Score, y = Cohort, fill = Cohort)) +
  geom_density_ridges(alpha = 0.7, scale = 1.2, rel_min_height = 0.01) +
  scale_fill_manual(values = c(dark_blue, warm_orange, light_green)) +
  labs(x = "Lipid Score", y = "", title = "B") +
  theme_bw(base_size = 10) +
  theme(
    legend.position = "none",
    plot.title = element_text(hjust = 0.5, face = "plain", size = 11, color = dark_blue,
                              family = "serif", margin = margin(b = 2)),
    axis.text.x = element_text(size = 9, color = "#333333", family = "serif"),
    axis.text.y = element_text(size = 9, color = "#333333", family = "serif"),
    axis.title.x = element_text(size = 10, color = "#333333", family = "serif"),
    panel.grid.minor = element_blank(),
    panel.border = element_rect(color = "gray70", linewidth = 0.5)
  )

# ============================================================
# 5. 图 C: 性别差异
# ============================================================

gender_df <- na.omit(data.frame(
  Gender = tcga_data$Gender,
  Score = tcga_scores$Lipid_Metabolism_Score
))

pC <- ggplot(gender_df, aes(x = Gender, y = Score, fill = Gender)) +
  geom_violin(alpha = 0.7, trim = FALSE, width = 0.8) +
  geom_boxplot(width = 0.15, fill = "white", alpha = 0.8, outlier.size = 0.5) +
  stat_compare_means(method = "wilcox.test", label = "p.signif", size = 4, label.x = 1.5) +
  scale_fill_manual(values = c(soft_pink, lake_blue)) +
  labs(x = "", y = "Lipid Score", title = "C") +
  theme_bw(base_size = 10) +
  theme(
    legend.position = "none",
    plot.title = element_text(hjust = 0.5, face = "plain", size = 11, color = dark_blue,
                              family = "serif", margin = margin(b = 2)),
    axis.text.x = element_text(size = 10, color = "#333333", face = "plain", family = "serif"),
    axis.text.y = element_text(size = 9, color = "#333333", family = "serif"),
    axis.title.y = element_text(size = 10, color = "#333333", family = "serif"),
    panel.grid.minor = element_blank(),
    panel.border = element_rect(color = "gray70", linewidth = 0.5)
  )

# ============================================================
# 6. 图 D: 生存状态差异
# ============================================================

surv_idx <- !is.na(clinical_valid$OS_time) & !is.na(clinical_valid$OS_status)
status_df <- data.frame(
  Status = factor(clinical_valid$OS_status[surv_idx], labels = c("Alive", "Dead")),
  Score = scores_valid$Lipid_Metabolism_Score[surv_idx]
)

pD <- ggplot(status_df, aes(x = Status, y = Score, fill = Status)) +
  geom_boxplot(alpha = 0.8, width = 0.6, outlier.size = 0.5) +
  stat_compare_means(method = "wilcox.test", label = "p.signif", size = 4, label.x = 1.5) +
  scale_fill_manual(values = c(light_green, bright_red)) +
  labs(x = "", y = "Lipid Score", title = "D") +
  theme_bw(base_size = 10) +
  theme(
    legend.position = "none",
    plot.title = element_text(hjust = 0.5, face = "plain", size = 11, color = dark_blue,
                              family = "serif", margin = margin(b = 2)),
    axis.text.x = element_text(size = 10, color = "#333333", face = "plain", family = "serif"),
    axis.text.y = element_text(size = 9, color = "#333333", family = "serif"),
    axis.title.y = element_text(size = 10, color = "#333333", family = "serif"),
    panel.grid.major.x = element_blank(),
    panel.grid.minor = element_blank(),
    panel.border = element_rect(color = "gray70", linewidth = 0.5)
  )

# ============================================================
# 7. 图 E: 评分分布
# ============================================================

median_score <- median(scores_valid$Lipid_Metabolism_Score, na.rm = TRUE)

pE <- ggplot(scores_valid, aes(x = Lipid_Metabolism_Score)) +
  geom_histogram(bins = 30, fill = indigo, alpha = 0.7, color = "white", linewidth = 0.2) +
  geom_density(aes(y = after_stat(count) * 0.6), color = bright_red, linewidth = 0.8) +
  geom_vline(xintercept = median_score, linetype = "dashed", color = bright_red, linewidth = 0.8) +
  annotate("text", x = median_score + 0.2, y = 8,
           label = paste0("Median = ", round(median_score, 2)),
           size = 3.2, color = dark_blue, fontface = "plain", family = "serif", hjust = 0) +
  labs(x = "Lipid Score", y = "Count", title = "E") +
  theme_bw(base_size = 10) +
  theme(
    plot.title = element_text(hjust = 0.5, face = "plain", size = 11, color = dark_blue,
                              family = "serif", margin = margin(b = 2)),
    axis.text.x = element_text(size = 9, color = "#333333", family = "serif"),
    axis.text.y = element_text(size = 9, color = "#333333", family = "serif"),
    axis.title.x = element_text(size = 10, color = "#333333", family = "serif"),
    axis.title.y = element_text(size = 10, color = "#333333", family = "serif"),
    panel.grid.minor = element_blank(),
    panel.border = element_rect(color = "gray70", linewidth = 0.5)
  )

# ============================================================
# 8. 图 F: 评分 vs 年龄
# ============================================================

pF <- ggplot(tcga_data, aes(x = Age, y = tcga_scores$Lipid_Metabolism_Score)) +
  geom_point(size = 2.5, alpha = 0.7, color = lake_blue) +
  geom_smooth(method = "lm", se = TRUE, color = bright_red, fill = "#F0C0C0", alpha = 0.3, linewidth = 0.8) +
  stat_cor(method = "spearman", size = 3.5, color = dark_blue, family = "serif") +
  labs(x = "Age (years)", y = "Lipid Score", title = "F") +
  theme_bw(base_size = 10) +
  theme(
    plot.title = element_text(hjust = 0.5, face = "plain", size = 11, color = dark_blue,
                              family = "serif", margin = margin(b = 2)),
    axis.text.x = element_text(size = 9, color = "#333333", family = "serif"),
    axis.text.y = element_text(size = 9, color = "#333333", family = "serif"),
    axis.title.x = element_text(size = 10, color = "#333333", family = "serif"),
    axis.title.y = element_text(size = 10, color = "#333333", family = "serif"),
    panel.grid.minor = element_blank(),
    panel.border = element_rect(color = "gray70", linewidth = 0.5)
  )

# ============================================================
# 9. 图 G: 脂质基因热图（修复下标出界错误）
# ============================================================

# 读取真实表达矩阵
log2cpm_raw <- readRDS("data/raw/TCGA_LPS_real_log2cpm.rds")

# 基因名转换
ensembl_ids <- gsub("\\.[0-9]+$", "", rownames(log2cpm_raw))
gene_symbols <- mapIds(org.Hs.eg.db, 
                       keys = ensembl_ids, 
                       column = "SYMBOL",
                       keytype = "ENSEMBL", 
                       multiVals = "first")

keep <- !is.na(gene_symbols)
log2cpm_filtered <- log2cpm_raw[keep, ]
rownames(log2cpm_filtered) <- gene_symbols[keep]

# 合并重复基因
log2cpm <- aggregate(log2cpm_filtered, by = list(rownames(log2cpm_filtered)), FUN = mean)
rownames(log2cpm) <- log2cpm$Group.1
log2cpm$Group.1 <- NULL
log2cpm <- as.matrix(log2cpm)

# 提取6个脂质基因的真实表达
lipid_genes <- c("FASN", "ACLY", "SCD", "SREBF1", "HMGCR", "PPARG")
available_genes <- lipid_genes[lipid_genes %in% rownames(log2cpm)]
lipid_expr_mat <- log2cpm[available_genes, ]

# ===== 关键修复：确保样本顺序一致 =====
# 获取共同样本（确保 scores_valid 和 lipid_expr_mat 的样本一致）
common_samples <- intersect(colnames(lipid_expr_mat), rownames(scores_valid))
lipid_expr_mat <- lipid_expr_mat[, common_samples, drop = FALSE]

# 按评分排序样本（使用匹配后的样本顺序）
score_order <- order(scores_valid[common_samples, "Lipid_Metabolism_Score"])
lipid_expr_mat <- lipid_expr_mat[, score_order, drop = FALSE]

# 转换为长格式用于ggplot热图
lipid_melt <- melt(lipid_expr_mat)
colnames(lipid_melt) <- c("Gene", "Sample", "Expression")

# 按评分高低添加分组标签用于注释
score_group <- ifelse(scores_valid[common_samples, "Lipid_Metabolism_Score"][score_order] > median_score, "High", "Low")
sample_annotation <- data.frame(
  Sample = colnames(lipid_expr_mat),
  Group = score_group
)

pG <- ggplot(lipid_melt, aes(x = Sample, y = Gene, fill = Expression)) +
  geom_tile(color = "white", linewidth = 0.1) +
  scale_fill_gradient2(low = dark_blue, high = bright_red, mid = "white", midpoint = 0,
                       name = "Expression") +
  labs(x = "", y = "", title = "G") +
  theme_bw(base_size = 8) +
  theme(
    plot.title = element_text(hjust = 0.5, face = "plain", size = 11, color = dark_blue,
                              family = "serif", margin = margin(b = 2)),
    axis.text.x = element_blank(),
    axis.ticks.x = element_blank(),
    axis.text.y = element_text(size = 8, color = "#333333", family = "serif", face = "italic"),
    panel.grid.major = element_blank(),
    panel.grid.minor = element_blank(),
    panel.border = element_rect(color = "gray70", linewidth = 0.5),
    legend.position = "bottom",
    legend.title = element_text(size = 7, family = "serif"),
    legend.text = element_text(size = 6, family = "serif"),
    legend.key.width = unit(0.8, "cm"),
    legend.key.height = unit(0.3, "cm")
  )

# ============================================================
# 10. 图 H: 瀑布图
# ============================================================

score_sorted <- data.frame(
  Rank = 1:nrow(scores_valid),
  Score = sort(scores_valid$Lipid_Metabolism_Score, decreasing = TRUE)
)

pH <- ggplot(score_sorted, aes(x = Rank, y = Score)) +
  geom_line(color = indigo, linewidth = 0.8) +
  geom_area(fill = indigo, alpha = 0.3) +
  labs(x = "Sample Rank", y = "Lipid Score", title = "H") +
  theme_bw(base_size = 10) +
  theme(
    plot.title = element_text(hjust = 0.5, face = "plain", size = 11, color = dark_blue,
                              family = "serif", margin = margin(b = 2)),
    axis.text.x = element_text(size = 9, color = "#333333", family = "serif"),
    axis.text.y = element_text(size = 9, color = "#333333", family = "serif"),
    axis.title.x = element_text(size = 10, color = "#333333", family = "serif"),
    axis.title.y = element_text(size = 10, color = "#333333", family = "serif"),
    panel.grid.minor = element_blank(),
    panel.border = element_rect(color = "gray70", linewidth = 0.5)
  )

# ============================================================
# 11. 组合
# ============================================================

row1 <- (pA | pB | pC | pD) + plot_layout(widths = c(1, 1, 0.85, 0.85))
row2 <- (pE | pF | pG | pH) + plot_layout(widths = c(1, 1, 0.9, 1))

figure2 <- (row1 / row2) +
  plot_layout(heights = c(1, 1)) +
  plot_annotation(
    title = "Figure 2: Lipid Metabolism Landscape",
    theme = theme(
      plot.title = element_text(hjust = 0.5, face = "bold", size = 14, color = dark_blue,
                                family = "serif", margin = margin(b = 8))
    )
  )

# ============================================================
# 12. 输出
# ============================================================

ggsave("figures/Figure2.pdf", figure2, width = 16, height = 9, units = "in", dpi = 600, limitsize = FALSE, device = cairo_pdf)
ggsave("figures/Figure2.tiff", figure2, width = 16, height = 9, units = "in", dpi = 600, compression = "lzw")
ggsave("figures/Figure2.png", figure2, width = 16, height = 9, units = "in", dpi = 300)

cat("\n============================================================\n")
cat("✅ Figure 2 已保存（出版级）\n")
cat("   PDF:  figures/Figure2.pdf\n")
cat("   TIFF: figures/Figure2.tiff\n")
cat("   PNG:  figures/Figure2.png\n")
cat("============================================================\n")