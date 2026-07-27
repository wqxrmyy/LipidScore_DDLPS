# ============================================================
# Figure 5: PPARG/CD36 Axis Validation and Clinical Relevance
# 出版级规范版本 - 与Figure 1-4风格统一
# 子图: A-H（8个子图）
# 布局: 2行 × 4列
# 字体: Times New Roman, 字号 9-11pt
# 输出: PDF（矢量）, TIFF（600 dpi LZW）, PNG（300 dpi）
# ============================================================

rm(list = ls())
setwd("C:/Users/Tao25/Documents/Liposarcoma/1")

# ============================================================
# 0. 清理图形设备和旧文件
# ============================================================

graphics.off()

if (file.exists("figures/Figure5.tiff")) {
  tryCatch({
    file.remove("figures/Figure5.tiff")
  }, error = function(e) {})
}

library(ggplot2)
library(pheatmap)
library(ggpubr)
library(survival)
library(survminer)
library(patchwork)
library(AnnotationDbi)
library(org.Hs.eg.db)
library(reshape2)
library(pROC)
library(png)
library(grid)

# ============================================================
# 1. 颜色系统（与Figure 1-4统一）
# ============================================================

dark_blue <- "#0F4C81"
bright_red <- "#C93312"
light_green <- "#8DAA7D"
lake_blue <- "#6BAED6"
indigo <- "#6053a0"
warm_orange <- "#E69F00"

# 确保输出目录存在
save_dir <- "figures"
if (!dir.exists(save_dir)) dir.create(save_dir, recursive = TRUE)

# ============================================================
# 2. 数据读取与预处理
# ============================================================

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

# 筛选共同样本
common <- intersect(colnames(log2cpm), rownames(scores))
common <- intersect(common, rownames(clinical))
log2cpm <- log2cpm[, common]
scores <- scores[common, ]
clinical <- clinical[common, ]

valid <- !is.na(clinical$OS_time) & !is.na(clinical$OS_status)
log2cpm <- log2cpm[, valid]
scores <- scores[valid, ]
clinical <- clinical[valid, ]

# 分组
median_score <- median(scores$Lipid_Metabolism_Score)
group <- ifelse(scores$Lipid_Metabolism_Score > median_score, "High", "Low")
group <- factor(group, levels = c("Low", "High"))

PPARG_expr <- as.numeric(log2cpm["PPARG", ])
CD36_expr <- as.numeric(log2cpm["CD36", ])

expr_data <- data.frame(
  Group = group,
  PPARG = PPARG_expr,
  CD36 = CD36_expr,
  Score = scores$Lipid_Metabolism_Score
)

# ============================================================
# 3. 图 A: PPARG 箱线图
# ============================================================

pA <- ggplot(expr_data, aes(x = Group, y = PPARG, fill = Group)) +
  geom_boxplot(alpha = 0.8, width = 0.6, outlier.size = 0.5) +
  geom_jitter(width = 0.15, size = 1, alpha = 0.5) +
  stat_compare_means(method = "wilcox.test", label = "p.signif", size = 4, label.x = 1.5) +
  scale_fill_manual(values = c(dark_blue, bright_red)) +
  labs(x = "", y = "PPARG Expression (log2CPM)", title = "A") +
  theme_bw(base_size = 10) +
  theme(
    legend.position = "none",
    plot.title = element_text(hjust = 0.5, face = "plain", size = 11, color = dark_blue,
                              family = "serif", margin = margin(b = 2)),
    axis.text.x = element_text(size = 9, color = "#333333", family = "serif"),
    axis.text.y = element_text(size = 9, color = "#333333", family = "serif"),
    axis.title.y = element_text(size = 10, color = "#333333", family = "serif"),
    panel.grid.major.x = element_blank(),
    panel.grid.minor = element_blank(),
    panel.border = element_rect(color = "gray70", fill = NA, linewidth = 0.5)
  )

# ============================================================
# 4. 图 B: CD36 箱线图
# ============================================================

pB <- ggplot(expr_data, aes(x = Group, y = CD36, fill = Group)) +
  geom_boxplot(alpha = 0.8, width = 0.6, outlier.size = 0.5) +
  geom_jitter(width = 0.15, size = 1, alpha = 0.5) +
  stat_compare_means(method = "wilcox.test", label = "p.signif", size = 4, label.x = 1.5) +
  scale_fill_manual(values = c(dark_blue, bright_red)) +
  labs(x = "", y = "CD36 Expression (log2CPM)", title = "B") +
  theme_bw(base_size = 10) +
  theme(
    legend.position = "none",
    plot.title = element_text(hjust = 0.5, face = "plain", size = 11, color = dark_blue,
                              family = "serif", margin = margin(b = 2)),
    axis.text.x = element_text(size = 9, color = "#333333", family = "serif"),
    axis.text.y = element_text(size = 9, color = "#333333", family = "serif"),
    axis.title.y = element_text(size = 10, color = "#333333", family = "serif"),
    panel.grid.major.x = element_blank(),
    panel.grid.minor = element_blank(),
    panel.border = element_rect(color = "gray70", fill = NA, linewidth = 0.5)
  )

# ============================================================
# 5. 图 C: PPARG vs Score 相关性
# ============================================================

pC <- ggplot(expr_data, aes(x = Score, y = PPARG)) +
  geom_point(size = 2.5, alpha = 0.7, color = lake_blue) +
  geom_smooth(method = "lm", se = TRUE, color = bright_red, linewidth = 0.8, 
              fill = light_green, alpha = 0.3) +
  stat_cor(method = "spearman", size = 3.5, color = dark_blue, family = "serif") +
  labs(x = "Lipid Metabolism Score", y = "PPARG Expression (log2CPM)", title = "C") +
  theme_bw(base_size = 10) +
  theme(
    plot.title = element_text(hjust = 0.5, face = "plain", size = 11, color = dark_blue,
                              family = "serif", margin = margin(b = 2)),
    axis.text.x = element_text(size = 9, color = "#333333", family = "serif"),
    axis.text.y = element_text(size = 9, color = "#333333", family = "serif"),
    axis.title.x = element_text(size = 10, color = "#333333", family = "serif"),
    axis.title.y = element_text(size = 10, color = "#333333", family = "serif"),
    panel.grid.minor = element_blank(),
    panel.border = element_rect(color = "gray70", fill = NA, linewidth = 0.5)
  )

# ============================================================
# 6. 图 D: CD36 vs Score 相关性
# ============================================================

pD <- ggplot(expr_data, aes(x = Score, y = CD36)) +
  geom_point(size = 2.5, alpha = 0.7, color = lake_blue) +
  geom_smooth(method = "lm", se = TRUE, color = bright_red, linewidth = 0.8,
              fill = light_green, alpha = 0.3) +
  stat_cor(method = "spearman", size = 3.5, color = dark_blue, family = "serif") +
  labs(x = "Lipid Metabolism Score", y = "CD36 Expression (log2CPM)", title = "D") +
  theme_bw(base_size = 10) +
  theme(
    plot.title = element_text(hjust = 0.5, face = "plain", size = 11, color = dark_blue,
                              family = "serif", margin = margin(b = 2)),
    axis.text.x = element_text(size = 9, color = "#333333", family = "serif"),
    axis.text.y = element_text(size = 9, color = "#333333", family = "serif"),
    axis.title.x = element_text(size = 10, color = "#333333", family = "serif"),
    axis.title.y = element_text(size = 10, color = "#333333", family = "serif"),
    panel.grid.minor = element_blank(),
    panel.border = element_rect(color = "gray70", fill = NA, linewidth = 0.5)
  )

# ============================================================
# 7. 图 E: PPARG 生存曲线
# ============================================================

median_pparg <- median(PPARG_expr)
pparg_group <- ifelse(PPARG_expr > median_pparg, "High", "Low")
surv_data <- data.frame(
  Time = clinical$OS_time / 365,
  Status = clinical$OS_status,
  Group = pparg_group
)
fit <- survfit(Surv(Time, Status) ~ Group, data = surv_data)

pE <- ggsurvplot(fit, data = surv_data, pval = TRUE, pval.method = TRUE,
                 conf.int = TRUE, risk.table = FALSE, palette = c(dark_blue, bright_red),
                 xlab = "Time (years)", ylab = "Overall Survival", title = "E",
                 legend.title = "PPARG", legend.labs = c("Low", "High"),
                 ggtheme = theme_bw(base_size = 10),
                 font.legend = 8, font.x = 10, font.y = 10)$plot +
  theme(
    plot.title = element_text(hjust = 0.5, face = "plain", size = 11, color = dark_blue,
                              family = "serif", margin = margin(b = 2)),
    axis.text.x = element_text(size = 9, color = "#333333", family = "serif"),
    axis.text.y = element_text(size = 9, color = "#333333", family = "serif"),
    panel.grid.minor = element_blank(),
    panel.border = element_rect(color = "gray70", fill = NA, linewidth = 0.5)
  )

# ============================================================
# 8. 图 F: 热图
# ============================================================

core_genes <- c("PPARG", "CD36", "FASN", "ACLY", "SCD", "SREBF1")
available_genes <- core_genes[core_genes %in% rownames(log2cpm)]
heat_mat <- log2cpm[available_genes, ]

# 按分组排序
order_idx <- order(group)
heat_mat <- heat_mat[, order_idx]
sample_groups <- group[order_idx]

annotation_col <- data.frame(Group = sample_groups)
rownames(annotation_col) <- colnames(heat_mat)
ann_colors <- list(Group = c(Low = dark_blue, High = bright_red))

# 生成热图并嵌入
temp_dir <- tempdir()
heatmap_file <- file.path(temp_dir, "heatmap.png")
png(heatmap_file, width = 600, height = 400, res = 100)
pheatmap(heat_mat, scale = "row",
         color = colorRampPalette(c(dark_blue, "white", bright_red))(100),
         annotation_col = annotation_col, annotation_colors = ann_colors,
         show_colnames = FALSE, fontsize = 8, border_color = NA,
         main = "F", treeheight_col = 10, treeheight_row = 10)
dev.off()

img_heat <- readPNG(heatmap_file)
pF <- ggplot() + 
  annotation_custom(rasterGrob(img_heat, interpolate = TRUE),
                    xmin = -Inf, xmax = Inf, ymin = -Inf, ymax = Inf) +
  coord_cartesian(xlim = c(0, 1), ylim = c(0, 1)) +
  theme_void() +
  labs(title = "F") +
  theme(plot.title = element_text(hjust = 0.5, face = "plain", size = 11,
                                  color = dark_blue, family = "serif", margin = margin(b = 2)))

# ============================================================
# 9. 图 G: 联合 ROC
# ============================================================

pparg_norm <- (PPARG_expr - min(PPARG_expr)) / (max(PPARG_expr) - min(PPARG_expr))
cd36_norm <- (CD36_expr - min(CD36_expr)) / (max(CD36_expr) - min(CD36_expr))
combined_score <- pparg_norm + cd36_norm

roc_obj <- roc(clinical$OS_status, combined_score)
roc_df <- data.frame(spec = 1 - roc_obj$specificities, sens = roc_obj$sensitivities)
auc_val <- round(auc(roc_obj), 3)

pG <- ggplot(roc_df, aes(x = spec, y = sens)) +
  geom_line(color = bright_red, linewidth = 1.2) +
  geom_abline(linetype = "dashed", color = "gray50", linewidth = 0.4) +
  labs(x = "1 - Specificity", y = "Sensitivity", title = "G") +
  annotate("text", x = 0.65, y = 0.20, label = paste("AUC =", auc_val),
           size = 3.5, color = dark_blue, family = "serif", hjust = 0) +
  theme_bw(base_size = 10) +
  theme(
    plot.title = element_text(hjust = 0.5, face = "plain", size = 11, color = dark_blue,
                              family = "serif", margin = margin(b = 2)),
    axis.text.x = element_text(size = 8, color = "#333333", family = "serif"),
    axis.text.y = element_text(size = 8, color = "#333333", family = "serif"),
    axis.title.x = element_text(size = 9, color = "#333333", family = "serif"),
    axis.title.y = element_text(size = 9, color = "#333333", family = "serif"),
    panel.grid.minor = element_blank(),
    aspect.ratio = 1,
    panel.border = element_rect(color = "gray70", fill = NA, linewidth = 0.5)
  )

# ============================================================
# 10. 图 H: PPARG vs CD36 相关性
# ============================================================

pair_df <- data.frame(PPARG = PPARG_expr, CD36 = CD36_expr, Group = group)

pH <- ggplot(pair_df, aes(x = PPARG, y = CD36, color = Group)) +
  geom_point(size = 2.5, alpha = 0.7) +
  geom_smooth(method = "lm", se = TRUE, color = "black", linewidth = 0.8, linetype = "dashed") +
  stat_cor(method = "spearman", size = 3.5, color = dark_blue,
           label.x.npc = 0.05, label.y.npc = 0.95, family = "serif") +
  scale_color_manual(values = c(dark_blue, bright_red)) +
  labs(x = "PPARG Expression (log2CPM)", y = "CD36 Expression (log2CPM)", title = "H") +
  theme_bw(base_size = 10) +
  theme(
    plot.title = element_text(hjust = 0.5, face = "plain", size = 11, color = dark_blue,
                              family = "serif", margin = margin(b = 2)),
    axis.text.x = element_text(size = 9, color = "#333333", family = "serif"),
    axis.text.y = element_text(size = 9, color = "#333333", family = "serif"),
    axis.title.x = element_text(size = 10, color = "#333333", family = "serif"),
    axis.title.y = element_text(size = 10, color = "#333333", family = "serif"),
    legend.position = "bottom",
    legend.title = element_text(size = 8, family = "serif"),
    legend.text = element_text(size = 7, family = "serif"),
    panel.grid.minor = element_blank(),
    panel.border = element_rect(color = "gray70", fill = NA, linewidth = 0.5)
  )

# ============================================================
# 11. 组合：2行 × 4列
# ============================================================

row1 <- (pA | pB | pC | pD) + plot_layout(widths = c(1, 1, 1, 1))
row2 <- (pE | pF | pG | pH) + plot_layout(widths = c(1, 1, 1, 1))

figure5 <- (row1 / row2) + 
  plot_layout(heights = c(1, 1.1)) +
  plot_annotation(
    title = "Figure 5: PPARG/CD36 Axis Validation and Clinical Relevance",
    subtitle = "A) PPARG expression | B) CD36 expression | C) PPARG vs Score | D) CD36 vs Score | E) PPARG survival | F) Core genes heatmap | G) Combined ROC | H) PPARG vs CD36 correlation",
    caption = paste(
      "Note: *p < 0.05, **p < 0.01, ***p < 0.001. Red: high-score group; Blue: low-score group.",
      "AUC: area under the curve.",
      sep = "\n"
    ),
    theme = theme(
      plot.title = element_text(hjust = 0.5, face = "bold", size = 14, color = dark_blue,
                                family = "serif", margin = margin(b = 6)),
      plot.subtitle = element_text(hjust = 0.5, face = "plain", size = 9, color = indigo,
                                   family = "serif", margin = margin(b = 4)),
      plot.caption = element_text(hjust = 0, size = 7, color = "gray50",
                                  family = "serif", margin = margin(t = 8), lineheight = 1.3)
    )
  )

# ============================================================
# 12. 出版级输出
# ============================================================

cat("\n开始保存图片...\n")

pdf("figures/Figure5.pdf", width = 16, height = 10, onefile = TRUE, paper = "special")
print(figure5)
dev.off()
cat("✅ PDF已保存: figures/Figure5.pdf\n")

tryCatch({
  tiff("figures/Figure5.tiff", width = 16, height = 10, units = "in", 
       res = 600, compression = "lzw", type = "windows")
  print(figure5)
  dev.off()
  cat("✅ TIFF已保存: figures/Figure5.tiff\n")
}, error = function(e) {
  cat("⚠️ TIFF 保存失败，使用 PNG 替代...\n")
  png("figures/Figure5.tiff.png", width = 16, height = 10, units = "in", res = 600)
  print(figure5)
  dev.off()
  cat("✅ 已保存为 PNG 格式（600 dpi）: figures/Figure5.tiff.png\n")
})

png("figures/Figure5.png", width = 16, height = 10, units = "in", res = 300)
print(figure5)
dev.off()
cat("✅ PNG已保存: figures/Figure5.png\n")

unlink(file.path(temp_dir, "*.png"))

cat("\n============================================================\n")
cat("✅ Figure 5 已保存（出版级）\n")
cat("   PDF:  figures/Figure5.pdf\n")
cat("   TIFF: figures/Figure5.tiff\n")
cat("   PNG:  figures/Figure5.png\n")
cat("\n布局: 2行 × 4列\n")
cat("  第一行: A | B | C | D\n")
cat("  第二行: E | F | G | H\n")
cat("\n字体规范: Times New Roman (serif)\n")
cat("字号规范: 标题 11pt, 轴标签 9-10pt\n")
cat("颜色规范: 与Figure 1-4统一\n")
cat("============================================================\n")