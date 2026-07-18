# ============================================================
# Figure 6: Prognostic Performance and TF Activity
# 出版级规范版本 - 与Figure 1-5风格统一
# 子图: A-G（7个子图）
# 布局: 3行（A-C第一行，D-E第二行，F-G第三行）
# 字体: Times New Roman, 字号 9-11pt
# 输出: PDF（矢量）, TIFF（600 dpi LZW）, PNG（300 dpi）
# ============================================================

rm(list = ls())
setwd("C:/Users/Tao25/Documents/Liposarcoma/1")

# ============================================================
# 0. 清理图形设备和旧文件
# ============================================================

graphics.off()

if (file.exists("figures/Figure6.tiff")) {
  tryCatch({
    file.remove("figures/Figure6.tiff")
  }, error = function(e) {})
}

library(ggplot2)
library(patchwork)
library(pROC)
library(ggpubr)
library(reshape2)
library(AnnotationDbi)
library(org.Hs.eg.db)
library(grid)

# ============================================================
# 1. 颜色系统（与Figure 1-5统一）
# ============================================================

dark_blue <- "#0F4C81"
bright_red <- "#C93312"
light_green <- "#8DAA7D"
lake_blue <- "#6BAED6"
indigo <- "#6053a0"
warm_orange <- "#E69F00"
gray_light <- "#E8E8E8"

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
common_samples <- intersect(colnames(log2cpm), rownames(scores))
common_samples <- intersect(common_samples, rownames(clinical))
log2cpm <- log2cpm[, common_samples]
scores <- scores[common_samples, ]
clinical <- clinical[common_samples, ]

# 去除生存数据缺失
valid_idx <- !is.na(clinical$OS_time) & !is.na(clinical$OS_status)
log2cpm <- log2cpm[, valid_idx]
scores <- scores[valid_idx, ]
clinical <- clinical[valid_idx, ]

# 分组
median_score <- median(scores$Lipid_Metabolism_Score)
group <- ifelse(scores$Lipid_Metabolism_Score > median_score, "High", "Low")
group_factor <- factor(group, levels = c("Low", "High"))

# 生存数据（年）
time <- clinical$OS_time / 365
status <- clinical$OS_status

# ============================================================
# 3. 图 A: ROC 曲线 (1/3/5年)
# ============================================================

years <- c(1, 3, 5)
roc_list <- list()
auc_values <- c()

for (i in seq_along(years)) {
  yr <- years[i]
  status_yr <- ifelse(time <= yr & status == 1, 1, ifelse(time > yr, 0, NA))
  valid_yr <- !is.na(status_yr)
  
  if (sum(valid_yr) > 10) {
    roc_obj <- roc(status_yr[valid_yr], scores$Lipid_Metabolism_Score[valid_yr], quiet = TRUE)
    auc_val <- round(auc(roc_obj), 3)
    auc_values <- c(auc_values, auc_val)
    roc_list[[i]] <- data.frame(spec = 1 - roc_obj$specificities, sens = roc_obj$sensitivities, year = yr, auc = auc_val)
  }
}

pA <- ggplot() +
  geom_abline(intercept = 0, slope = 1, linetype = "dashed", color = "gray50", linewidth = 0.4) +
  scale_x_reverse(expand = c(0, 0), limits = c(1, 0)) +
  scale_y_continuous(expand = c(0, 0), limits = c(0, 1)) +
  xlab("1 - Specificity") + ylab("Sensitivity") +
  theme_bw(base_size = 10) +
  theme(
    plot.title = element_text(hjust = 0.5, face = "plain", size = 11, color = dark_blue,
                              family = "serif", margin = margin(b = 2)),
    panel.grid.minor = element_blank(),
    axis.text.x = element_text(size = 8, color = "#333333", family = "serif"),
    axis.text.y = element_text(size = 8, color = "#333333", family = "serif"),
    axis.title.x = element_text(size = 9, color = "#333333", family = "serif"),
    axis.title.y = element_text(size = 9, color = "#333333", family = "serif"),
    aspect.ratio = 1,
    panel.border = element_rect(color = "gray70", fill = NA, linewidth = 0.5)
  )

colors_roc <- c(bright_red, lake_blue, warm_orange)
for (i in seq_along(roc_list)) {
  pA <- pA + geom_line(data = roc_list[[i]], aes(x = spec, y = sens, color = factor(year)), linewidth = 1)
}

pA <- pA + 
  scale_color_manual(
    values = colors_roc[1:length(roc_list)],
    name = "Year",
    labels = paste0(years[1:length(roc_list)], "-year (AUC=", auc_values, ")")
  ) +
  labs(title = "A") +
  theme(
    legend.position = "bottom",
    legend.title = element_text(size = 8, family = "serif"),
    legend.text = element_text(size = 7, family = "serif")
  )

# ============================================================
# 4. 图 B: 时间依赖性 AUC
# ============================================================

time_points <- seq(0.5, 8, by = 0.5)
auc_time <- data.frame()

for (t in time_points) {
  if (t <= max(time)) {
    status_t <- ifelse(time <= t & status == 1, 1, 0)
    if (sum(status_t) > 2 && sum(status_t) < length(status_t)) {
      roc_t <- tryCatch({
        roc(status_t, scores$Lipid_Metabolism_Score, quiet = TRUE)
      }, error = function(e) NULL)
      if (!is.null(roc_t)) {
        auc_time <- rbind(auc_time, data.frame(Time = t, AUC = as.numeric(auc(roc_t))))
      }
    }
  }
}

pB <- ggplot(auc_time, aes(x = Time, y = AUC)) +
  geom_line(color = bright_red, linewidth = 1.2) +
  geom_point(color = dark_blue, size = 2.5) +
  geom_hline(yintercept = 0.5, linetype = "dashed", color = "gray50", linewidth = 0.4) +
  ylim(0, 1) + xlim(0, 8) +
  labs(x = "Time (years)", y = "AUC", title = "B") +
  theme_bw(base_size = 10) +
  theme(
    plot.title = element_text(hjust = 0.5, face = "plain", size = 11, color = dark_blue,
                              family = "serif", margin = margin(b = 2)),
    axis.text.x = element_text(size = 8, color = "#333333", family = "serif"),
    axis.text.y = element_text(size = 8, color = "#333333", family = "serif"),
    axis.title.x = element_text(size = 9, color = "#333333", family = "serif"),
    axis.title.y = element_text(size = 9, color = "#333333", family = "serif"),
    panel.grid.minor = element_blank(),
    panel.border = element_rect(color = "gray70", fill = NA, linewidth = 0.5)
  )

# ============================================================
# 5. 图 C: 决策曲线分析 (DCA)
# ============================================================

calculate_net_benefit <- function(predictions, outcomes) {
  thresholds <- seq(0.05, 0.95, by = 0.05)
  net_benefit <- numeric(length(thresholds))
  prob <- scale(predictions)
  prob <- (prob - min(prob)) / (max(prob) - min(prob))
  for (i in seq_along(thresholds)) {
    pt <- thresholds[i]
    high_risk <- prob > pt
    TP <- sum(high_risk & outcomes == 1)
    FP <- sum(high_risk & outcomes == 0)
    n <- length(outcomes)
    net_benefit[i] <- (TP / n) - (FP / n) * (pt / (1 - pt))
  }
  return(data.frame(Threshold = thresholds, NetBenefit = net_benefit))
}

status_3yr <- ifelse(time <= 3 & status == 1, 1, 0)
status_3yr[time > 3] <- 0
valid_3yr <- !is.na(status_3yr)
dca_result <- calculate_net_benefit(scores$Lipid_Metabolism_Score[valid_3yr], status_3yr[valid_3yr])

pC <- ggplot(dca_result, aes(x = Threshold, y = NetBenefit)) +
  geom_line(color = bright_red, linewidth = 1.2) +
  geom_hline(yintercept = 0, linetype = "dashed", color = "gray50", linewidth = 0.4) +
  labs(x = "Threshold Probability", y = "Net Benefit", title = "C") +
  theme_bw(base_size = 10) +
  theme(
    plot.title = element_text(hjust = 0.5, face = "plain", size = 11, color = dark_blue,
                              family = "serif", margin = margin(b = 2)),
    axis.text.x = element_text(size = 8, color = "#333333", family = "serif"),
    axis.text.y = element_text(size = 8, color = "#333333", family = "serif"),
    axis.title.x = element_text(size = 9, color = "#333333", family = "serif"),
    axis.title.y = element_text(size = 9, color = "#333333", family = "serif"),
    panel.grid.minor = element_blank(),
    panel.border = element_rect(color = "gray70", fill = NA, linewidth = 0.5)
  )

# ============================================================
# 6. 图 D: 核心基因表达箱线图
# ============================================================

core_genes <- c("FASN", "ACLY", "SCD", "SREBF1")
available_genes <- core_genes[core_genes %in% rownames(log2cpm)]

expr_data <- data.frame(t(log2cpm[available_genes, ]))
expr_data$Group <- group_factor
expr_long <- melt(expr_data, id.vars = "Group", variable.name = "Gene", value.name = "Expression")

pD <- ggplot(expr_long, aes(x = Group, y = Expression, fill = Group)) +
  geom_boxplot(alpha = 0.8, width = 0.6, outlier.size = 0.5) +
  geom_jitter(width = 0.15, size = 0.8, alpha = 0.5) +
  stat_compare_means(method = "wilcox.test", label = "p.signif", size = 3.5, label.x = 1.5) +
  facet_wrap(~Gene, scales = "free_y", nrow = 2) +
  scale_fill_manual(values = c(dark_blue, bright_red)) +
  labs(x = "", y = "log2CPM", title = "D") +
  theme_bw(base_size = 9) +
  theme(
    legend.position = "none",
    plot.title = element_text(hjust = 0.5, face = "plain", size = 11, color = dark_blue,
                              family = "serif", margin = margin(b = 2)),
    strip.text = element_text(face = "plain", size = 9, family = "serif", color = dark_blue),
    axis.text.x = element_text(angle = 45, hjust = 1, size = 8, color = "#333333", family = "serif"),
    axis.text.y = element_text(size = 8, color = "#333333", family = "serif"),
    axis.title.y = element_text(size = 9, color = "#333333", family = "serif"),
    panel.grid.major.x = element_blank(),
    panel.grid.minor = element_blank(),
    panel.border = element_rect(color = "gray70", fill = NA, linewidth = 0.5)
  )

# ============================================================
# 7. 图 E: 基因相关性网络（热图）
# ============================================================

network_genes <- c("FASN", "ACLY", "SCD", "SREBF1", "HMGCR", "PPARG")
available_network <- network_genes[network_genes %in% rownames(log2cpm)]

if (length(available_network) >= 3) {
  expr_network <- t(log2cpm[available_network, ])
  cor_matrix <- cor(expr_network, method = "spearman")
  cor_melt <- melt(cor_matrix)
  colnames(cor_melt) <- c("Gene1", "Gene2", "Correlation")
  
  pE <- ggplot(cor_melt, aes(x = Gene1, y = Gene2, fill = Correlation)) +
    geom_tile(color = "white", linewidth = 0.2) +
    scale_fill_gradient2(low = dark_blue, high = bright_red, mid = "white",
                         midpoint = 0, limit = c(-1, 1), name = "Spearman\ncorrelation") +
    geom_text(aes(label = round(Correlation, 2)), size = 2.5, color = "gray30", family = "serif") +
    labs(x = "", y = "", title = "E") +
    theme_bw(base_size = 9) +
    theme(
      plot.title = element_text(hjust = 0.5, face = "plain", size = 11, color = dark_blue,
                                family = "serif", margin = margin(b = 2)),
      axis.text.x = element_text(angle = 45, hjust = 1, size = 8, color = "#333333", family = "serif", face = "italic"),
      axis.text.y = element_text(size = 8, color = "#333333", family = "serif", face = "italic"),
      legend.position = "bottom",
      legend.title = element_text(size = 7, family = "serif"),
      legend.text = element_text(size = 6, family = "serif"),
      legend.key.width = unit(0.5, "cm"),
      legend.key.height = unit(0.3, "cm"),
      panel.grid.major = element_blank(),
      panel.grid.minor = element_blank(),
      panel.border = element_rect(color = "gray70", fill = NA, linewidth = 0.5)
    )
} else {
  pE <- ggplot() + 
    annotate("text", x = 0.5, y = 0.5, label = "Genes not available", size = 4, family = "serif") +
    theme_void() +
    labs(title = "E") +
    theme(plot.title = element_text(hjust = 0.5, face = "plain", size = 11, color = dark_blue,
                                    family = "serif", margin = margin(b = 2)))
}

# ============================================================
# 8. 图 F: 转录因子活性条形图
# ============================================================

tf_targets <- list(
  "PPARG" = c("CD36", "LPL", "FABP4", "ADIPOQ", "PLIN1", "FASN"),
  "SREBF1" = c("FASN", "ACLY", "SCD", "ACACA", "SREBF1"),
  "E2F1" = c("CCNA2", "CCNE1", "CCNB1", "CDK1", "PCNA", "MCM2"),
  "MYC" = c("CDK4", "CCND2", "E2F2", "G6PD", "LDHA", "PKM2")
)

tf_scores <- matrix(NA, nrow = ncol(log2cpm), ncol = length(tf_targets))
colnames(tf_scores) <- names(tf_targets)
rownames(tf_scores) <- colnames(log2cpm)

for (i in seq_along(tf_targets)) {
  genes <- tf_targets[[i]]
  available_genes <- genes[genes %in% rownames(log2cpm)]
  if (length(available_genes) >= 2) {
    tf_scores[, i] <- colMeans(log2cpm[available_genes, , drop = FALSE], na.rm = TRUE)
  } else {
    tf_scores[, i] <- NA
  }
}

high_idx <- group_factor == "High"
low_idx <- group_factor == "Low"

tf_logfc <- colMeans(tf_scores[high_idx, , drop = FALSE], na.rm = TRUE) -
  colMeans(tf_scores[low_idx, , drop = FALSE], na.rm = TRUE)

tf_df <- data.frame(TF = names(tf_targets), LogFC = tf_logfc)
tf_df$Direction <- ifelse(tf_df$LogFC > 0, "High-score activated", "High-score suppressed")
tf_df$LogFC_rounded <- round(tf_df$LogFC, 2)

pF <- ggplot(tf_df, aes(x = reorder(TF, LogFC), y = LogFC, fill = Direction)) +
  geom_bar(stat = "identity", width = 0.7) +
  coord_flip() +
  scale_fill_manual(values = c("High-score activated" = bright_red,
                               "High-score suppressed" = dark_blue)) +
  geom_text(aes(label = LogFC_rounded),
            hjust = ifelse(tf_df$LogFC > 0, -0.2, 1.2),
            size = 3.5, family = "serif", color = "gray30") +
  labs(x = "", y = "Log2 Fold Change (High/Low)", title = "F") +
  theme_bw(base_size = 10) +
  theme(
    plot.title = element_text(hjust = 0.5, face = "plain", size = 11, color = dark_blue,
                              family = "serif", margin = margin(b = 2)),
    axis.text.y = element_text(size = 9, color = "#333333", family = "serif"),
    axis.text.x = element_text(size = 8, color = "#333333", family = "serif"),
    axis.title.x = element_text(size = 9, color = "#333333", family = "serif"),
    legend.position = "bottom",
    legend.title = element_blank(),
    legend.text = element_text(size = 8, family = "serif"),
    panel.grid.major.y = element_blank(),
    panel.grid.minor = element_blank(),
    panel.border = element_rect(color = "gray70", fill = NA, linewidth = 0.5)
  )

# ============================================================
# 9. 图 G: 免疫评分 vs 脂质评分相关性
# ============================================================

immune_genes <- list(
  "CD8+ T" = c("CD8A", "CD8B"),
  "CD4+ T" = c("CD4", "CD40LG"),
  "Treg" = c("FOXP3", "IL2RA"),
  "B" = c("CD19", "MS4A1"),
  "NK" = c("NCAM1", "KLRK1"),
  "Macrophage" = c("CD68", "CD163")
)

immune_scores <- matrix(NA, nrow = ncol(log2cpm), ncol = length(immune_genes))
colnames(immune_scores) <- names(immune_genes)
rownames(immune_scores) <- colnames(log2cpm)

for (i in seq_along(immune_genes)) {
  genes <- immune_genes[[i]]
  available_genes <- genes[genes %in% rownames(log2cpm)]
  if (length(available_genes) >= 1) {
    immune_scores[, i] <- colMeans(log2cpm[available_genes, , drop = FALSE], na.rm = TRUE)
  }
}

cor_immune <- sapply(1:ncol(immune_scores), function(i) {
  cor(immune_scores[, i], scores$Lipid_Metabolism_Score, method = "spearman", use = "complete.obs")
})
cor_df <- data.frame(Immune_Type = names(immune_genes), Correlation = cor_immune)

pG <- ggplot(cor_df, aes(x = reorder(Immune_Type, Correlation), y = Correlation, fill = Correlation > 0)) +
  geom_bar(stat = "identity", width = 0.7) +
  coord_flip() +
  scale_fill_manual(values = c("TRUE" = bright_red, "FALSE" = dark_blue)) +
  geom_text(aes(label = round(Correlation, 2)),
            hjust = ifelse(cor_df$Correlation > 0, -0.2, 1.2),
            size = 3.5, family = "serif", color = "gray30") +
  labs(x = "", y = "Spearman Correlation with Lipid Score", title = "G") +
  theme_bw(base_size = 10) +
  theme(
    plot.title = element_text(hjust = 0.5, face = "plain", size = 11, color = dark_blue,
                              family = "serif", margin = margin(b = 2)),
    axis.text.y = element_text(size = 8, color = "#333333", family = "serif"),
    axis.text.x = element_text(size = 8, color = "#333333", family = "serif"),
    axis.title.x = element_text(size = 9, color = "#333333", family = "serif"),
    legend.position = "none",
    panel.grid.major.y = element_blank(),
    panel.grid.minor = element_blank(),
    panel.border = element_rect(color = "gray70", fill = NA, linewidth = 0.5)
  )

# ============================================================
# 10. 组合
# ============================================================

row1 <- (pA | pB | pC) + plot_layout(widths = c(1.2, 1, 1))
row2 <- (pD | pE) + plot_layout(widths = c(1.2, 1))
row3 <- (pF | pG) + plot_layout(widths = c(1, 0.9))

figure6 <- (row1 / row2 / row3) +
  plot_layout(heights = c(1, 1, 0.8)) +
  plot_annotation(
    title = "Figure 6: Prognostic Performance and TF Activity",
    subtitle = "A) ROC curves | B) Time-dependent AUC | C) Decision curve | D) Core gene expression | E) Gene correlation network | F) TF activity | G) Immune correlation",
    caption = paste(
      "Note: *p < 0.05, **p < 0.01, ***p < 0.001.",
      "TF: transcription factor; AUC: area under the curve.",
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
# 11. 出版级输出
# ============================================================

cat("\n开始保存图片...\n")

pdf("figures/Figure6.pdf", width = 15, height = 13, onefile = TRUE, paper = "special")
print(figure6)
dev.off()
cat("✅ PDF已保存: figures/Figure6.pdf\n")

tryCatch({
  tiff("figures/Figure6.tiff", width = 15, height = 13, units = "in",
       res = 600, compression = "lzw", type = "windows")
  print(figure6)
  dev.off()
  cat("✅ TIFF已保存: figures/Figure6.tiff\n")
}, error = function(e) {
  cat("⚠️ TIFF 保存失败，使用 PNG 替代...\n")
  png("figures/Figure6.tiff.png", width = 15, height = 13, units = "in", res = 600)
  print(figure6)
  dev.off()
  cat("✅ 已保存为 PNG 格式（600 dpi）: figures/Figure6.tiff.png\n")
})

png("figures/Figure6.png", width = 15, height = 13, units = "in", res = 300)
print(figure6)
dev.off()
cat("✅ PNG已保存: figures/Figure6.png\n")

cat("\n============================================================\n")
cat("✅ Figure 6 已保存（出版级）\n")
cat("   PDF:  figures/Figure6.pdf\n")
cat("   TIFF: figures/Figure6.tiff\n")
cat("   PNG:  figures/Figure6.png\n")
cat("\n布局: 3行\n")
cat("  第一行: A | B | C (ROC, AUC, DCA)\n")
cat("  第二行: D | E (基因表达, 相关性网络)\n")
cat("  第三行: F | G (TF活性, 免疫相关性)\n")
cat("\n字体规范: Times New Roman (serif)\n")
cat("字号规范: 标题 11pt, 轴标签 8-9pt\n")
cat("颜色规范: 与Figure 1-5统一\n")
cat("============================================================\n")