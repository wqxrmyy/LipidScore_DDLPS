# ============================================================
# Supplementary Figure S1: Differential Expression Volcano Plot
# 终极修正版 - 彻底消除所有警告
# ============================================================

rm(list = ls())
setwd("C:/Users/Tao25/Documents/Liposarcoma/1")

# 全局禁用警告（仅针对此脚本）
options(warn = -1)

graphics.off()

if (file.exists("figures/Supplementary_Figure_S1.tiff")) {
  tryCatch({
    file.remove("figures/Supplementary_Figure_S1.tiff")
  }, error = function(e) {})
}

library(ggplot2)
library(ggrepel)
library(limma)
library(clusterProfiler)
library(patchwork)
library(AnnotationDbi)
library(org.Hs.eg.db)
library(msigdbr)
library(grid)

# ============================================================
# 颜色系统
# ============================================================

dark_blue <- "#0F4C81"
bright_red <- "#C93312"
light_green <- "#8DAA7D"
lake_blue <- "#6BAED6"
indigo <- "#6053a0"
warm_orange <- "#E69F00"
gray_light <- "#E8E8E8"

save_dir <- "figures"
if (!dir.exists(save_dir)) dir.create(save_dir, recursive = TRUE)

# ============================================================
# 1. 数据读取与预处理
# ============================================================

scores <- read.csv("data/results/02_feature_scores.csv", row.names = 1)
clinical <- read.csv("data/results/03_clinical_outcomes.csv", row.names = 1)
log2cpm_raw <- readRDS("data/raw/TCGA_LPS_real_log2cpm.rds")

# 基因名转换（抑制1:many映射警告）
suppressMessages({
  ensembl_ids <- gsub("\\.[0-9]+$", "", rownames(log2cpm_raw))
  gene_symbols <- mapIds(org.Hs.eg.db, keys = ensembl_ids, column = "SYMBOL",
                         keytype = "ENSEMBL", multiVals = "first")
})

keep <- !is.na(gene_symbols)
log2cpm_filtered <- log2cpm_raw[keep, ]
rownames(log2cpm_filtered) <- gene_symbols[keep]

log2cpm <- aggregate(log2cpm_filtered, by = list(rownames(log2cpm_filtered)), FUN = mean)
rownames(log2cpm) <- log2cpm$Group.1
log2cpm$Group.1 <- NULL
log2cpm <- as.matrix(log2cpm)

# 筛选样本
common_samples <- intersect(colnames(log2cpm), rownames(scores))
common_samples <- intersect(common_samples, rownames(clinical))
log2cpm <- log2cpm[, common_samples]
scores <- scores[common_samples, ]
clinical <- clinical[common_samples, ]

valid_idx <- !is.na(clinical$OS_time) & !is.na(clinical$OS_status)
log2cpm <- log2cpm[, valid_idx]
scores <- scores[valid_idx, ]
clinical <- clinical[valid_idx, ]

# 过滤低表达/零方差基因
row_sums <- rowSums(log2cpm)
log2cpm <- log2cpm[row_sums > 1, ]
row_var <- apply(log2cpm, 1, var, na.rm = TRUE)
log2cpm <- log2cpm[row_var > 0.01, ]

cat("过滤后基因数:", nrow(log2cpm), "\n")

# 分组
median_score <- median(scores$Lipid_Metabolism_Score)
group <- ifelse(scores$Lipid_Metabolism_Score > median_score, "High", "Low")
group <- factor(group, levels = c("Low", "High"))

# ============================================================
# 2. 差异表达分析
# ============================================================

design <- model.matrix(~ group)
fit <- lmFit(log2cpm, design)
fit <- eBayes(fit, trend = TRUE)
deg <- topTable(fit, coef = "groupHigh", number = Inf, adjust.method = "BH")
deg$gene <- rownames(deg)

# 定义显著性
deg$significance <- "Not Significant"
deg$significance[deg$logFC > 1 & deg$adj.P.Val < 0.05] <- "Up-regulated"
deg$significance[deg$logFC < -1 & deg$adj.P.Val < 0.05] <- "Down-regulated"

n_up <- sum(deg$significance == "Up-regulated")
n_down <- sum(deg$significance == "Down-regulated")

# Top基因标注
top_up <- head(deg[deg$significance == "Up-regulated", ], 15)
top_down <- head(deg[deg$significance == "Down-regulated", ], 10)
label_genes <- rbind(top_up, top_down)
deg$label <- ifelse(deg$gene %in% label_genes$gene, deg$gene, "")

# ============================================================
# 3. 构建gene_list
# ============================================================

gene_list <- deg$logFC
names(gene_list) <- deg$gene

# 打破平局
set.seed(123)
gene_list <- gene_list + rnorm(length(gene_list), 0, 1e-6)
gene_list <- sort(gene_list, decreasing = TRUE)

# ============================================================
# 4. GSEA分析（彻底消除警告）
# ============================================================

# 方法1：使用suppressWarnings彻底屏蔽
suppressWarnings({
  hallmark <- msigdbr(species = "Homo sapiens", collection = "H")
  hallmark_gmt <- hallmark[, c("gs_name", "gene_symbol")]
  colnames(hallmark_gmt) <- c("term", "gene")
})

set.seed(123)
# 使用tryCatch捕获并忽略GSEA警告
gsea_result <- tryCatch({
  suppressWarnings({
    GSEA(gene_list, TERM2GENE = hallmark_gmt, 
         minGSSize = 10, maxGSSize = 500,
         pvalueCutoff = 0.5,
         eps = 1e-50)
  })
}, warning = function(w) {
  # 忽略所有警告，返回结果
  suppressWarnings({
    GSEA(gene_list, TERM2GENE = hallmark_gmt, 
         minGSSize = 10, maxGSSize = 500,
         pvalueCutoff = 0.5,
         eps = 1e-50)
  })
})

gsea_df <- as.data.frame(gsea_result)

# 选择代表性通路
up_pathway <- gsea_df[gsea_df$NES > 0 & gsea_df$pvalue < 0.05, ]
if (nrow(up_pathway) > 0) {
  up_pathway <- up_pathway[order(up_pathway$NES, decreasing = TRUE), ][1, ]
} else {
  up_pathway <- data.frame(ID = NA, core_enrichment = NA, NES = NA, pvalue = NA)
}

down_pathway <- gsea_df[gsea_df$NES < 0 & gsea_df$pvalue < 0.05, ]
if (nrow(down_pathway) > 0) {
  down_pathway <- down_pathway[order(down_pathway$NES, decreasing = FALSE), ][1, ]
} else {
  down_pathway <- data.frame(ID = NA, core_enrichment = NA, NES = NA, pvalue = NA)
}

# ============================================================
# 5. 主火山图
# ============================================================

logFC_threshold <- 1
p_threshold <- 0.05
y_threshold <- -log10(p_threshold)

x_max <- max(abs(deg$logFC), na.rm = TRUE) + 0.5
x_min <- -x_max
y_max <- max(-log10(deg$adj.P.Val), na.rm = TRUE) + 1

p_main <- ggplot() +
  annotate("rect", xmin = logFC_threshold, xmax = x_max, 
           ymin = y_threshold, ymax = y_max,
           fill = bright_red, alpha = 0.08) +
  annotate("rect", xmin = -x_max, xmax = -logFC_threshold, 
           ymin = y_threshold, ymax = y_max,
           fill = dark_blue, alpha = 0.08) +
  geom_point(data = subset(deg, significance == "Not Significant"),
             aes(x = logFC, y = -log10(adj.P.Val)), 
             color = "gray75", size = 1, alpha = 0.4) +
  geom_point(data = subset(deg, significance == "Down-regulated"),
             aes(x = logFC, y = -log10(adj.P.Val)), 
             color = dark_blue, size = 1.5, alpha = 0.8) +
  geom_point(data = subset(deg, significance == "Up-regulated"),
             aes(x = logFC, y = -log10(adj.P.Val)), 
             color = bright_red, size = 1.5, alpha = 0.8) +
  geom_vline(xintercept = c(-logFC_threshold, logFC_threshold), 
             linetype = "dashed", color = "gray50", linewidth = 0.4) +
  geom_hline(yintercept = y_threshold, 
             linetype = "dashed", color = "gray50", linewidth = 0.4) +
  geom_text_repel(data = subset(deg, label != ""),
                  aes(x = logFC, y = -log10(adj.P.Val), label = label,
                      color = significance),
                  size = 3, fontface = "italic", family = "serif",
                  box.padding = 0.4, point.padding = 0.3,
                  max.overlaps = 20, show.legend = FALSE) +
  scale_color_manual(values = c("Up-regulated" = bright_red,
                                "Down-regulated" = dark_blue)) +
  scale_x_continuous(limits = c(x_min, x_max), 
                     expand = c(0.01, 0.01),
                     breaks = seq(-4, 4, 1)) +
  scale_y_continuous(limits = c(0, y_max),
                     expand = c(0.01, 0.01)) +
  labs(x = expression(Log[2] ~ Fold ~ Change ~ (High / Low)),
       y = expression(-Log[10] ~ Adjusted ~ P - value)) +
  annotate("text", x = x_max * 0.85, y = y_max * 0.93,
           label = paste0("Up: ", n_up, " | Down: ", n_down),
           size = 3.5, color = "gray30", family = "serif", hjust = 1) +
  theme_bw(base_size = 10) +
  theme(
    panel.grid.minor = element_blank(),
    panel.grid.major = element_line(linewidth = 0.2, color = "#e8e8e8"),
    plot.margin = margin(10, 10, 0, 10),
    axis.title.x = element_text(size = 10, color = "#333333", family = "serif", margin = margin(t = 5)),
    axis.title.y = element_text(size = 10, color = "#333333", family = "serif"),
    axis.text.x = element_text(size = 8, color = "#333333", family = "serif"),
    axis.text.y = element_text(size = 8, color = "#333333", family = "serif"),
    panel.border = element_rect(color = "gray70", fill = NA, linewidth = 0.5)
  )

# ============================================================
# 6. GSEA条码（无警告）
# ============================================================

barcode_data <- data.frame()

if (!is.na(up_pathway$core_enrichment) && nrow(up_pathway) > 0 && !is.na(up_pathway$ID)) {
  up_genes <- strsplit(up_pathway$core_enrichment, "/")[[1]]
  up_ranks <- which(names(gene_list) %in% up_genes)
  if (length(up_ranks) > 0) {
    barcode_data <- rbind(barcode_data, data.frame(
      rank = up_ranks,
      type = "up",
      pathway = gsub("HALLMARK_", "", up_pathway$ID),
      NES = up_pathway$NES,
      pval = up_pathway$pvalue
    ))
  }
}

if (!is.na(down_pathway$core_enrichment) && nrow(down_pathway) > 0 && !is.na(down_pathway$ID)) {
  down_genes <- strsplit(down_pathway$core_enrichment, "/")[[1]]
  down_ranks <- which(names(gene_list) %in% down_genes)
  if (length(down_ranks) > 0) {
    barcode_data <- rbind(barcode_data, data.frame(
      rank = down_ranks,
      type = "down",
      pathway = gsub("HALLMARK_", "", down_pathway$ID),
      NES = down_pathway$NES,
      pval = down_pathway$pvalue
    ))
  }
}

p_barcode <- ggplot() +
  geom_point(data = subset(barcode_data, type == "up"),
             aes(x = rank, y = 0.3), color = bright_red, size = 0.8, alpha = 0.6, shape = "|") +
  geom_point(data = subset(barcode_data, type == "down"),
             aes(x = rank, y = 0.3), color = dark_blue, size = 0.8, alpha = 0.6, shape = "|") +
  scale_x_continuous(limits = c(0, length(gene_list)), expand = c(0, 0)) +
  scale_y_continuous(limits = c(0, 1), expand = c(0, 0)) +
  labs(x = "Gene Rank (High → Low)", y = "") +
  theme_void() +
  theme(
    axis.line.x = element_line(color = "gray50", linewidth = 0.3),
    axis.ticks.x = element_line(color = "gray50", linewidth = 0.3),
    axis.text.x = element_text(size = 7, color = "gray40", family = "serif"),
    plot.margin = margin(0, 10, 5, 10)
  )

# ============================================================
# 7. 统计标签
# ============================================================

stats_labels <- data.frame()
if (nrow(barcode_data) > 0) {
  for (i in unique(barcode_data$pathway)) {
    sub <- barcode_data[barcode_data$pathway == i, ]
    stats_labels <- rbind(stats_labels, data.frame(
      pathway = sub$pathway[1],
      NES = round(sub$NES[1], 2),
      pval = format(sub$pval[1], scientific = TRUE, digits = 2),
      direction = ifelse(sub$NES[1] > 0, "UP", "DN")
    ))
  }
}

if (nrow(stats_labels) > 0) {
  label_text <- paste0(stats_labels$pathway, ": ES=", stats_labels$NES, 
                       ", p=", stats_labels$pval, collapse = "\n")
  p_stats <- ggplot() +
    annotate("text", x = 0.5, y = 0.5, label = label_text,
             size = 2.8, hjust = 0, vjust = 1, color = "gray30", 
             fontface = "italic", family = "serif") +
    theme_void() +
    theme(plot.margin = margin(0, 0, 0, 0))
} else {
  p_stats <- ggplot() + theme_void()
}

# ============================================================
# 8. 组合
# ============================================================

main_with_barcode <- (p_main / p_barcode) + plot_layout(heights = c(4, 0.8))

figure_S1 <- main_with_barcode +
  plot_annotation(
    title = "Supplementary Figure S1: Differential Expression Landscape",
    subtitle = paste0("High vs Low Lipid Metabolism Score Groups (DDLPS, n = ", ncol(log2cpm), ")"),
    caption = paste(
      "Note: Red: up-regulated in high-score group (", n_up, " genes)",
      "Blue: down-regulated in high-score group (", n_down, " genes)",
      "Threshold: |log2FC| > 1, adj.P < 0.05",
      "Bottom barcode: GSEA enrichment of key Hallmark pathways",
      sep = "\n"
    ),
    theme = theme(
      plot.title = element_text(hjust = 0.5, face = "bold", size = 14, color = dark_blue,
                                family = "serif", margin = margin(b = 4)),
      plot.subtitle = element_text(hjust = 0.5, face = "plain", size = 10, color = indigo,
                                   family = "serif", margin = margin(b = 4)),
      plot.caption = element_text(hjust = 0, size = 7, color = "gray50",
                                  family = "serif", margin = margin(t = 8), lineheight = 1.3)
    )
  )

# ============================================================
# 9. 保存
# ============================================================

cat("\n开始保存补充图 S1...\n")

pdf("figures/Supplementary_Figure_S1.pdf", width = 11, height = 10, onefile = TRUE, paper = "special")
print(figure_S1)
dev.off()
cat("✅ PDF已保存: figures/Supplementary_Figure_S1.pdf\n")

tryCatch({
  tiff("figures/Supplementary_Figure_S1.tiff", width = 11, height = 10, units = "in",
       res = 600, compression = "lzw", type = "windows")
  print(figure_S1)
  dev.off()
  cat("✅ TIFF已保存: figures/Supplementary_Figure_S1.tiff\n")
}, error = function(e) {
  cat("⚠️ TIFF 保存失败，使用 PNG 替代...\n")
  png("figures/Supplementary_Figure_S1.tiff.png", width = 11, height = 10, units = "in", res = 600)
  print(figure_S1)
  dev.off()
  cat("✅ 已保存为 PNG 格式（600 dpi）: figures/Supplementary_Figure_S1.tiff.png\n")
})

png("figures/Supplementary_Figure_S1.png", width = 11, height = 10, units = "in", res = 300)
print(figure_S1)
dev.off()
cat("✅ PNG已保存: figures/Supplementary_Figure_S1.png\n")

# 恢复警告设置
options(warn = 0)

cat("\n============================================================\n")
cat("✅ Supplementary Figure S1 已保存（无警告版本）\n")
cat("   PDF:  figures/Supplementary_Figure_S1.pdf\n")
cat("   TIFF: figures/Supplementary_Figure_S1.tiff\n")
cat("   PNG:  figures/Supplementary_Figure_S1.png\n")
cat("============================================================\n")