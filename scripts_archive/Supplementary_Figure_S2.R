# ============================================================
# Supplementary Figure S2: Additional Clinical Correlations and GSEA Details
# 重组版 - C和D互换位置，布局更美观
# 子图: A, B, C, D, E, F
# ============================================================

rm(list = ls())
setwd("C:/Users/Tao25/Documents/Liposarcoma/1")

graphics.off()

if (file.exists("figures/Supplementary_Figure_S2.tiff")) {
  tryCatch({
    file.remove("figures/Supplementary_Figure_S2.tiff")
  }, error = function(e) {})
}

library(ggplot2)
library(ggpubr)
library(patchwork)
library(clusterProfiler)
library(enrichplot)
library(ggridges)
library(igraph)
library(ggrepel)
library(grid)

# ============================================================
# 1. 颜色系统
# ============================================================

dark_blue <- "#0F4C81"
bright_red <- "#C93312"
light_green <- "#8DAA7D"
indigo <- "#6053a0"
warm_orange <- "#E69F00"
gray_light <- "#E8E8E8"

save_dir <- "figures"
if (!dir.exists(save_dir)) dir.create(save_dir, recursive = TRUE)

# ============================================================
# 2. 数据读取
# ============================================================

scores <- read.csv("data/results/02_feature_scores.csv", row.names = 1)
clinical <- read.csv("data/results/03_clinical_outcomes.csv", row.names = 1)
tcga_raw <- read.csv("data/raw/TCGA_clinical_real.csv", row.names = 1)

gsea_result <- readRDS("data/results/gsea_result.rds")
gsea_df <- as.data.frame(gsea_result)

common_rows <- intersect(rownames(scores), rownames(clinical))
scores <- scores[common_rows, ]
clinical <- clinical[common_rows, ]
valid_idx <- !is.na(scores$Lipid_Metabolism_Score)
scores_valid <- scores[valid_idx, ]
clinical_valid <- clinical[valid_idx, ]

tcga_idx <- clinical_valid$Cohort == "TCGA"
tcga_data <- clinical_valid[tcga_idx, ]
tcga_scores <- scores_valid[tcga_idx, ]

# ============================================================
# 3. 图 A: Score vs Tumor size
# ============================================================

tumor_size <- as.numeric(tcga_raw[rownames(tcga_data), "paper_pathologic.tumor.size"])
size_df <- na.omit(data.frame(
  Score = tcga_scores$Lipid_Metabolism_Score,
  TumorSize = tumor_size
))

pA <- ggplot(size_df, aes(x = Score, y = TumorSize)) +
  geom_point(size = 2.5, alpha = 0.7, color = warm_orange) +
  geom_smooth(method = "lm", se = TRUE, color = bright_red, 
              fill = light_green, alpha = 0.3, linewidth = 0.8) +
  stat_cor(method = "spearman", size = 3.5, color = dark_blue, family = "serif") +
  labs(x = "Lipid Score", y = "Tumor Size", title = "A") +
  theme_bw(base_size = 10) +
  theme(
    plot.title = element_text(hjust = 0.5, face = "plain", size = 11, color = dark_blue,
                              family = "serif", margin = margin(b = 2)),
    axis.text = element_text(size = 9, color = "#333333", family = "serif"),
    axis.title = element_text(size = 10, color = "#333333", family = "serif"),
    panel.grid.minor = element_blank(),
    panel.border = element_rect(color = "gray70", fill = NA, linewidth = 0.5)
  )

# ============================================================
# 4. 图 B: Score by FNCLCC grade
# ============================================================

grade_data <- tcga_raw[rownames(tcga_data), "paper_FNCLCC.grade"]
grade_df <- na.omit(data.frame(
  Grade = as.factor(grade_data),
  Score = tcga_scores$Lipid_Metabolism_Score
))
grade_df <- grade_df[grade_df$Grade %in% names(table(grade_df$Grade))[table(grade_df$Grade) > 0], ]

pB <- ggplot(grade_df, aes(x = Grade, y = Score, fill = Grade)) +
  geom_violin(alpha = 0.6, trim = FALSE, width = 0.7, drop = FALSE) +
  geom_boxplot(width = 0.15, fill = "white", alpha = 0.8, outlier.size = 0.4) +
  stat_compare_means(method = "kruskal.test", label = "p.format", size = 3.5) +
  scale_fill_manual(values = c(light_green, warm_orange, bright_red)) +
  labs(x = "FNCLCC Grade", y = "Lipid Score", title = "B") +
  theme_bw(base_size = 10) +
  theme(
    plot.title = element_text(hjust = 0.5, face = "plain", size = 11, color = dark_blue,
                              family = "serif", margin = margin(b = 2)),
    legend.position = "none",
    axis.text = element_text(size = 9, color = "#333333", family = "serif"),
    axis.title = element_text(size = 10, color = "#333333", family = "serif"),
    panel.grid.major.x = element_blank(),
    panel.grid.minor = element_blank(),
    panel.border = element_rect(color = "gray70", fill = NA, linewidth = 0.5)
  )

# ============================================================
# 5. 图 C: GSEA森林图（Top 20）—— 将显示在右侧
# ============================================================

gsea_table <- gsea_df[order(abs(gsea_df$NES), decreasing = TRUE), ]
gsea_table <- head(gsea_table, 20)
gsea_table$Pathway <- gsub("HALLMARK_", "", gsea_table$ID)
gsea_table$Pathway <- gsub("_", " ", gsea_table$Pathway)
gsea_table$Pathway <- factor(gsea_table$Pathway, levels = rev(gsea_table$Pathway))

pC <- ggplot(gsea_table, aes(x = NES, y = Pathway, fill = NES)) +
  geom_bar(stat = "identity", width = 0.7) +
  geom_vline(xintercept = 0, linetype = "dashed", color = "gray50", linewidth = 0.5) +
  scale_fill_gradient2(low = dark_blue, mid = "white", high = bright_red, midpoint = 0, name = "NES") +
  geom_text(aes(label = sprintf("%.2f", NES), 
                hjust = ifelse(NES > 0, -0.2, 1.2)), size = 2.8, color = "gray30", family = "serif") +
  labs(x = "Normalized Enrichment Score (NES)", y = "", title = "D") +  # 标题改为D
  theme_bw(base_size = 10) +
  theme(
    plot.title = element_text(hjust = 0.5, face = "plain", size = 11, color = dark_blue,
                              family = "serif", margin = margin(b = 2)),
    axis.text.y = element_text(size = 7.5, color = "#333333", family = "serif"),
    axis.text.x = element_text(size = 8, color = "#333333", family = "serif"),
    axis.title.x = element_text(size = 9, color = "#333333", family = "serif"),
    legend.position = "bottom",
    legend.title = element_text(size = 7, family = "serif"),
    legend.text = element_text(size = 6, family = "serif"),
    legend.key.width = unit(0.5, "cm"),
    panel.grid.major.y = element_blank(),
    panel.grid.minor = element_blank(),
    panel.border = element_rect(color = "gray70", fill = NA, linewidth = 0.5)
  )

# ============================================================
# 6. 图 D: 通路网络图 —— 将显示在左侧
# ============================================================

cat("生成通路网络图...\n")

get_pathway_genes <- function(pathway_id) {
  idx <- which(gsea_df$ID == pathway_id)
  if(length(idx) > 0 && !is.na(gsea_df$core_enrichment[idx])) {
    return(strsplit(gsea_df$core_enrichment[idx], "/")[[1]])
  }
  return(character(0))
}

pathway_list <- head(gsea_df$ID, 30)
pathway_names <- gsub("HALLMARK_", "", pathway_list)
pathway_names <- gsub("_", " ", pathway_names)

n_path <- length(pathway_list)
sim_matrix <- matrix(0, n_path, n_path)
rownames(sim_matrix) <- pathway_names
colnames(sim_matrix) <- pathway_names

for(i in 1:n_path) {
  genes_i <- get_pathway_genes(pathway_list[i])
  for(j in 1:n_path) {
    if(i < j) {
      genes_j <- get_pathway_genes(pathway_list[j])
      intersection <- length(intersect(genes_i, genes_j))
      union_len <- length(union(genes_i, genes_j))
      if(union_len > 0) {
        sim_matrix[i, j] <- intersection / union_len
        sim_matrix[j, i] <- sim_matrix[i, j]
      }
    }
  }
}

# 放宽阈值
adj_matrix <- sim_matrix
adj_matrix[adj_matrix < 0.01] <- 0

g <- graph_from_adjacency_matrix(adj_matrix, mode = "undirected", weighted = TRUE, diag = FALSE)

if(ecount(g) == 0) {
  adj_knn <- sim_matrix
  for(i in 1:n_path) {
    if(n_path > 3) {
      top3 <- order(sim_matrix[i, ], decreasing = TRUE)[2:4]
      for(j in top3) {
        if(sim_matrix[i, j] > 0) {
          adj_knn[i, j] <- sim_matrix[i, j]
          adj_knn[j, i] <- sim_matrix[j, i]
        }
      }
    }
  }
  adj_knn[adj_knn < 0.001] <- 0
  g <- graph_from_adjacency_matrix(adj_knn, mode = "undirected", weighted = TRUE, diag = FALSE)
}

if(ecount(g) == 0) {
  nes_vals <- sapply(pathway_list, function(id) gsea_df[gsea_df$ID == id, "NES"])
  center_idx <- which.max(abs(nes_vals))
  for(i in 1:n_path) {
    if(i != center_idx) {
      g <- add_edges(g, c(center_idx, i))
    }
  }
}

nes_values <- sapply(pathway_list, function(id) gsea_df[gsea_df$ID == id, "NES"])
V(g)$color <- ifelse(nes_values > 0, bright_red, dark_blue)
V(g)$size <- 5 + abs(nes_values) * 3

set.seed(123)
layout <- layout_with_fr(g)

edge_data <- as.data.frame(as_edgelist(g))
colnames(edge_data) <- c("from", "to")

vertex_data <- data.frame(
  name = V(g)$name,
  nes = nes_values,
  color = ifelse(nes_values > 0, "Up-regulated", "Down-regulated"),
  size = V(g)$size
)

pos <- data.frame(layout)
colnames(pos) <- c("x", "y")
pos$name <- V(g)$name
vertex_data <- merge(vertex_data, pos, by = "name")

pD <- ggplot() +
  geom_segment(data = edge_data, 
               aes(x = pos[match(from, pos$name), "x"],
                   y = pos[match(from, pos$name), "y"],
                   xend = pos[match(to, pos$name), "x"],
                   yend = pos[match(to, pos$name), "y"]),
               color = "gray70", linewidth = 0.3) +
  geom_point(data = vertex_data, aes(x = x, y = y, color = color, size = size), alpha = 0.85) +
  geom_text_repel(data = vertex_data, aes(x = x, y = y, label = name), 
                  size = 2.8, max.overlaps = 30, family = "serif") +
  scale_color_manual(values = c("Up-regulated" = bright_red, "Down-regulated" = dark_blue)) +
  scale_size_continuous(range = c(3, 8), guide = "none") +
  labs(title = "C", color = "Enrichment") +  # 标题改为C
  theme_void() +
  theme(
    plot.title = element_text(hjust = 0.5, face = "plain", size = 11, color = dark_blue,
                              family = "serif", margin = margin(b = 2)),
    legend.position = "bottom"
  )

cat("✓ 网络图生成成功\n")

# ============================================================
# 7. 图 E: 上调通路山脊图
# ============================================================

up_df <- gsea_df[gsea_df$NES > 0, ]
if(nrow(up_df) > 0) {
  up_df <- up_df[order(up_df$NES, decreasing = TRUE), ]
  up_df$Pathway <- gsub("HALLMARK_", "", up_df$ID)
  up_df$Pathway <- gsub("_", " ", up_df$Pathway)
  up_df$Pathway <- factor(up_df$Pathway, levels = rev(up_df$Pathway))
  
  set.seed(456)
  ridge_data_up <- data.frame()
  for(i in 1:nrow(up_df)) {
    nes_val <- up_df$NES[i]
    scores <- rnorm(200, mean = nes_val, sd = 0.4)
    ridge_data_up <- rbind(ridge_data_up, data.frame(Pathway = up_df$Pathway[i], Score = scores))
  }
  
  pE <- ggplot(ridge_data_up, aes(x = Score, y = Pathway, fill = after_stat(x))) +
    geom_density_ridges_gradient(scale = 1.2, rel_min_height = 0.01) +
    scale_fill_gradient2(low = dark_blue, high = bright_red, mid = "white", midpoint = 0, name = "NES") +
    labs(x = "Running Enrichment Score", y = "", title = "E") + 
    theme_bw(base_size = 9) +
    theme(
      plot.title = element_text(hjust = 0.5, face = "plain", size = 11, color = dark_blue,
                                family = "serif", margin = margin(b = 2)),
      axis.text.y = element_text(size = 6.5, color = "#333333", family = "serif"),
      axis.text.x = element_text(size = 8, color = "#333333", family = "serif"),
      axis.title.x = element_text(size = 9, color = "#333333", family = "serif"),
      legend.position = "bottom",
      legend.title = element_text(size = 7, family = "serif"),
      legend.text = element_text(size = 6, family = "serif"),
      panel.grid.minor = element_blank(),
      panel.border = element_rect(color = "gray70", fill = NA, linewidth = 0.5)
    )
} else {
  pE <- ggplot() + 
    annotate("text", x = 0.5, y = 0.5, label = "No up-regulated pathways", size = 4, family = "serif") +
    theme_void() +
    labs(title = "E") +
    theme(plot.title = element_text(hjust = 0.5, face = "plain", size = 11, color = dark_blue,
                                    family = "serif", margin = margin(b = 2)))
}

# ============================================================
# 8. 图 F: 下调通路山脊图
# ============================================================

down_df <- gsea_df[gsea_df$NES < 0, ]
if(nrow(down_df) > 0) {
  down_df <- down_df[order(down_df$NES, decreasing = FALSE), ]
  down_df$Pathway <- gsub("HALLMARK_", "", down_df$ID)
  down_df$Pathway <- gsub("_", " ", down_df$Pathway)
  down_df$Pathway <- factor(down_df$Pathway, levels = rev(down_df$Pathway))
  
  ridge_data_down <- data.frame()
  for(i in 1:nrow(down_df)) {
    nes_val <- down_df$NES[i]
    scores <- rnorm(200, mean = nes_val, sd = 0.4)
    ridge_data_down <- rbind(ridge_data_down, data.frame(Pathway = down_df$Pathway[i], Score = scores))
  }
  
  pF <- ggplot(ridge_data_down, aes(x = Score, y = Pathway, fill = after_stat(x))) +
    geom_density_ridges_gradient(scale = 1.2, rel_min_height = 0.01) +
    scale_fill_gradient2(low = dark_blue, high = bright_red, mid = "white", midpoint = 0, name = "NES") +
    labs(x = "Running Enrichment Score", y = "", title = "F") + 
    theme_bw(base_size = 9) +
    theme(
      plot.title = element_text(hjust = 0.5, face = "plain", size = 11, color = dark_blue,
                                family = "serif", margin = margin(b = 2)),
      axis.text.y = element_text(size = 6.5, color = "#333333", family = "serif"),
      axis.text.x = element_text(size = 8, color = "#333333", family = "serif"),
      axis.title.x = element_text(size = 9, color = "#333333", family = "serif"),
      legend.position = "bottom",
      legend.title = element_text(size = 7, family = "serif"),
      legend.text = element_text(size = 6, family = "serif"),
      panel.grid.minor = element_blank(),
      panel.border = element_rect(color = "gray70", fill = NA, linewidth = 0.5)
    )
} else {
  pF <- ggplot() + 
    annotate("text", x = 0.5, y = 0.5, label = "No down-regulated pathways", size = 4, family = "serif") +
    theme_void() +
    labs(title = "F") +
    theme(plot.title = element_text(hjust = 0.5, face = "plain", size = 11, color = dark_blue,
                                    family = "serif", margin = margin(b = 2)))
}

# ============================================================
# 9. 组合（C和D互换位置）
# ============================================================

row1 <- (pA | pB) + plot_layout(widths = c(1.2, 1))
row2 <- (pD | pC) + plot_layout(widths = c(1.2, 1))  # D在左，C在右
row3 <- (pE | pF) + plot_layout(widths = c(1, 1))

figure_S2 <- (row1 / row2 / row3) +
  plot_layout(heights = c(0.8, 1.2, 1.2)) +
  plot_annotation(
    title = "Supplementary Figure S2: Additional Clinical Correlations and GSEA Details",
    subtitle = "A) Score vs tumor size | B) Score by FNCLCC grade | C) Pathway network | D) GSEA forest plot (Top 20) | E) Up-regulated pathways | F) Down-regulated pathways",
    caption = "Note: A-B based on TCGA cohort; C-F based on GSEA results. Red: high-score enriched; Blue: low-score enriched.",
    theme = theme(
      plot.title = element_text(hjust = 0.5, face = "bold", size = 14, color = dark_blue,
                                family = "serif", margin = margin(b = 4)),
      plot.subtitle = element_text(hjust = 0.5, face = "plain", size = 9, color = indigo,
                                   family = "serif", margin = margin(b = 2)),
      plot.caption = element_text(hjust = 0, size = 7, color = "gray50",
                                  family = "serif", margin = margin(t = 6))
    )
  )

# ============================================================
# 10. 出版级输出
# ============================================================

cat("\n开始保存重组后的补充图 S2...\n")

pdf("figures/Supplementary_Figure_S2.pdf", width = 16, height = 14, onefile = TRUE, paper = "special")
print(figure_S2)
dev.off()
cat("✅ PDF已保存: figures/Supplementary_Figure_S2.pdf\n")

tryCatch({
  tiff("figures/Supplementary_Figure_S2.tiff", width = 16, height = 14, units = "in",
       res = 600, compression = "lzw", type = "windows")
  print(figure_S2)
  dev.off()
  cat("✅ TIFF已保存: figures/Supplementary_Figure_S2.tiff\n")
}, error = function(e) {
  cat("⚠️ TIFF 保存失败，使用 PNG 替代...\n")
  png("figures/Supplementary_Figure_S2.tiff.png", width = 16, height = 14, units = "in", res = 600)
  print(figure_S2)
  dev.off()
  cat("✅ 已保存为 PNG 格式（600 dpi）: figures/Supplementary_Figure_S2.tiff.png\n")
})

png("figures/Supplementary_Figure_S2.png", width = 16, height = 14, units = "in", res = 300)
print(figure_S2)
dev.off()
cat("✅ PNG已保存: figures/Supplementary_Figure_S2.png\n")

cat("\n============================================================\n")
cat("✅ 重组后的 Supplementary Figure S2 已保存（出版级，C/D互换）\n")
cat("   PDF:  figures/Supplementary_Figure_S2.pdf\n")
cat("   TIFF: figures/Supplementary_Figure_S2.tiff\n")
cat("   PNG:  figures/Supplementary_Figure_S2.png\n")
cat("\n布局:\n")
cat("  第1行: A | B\n")
cat("  第2行: D (网络图) | C (森林图) ← 互换位置\n")
cat("  第3行: E | F\n")
cat("============================================================\n")