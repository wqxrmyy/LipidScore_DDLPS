# ============================================================
# Figure 4: GSEA Dual Phenotype
# 出版级规范版本 - 修复 TIFF 文件无法写入问题
# 子图: A-H（8个子图），均含字母编号
# 字体: Times New Roman, 字号 11pt
# 输出: PDF（矢量）, TIFF（600 dpi LZW）, PNG（300 dpi）
# ============================================================

rm(list = ls())
setwd("C:/Users/Tao25/Documents/Liposarcoma/1")

# ============================================================
# 0. 清理图形设备和旧文件（修复 TIFF 错误）
# ============================================================

graphics.off()  # 关闭所有图形设备

# 删除可能被占用的旧文件
if (file.exists("figures/Figure4.tiff")) {
  tryCatch({
    file.remove("figures/Figure4.tiff")
    cat("✅ 已删除旧的 Figure4.tiff\n")
  }, error = function(e) {
    cat("⚠️ 无法删除旧文件，请手动关闭后重试\n")
  })
}

library(ggplot2)
library(clusterProfiler)
library(enrichplot)
library(patchwork)
library(grid)
library(png)

# ============================================================
# 1. 颜色系统
# ============================================================

dark_blue <- "#0F4C81"
bright_red <- "#C93312"
light_green <- "#8DAA7D"
warm_orange <- "#E69F00"
indigo <- "#6053a0"
gray_light <- "#E8E8E8"

# 确保输出目录存在
if (!dir.exists("figures")) {
  dir.create("figures", recursive = TRUE)
}
if (!dir.exists("figures/Figure4_temp")) {
  dir.create("figures/Figure4_temp", recursive = TRUE)
}

# ============================================================
# 2. 数据读取
# ============================================================

gsea_result <- readRDS("data/results/gsea_result.rds")
gsea_df <- as.data.frame(gsea_result)

cat("GSEA result loaded, pathways:", nrow(gsea_df), "\n")

# ============================================================
# 3. 简化通路名称
# ============================================================

shorten_pathway <- function(name) {
  name <- gsub("HALLMARK_", "", name)
  name <- gsub("_", " ", name)
  name <- gsub("INTERFERON GAMMA RESPONSE", "IFN-gamma", name)
  name <- gsub("INTERFERON ALPHA RESPONSE", "IFN-alpha", name)
  name <- gsub("EPITHELIAL MESENCHYMAL TRANSITION", "EMT", name)
  name <- gsub("INFLAMMATORY RESPONSE", "Inflammation", name)
  name <- gsub("MITOTIC SPINDLE", "Mitotic Spindle", name)
  name <- gsub("G2M CHECKPOINT", "G2M Checkpoint", name)
  name <- gsub("E2F TARGETS", "E2F Targets", name)
  name <- gsub("TNFA SIGNALING VIA NFKB", "TNFa/NF-kB", name)
  name <- gsub("IL6 JAK STAT3 SIGNALING", "IL-6/JAK/STAT3", name)
  name <- gsub("ALLOGRAFT REJECTION", "Allograft Rejection", name)
  name <- gsub("OXIDATIVE PHOSPHORYLATION", "OxPhos", name)
  name <- gsub("MTORC1 SIGNALING", "mTORC1", name)
  return(name)
}

# ============================================================
# 4. 图 A: NES条形图（Top 15）
# ============================================================

plot_df <- gsea_df[order(abs(gsea_df$NES), decreasing = TRUE), ]
plot_df <- head(plot_df, 15)
plot_df$Pathway <- shorten_pathway(plot_df$ID)
plot_df$Pathway <- factor(plot_df$Pathway, levels = rev(plot_df$Pathway))

pA <- ggplot(plot_df, aes(x = NES, y = Pathway, fill = NES > 0)) +
  geom_bar(stat = "identity", width = 0.7) +
  scale_fill_manual(values = c(dark_blue, bright_red),
                    labels = c("Low Score", "High Score")) +
  labs(x = "Normalized Enrichment Score (NES)", y = "", title = "A") +
  theme_bw(base_size = 10) +
  theme(
    plot.title = element_text(hjust = 0.5, face = "plain", size = 11, color = dark_blue,
                              family = "serif", margin = margin(b = 2)),
    legend.position = "bottom",
    legend.text = element_text(size = 8, family = "serif"),
    axis.text.y = element_text(size = 8, color = "#333333", family = "serif"),
    axis.text.x = element_text(size = 9, color = "#333333", family = "serif"),
    axis.title.x = element_text(size = 10, color = "#333333", family = "serif"),
    panel.grid.major.y = element_blank(),
    panel.grid.minor = element_blank(),
    panel.border = element_rect(color = "gray70", fill = NA, linewidth = 0.5)
  )

# ============================================================
# 5. 图 B-E: GSEA曲线（使用 labs(title) 统一大小）
# ============================================================

top_up <- head(gsea_df[gsea_df$NES > 0, "ID"], 2)
top_down <- head(gsea_df[gsea_df$NES < 0, "ID"], 2)

temp_dir <- "figures/Figure4_temp"

for (i in seq_along(top_up)) {
  p <- gseaplot2(gsea_result, top_up[i], title = "",
                 color = bright_red, base_size = 9, rel_heights = c(1.5, 0.5, 1))
  ggsave(file.path(temp_dir, paste0("up_", i, ".png")), p, width = 5, height = 4.5, dpi = 150)
}
for (i in seq_along(top_down)) {
  p <- gseaplot2(gsea_result, top_down[i], title = "",
                 color = dark_blue, base_size = 9, rel_heights = c(1.5, 0.5, 1))
  ggsave(file.path(temp_dir, paste0("down_", i, ".png")), p, width = 5, height = 4.5, dpi = 150)
}

embed_plot_with_label <- function(file, label) {
  if (!file.exists(file)) {
    return(ggplot() + 
             annotate("text", x = 0.5, y = 0.5, label = paste(label, "not available"), 
                      size = 4, family = "serif") +
             theme_void())
  }
  img <- readPNG(file)
  ggplot() +
    annotation_custom(rasterGrob(img, interpolate = TRUE),
                      xmin = -Inf, xmax = Inf, ymin = -Inf, ymax = Inf) +
    coord_cartesian(xlim = c(0, 1), ylim = c(0, 1)) +
    labs(title = label) +
    theme_void() +
    theme(
      plot.title = element_text(
        hjust = 0.02,
        vjust = 0.99,
        size = 11,
        face = "plain",
        family = "serif",
        color = dark_blue,
        margin = margin(0, 0, 0, 0)
      )
    )
}

pB <- embed_plot_with_label(file.path(temp_dir, "up_1.png"), "B")
pC <- embed_plot_with_label(file.path(temp_dir, "up_2.png"), "C")
pD <- embed_plot_with_label(file.path(temp_dir, "down_1.png"), "D")
pE <- embed_plot_with_label(file.path(temp_dir, "down_2.png"), "E")

# ============================================================
# 6. 图 F: 上调通路TOP10
# ============================================================

up_top10 <- head(gsea_df[gsea_df$NES > 0, ], 10)
up_top10$Pathway <- shorten_pathway(up_top10$ID)
up_top10$Pathway <- factor(up_top10$Pathway, levels = rev(up_top10$Pathway))

pF <- ggplot(up_top10, aes(x = NES, y = Pathway, fill = NES)) +
  geom_bar(stat = "identity", width = 0.7) +
  scale_fill_gradient2(low = light_green, mid = warm_orange, high = bright_red,
                       midpoint = 0.5, name = "NES") +
  labs(x = "NES", y = "", title = "F") +
  theme_bw(base_size = 9) +
  theme(
    plot.title = element_text(hjust = 0.5, face = "plain", size = 11, color = dark_blue,
                              family = "serif", margin = margin(b = 2)),
    axis.text.y = element_text(size = 8, color = "#333333", family = "serif"),
    axis.text.x = element_text(size = 8, color = "#333333", family = "serif"),
    axis.title.x = element_text(size = 9, color = "#333333", family = "serif"),
    legend.position = "bottom",
    legend.title = element_text(size = 7, family = "serif"),
    legend.text = element_text(size = 6, family = "serif"),
    panel.grid.major.y = element_blank(),
    panel.grid.minor = element_blank(),
    panel.border = element_rect(color = "gray70", fill = NA, linewidth = 0.5)
  )

# ============================================================
# 7. 图 G: 森林图
# ============================================================

sig_df <- gsea_df[gsea_df$p.adjust < 0.05, ]
sig_df <- sig_df[order(abs(sig_df$NES), decreasing = TRUE), ]
sig_df <- head(sig_df, 25)
sig_df$Pathway <- shorten_pathway(sig_df$ID)
sig_df$Pathway <- factor(sig_df$Pathway, levels = rev(sig_df$Pathway))
sig_df$Category <- ifelse(sig_df$NES > 0, "High-score enriched", "Low-score enriched")
sig_df$logP <- -log10(sig_df$p.adjust)
sig_df$logP <- pmin(sig_df$logP, 8)

pG <- ggplot(sig_df, aes(x = NES, y = Pathway, color = Category, size = logP)) +
  geom_point(alpha = 0.85) +
  geom_vline(xintercept = 0, linetype = "dashed", color = "gray50", linewidth = 0.5) +
  scale_color_manual(values = c("High-score enriched" = bright_red,
                                "Low-score enriched" = dark_blue)) +
  scale_size_continuous(range = c(2.5, 6.5), name = "-log10(FDR)") +
  labs(x = "Normalized Enrichment Score (NES)", y = "", title = "G", color = NULL) +
  theme_bw(base_size = 9) +
  theme(
    plot.title = element_text(hjust = 0.5, face = "plain", size = 11, color = dark_blue,
                              family = "serif", margin = margin(b = 2)),
    axis.text.y = element_text(size = 7.5, color = "#333333", family = "serif"),
    axis.text.x = element_text(size = 8, color = "#333333", family = "serif"),
    axis.title.x = element_text(size = 9, color = "#333333", family = "serif"),
    legend.position = "bottom",
    legend.box = "horizontal",
    legend.title = element_text(size = 7, family = "serif"),
    legend.text = element_text(size = 7, family = "serif"),
    panel.grid.major.y = element_blank(),
    panel.grid.minor = element_blank(),
    panel.border = element_rect(color = "gray70", fill = NA, linewidth = 0.5)
  )

# ============================================================
# 8. 图 H: 下调通路TOP10
# ============================================================

down_top10 <- head(gsea_df[gsea_df$NES < 0, ], 10)
down_top10$Pathway <- shorten_pathway(down_top10$ID)
down_top10$Pathway <- factor(down_top10$Pathway, levels = rev(down_top10$Pathway))

pH <- ggplot(down_top10, aes(x = NES, y = Pathway, fill = NES)) +
  geom_bar(stat = "identity", width = 0.7) +
  scale_fill_gradient2(low = dark_blue, mid = light_green, high = warm_orange,
                       midpoint = -0.5, name = "NES") +
  labs(x = "NES", y = "", title = "H") +
  theme_bw(base_size = 9) +
  theme(
    plot.title = element_text(hjust = 0.5, face = "plain", size = 11, color = dark_blue,
                              family = "serif", margin = margin(b = 2)),
    axis.text.y = element_text(size = 8, color = "#333333", family = "serif"),
    axis.text.x = element_text(size = 8, color = "#333333", family = "serif"),
    axis.title.x = element_text(size = 9, color = "#333333", family = "serif"),
    legend.position = "bottom",
    legend.title = element_text(size = 7, family = "serif"),
    legend.text = element_text(size = 6, family = "serif"),
    panel.grid.major.y = element_blank(),
    panel.grid.minor = element_blank(),
    panel.border = element_rect(color = "gray70", fill = NA, linewidth = 0.5)
  )

# ============================================================
# 9. 组合
# ============================================================

spacer <- ggplot() + theme_void()

row1 <- pA
row2 <- (pB | pC | pD | pE) + plot_layout(widths = c(1, 1, 1, 1))
row3 <- (pF | pH | pG) + plot_layout(widths = c(0.9, 0.9, 1.2))

figure4 <- (row1 / spacer / row2 / spacer / row3) +
  plot_layout(heights = c(0.8, 0.08, 1.0, 0.08, 1.2)) +
  plot_annotation(
    title = "Figure 4: GSEA Dual Phenotype",
    subtitle = "A) NES bar plot (Top 15) | B-E) Representative GSEA curves | F) Up-regulated pathways | H) Down-regulated pathways | G) Forest plot of significant pathways (FDR < 0.05)",
    caption = paste(
      "Note: All plots based on real GSEA results from clusterProfiler.",
      "Red: enriched in high-score group; Blue: enriched in low-score group.",
      "Forest plot point size represents -log10(FDR).",
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
# 10. 出版级输出
# ============================================================

cat("\n开始保存图片...\n")

# PDF
pdf("figures/Figure4.pdf", width = 15, height = 12, onefile = TRUE, paper = "special")
print(figure4)
dev.off()
cat("✅ PDF已保存: figures/Figure4.pdf\n")

# TIFF（使用更稳定的参数）
tryCatch({
  tiff("figures/Figure4.tiff", width = 15, height = 12, units = "in", 
       res = 600, compression = "lzw", type = "windows")
  print(figure4)
  dev.off()
  cat("✅ TIFF已保存: figures/Figure4.tiff\n")
}, error = function(e) {
  cat("⚠️ TIFF 保存失败，尝试使用 png 替代...\n")
  # 如果 TIFF 失败，生成高质量 PNG 作为备选
  png("figures/Figure4.tiff.png", width = 15, height = 12, units = "in", res = 600)
  print(figure4)
  dev.off()
  cat("✅ 已保存为 PNG 格式（600 dpi）: figures/Figure4.tiff.png\n")
})

# PNG
png("figures/Figure4.png", width = 15, height = 12, units = "in", res = 300)
print(figure4)
dev.off()
cat("✅ PNG已保存: figures/Figure4.png\n")

# 清理临时文件
unlink(file.path(temp_dir, "*.png"))

cat("\n============================================================\n")
cat("✅ Figure 4 已保存（出版级）\n")
cat("   PDF:  figures/Figure4.pdf\n")
cat("   TIFF: figures/Figure4.tiff\n")
cat("   PNG:  figures/Figure4.png\n")
cat("\n字母统一:\n")
cat("  A-H 全部使用 plot.title，size = 11pt\n")
cat("  B-E: hjust=0.02, vjust=0.99（左上角贴近边缘）\n")
cat("============================================================\n")