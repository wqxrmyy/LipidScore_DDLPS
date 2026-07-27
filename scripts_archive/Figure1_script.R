# ============================================================
# Figure 1: Study Design and Cohort Characteristics
# 出版级规范版本（Biomarker Research / Cancer Science 通用）
# 子图: A, B, C, D, E, F
# 字体: Times New Roman, 字号 9-11 pt（SCI标准）
# 输出: PDF（矢量）, TIFF（600 dpi LZW）, PNG（300 dpi）
# ============================================================

rm(list = ls())
setwd("C:/Users/Tao25/Documents/Liposarcoma/1")

# 加载包（无需额外安装）
library(ggplot2)
library(patchwork)
library(ggpubr)

# ============================================================
# 1. 颜色系统（色盲友好 + 灰度兼容 + 高对比度）
# ============================================================

dark_blue <- "#0F4C81"      # 主色（流程图、TCGA）
soft_pink <- "#E69F8C"       # 女性
lake_blue <- "#6BAED6"       # 男性
light_green <- "#8DAA7D"     # 存活
bright_red <- "#C93312"      # 死亡
warm_orange <- "#E69F00"     # 队列/亚型
gray_light <- "#E8E8E8"      # 流程图中间节点
indigo <- "#6053a0"

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
# 3. 图 A: 流程图（纯 annotate，无额外依赖）
# 字体: 节点标题 11pt，描述 10pt
# ============================================================

pA <- ggplot() +
  # 节点1
  annotate("rect", xmin = 0.08, xmax = 0.92, ymin = 0.72, ymax = 0.92,
           fill = dark_blue, color = "white", alpha = 0.9, linewidth = 0.5) +
  annotate("text", x = 0.5, y = 0.82, label = "TCGA / GEO\nDatasets", 
           size = 4.5, color = "white", fontface = "bold", family = "serif", lineheight = 1.2) +
  # 箭头1
  annotate("segment", x = 0.5, xend = 0.5, y = 0.72, yend = 0.58,
           arrow = arrow(length = unit(0.2, "cm"), type = "closed"),
           linewidth = 0.8, color = "gray40") +
  # 节点2
  annotate("rect", xmin = 0.08, xmax = 0.92, ymin = 0.42, ymax = 0.58,
           fill = gray_light, color = "gray50", linewidth = 0.5) +
  annotate("text", x = 0.5, y = 0.50, label = "Pathology Review &\nData Filtering",
           size = 4.0, color = "#222222", family = "serif", lineheight = 1.2) +
  # 箭头2
  annotate("segment", x = 0.5, xend = 0.5, y = 0.42, yend = 0.28,
           arrow = arrow(length = unit(0.2, "cm"), type = "closed"),
           linewidth = 0.8, color = "gray40") +
  # 节点3
  annotate("rect", xmin = 0.08, xmax = 0.92, ymin = 0.12, ymax = 0.28,
           fill = dark_blue, color = "white", alpha = 0.9, linewidth = 0.5) +
  annotate("text", x = 0.5, y = 0.20, label = "Final Cohort\n(n = 253)",
           size = 4.5, color = "white", fontface = "bold", family = "serif", lineheight = 1.2) +
  xlim(0, 1) + ylim(0, 1) +
  theme_void() +
  labs(title = "A") +
  theme(
    plot.title = element_text(hjust = 0.5, face = "plain", size = 11, color = dark_blue,
                              family = "serif", margin = margin(b = 2))
  )

# ============================================================
# 4. 图 B: 样本来源
# ============================================================

source_data <- data.frame(
  Cohort = factor(c("TCGA", "GSE21122", "GSE159659"), 
                  levels = c("TCGA", "GSE21122", "GSE159659")),
  Count = c(50, 158, 45)
)
source_data$Percent <- round(source_data$Count / sum(source_data$Count) * 100, 1)

pB <- ggplot(source_data, aes(x = Cohort, y = Count, fill = Cohort)) +
  geom_bar(stat = "identity", width = 0.7, color = "white", linewidth = 0.3) +
  scale_fill_manual(values = c(dark_blue, warm_orange, light_green)) +
  geom_text(aes(label = paste0(Count, "\n(", Percent, "%)")),
            vjust = -0.8, size = 3.5, fontface = "plain", family = "serif", lineheight = 1.2) +
  labs(x = "", y = "Number of samples", title = "B") +
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
    panel.border = element_rect(color = "gray70", linewidth = 0.5)
  )

# ============================================================
# 5. 图 C: 亚型分布
# ============================================================

subtype_data <- data.frame(
  Subtype = factor(c("DDLPS", "WDLPS", "Myxoid", "Pleomorphic", "Other"),
                   levels = c("DDLPS", "WDLPS", "Myxoid", "Pleomorphic", "Other")),
  Count = c(61, 50, 15, 20, 23)
)
subtype_data$Percent <- round(subtype_data$Count / sum(subtype_data$Count) * 100, 1)

pC <- ggplot(subtype_data, aes(x = Subtype, y = Count, fill = Subtype)) +
  geom_bar(stat = "identity", width = 0.7, color = "white", linewidth = 0.3) +
  scale_fill_manual(values = c(dark_blue, light_green, warm_orange, bright_red, "gray65")) +
  geom_text(aes(label = paste0(Count, "\n(", Percent, "%)")),
            vjust = -0.8, size = 3.2, fontface = "plain", family = "serif", lineheight = 1.1) +
  labs(x = "", y = "Number of patients", title = "C") +
  theme_bw(base_size = 10) +
  theme(
    legend.position = "none",
    plot.title = element_text(hjust = 0.5, face = "plain", size = 11, color = dark_blue,
                              family = "serif", margin = margin(b = 2)),
    axis.text.x = element_text(angle = 45, hjust = 1, size = 9, color = "#333333", family = "serif"),
    axis.text.y = element_text(size = 9, color = "#333333", family = "serif"),
    axis.title.y = element_text(size = 10, color = "#333333", family = "serif"),
    panel.grid.major.x = element_blank(),
    panel.grid.minor = element_blank(),
    panel.border = element_rect(color = "gray70", linewidth = 0.5)
  )

# ============================================================
# 6. 图 D: 年龄分布
# ============================================================

median_age <- median(tcga_data$Age, na.rm = TRUE)

pD <- ggplot(tcga_data, aes(x = Age)) +
  geom_histogram(bins = 15, fill = lake_blue, alpha = 0.8, color = "white", linewidth = 0.2) +
  geom_vline(xintercept = median_age, linetype = "dashed", color = bright_red, linewidth = 0.8) +
  annotate("text", x = median_age + 6, y = 6, 
           label = paste0("Median = ", round(median_age, 1), " years"),
           size = 3.2, color = dark_blue, fontface = "plain", family = "serif", hjust = 0) +
  labs(x = "Age (years)", y = "Count", title = "D") +
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
# 7. 图 E: 性别分布
# ============================================================

gender_counts <- table(tcga_data$Gender)
gender_data <- data.frame(
  Gender = names(gender_counts),
  Count = as.numeric(gender_counts),
  Percent = round(as.numeric(gender_counts) / sum(gender_counts) * 100, 1)
)

pE <- ggplot(gender_data, aes(x = Gender, y = Count, fill = Gender)) +
  geom_bar(stat = "identity", width = 0.6, color = "white", linewidth = 0.3) +
  scale_fill_manual(values = c(soft_pink, lake_blue)) +
  geom_text(aes(label = paste0(Count, "\n(", Percent, "%)")),
            vjust = -0.8, size = 3.5, fontface = "plain", family = "serif", lineheight = 1.2) +
  labs(x = "", y = "Number of patients", title = "E") +
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
# 8. 图 F: 生存状态
# ============================================================

survival_counts <- table(tcga_data$OS_status)
survival_data <- data.frame(
  Status = c("Alive", "Dead"),
  Count = as.numeric(survival_counts),
  Percent = round(as.numeric(survival_counts) / sum(survival_counts) * 100, 1)
)

pF <- ggplot(survival_data, aes(x = Status, y = Count, fill = Status)) +
  geom_bar(stat = "identity", width = 0.6, color = "white", linewidth = 0.3) +
  scale_fill_manual(values = c(light_green, bright_red)) +
  geom_text(aes(label = paste0(Count, "\n(", Percent, "%)")),
            vjust = -0.8, size = 3.5, fontface = "plain", family = "serif", lineheight = 1.2) +
  labs(x = "", y = "Number of patients", title = "F") +
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
# 9. 组合（优化布局 + 统一对齐）
# ============================================================

row1 <- (pA | pB | pC) + plot_layout(widths = c(1.1, 1, 1.1))
row2 <- (pD | pE | pF) + plot_layout(widths = c(1.1, 0.9, 0.9))

figure1 <- (row1 / row2) +
  plot_layout(heights = c(1, 0.9)) +
  plot_annotation(
    title = "Figure 1: Study Design and Cohort Characteristics",
    theme = theme(
      plot.title = element_text(hjust = 0.5, face = "bold", size = 14, color = dark_blue,
                                family = "serif", margin = margin(b = 8))
    )
  )

# ============================================================
# 10. 出版级输出（多格式 + 高分辨率）
# ============================================================

# PDF（矢量格式，期刊首选）
ggsave("figures/Figure1.pdf",
       figure1,
       width = 14,
       height = 9,
       units = "in",
       dpi = 600,
       limitsize = FALSE,
       device = cairo_pdf)

# TIFF（高分辨率位图，LZW压缩）
ggsave("figures/Figure1.tiff",
       figure1,
       width = 14,
       height = 9,
       units = "in",
       dpi = 600,
       compression = "lzw")

# PNG（预览用）
ggsave("figures/Figure1.png",
       figure1,
       width = 14,
       height = 9,
       units = "in",
       dpi = 300)

cat("\n============================================================\n")
cat("✅ Figure 1 已保存（出版级）\n")
cat("   PDF:  figures/Figure1.pdf  (矢量, 600 dpi)\n")
cat("   TIFF: figures/Figure1.tiff (LZW压缩, 600 dpi)\n")
cat("   PNG:  figures/Figure1.png  (300 dpi)\n")
cat("\n字体规范: Times New Roman (serif)\n")
cat("字号规范: 标题 11pt, 轴标签 9-10pt, 数值标注 9-10pt\n")
cat("颜色规范: 色盲友好 + 灰度打印兼容\n")
cat("============================================================\n")