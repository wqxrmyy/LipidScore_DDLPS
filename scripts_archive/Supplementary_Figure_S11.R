# ============================================================
# Figure SX: 脂评分分布小提琴图（四队列）
# ============================================================
library(ggplot2)

# 读取各队列评分
tcga_score <- read.csv("raw/TCGA_Score.csv")
gse21122_score <- read.csv("raw/GSE21122_Score.csv")
gse159659_score <- read.csv("raw/GSE159659_Score.csv")
gse30929_score <- read.csv("raw/bulk/GSE30929_Score_all.csv")

# 统一格式并合并
tcga_score$Cohort <- "TCGA-SARC\n(n=50)"
gse21122_score$Cohort <- "GSE21122\n(n=158)"
gse159659_score$Cohort <- "GSE159659\n(n=45)"
gse30929_score$Cohort <- "GSE30929\n(n=140)"

# 提取脂评分列
tcga_df <- data.frame(Cohort = tcga_score$Cohort, Score = tcga_score$Score)
gse21122_df <- data.frame(Cohort = gse21122_score$Cohort, Score = gse21122_score$Score)
gse159659_df <- data.frame(Cohort = gse159659_score$Cohort, Score = gse159659_score$Score)
gse30929_df <- data.frame(Cohort = gse30929_score$Cohort, Score = gse30929_score$LipidScore)

all_scores <- rbind(tcga_df, gse21122_df, gse159659_df, gse30929_df)
all_scores$Cohort <- factor(all_scores$Cohort, 
                            levels = c("TCGA-SARC\n(n=50)", "GSE21122\n(n=158)", 
                                       "GSE159659\n(n=45)", "GSE30929\n(n=140)"))

# 小提琴图 + 箱线图
p_violin <- ggplot(all_scores, aes(x = Cohort, y = Score, fill = Cohort)) +
  geom_violin(alpha = 0.7, trim = FALSE, color = "grey30", linewidth = 0.3) +
  geom_boxplot(width = 0.12, outlier.shape = NA, fill = "white", alpha = 0.6) +
  scale_fill_manual(values = c("#377EB8", "#4DAF4A", "#984EA3", "#E41A1C")) +
  theme_bw(base_size = 13) +
  theme(
    legend.position = "none",
    panel.grid.minor = element_blank(),
    axis.text.x = element_text(size = 11, face = "bold"),
    axis.text.y = element_text(size = 11),
    axis.title.y = element_text(size = 13),
    plot.title = element_text(face = "bold", size = 14, hjust = 0.5)
  ) +
  labs(
    title = "Lipid Metabolism Score Distribution Across Cohorts",
    x = "",
    y = "Lipid Score (PCA first principal component)"
  )

cairo_pdf("results/figures/extended/FigSX_score_violin.pdf", 
          width = 10, height = 6, family = "Arial")
print(p_violin)
dev.off()

ggsave("results/figures/extended/FigSX_score_violin.png", 
       p_violin, width = 10, height = 6, dpi = 600)

cat("脂评分小提琴图已保存\n")