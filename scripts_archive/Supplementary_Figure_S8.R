# ============================================================
# Figure S8: C-index（仅生存验证队列）
# ============================================================
library(ggplot2)

cindex_data <- data.frame(
  Cohort = c("TCGA-SARC\n(OS, n=50)", "GSE30929 all\n(DFS, n=140)", "GSE30929 DDLPS\n(DFS, n=40)"),
  C_index = c(0.72, 0.629, 0.58),
  Endpoint = c("OS", "DFS", "DFS"),
  stringsAsFactors = FALSE
)

cindex_data$Cohort <- factor(cindex_data$Cohort, levels = cindex_data$Cohort)

p_cindex <- ggplot(cindex_data, aes(x = Cohort, y = C_index, fill = Endpoint)) +
  geom_bar(stat = "identity", width = 0.55, color = "grey30", linewidth = 0.3) +
  geom_text(aes(label = sprintf("%.3f", C_index)), vjust = -0.6, size = 5, 
            fontface = "bold", color = "grey20") +
  scale_fill_manual(values = c("OS" = "#377EB8", "DFS" = "#E41A1C")) +
  ylim(0, 0.85) +
  theme_bw(base_size = 13) +
  theme(
    legend.position = c(0.88, 0.88),
    legend.background = element_rect(fill = "white", color = "grey80", linewidth = 0.3),
    legend.title = element_text(size = 10),
    legend.text = element_text(size = 9),
    panel.grid.major.x = element_blank(),
    panel.grid.minor = element_blank(),
    axis.text.x = element_text(size = 11, face = "bold"),
    axis.text.y = element_text(size = 11),
    axis.title.y = element_text(size = 13),
    plot.title = element_text(face = "bold", size = 14, hjust = 0.5),
    plot.subtitle = element_text(size = 9, hjust = 0.5, color = "grey40")
  ) +
  labs(
    title = "Model Discrimination Across Survival Cohorts",
    subtitle = "Harrell's C-index by endpoint | Additional expression validation in GSE21122 (n=158) and GSE159659 (n=45)",
    x = "",
    y = "C-index",
    fill = "Endpoint"
  )

# 保存
cairo_pdf("results/figures/extended/FigS8_cindex_barplot_final.pdf", 
          width = 8, height = 5.5, family = "Arial")
print(p_cindex)
dev.off()

ggsave("results/figures/extended/FigS8_cindex_barplot_final.png", 
       p_cindex, width = 8, height = 5.5, dpi = 600)

file_size_pdf <- file.info("results/figures/extended/FigS8_cindex_barplot_final.pdf")$size / 1024
cat(sprintf("Figure S8 PDF: %.1f KB\n", file_size_pdf))
cat("Figure S8 修复版已保存\n")