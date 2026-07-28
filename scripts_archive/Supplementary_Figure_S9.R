setwd("C:/Users/Tao25/Documents/Liposarcoma/2")
library(ggplot2)
library(survival)
library(survminer)

# ==================== Figure 7C：GSE30929按组织学亚型分层KM ====================
clinical_all <- readRDS("raw/bulk/GSE30929_clinical_all.rds")

# 高低级别分组
clinical_all$grade <- ifelse(
  clinical_all$subtype %in% c("dedifferentiated", "pleomorphic"),
  "High-grade", "Low-grade"
)

fit_grade <- survfit(Surv(OS_time, OS_status) ~ ScoreGroup + grade, data = clinical_all)

p_grade <- ggsurvplot(fit_grade, data = clinical_all,
                      pval = TRUE,
                      palette = c("#377EB8", "#E41A1C", "#4DAF4A", "#984EA3"),
                      title = "GSE30929 stratified by histologic grade",
                      xlab = "Time (months)", ylab = "DFS")

ggsave("results/figures/extended/KM_GSE30929_by_grade.pdf", p_grade$plot, width = 10, height = 7)
ggsave("results/figures/extended/KM_GSE30929_by_grade.png", p_grade$plot, width = 10, height = 7, dpi = 300)
cat("Figure 7C 已保存\n")


# ==================== Figure S8：四队列C-index对比柱状图 ====================
cindex_data <- data.frame(
  Cohort = c("TCGA-SARC", "GSE21122", "GSE159659", "GSE30929\n(all subtypes)", "GSE30929\n(DDLPS)"),
  C_index = c(0.72, NA, NA, 0.629, 0.58),
  N = c(50, 158, 45, 140, 40),
  stringsAsFactors = FALSE
)

p_cindex <- ggplot(cindex_data, aes(x = Cohort, y = C_index, fill = Cohort)) +
  geom_bar(stat = "identity", width = 0.6) +
  geom_text(aes(label = ifelse(is.na(C_index), "No survival\ndata", round(C_index, 3))), 
            vjust = -0.5, size = 3.5) +
  scale_fill_manual(values = c("#377EB8", "#4DAF4A", "#984EA3", "#E41A1C", "#FF7F00")) +
  ylim(0, 0.85) +
  theme_bw(base_size = 12) +
  theme(legend.position = "none") +
  labs(title = "C-index across validation cohorts", y = "Harrell's C-index", x = "")

ggsave("results/figures/extended/FigS8_cindex_barplot.pdf", p_cindex, width = 8, height = 5)
ggsave("results/figures/extended/FigS8_cindex_barplot.png", p_cindex, width = 8, height = 5, dpi = 300)
cat("Figure S8 已保存\n")


# ==================== Figure S9：四队列森林图 ====================
forest_data <- data.frame(
  Cohort = c("TCGA-SARC\n(male only)", "GSE30929\n(all subtypes)", "GSE30929\n(DDLPS subset)", "GSE30929\n(High-grade)"),
  HR = c(5.83, 3.11, 0.52, 2.85),
  lower = c(1.55, 1.64, 0.23, 1.47),
  upper = c(21.95, 5.88, 1.22, 5.52),
  P = c("<0.01", "0.00049", "0.13", "0.002"),
  stringsAsFactors = FALSE
)

forest_data$Cohort <- factor(forest_data$Cohort, levels = rev(forest_data$Cohort))

p_forest <- ggplot(forest_data, aes(x = HR, y = Cohort)) +
  geom_vline(xintercept = 1, linetype = "dashed", color = "grey50") +
  geom_point(aes(color = P < 0.05), size = 3) +
  geom_errorbarh(aes(xmin = lower, xmax = upper, color = P < 0.05), height = 0.2, linewidth = 1) +
  geom_text(aes(x = max(upper) * 1.3, label = paste0("HR=", round(HR, 2), "\nP=", P))),
size = 3.5, hjust = 0) +
  scale_color_manual(values = c("TRUE" = "#E41A1C", "FALSE" = "grey60")) +
  xlim(0.1, max(forest_data$upper) * 2) +
  theme_bw(base_size = 12) +
  theme(legend.position = "none") +
  labs(title = "Prognostic value of lipid score across cohorts", 
       x = "Hazard Ratio (95% CI)", y = "")

ggsave("results/figures/extended/FigS9_forest_plot.pdf", p_forest, width = 10, height = 5)
ggsave("results/figures/extended/FigS9_forest_plot.png", p_forest, width = 10, height = 5, dpi = 300)
cat("Figure S9 已保存\n")