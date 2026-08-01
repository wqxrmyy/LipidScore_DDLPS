# ============================================================
# Fig. 7 完整代码：GSE30929 外部验证（5子图，无警告）
# ============================================================
setwd("C:/Users/Tao25/Documents/Liposarcoma/2")
library(survival)
library(survminer)
library(ggplot2)
library(ggpubr)
library(patchwork)

# 读取数据
clinical_all <- readRDS("raw/bulk/GSE30929_clinical_all.rds")
expr_6genes <- readRDS("raw/bulk/GSE30929_6genes_expr.rds")

clinical_all$ScoreGroup <- ifelse(
  clinical_all$LipidScore > median(clinical_all$LipidScore), "High", "Low"
)
clinical_all$ScoreGroup <- factor(clinical_all$ScoreGroup, levels = c("Low", "High"))

# ---- Panel A: 全部140例 DFS ----
fit_all <- survfit(Surv(OS_time, OS_status) ~ ScoreGroup, data = clinical_all)
pA <- ggsurvplot(fit_all, data = clinical_all, pval = TRUE,
                 legend.title = "Lipid Score", legend.labs = c("Low", "High"),
                 title = "All patients (n=140)",
                 xlab = "Time (months)", ylab = "Disease-Free Survival",
                 risk.table = TRUE, risk.table.height = 0.2,
                 ggtheme = theme_bw())
pA$plot <- pA$plot + scale_color_manual(values = c("Low" = "#377EB8", "High" = "#E41A1C"))

# ---- Panel B: DDLPS 亚组 (n=40) ----
clinical_ddlps <- subset(clinical_all, subtype == "dedifferentiated")
clinical_ddlps$ScoreGroup <- ifelse(
  clinical_ddlps$LipidScore > median(clinical_ddlps$LipidScore), "High", "Low"
)
clinical_ddlps$ScoreGroup <- factor(clinical_ddlps$ScoreGroup, levels = c("Low", "High"))
fit_ddlps <- survfit(Surv(OS_time, OS_status) ~ ScoreGroup, data = clinical_ddlps)
pB <- ggsurvplot(fit_ddlps, data = clinical_ddlps, pval = TRUE,
                 legend.title = "Lipid Score", legend.labs = c("Low", "High"),
                 title = "DDLPS subset (n=40)",
                 xlab = "Time (months)", ylab = "Disease-Free Survival",
                 risk.table = TRUE, risk.table.height = 0.2,
                 ggtheme = theme_bw())
pB$plot <- pB$plot + scale_color_manual(values = c("Low" = "#377EB8", "High" = "#E41A1C"))

# ---- Panel C: 组织学分层 ----
clinical_all$grade <- ifelse(
  clinical_all$subtype %in% c("dedifferentiated", "pleomorphic"), "High-grade", "Low-grade"
)
fit_grade <- survfit(Surv(OS_time, OS_status) ~ ScoreGroup + grade, data = clinical_all)
pC <- ggsurvplot(fit_grade, data = clinical_all, pval = TRUE,
                 legend.title = "Group",
                 legend.labs = c("Low, Low-grade", "Low, High-grade", 
                                 "High, Low-grade", "High, High-grade"),
                 title = "Stratified by histologic grade",
                 xlab = "Time (months)", ylab = "Disease-Free Survival",
                 risk.table = TRUE, risk.table.height = 0.25,
                 ggtheme = theme_bw())
pC$plot <- pC$plot + scale_color_manual(
  values = c("Low, Low-grade" = "#377EB8", "Low, High-grade" = "#4DAF4A",
             "High, Low-grade" = "#E41A1C", "High, High-grade" = "#984EA3")
)

# ---- Panel D: 三分位 KM ----
tertiles <- quantile(clinical_all$LipidScore, probs = c(0, 1/3, 2/3, 1))
clinical_all$ScoreTertile <- cut(clinical_all$LipidScore, 
                                 breaks = tertiles, 
                                 labels = c("Low", "Mid", "High"),
                                 include.lowest = TRUE)
clinical_all$ScoreTertile <- factor(clinical_all$ScoreTertile, levels = c("Low", "Mid", "High"))
fit_tertile <- survfit(Surv(OS_time, OS_status) ~ ScoreTertile, data = clinical_all)
pD <- ggsurvplot(fit_tertile, data = clinical_all, pval = TRUE,
                 legend.title = "Score Tertile", legend.labs = c("Low", "Mid", "High"),
                 title = "Tertile stratification (n=140)",
                 xlab = "Time (months)", ylab = "Disease-Free Survival",
                 risk.table = TRUE, risk.table.height = 0.25,
                 ggtheme = theme_bw())
pD$plot <- pD$plot + scale_color_manual(
  values = c("Low" = "#377EB8", "Mid" = "#4DAF4A", "High" = "#E41A1C")
)

# ---- Panel E: 六基因表达箱线图 ----
expr_long <- data.frame(
  Gene = rep(rownames(expr_6genes), each = ncol(expr_6genes)),
  Expression = as.vector(t(expr_6genes)),
  Sample = rep(colnames(expr_6genes), nrow(expr_6genes)),
  ScoreGroup = clinical_all$ScoreGroup[match(colnames(expr_6genes), clinical_all$sample)]
)
expr_long$ScoreGroup <- factor(expr_long$ScoreGroup, levels = c("Low", "High"))

pE <- ggplot(expr_long, aes(x = Gene, y = Expression, fill = ScoreGroup)) +
  geom_boxplot(outlier.shape = NA, alpha = 0.7, position = position_dodge(0.8)) +
  geom_jitter(aes(color = ScoreGroup), size = 0.5, alpha = 0.4,
              position = position_jitterdodge(jitter.width = 0.15, dodge.width = 0.8)) +
  scale_fill_manual(values = c("Low" = "#377EB8", "High" = "#E41A1C")) +
  scale_color_manual(values = c("Low" = "#377EB8", "High" = "#E41A1C")) +
  stat_compare_means(aes(group = ScoreGroup), label = "p.signif", 
                     method = "wilcox.test", hide.ns = TRUE) +
  theme_bw(base_size = 11) +
  theme(legend.position = "top",
        axis.text.x = element_text(angle = 30, hjust = 1, face = "italic")) +
  labs(title = "Six-gene expression by score group (GSE30929)", x = "", y = "Expression (log2)")

# ---- 组合 Fig 7 ----
layout <- "
AABB
AABB
CCDD
CCDD
EEEE
EEEE
"

p_combined <- wrap_plots(
  A = pA$plot, B = pB$plot, C = pC$plot, D = pD$plot, E = pE,
  design = layout
) + plot_annotation(
  title = "Fig. 7. External Validation of the Lipid Metabolism Score in GSE30929",
  theme = theme(plot.title = element_text(face = "bold", size = 14, hjust = 0.5))
)

# ---- 保存 ----
cairo_pdf("results/figures/Fig7_GSE30929_validation.pdf", width = 16, height = 20)
print(p_combined)
dev.off()

ggsave("results/figures/Fig7_GSE30929_validation.png", p_combined, width = 16, height = 20, dpi = 600)

file.copy("results/figures/Fig7_GSE30929_validation.pdf",
          "C:/Users/Tao25/Documents/Liposarcoma/2/results/figures/Figure7.pdf", overwrite = TRUE)

cat("Fig. 7 完成（无警告版本）\n")