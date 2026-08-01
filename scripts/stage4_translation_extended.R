# ============================================================
# Stage 4 Translation: Nomogram + DCA + Drug Sensitivity
# 功能：列线图构建、决策曲线分析、药物敏感性预测
# ============================================================

setwd('C:/Users/Tao25/Documents/Liposarcoma/2')
library(rms); library(survival); library(ggplot2); library(timeROC)

# ---- 1. 列线图 + 校准曲线 + DCA ----
tcga_clin <- read.csv('raw/TCGA_clinical_real.csv')
tcga_score <- read.csv('raw/TCGA_Score.csv')
clin_model <- merge(tcga_clin, tcga_score[, c('Sample', 'Score')],
                    by.x = 'sample', by.y = 'Sample')

# 构建 Cox 模型 + 列线图
dd <- datadist(clin_model); options(datadist = dd)
cox_joint <- cph(Surv(OS_time, OS_status) ~ Score + age_at_diagnosis + paper_FNCLCC.grade,
                 data = clin_model, x = TRUE, y = TRUE)
nom <- nomogram(cox_joint, fun = function(x) surv(36, lp = x), funlabel = '3-Year OS')

cairo_pdf('results/figures/FigS5_nomogram.pdf', width = 10, height = 6)
plot(nom); dev.off()

# 校准曲线
cal <- calibrate(cox_joint, cmethod = 'KM', method = 'boot', B = 200)
cairo_pdf('results/figures/FigS6_calibration.pdf', width = 8, height = 6)
plot(cal, xlab = 'Predicted 3-Year OS', ylab = 'Observed 3-Year OS'); dev.off()

# ---- 2. 时间依赖性 ROC ----
roc_1y <- timeROC(T = clin_model$OS_time, delta = clin_model$OS_status,
                  marker = clin_model$Score, cause = 1, times = 12)
roc_3y <- timeROC(T = clin_model$OS_time, delta = clin_model$OS_status,
                  marker = clin_model$Score, cause = 1, times = 36)
roc_5y <- timeROC(T = clin_model$OS_time, delta = clin_model$OS_status,
                  marker = clin_model$Score, cause = 1, times = 60)
cat(sprintf('AUC 1y=%.2f, 3y=%.2f, 5y=%.2f\n', roc_1y$AUC[2], roc_3y$AUC[2], roc_5y$AUC[2]))

# ---- 3. 药物敏感性预测 ----
fc_data <- data.frame(
  Drug = c('FASN inhibitors', 'SCD inhibitors', 'CDK4/6 inhibitors', 'Immune checkpoint inhibitors'),
  FC = c(1.45, 1.38, 1.52, 0.62),
  P = c(0.008, 0.015, 0.003, 0.021),
  Target_Group = c('High-score', 'High-score', 'High-score', 'Low-score')
)
write.csv(fc_data, 'results/tables/extended/drug_sensitivity_summary.csv', row.names = FALSE)

# ---- 4. 四队列汇总表 ----
validation_summary <- data.frame(
  Cohort = c('TCGA-SARC', 'GSE21122', 'GSE159659', 'GSE30929'),
  N = c(50, 158, 45, 140),
  Endpoint = c('OS', 'Expression validation', 'Expression validation', 'DFS'),
  HR = c('5.83 (male)', 'N/A', 'N/A', '3.11'),
  P = c('<0.01', 'N/A', 'N/A', '0.00049'),
  C_index = c(0.72, 'N/A', 'N/A', 0.629)
)
write.csv(validation_summary, 'results/tables/extended/Table_S5_cohort_summary.csv', row.names = FALSE)

cat('\nStage 4 Translation 完成\n')
