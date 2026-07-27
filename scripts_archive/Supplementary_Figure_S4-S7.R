# =============================================
# 最终精简版：补充图 S4-S7
# S4 = KM 曲线 | S5 = 列线图 | S6 = 校准曲线 | S7 = 亚型验证合并图
# =============================================

library(survival)
library(survminer)
library(ggplot2)
library(rms)
library(cowplot)
library(patchwork)
library(ggpubr)

out_dir <- "C:/Users/Tao25/Documents/Liposarcoma/1/data/raw/"

# ---------- 确保 clinical_merged 存在 ----------
if (!exists("clinical_merged")) {
  clinical <- read.csv(paste0(out_dir, "TCGA_clinical_real.csv"),
                       stringsAsFactors = FALSE, check.names = FALSE)
  clinical_sub <- clinical[, c("barcode", "OS_time", "OS_status", 
                               "age_at_diagnosis", "paper_FNCLCC.grade")]
  colnames(clinical_sub) <- c("SampleID", "OS_time", "OS_status", "Age", "Grade")
  clinical_sub$OS_status <- as.numeric(clinical_sub$OS_status)
  clinical_sub$OS_time <- as.numeric(clinical_sub$OS_time)
  clinical_sub$Grade <- as.numeric(clinical_sub$Grade)
  clinical_sub <- clinical_sub[!is.na(clinical_sub$OS_time) & !is.na(clinical_sub$OS_status), ]
  score <- readRDS(paste0(out_dir, "TCGA_Lipid_Score.rds"))
  score_df <- data.frame(SampleID = names(score), Lipid_Score = score)
  clinical_merged <- merge(clinical_sub, score_df, by = "SampleID", all.x = TRUE)
  clinical_merged <- clinical_merged[!is.na(clinical_merged$Lipid_Score), ]
}

# ---------- Figure S4: KM 曲线 ----------
cutpoint <- surv_cutpoint(clinical_merged, time = "OS_time", event = "OS_status", variables = "Lipid_Score")
optimal_cut <- cutpoint$cutpoint$cutpoint
clinical_merged$RiskGroup <- ifelse(clinical_merged$Lipid_Score > optimal_cut, "High Risk", "Low Risk")
clinical_merged$RiskGroup <- factor(clinical_merged$RiskGroup, levels = c("Low Risk", "High Risk"))

fit_km <- survfit(Surv(OS_time, OS_status) ~ RiskGroup, data = clinical_merged)
p_km <- ggsurvplot(fit_km, data = clinical_merged,
                   pval = TRUE, pval.method = TRUE, conf.int = TRUE,
                   risk.table = TRUE, risk.table.col = "strata",
                   palette = c("#377EB8", "#E41A1C"),
                   xlab = "Time (days)", ylab = "Overall Survival Probability",
                   title = "TCGA-LPS: Overall Survival by Lipid Score",
                   legend.title = "Risk Group", ggtheme = theme_classic2())

# 拼合 KM 曲线 + 风险表
km_plot <- p_km$plot + theme(legend.position = "top")
km_table <- p_km$table + theme(legend.position = "none")
combined_km <- plot_grid(km_plot, km_table, ncol = 1, rel_heights = c(3, 1))

ggsave(paste0(out_dir, "Figure_S4.pdf"), plot = combined_km, width = 8, height = 6)
ggsave(paste0(out_dir, "Figure_S4.png"), plot = combined_km, width = 8, height = 6, dpi = 300)
ggsave(paste0(out_dir, "Figure_S4.tiff"), plot = combined_km, width = 8, height = 6, dpi = 600, compression = "lzw")

# ---------- Figure S5: 列线图 ----------
dd <- datadist(clinical_merged)
options(datadist = "dd")
cox_nomo <- cph(Surv(OS_time, OS_status) ~ Age + Grade + Lipid_Score,
                data = clinical_merged, x = TRUE, y = TRUE, surv = TRUE)
sf <- survfit(cox_nomo)
S0_365 <- summary(sf, times = 365)$surv[1]
S0_1095 <- summary(sf, times = 1095)$surv[1]
S0_1825 <- summary(sf, times = 1825)$surv[1]

nom <- nomogram(cox_nomo, fun = list(function(x) S0_365^exp(x),
                                     function(x) S0_1095^exp(x),
                                     function(x) S0_1825^exp(x)),
                funlabel = c("1-year OS", "3-year OS", "5-year OS"), lp = FALSE)

pdf(paste0(out_dir, "Figure_S5.pdf"), width = 10, height = 8)
plot(nom)
dev.off()
png(paste0(out_dir, "Figure_S5.png"), width = 1000, height = 800, res = 300)
plot(nom)
dev.off()
tiff(paste0(out_dir, "Figure_S5.tiff"), width = 2000, height = 1600, res = 600, compression = "lzw", type = "windows")
plot(nom)
dev.off()

# ---------- Figure S6: 校准曲线 ----------
pdf(paste0(out_dir, "Figure_S6.pdf"), width = 12, height = 4)
par(mfrow = c(1, 3))
for (t in c(365, 1095, 1825)) {
  cal <- calibrate(cox_nomo, method = "boot", u = t, B = 200)
  plot(cal, xlab = "Predicted Probability", ylab = "Observed Probability",
       main = paste(t/365, "-year", sep = ""), col = "blue", lwd = 2)
  abline(0, 1, col = "red", lty = 2)
}
dev.off()

png(paste0(out_dir, "Figure_S6.png"), width = 1200, height = 400, res = 300)
par(mfrow = c(1, 3))
for (t in c(365, 1095, 1825)) {
  cal <- calibrate(cox_nomo, method = "boot", u = t, B = 200)
  plot(cal, xlab = "Predicted Probability", ylab = "Observed Probability",
       main = paste(t/365, "-year", sep = ""), col = "blue", lwd = 2)
  abline(0, 1, col = "red", lty = 2)
}
dev.off()

tiff(paste0(out_dir, "Figure_S6.tiff"), width = 2400, height = 800, res = 600, compression = "lzw", type = "windows")
par(mfrow = c(1, 3))
for (t in c(365, 1095, 1825)) {
  cal <- calibrate(cox_nomo, method = "boot", u = t, B = 200)
  plot(cal, xlab = "Predicted Probability", ylab = "Observed Probability",
       main = paste(t/365, "-year", sep = ""), col = "blue", lwd = 2)
  abline(0, 1, col = "red", lty = 2)
}
dev.off()

# ---------- Figure S7: 亚型验证合并图 (GSE159659 + GSE21122) ----------
# ---- S7a: GSE159659 ----
con <- gzfile(paste0(out_dir, "GSE159659/GSE159659_series_matrix.txt.gz"), "r")
lines <- readLines(con, n = 200)
close(con)
title_line <- grep("!Sample_title", lines, value = TRUE)
title_parts <- strsplit(title_line, "\t")[[1]]
titles <- title_parts[-1]
subtype_raw <- sapply(strsplit(titles, ","), function(x) trimws(x[1]))
subtype_clean <- gsub('"', '', subtype_raw)
subtype_map159 <- rep(NA, length(subtype_clean))
subtype_map159[grepl("dedifferentiated", subtype_clean, ignore.case = TRUE)] <- "DDLPS"
subtype_map159[grepl("well differentiated", subtype_clean, ignore.case = TRUE)] <- "WDLPS"
subtype_map159[grepl("adipose", subtype_clean, ignore.case = TRUE)] <- "Normal"
gsm_line <- grep("!Sample_geo_accession", lines, value = TRUE)
gsm_parts <- strsplit(gsm_line, "\t")[[1]]
gsm_ids <- gsub('"', '', gsm_parts[-1])
names(subtype_map159) <- gsm_ids

expr_raw <- read.table(paste0(out_dir, "GSE159659/GSE159659_series_matrix.txt.gz"),
                       header = TRUE, sep = "\t", row.names = 1,
                       comment.char = "!", check.names = FALSE, stringsAsFactors = FALSE)
score_159659 <- readRDS(paste0(out_dir, "GSE159659_Lipid_Score.rds"))
score_159659 <- score_159659[colnames(expr_raw)]
subtype_vec159 <- subtype_map159[colnames(expr_raw)]
df159 <- data.frame(Score = score_159659, Subtype = subtype_vec159, stringsAsFactors = FALSE)
df159 <- df159[!is.na(df159$Subtype), ]
df159$Subtype <- factor(df159$Subtype, levels = c("DDLPS", "WDLPS", "Normal"))

p7a <- ggplot(df159, aes(x = Subtype, y = Score, fill = Subtype)) +
  geom_boxplot(outlier.shape = NA, alpha = 0.6) +
  geom_jitter(width = 0.2, size = 1.5, alpha = 0.7) +
  stat_compare_means(comparisons = combn(unique(df159$Subtype), 2, simplify = FALSE),
                     method = "t.test", label = "p.signif") +
  labs(title = "GSE159659", x = "", y = "Lipid Score (z-score)") +
  theme_classic() +
  theme(legend.position = "none", axis.text.x = element_text(angle = 0, hjust = 0.5))

# ---- S7b: GSE21122 ----
con <- gzfile(paste0(out_dir, "GSE21122/GSE21122_series_matrix.txt.gz"), "r")
lines <- readLines(con, n = 500)
close(con)
gsm_line <- grep("!Sample_geo_accession", lines, value = TRUE)
gsm_parts <- strsplit(gsm_line, "\t")[[1]]
gsm_ids <- gsub('"', '', gsm_parts[-1])
title_line <- grep("!Sample_title", lines, value = TRUE)
title_parts <- strsplit(title_line, "\t")[[1]]
titles <- gsub('"', '', title_parts[-1])
subtype_code <- gsub(".*PT[0-9]+(.*)", "\\1", titles)
subtype_code <- trimws(subtype_code)
subtype_code[subtype_code == ""] <- "Normal"
subtype_map211 <- setNames(subtype_code, gsm_ids)

score_21122 <- readRDS(paste0(out_dir, "GSE21122_Lipid_Score.rds"))
sample_names <- names(score_21122)
sample_names_clean <- gsub("\\..*$", "", sample_names)
subtype_vec211 <- subtype_map211[sample_names_clean]
df211 <- data.frame(Score = score_21122, Subtype = subtype_vec211, stringsAsFactors = FALSE)
df211 <- df211[!is.na(df211$Subtype) & df211$Subtype != "", ]
main_subtypes <- names(which(table(df211$Subtype) >= 5))
df211 <- df211[df211$Subtype %in% main_subtypes, ]
df211$Subtype <- factor(df211$Subtype)

p7b <- ggplot(df211, aes(x = Subtype, y = Score, fill = Subtype)) +
  geom_boxplot(outlier.shape = NA, alpha = 0.6) +
  geom_jitter(width = 0.2, size = 1.5, alpha = 0.7) +
  stat_compare_means(method = "anova", label = "p.format") +
  labs(title = "GSE21122", x = "", y = "Lipid Score (z-score)") +
  theme_classic() +
  theme(legend.position = "none", axis.text.x = element_text(angle = 45, hjust = 1))

# ---- 合并 S7a + S7b ----
p7 <- p7a + p7b + 
  plot_annotation(title = "Lipid Score Distribution by Subtype in Independent Cohorts",
                  theme = theme(plot.title = element_text(hjust = 0.5, size = 14)))

ggsave(paste0(out_dir, "Figure_S7.pdf"), plot = p7, width = 12, height = 5, device = "pdf")
ggsave(paste0(out_dir, "Figure_S7.png"), plot = p7, width = 12, height = 5, dpi = 300)
ggsave(paste0(out_dir, "Figure_S7.tiff"), plot = p7, width = 12, height = 5, dpi = 600, compression = "lzw")

cat("✅ 补充图 S4-S7 已全部生成（PDF/PNG/TIFF）。\n")
cat("文件列表：Figure_S4.pdf/png/tiff, Figure_S5.pdf/png/tiff, Figure_S6.pdf/png/tiff, Figure_S7.pdf/png/tiff\n")