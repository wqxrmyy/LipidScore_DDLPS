# ============================================================
# Stage 2 Advanced: GSE30929 External Validation + Mediation + MDM2
# 功能：外部独立队列验证、Bootstrap中介分析、MDM2拷贝数关联
# ============================================================

setwd('C:/Users/Tao25/Documents/Liposarcoma/2')
library(survival); library(survminer); library(ggplot2); library(GEOquery)

# ---- 1. GSE30929 外部验证 ----
gse30929 <- getGEO(filename = 'raw/bulk/GSE30929/GSE30929_series_matrix.txt.gz',
                   GSEMatrix = TRUE, getGPL = FALSE)
expr_30929 <- exprs(gse30929)

# 提取六基因
gpl_table <- read.delim('raw/bulk/GPL570-55999.txt', header = TRUE, sep = '	',
                        stringsAsFactors = FALSE, comment.char = '#')
probe_to_gene <- gpl_table[, c('ID', 'Gene.Symbol')]
colnames(probe_to_gene) <- c('PROBEID', 'GENE')
probe_to_gene <- probe_to_gene[probe_to_gene$GENE != '' & !grepl('///', probe_to_gene$GENE), ]

core_genes <- c('FASN', 'ACLY', 'SCD', 'SREBF1', 'HMGCR', 'PPARG')
expr_6genes <- matrix(NA, nrow = length(core_genes), ncol = ncol(expr_30929))
rownames(expr_6genes) <- core_genes; colnames(expr_6genes) <- colnames(expr_30929)
for (i in seq_along(core_genes)) {
  gene <- core_genes[i]
  probes <- probe_to_gene$PROBEID[probe_to_gene$GENE == gene]
  probes <- intersect(probes, rownames(expr_30929))
  if (length(probes) == 0) next
  if (length(probes) > 1) {
    best_probe <- names(which.max(rowMeans(expr_30929[probes, , drop = FALSE])))
  } else { best_probe <- probes }
  expr_6genes[gene, ] <- expr_30929[best_probe, ]
}
saveRDS(expr_6genes, 'raw/bulk/GSE30929_6genes_expr.rds')

# 构建临床信息 + 脂评分
pheno_30929 <- pData(gse30929)
clinical_all <- data.frame(
  sample = pheno_30929$geo_accession,
  subtype = pheno_30929$`subtype:ch1`,
  OS_time = as.numeric(gsub('tt.drfs: ', '', pheno_30929$characteristics_ch1.5)),
  OS_status = ifelse(grepl('TRUE', pheno_30929$characteristics_ch1.4), 1, 0),
  stringsAsFactors = FALSE
)
pca_res <- prcomp(t(expr_6genes), center = TRUE, scale. = TRUE)
clinical_all$LipidScore <- pca_res$x[match(clinical_all$sample, colnames(expr_6genes)), 1]
clinical_all$ScoreGroup <- ifelse(clinical_all$LipidScore > median(clinical_all$LipidScore), 'High', 'Low')
saveRDS(clinical_all, 'raw/bulk/GSE30929_clinical_all.rds')

# Cox 回归
cox_all <- coxph(Surv(OS_time, OS_status) ~ ScoreGroup, data = clinical_all)
summary(cox_all)
saveRDS(cox_all, 'raw/bulk/GSE30929_cox_all.rds')

# ---- 2. Bootstrap 中介分析 ----
tcga_expr <- readRDS('raw/TCGA_6genes_expr.rds')
tcga_full <- readRDS('raw/TCGA_LPS_real_log2cpm.rds')
tcga_score <- read.csv('raw/TCGA_Score.csv')
tcga_clin <- read.csv('raw/TCGA_clinical_real.csv')

med_df <- data.frame(
  Sample = colnames(tcga_expr),
  LipidScore = tcga_score$Score[match(colnames(tcga_expr), tcga_score$Sample)],
  PPARG = as.numeric(tcga_expr['PPARG', ]),
  CD36 = as.numeric(tcga_full['ENSG00000135218.19', colnames(tcga_expr)]),
  stringsAsFactors = FALSE
)
med_df$Sample_short <- substr(med_df$Sample, 1, 15)
tcga_clin$Sample_short <- substr(tcga_clin$sample, 1, 15)
med_df <- merge(med_df, tcga_clin[, c('Sample_short', 'OS_time', 'OS_status')], by = 'Sample_short')
med_df <- med_df[!duplicated(med_df$Sample_short), ]

bootstrap_mediation <- function(data, mediator_name, n_boot = 1000) {
  set.seed(123); n <- nrow(data)
  a_model <- lm(as.formula(paste(mediator_name, '~ LipidScore')), data = data)
  a_orig <- coef(a_model)[2]
  cprime_model <- coxph(Surv(OS_time, OS_status) ~ LipidScore + get(mediator_name), data = data)
  b_orig <- coef(cprime_model)[2]
  indirect_orig <- a_orig * b_orig
  boot_indirect <- numeric(n_boot)
  for (i in 1:n_boot) {
    idx <- sample(1:n, n, replace = TRUE); boot_data <- data[idx, ]
    a_boot <- tryCatch(coef(lm(as.formula(paste(mediator_name, '~ LipidScore')), data = boot_data))[2], error = function(e) NA)
    b_boot <- tryCatch(coef(coxph(Surv(OS_time, OS_status) ~ LipidScore + get(mediator_name), data = boot_data))[2], error = function(e) NA)
    if (!is.na(a_boot) && !is.na(b_boot)) boot_indirect[i] <- a_boot * b_boot else boot_indirect[i] <- NA
  }
  boot_indirect <- boot_indirect[!is.na(boot_indirect)]
  ci_low <- quantile(boot_indirect, 0.025); ci_high <- quantile(boot_indirect, 0.975)
  p_boot <- ifelse(ci_low * ci_high > 0, 0.001, min(mean(boot_indirect <= 0), mean(boot_indirect >= 0)) * 2)
  list(indirect = indirect_orig, ci = c(ci_low, ci_high), p = p_boot)
}

pparg_boot <- bootstrap_mediation(med_df, 'PPARG', 1000)
cd36_boot <- bootstrap_mediation(med_df, 'CD36', 1000)

sink('results/tables/extended/bootstrap_mediation_results.txt')
cat('Bootstrap Mediation Analysis (1000 resamples)\n')
cat(sprintf('PPARG: Indirect=%.4f, 95%%CI=[%.4f, %.4f], P=%.4f\n', pparg_boot$indirect, pparg_boot$ci[1], pparg_boot$ci[2], pparg_boot$p))
cat(sprintf('CD36: Indirect=%.4f, 95%%CI=[%.4f, %.4f], P=%.4f\n', cd36_boot$indirect, cd36_boot$ci[1], cd36_boot$ci[2], cd36_boot$p))
sink()

# ---- 3. MDM2 拷贝数分析 ----
cna <- read.delim('raw/cnv_mutation/sarc_tcga/sarc_tcga_pub/data_cna.txt', check.names = FALSE)
cna_t <- as.data.frame(t(cna[, -c(1,2)]))
colnames(cna_t) <- cna[, 1]
cna_t$Sample <- gsub('\\.', '-', rownames(cna_t))
tcga_score$Sample_short3 <- sapply(strsplit(tcga_score$Sample, '-'), function(x) paste(x[1:3], collapse = '-'))
cna_t$Sample_short3 <- sapply(strsplit(cna_t$Sample, '-'), function(x) paste(x[1:3], collapse = '-'))
cna_merge <- merge(tcga_score, cna_t[, c('Sample_short3', 'MDM2', 'CDK4')], by = 'Sample_short3')
cor_mdm2 <- cor.test(cna_merge$Score, cna_merge$MDM2, method = 'spearman')
cat(sprintf('\nMDM2 copy-number vs LipidScore: R=%.3f, P=%.4f\n', cor_mdm2$estimate, cor_mdm2$p.value))
write.csv(cna_merge, 'results/tables/extended/MDM2_LipidScore_association.csv', row.names = FALSE)

cat('\nStage 2 Advanced 完成\n')
