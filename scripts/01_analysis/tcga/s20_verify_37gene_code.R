# S20：验证用户提供的 37 基因 PCA 代码，复现的是 B 还是 C
sink(file.path(Sys.getenv("LPS_WORK_ROOT", "."), "_work", "s20_verify_out.txt"))

setwd(file.path(Sys.getenv("LPS_DATA_ROOT", "."), "1"))

lipid_genes_37 <- c(
  "ACSL1", "ACSL3", "ACSL4", "ACSL5", "ACACA", "ACACB", "ACADM",
  "ACADS", "ACADVL", "ACLY", "ACOX1", "ACSS2", "ALDH1B1", "ALDH2",
  "ALDH3A2", "ALDH9A1", "ANGPTL4", "APOA1", "APOA2", "APOA4",
  "APOB", "APOC1", "APOC2", "APOC3", "APOC4", "APOE", "CETP",
  "CPT1A", "CPT2", "FABP1", "FABP2", "FABP3", "FABP4", "FABP5",
  "FABP6", "FABP7", "FASN"
)
cat("37 基因列表长度:", length(lipid_genes_37), "\n")

log2cpm <- readRDS("data/raw/TCGA_LPS_real_log2cpm.rds")
cat("log2cpm dim:", dim(log2cpm), "\n")

suppressMessages({
  library(AnnotationDbi)
  library(org.Hs.eg.db)
})

ensembl_ids <- gsub("\\.[0-9]+$", "", rownames(log2cpm))
gene_symbols <- mapIds(org.Hs.eg.db, keys = ensembl_ids,
                       column = "SYMBOL", keytype = "ENSEMBL",
                       multiVals = "first")

keep <- !is.na(gene_symbols)
log2cpm_filtered <- log2cpm[keep, ]
rownames(log2cpm_filtered) <- gene_symbols[keep]

log2cpm_clean <- aggregate(log2cpm_filtered,
                           by = list(rownames(log2cpm_filtered)),
                           FUN = mean)
rownames(log2cpm_clean) <- log2cpm_clean$Group.1
log2cpm_clean$Group.1 <- NULL
log2cpm_clean <- as.matrix(log2cpm_clean)

available_genes <- lipid_genes_37[lipid_genes_37 %in% rownames(log2cpm_clean)]
cat("实际可用基因数:", length(available_genes), "/ 37\n")
cat("缺失基因:", setdiff(lipid_genes_37, rownames(log2cpm_clean)), "\n")

expr_matrix <- t(log2cpm_clean[available_genes, ])
cat("expr_matrix dim (样本 x 基因):", dim(expr_matrix), "\n")

pca_result <- prcomp(expr_matrix, scale. = TRUE, center = TRUE)
cat("PC1 方差占比:", round(summary(pca_result)$importance[2,1], 4), "\n")
lipid_score_new <- pca_result$x[, 1]

# 与 B 评分（37genes csv）比对
b_file <- "data/results/02_feature_scores_37genes.csv"
b <- read.csv(b_file)
cat("\n=== B 文件（02_feature_scores_37genes.csv）===\n")
cat("B 行数:", nrow(b), "| B 列名:", paste(colnames(b), collapse=", "), "\n")
cat("B 值域:", range(b[[2]], na.rm=TRUE), "\n")
cat("B 前3行 sample:", head(b[[1]], 3), "\n")

# 对齐样本名
# 新评分 rownames 是 TCGA 样本（可能带 .01A 等）；B 的 sample 列需核对
cat("\n新评分 rownames 前3:", head(rownames(expr_matrix), 3), "\n")
cat("B 的 sample 前3:", head(b[[1]], 3), "\n")

# 尝试匹配：按 rownames 或 sample
new_score_vec <- setNames(lipid_score_new, rownames(expr_matrix))
b_vec <- setNames(b[[2]], b[[1]])

common <- intersect(names(new_score_vec), names(b_vec))
cat("\n精确名匹配样本数:", length(common), "\n")
if (length(common) == 0) {
  # 尝试截断到 12 位 TCGA barcode
  short_new <- substr(names(new_score_vec), 1, 12)
  short_b <- substr(names(b_vec), 1, 12)
  common <- intersect(short_new, short_b)
  cat("截断12位匹配样本数:", length(common), "\n")
  if (length(common) > 0) {
    m_new <- new_score_vec[match(common, short_new)]
    m_b <- b_vec[match(common, short_b)]
  }
} else {
  m_new <- new_score_vec[common]
  m_b <- b_vec[common]
}

if (length(common) > 1) {
  r <- cor(m_new, m_b, method="pearson")
  cat("新评分 vs B 的 Pearson r:", round(r, 6), "\n")
  # 线性回归
  fit <- lm(m_b ~ m_new)
  cat("B = ", round(coef(fit)[2],6), " * new + ", round(coef(fit)[1],6), "\n")
  cat("R^2:", round(summary(fit)$r.squared, 6), "\n")
  cat("残差 max abs:", round(max(abs(resid(fit))), 6), "\n")
}

# 与 C 评分比对
c_file <- "data/results/02_feature_scores.csv"
c_rows <- read.csv(c_file, check.names=FALSE)
cat("\n=== C 文件（02_feature_scores.csv）===\n")
cat("C 行数:", nrow(c_rows), "| C 列名:", paste(colnames(c_rows), collapse=", "), "\n")
c_vec <- setNames(c_rows[["Lipid_Metabolism_Score"]], c_rows[[1]])
cat("C 值域:", range(c_vec, na.rm=TRUE), "\n")

# 新评分只在 TCGA 50 样本上；C 里 TCGA 样本也是这些
common_c <- intersect(names(new_score_vec), names(c_vec))
cat("新评分 vs C 的精确名匹配:", length(common_c), "\n")
if (length(common_c) == 0) {
  short_new <- substr(names(new_score_vec), 1, 12)
  short_c <- substr(names(c_vec), 1, 12)
  common_c <- intersect(short_new, short_c)
  cat("截断12位匹配:", length(common_c), "\n")
  if (length(common_c) > 0) {
    m_new2 <- new_score_vec[match(common_c, short_new)]
    m_c <- c_vec[match(common_c, short_c)]
  }
} else {
  m_new2 <- new_score_vec[common_c]
  m_c <- c_vec[common_c]
}
if (length(common_c) > 1) {
  r2 <- cor(m_new2, m_c, method="pearson")
  cat("新评分 vs C 的 Pearson r:", round(r2, 6), "\n")
  fit2 <- lm(m_c ~ m_new2)
  cat("C = ", round(coef(fit2)[2],6), " * new + ", round(coef(fit2)[1],6), "\n")
  cat("R^2:", round(summary(fit2)$r.squared, 6), "\n")
  cat("残差 max abs:", round(max(abs(resid(fit2))), 6), "\n")
}

sink()
cat("DONE\n")
