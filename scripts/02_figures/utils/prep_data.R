# ============================================================
# P5 出图 · 数据准备（统一加载，供所有子图脚本 source）
# 输出：一个干净的 RData 到 _work/p5_figures/p5_data.RData
# 基准：A 评分（TCGA_Score.csv），探针 ID 已经 S18 校验
# ============================================================
DATA_RAW <- file.path(Sys.getenv("LPS_DATA_ROOT", "."), "1", "data", "raw")
DATA_RES <- file.path(Sys.getenv("LPS_DATA_ROOT", "."), "1", "data", "results")
OUT <- file.path(Sys.getenv("LPS_FIG_ROOT", "."))

# ---- 1. TCGA 评分 A（6 基因 z 分均值）----
score_A <- read.csv(file.path(DATA_RAW, "TCGA_Score.csv"), stringsAsFactors = FALSE)
rownames(score_A) <- score_A$Sample

# ---- 2. 临床（含 OS + 拷贝数变量）----
clin <- read.csv(file.path(DATA_RES, "03_clinical_outcomes.csv"), stringsAsFactors = FALSE)
rownames(clin) <- clin$Sample_ID

clin_real <- read.csv(file.path(DATA_RAW, "TCGA_clinical_real.csv"), stringsAsFactors = FALSE,
                      check.names = FALSE)
rownames(clin_real) <- clin_real$barcode

# 拷贝数簇等变量（从 clin_real 的 paper_* 列取）
cn_vars <- c("paper_CN.cluster", "paper_CIN", "paper_ploidy", "paper_Genome.doublings",
             "paper_Cancer.DNA.fraction", "paper_Subclona.genome.fraction",
             "paper_Mutation.Load", "paper_avg.CN.12q13.15", "paper_purity",
             "paper_leukocyte.fraction", "paper_mitotic.rate", "paper_FNCLCC.grade")

# ---- 3. 6 基因表达（log2CPM）----
# 从全转录组 log2cpm 提取 6 基因（用正确 Ensembl ID）
log2cpm <- readRDS(file.path(DATA_RAW, "TCGA_LPS_real_log2cpm.rds"))

idx <- setNames(rownames(log2cpm), sub("\\.[0-9]+$", "", rownames(log2cpm)))
ENS6 <- c(SREBF1 = "ENSG00000072310", SCD = "ENSG00000099194", HMGCR = "ENSG00000113161",
          ACLY = "ENSG00000131473", PPARG = "ENSG00000132170", FASN = "ENSG00000169710")
expr6 <- matrix(NA, nrow = 6, ncol = ncol(log2cpm),
                dimnames = list(names(ENS6), colnames(log2cpm)))
for (g in names(ENS6)) {
  row <- idx[ENS6[g]]
  if (!is.na(row)) expr6[g, ] <- as.numeric(log2cpm[row, ])
}

# ---- 4. 其他队列评分（GSE21122 / GSE159659）----
score_GSE21122 <- read.csv(file.path(DATA_RAW, "GSE21122_Score.csv"), stringsAsFactors = FALSE)
score_GSE21122$Sample <- sub("\\.CEL\\.gz$", "", score_GSE21122$Sample)
score_GSE159659 <- read.csv(file.path(DATA_RAW, "GSE159659_Score.csv"), stringsAsFactors = FALSE)

# ---- 5. GSE221492（支柱 1 核心，bulk log2 FPKM）----
gse221492_dir <- file.path(Sys.getenv("LPS_WORK_ROOT", "."), "_work", "gse221492")
# 读取表达矩阵（tsv.gz）
gse_expr <- NULL
gse_meta <- NULL
tryCatch({
  gse_expr <- read.delim(file.path(gse221492_dir, "GSE221492_bulk_LPS_atlas_log2_FPKM_counts_matrix.tsv.gz"),
                         check.names = FALSE, row.names = 1)
}, error = function(e) message("GSE221492 表达矩阵读取失败: ", e$message))

# ---- 保存 ----
save(score_A, clin, clin_real, cn_vars, expr6, ENS6, score_GSE21122, score_GSE159659,
     gse_expr, gse_meta, file = file.path(OUT, "p5_data.RData"))
cat("数据准备完成，已存", file.path(OUT, "p5_data.RData"), "\n")
cat("score_A:", nrow(score_A), "样本\n")
cat("clin:", nrow(clin), "样本\n")
cat("expr6 维度:", dim(expr6), "\n")
if (!is.null(gse_expr)) cat("gse_expr 维度:", dim(gse_expr), "\n") else cat("gse_expr: NULL\n")
