# ============================================================
# Fig 3 数据准备：GSE221492 样本标签解析 + 6 基因 + 身份基因提取
# 产出：gse2_df（样本×{score, label, patient, group}）、gse_genes（基因表达矩阵）
# ============================================================
suppressWarnings(suppressMessages(library(jsonlite)))
GSE_DIR <- file.path(Sys.getenv("LPS_WORK_ROOT", "."), "_work", "gse221492")
OUT <- file.path(Sys.getenv("LPS_FIG_ROOT", "."))

# 表达矩阵（log2 FPKM）
gse_expr <- read.delim(file.path(GSE_DIR, "GSE221492_bulk_LPS_atlas_log2_FPKM_counts_matrix.tsv.gz"),
                       check.names = FALSE, row.names = 1)
# 列名形如 P1_DDLPS_DD（三段）或 P5_WDLPS（两段）
samples <- colnames(gse_expr)
gse_meta <- data.frame(Sample = samples, stringsAsFactors = FALSE)
# 病人 ID = 第一段；标签 = 其余段合并
gse_meta$Patient <- sub("_.*$", "", samples)
gse_meta$Label <- sub("^[^_]+_", "", samples)
# 分组归类（把 DDLPS_ 前缀归一化）
lab_clean <- gsub("^DDLPS_", "", gse_meta$Label)
lab_map <- c(DD = "DD", WD = "WD", sWD = "sWD",
             WDLPS = "WDLPS", lipoma = "Lipoma", peritumoral_fat = "Fat", Hibernoma = "Hibernoma")
gse_meta$Group <- unname(lab_map[lab_clean])
stopifnot(!any(is.na(gse_meta$Group)))

# 基因名（行名是 Ensembl 或 symbol？检查）
rn <- rownames(gse_expr)
if (all(startsWith(head(rn, 5), "ENSG"))) {
  # 需要 symbol 映射：用与 TCGA 相同的 Ensembl ID（无版本号）
  ENS_MAP <- c(SREBF1 = "ENSG00000072310", SCD = "ENSG00000099194", HMGCR = "ENSG00000113161",
               ACLY = "ENSG00000131473", PPARG = "ENSG00000132170", FASN = "ENSG00000169710",
               CD36 = "ENSG00000135218", FABP4 = "ENSG00000091131", PLIN1 = "ENSG00000166819",
               LPL = "ENSG00000175592", ADIPOQ = "ENSG00000181092", CPT1A = "ENSG00000110090",
               CD8A = "ENSG00000153563", MDM2 = "ENSG00000135679", CDK4 = "ENSG00000135446")
  rn_clean <- sub("\\.[0-9]+$", "", rn)
  idx <- match(ENS_MAP, rn_clean)
  hit <- !is.na(idx)
  gse_genes <- gse_expr[idx[hit], , drop = FALSE]
  rownames(gse_genes) <- names(ENS_MAP)[hit]
} else {
  # 已经是 symbol
  want <- c("SREBF1", "SCD", "HMGCR", "ACLY", "PPARG", "FASN", "CD36", "FABP4",
            "PLIN1", "LPL", "ADIPOQ", "CPT1A", "CD8A", "MDM2", "CDK4")
  hit <- want %in% rn
  gse_genes <- gse_expr[want[hit], , drop = FALSE]
}

# 评分（6 基因 z 分均值，与 A 评分同一定义）
six <- c("SREBF1", "SCD", "HMGCR", "ACLY", "PPARG", "FASN")
six_hit <- six %in% rownames(gse_genes)
z <- t(scale(t(gse_genes[six_hit, , drop = FALSE])))
gse_meta$Score <- as.numeric(colMeans(z))

save(gse_expr, gse_meta, gse_genes, file = file.path(OUT, "p5_gse221492.RData"))
cat("GSE221492 数据准备完成\n")
cat("样本:", nrow(gse_meta), "| 分组:", paste(names(table(gse_meta$Group)), table(gse_meta$Group), sep = "=", collapse = ", "), "\n")
cat("基因矩阵:", nrow(gse_genes), "基因\n")
cat("配对患者数（同时有 DD 和 WD）:", sum(table(gse_meta$Patient[gse_meta$Group %in% c("DD", "WD")]) == 2), "\n")
