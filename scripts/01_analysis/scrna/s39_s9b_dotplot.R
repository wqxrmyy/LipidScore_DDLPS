# ============================================================
# Step 9b · Fig5 面板 C（点图）的表达层准备
#   来源：s39_seurat.s6.rds（归一化 data 层；scale.data 已在 Step 6 释放）
#   产出：fig5_data/dotplot_genes.csv —— L1 细胞类型 × 6 分量基因的
#         pct_exp（表达细胞比例，data > 0）与 mean_exp（expm1 均值，即非对数尺度平均表达）
#   ⚠️ 只用 6 个 bulk 同款基因（SREBF1/SCD/HMGCR/ACLY/PPARG/FASN），不额外增删
# ============================================================
suppressWarnings(suppressMessages({library(Seurat); library(Matrix); library(jsonlite)}))
W    <- file.path(Sys.getenv("LPS_WORK_ROOT", "."), "_work")
OUTD <- file.path(W, "p5_figures", "fig5_data")
dir.create(OUTD, recursive = TRUE, showWarnings = FALSE)

GENES <- c("SREBF1", "SCD", "HMGCR", "ACLY", "PPARG", "FASN")
REP <- character(); A <- function(...) { l <- paste(...); cat(l, "\n"); REP <<- c(REP, l) }

A("=== Step 9b · 点图表达层 ===")
obj <- readRDS(file.path(W, "s39_seurat.s6.rds"))
A("细胞数:", ncol(obj), "| 基因数:", nrow(obj))

miss <- setdiff(GENES, rownames(obj))
stopifnot(length(miss) == 0)
A("6 基因全部在库 ✓")

# L1 细胞类型列名探测（防凭记忆写错列名）
md <- obj@meta.data
cand <- c("L1", "celltype", "cell_type", "L1_annot")
hit <- cand[cand %in% colnames(md)]
A("meta 中可用的 L1 列:", paste(hit, collapse = ", "))
stopifnot(length(hit) >= 1)
l1col <- hit[1]
grp <- as.character(md[[l1col]])
A("L1 类数:", length(unique(grp)))

# 归一化表达（data 层）
mat <- GetAssayData(obj, assay = "RNA", layer = "data")[GENES, , drop = FALSE]
mat <- as(mat, "dgCMatrix")
A("抽取矩阵:", nrow(mat), "×", ncol(mat), "| 非零元:", length(mat@x))

rows <- list()
for (g in GENES) {
  v <- mat[g, ]
  for (cl in sort(unique(grp))) {
    sel <- grp == cl
    x <- v[sel]
    rows[[length(rows) + 1]] <- data.frame(
      gene = g, celltype = cl,
      n_cells = sum(sel),
      pct_exp = 100 * mean(x > 0),
      mean_exp = mean(expm1(x)),          # 非对数尺度平均表达（Seurat DotPlot 口径）
      mean_log = mean(x),                 # 对数尺度均值（备查）
      stringsAsFactors = FALSE)
  }
}
df <- do.call(rbind, rows)

# 与逐细胞评分表做一次类数交叉核对（防 L1 口径与 Python 侧漂移）
cs <- read.csv(file.path(W, "s6", "cell_scores.csv.gz"), stringsAsFactors = FALSE)
stopifnot(nrow(cs) == ncol(obj))
l1_py <- sort(unique(cs$L1)); l1_r <- sort(unique(grp))
A("L1 类别（R）:", paste(l1_r, collapse = " | "))
A("与 Python 侧 cell_scores$L1 是否同一集合:", identical(l1_py, l1_r))
stopifnot(identical(l1_py, l1_r))
stopifnot(identical(sort(cs$barcode), sort(colnames(obj))))

write.csv(df, file.path(OUTD, "dotplot_genes.csv"), row.names = FALSE,
          fileEncoding = "utf-8")
A("→ dotplot_genes.csv  ", nrow(df), "行 =", length(GENES), "基因 ×", length(l1_r), "类")

# 关键读数（供图注核对）
sub <- df[df$celltype == "Normal adipocytes", ]
A("Normal adipocytes 的 6 基因 pct_exp:")
for (i in seq_len(nrow(sub))) A(sprintf("   %-7s %5.1f%%  mean_exp=%.3f", sub$gene[i], sub$pct_exp[i], sub$mean_exp[i]))
ord <- sapply(GENES, function(g) df$pct_exp[df$gene == g][which.max(df$pct_exp[df$gene == g])])
A("各基因 pct_exp 最高的细胞类型:")
for (i in seq_along(GENES)) A(sprintf("   %-7s -> %s", GENES[i],
   df$celltype[df$gene == GENES[i]][which.max(df$pct_exp[df$gene == GENES[i]])]))

writeLines(REP, file.path(W, "s39_s9b_report.txt"))
A("DONE")
