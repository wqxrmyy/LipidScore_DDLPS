# ============================================================
# Step 9c · L2 层基因表达（面板 C/D 的判读依据 + 关键核查）
#   目的：看清 「评分高分类」与「注释为 Adipo diff 的类」到底表达什么
#   抽 6 个评分基因 + 3 个脂肪身份锚基因（Step 2 已定：FABP4 / CD36 / PPARG）
#   ⚠️ 本脚本只做描述，不做任何过滤/取舍；输出供人工判读
# ============================================================
suppressWarnings(suppressMessages({library(Seurat); library(Matrix)}))
W <- file.path(Sys.getenv("LPS_WORK_ROOT", "."), "_work")
OUTD <- file.path(W, "p5_figures", "fig5_data")
GENES <- c("SREBF1", "SCD", "HMGCR", "ACLY", "PPARG", "FASN",   # 评分 6 基因
           "FABP4", "CD36", "ADIPOQ", "PLIN1", "LPL", "CEBPA",  # 脂肪身份 / 分化
           "COL1A1", "MDM2", "CDK4", "PTPRC", "PECAM1", "ACTA2")# 谱系反证
REP <- character(); A <- function(...) { l <- paste(...); cat(l, "\n"); REP <<- c(REP, l) }

A("=== Step 9c · L2 层基因表达 ===")
obj <- readRDS(file.path(W, "s39_seurat.s6.rds"))
md  <- obj@meta.data
stopifnot("L2" %in% colnames(md))
grp <- as.character(md$L2)
miss <- setdiff(GENES, rownames(obj)); A("缺失基因:", if (length(miss)) paste(miss, collapse=",") else "无")

mat <- as(GetAssayData(obj, assay = "RNA", layer = "data")[GENES, , drop = FALSE], "dgCMatrix")
A("矩阵:", nrow(mat), "×", ncol(mat))

rows <- list()
for (g in GENES) for (cl in sort(unique(grp))) {
  x <- mat[g, grp == cl]
  rows[[length(rows) + 1]] <- data.frame(
    gene = g, L2 = cl, n_cells = length(x),
    pct_exp = 100 * mean(x > 0), mean_exp = mean(expm1(x)),
    stringsAsFactors = FALSE)
}
df <- do.call(rbind, rows)
write.csv(df, file.path(OUTD, "L2_genes.csv"), row.names = FALSE, fileEncoding = "utf-8")

# ---- 判读表 1：各 L2 类在「6 评分基因」上的平均 z（跨类标准化）----------
pv <- tapply(df$mean_exp[df$gene %in% GENES[1:6]],
             list(df$L2[df$gene %in% GENES[1:6]], df$gene[df$gene %in% GENES[1:6]]),
             identity)
pv <- pv[, GENES[1:6]]
z  <- scale(pv)
A("")
A("--- 6 评分基因在 L2 各类的类内 z（跨类标准化）---")
ord <- order(-rowMeans(z[, c("SCD", "FASN", "PPARG"), drop = FALSE]))
print(round(z[ord, ], 2))

# ---- 判读表 2：脂肪身份锚基因（pct_exp%）-------------------------------
pct <- tapply(df$pct_exp, list(df$L2, df$gene), identity)
A("")
A("--- 脂肪身份锚基因 pct_exp(%) ---")
print(round(pct[ord, c("FABP4", "CD36", "PPARG", "ADIPOQ", "PLIN1", "LPL", "SCD", "FASN")], 1))

# ---- 判读表 3：谱系反证基因（pct_exp%）---------------------------------
A("")
A("--- 谱系反证基因 pct_exp(%) ---")
print(round(pct[ord, c("PTPRC", "PECAM1", "ACTA2", "COL1A1", "MDM2", "CDK4")], 1))

# ---- 关键对比：Adipo diff vs Normal adipocytes -------------------------
A("")
A("--- ★ Adipo diff vs Normal adipocytes（mean_exp / pct_exp）---")
sub <- df[df$L2 %in% c("Adipo diff", "Normal adipocytes"), ]
for (g in GENES) {
  a <- sub[sub$gene == g & sub$L2 == "Adipo diff", ]
  b <- sub[sub$gene == g & sub$L2 == "Normal adipocytes", ]
  if (nrow(a) && nrow(b))
    A(sprintf("  %-8s Adipo diff %6.3f (%5.1f%%) | Normal adipocytes %6.3f (%5.1f%%)",
              g, a$mean_exp, a$pct_exp, b$mean_exp, b$pct_exp))
}
writeLines(REP, file.path(W, "s39_s9c_report.txt"))
A("DONE")
