# ============================================================
# Figure SX: GSE30929 六基因表达热图（按评分分组）
# ============================================================
library(pheatmap)
library(RColorBrewer)

setwd("C:/Users/Tao25/Documents/Liposarcoma/2")

# 读取数据
expr_6genes <- readRDS("raw/bulk/GSE30929_6genes_expr.rds")
clinical_all <- readRDS("raw/bulk/GSE30929_clinical_all.rds")

# 按脂评分排序
clinical_all <- clinical_all[order(clinical_all$LipidScore), ]
expr_ordered <- expr_6genes[, clinical_all$sample]

# 注释条
annotation_col <- data.frame(
  ScoreGroup = clinical_all$ScoreGroup,
  Subtype = clinical_all$subtype,
  row.names = clinical_all$sample
)

ann_colors <- list(
  ScoreGroup = c(High = "#E41A1C", Low = "#377EB8"),
  Subtype = c(
    dedifferentiated = "#E41A1C",
    pleomorphic = "#FF7F00",
    `well-differentiated` = "#4DAF4A",
    `myxoid/round cell` = "#984EA3",
    myxoid = "#377EB8"
  )
)

# 绘制热图
cairo_pdf("results/figures/extended/FigSX_GSE30929_heatmap.pdf", 
          width = 16, height = 6, family = "Arial")
pheatmap(expr_ordered,
         scale = "row",
         cluster_cols = FALSE,
         cluster_rows = TRUE,
         show_colnames = FALSE,
         annotation_col = annotation_col,
         annotation_colors = ann_colors,
         color = colorRampPalette(rev(brewer.pal(11, "RdBu")))(100),
         main = "GSE30929: Six-Gene Expression Heatmap (ordered by LipidScore)",
         fontsize = 11,
         border_color = NA)
dev.off()

cat("GSE30929 六基因热图已保存\n")

