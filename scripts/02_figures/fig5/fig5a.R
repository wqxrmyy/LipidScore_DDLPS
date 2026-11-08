# ============================================================
# Fig 5A · UMAP（官方嵌入坐标）按细胞类型着色
#   坐标来源：GSE221493 官方 normalized barcodes 第 9/10 列（非自算 UMAP）
#     依据：s3b2 实测「官方坐标 L3 纯度 0.8514 > 自算 0.6517」→ Fig5 底图用官方坐标
#   细胞类型：L1（10 类，图例/点图/箱线共用同一顺序与同一配色）
#   ★ 2026-09-22（PI 裁定，见 S25 §9）：标题改 "Cell types (L1)" ——
#     显式标注注释层级，避免读者把 L1 的 `Cycling cells`（细胞周期**状态**，
#     1,989 细胞，与 `Tumor cells` 在 L1 并列）误读为与 T cells 平级的细胞类型。
#     该簇性质须在 caption 注明（增殖状态簇）。
#   10 万点渲染：scattermore 栅格化（矢量 PDF 会达数百 MB）
# ============================================================
suppressWarnings(suppressMessages({library(ggplot2); library(scattermore)}))
source(file.path(Sys.getenv("LPS_FIG_ROOT", "."), "utils", "colors.R"))
source(file.path(Sys.getenv("LPS_FIG_ROOT", "."), "utils", "themes.R"))
source(file.path(Sys.getenv("LPS_FIG_ROOT", "."), "utils", "fig5_common.R"))

um <- read.csv(file.path(FIG5_DATA, "umap_cells.csv.gz"), stringsAsFactors = FALSE)
names(um)[names(um) == "celltype"] <- "L1"
stopifnot(nrow(um) == 101890)

# 画序与图例序分离：
#   图例序 = 类平均评分降序（与 5C/5D 的排序原则一致）
#   绘制序 = 细胞数降序（大类垫底、稀有类浮在最上，否则小簇被大簇盖住）
ag  <- aggregate(score_ams ~ L1, data = um, FUN = mean)
ord <- ag$L1[order(-ag$score_ams)]
cnt <- sort(table(um$L1), decreasing = TRUE)
um$L1 <- factor(um$L1, levels = ord)
um <- um[order(match(um$L1, names(cnt))), ]

XL <- range(um$umap1); YL <- range(um$umap2)
stopifnot(all(names(L1_COL) %in% levels(um$L1)))

p <- ggplot(um, aes(umap1, umap2, colour = L1)) +
  geom_scattermore(pointsize = 0.78, alpha = 0.92, pixels = c(1400, 900)) +
  scale_colour_manual(values = L1_COL, name = NULL) +
  guides(colour = guide_legend(override.aes = list(size = 1.6, alpha = 1),
                               keyheight = unit(8.5, "pt"), keywidth = unit(8.5, "pt"))) +
  labs(x = "UMAP 1", y = "UMAP 2", title = "Cell types (L1)") +
  coord_fixed(xlim = XL, ylim = YL, expand = FALSE) +
  theme_publication() +
  theme(legend.position = "right",
        legend.text = element_text(size = 6.2),
        legend.key.spacing.y = unit(0.5, "pt"),
        plot.margin = margin(3, 2, 3, 3))

export_fig(p, file.path(FIG5_OUT, "fig5a.pdf"), w = 3.2, h = 2.1)
cat(sprintf("fig5a: %d 细胞 | %d 类 | 坐标域 x[%.2f,%.2f] y[%.2f,%.2f]\n",
            nrow(um), nlevels(um$L1), XL[1], XL[2], YL[1], YL[2]))
