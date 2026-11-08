# ============================================================
# Fig S3G · 可塑性：肿瘤细胞室拟时序（slingshot，WD -> DD；原 Fig 6E 上栏）
#   V1 嵌入（tumour cells 36,592）按主谱系伪时序着色
#   ★ 叠加 WD / DD 两组的质心标注
#   ★ 不承载评分方向主张（患者级 5/5 不显著，MDE|rho| = 0.634 -> 不可判定）
#   ★ V3（PCA 空间）符号翻转 -> 任何「评分/SCD 沿轨迹」的表述一律禁引
#   数据：figS6_data/slingshot_pseudotime_V1.csv
# ============================================================
suppressWarnings(suppressMessages({library(ggplot2); library(patchwork); library(scattermore); library(scales)}))
source(file.path(file.path(Sys.getenv("LPS_FIG_ROOT", "."), "utils"), "themes.R"))
source(file.path(file.path(Sys.getenv("LPS_FIG_ROOT", "."), "utils"), "figS6_common.R"))

pt <- read.csv(file.path(FIG6_DATA, "slingshot_pseudotime_V1.csv"),
               stringsAsFactors = FALSE, fileEncoding = "UTF-8-BOM")
pt <- pt[!is.na(pt$pt_L1) & pt$pt_L1 != "NA", ]
pt$pt_L1 <- as.numeric(pt$pt_L1)
stopifnot(nrow(pt) == 22892)
cat(sprintf("S3G: cells=%d | pt range [%.2f, %.2f]\n", nrow(pt), min(pt$pt_L1), max(pt$pt_L1)))
pt <- pt[order(pt$pt_L1), ]

XL <- range(pt$umap1); YL <- range(pt$umap2)
ctr <- aggregate(cbind(umap1, umap2) ~ tissue, data = pt, FUN = mean)
cw <- ctr[ctr$tissue == "DDLPS-WD", ]
cd <- ctr[ctr$tissue == "DDLPS-DD", ]
stopifnot(nrow(cw) == 1, nrow(cd) == 1)
cat(sprintf("   centroid WD=(%.2f, %.2f)  DD=(%.2f, %.2f)\n",
            cw$umap1, cw$umap2, cd$umap1, cd$umap2))

pE1 <- ggplot(pt, aes(umap1, umap2, colour = pt_L1)) +
  geom_point(size = 0.24, alpha = 0.55, shape = 16, stroke = 0) +
  scale_colour_gradientn(colours = SCORE_GRAD, name = "Pseudotime\n(L1)",
                         guide = guide_colourbar(barwidth = unit(5.5, "pt"),
                                                 barheight = unit(52, "pt"),
                                                 ticks.colour = "black",
                                                 frame.colour = "black",
                                                 title.position = "top")) +
  annotate("point", x = cw$umap1, y = cw$umap2, shape = 21, size = 2.4,
           fill = "white", colour = "grey10", stroke = 0.6) +
  annotate("point", x = cd$umap1, y = cd$umap2, shape = 21, size = 2.4,
           fill = "grey20", colour = "grey10", stroke = 0.6) +
  annotate("text", x = cw$umap1, y = cw$umap2 + 0.75, label = "WD",
           size = 2.4, fontface = "bold", colour = "grey10") +
  annotate("text", x = cd$umap1, y = cd$umap2 - 0.95, label = "DD",
           size = 2.4, fontface = "bold", colour = "grey10") +
  labs(x = "UMAP 1", y = "UMAP 2", title = "Pseudotime (lineage 1)") +
  coord_fixed(xlim = XL, ylim = YL, expand = FALSE) +
  theme_publication(base_size = 7.5) +
  theme(legend.position = "right",
        legend.title = element_text(size = 6.2),
        legend.text = element_text(size = 6.2),
        plot.title = element_text(size = 7.2),
        plot.margin = margin(3, 2, 2, 3))

save6(pE1, "figS3g")
