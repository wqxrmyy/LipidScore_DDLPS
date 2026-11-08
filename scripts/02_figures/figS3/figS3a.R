# ============================================================
# Fig S3A · 评分 A 在肿瘤细胞的 UMAP 投影（原 Fig 6A 上栏，2026-09-28 重构）
#   灰底 = 非肿瘤细胞，彩色 = 肿瘤细胞按评分 A
#   ★ 与 Fig5B（全细胞评分投影）的区别：此处只上色肿瘤细胞室，
#     用于回答「评分在同一细胞室内是否有异质性」，不是全图谱重画
#   ★ 2026-09-24（B 方案）防「与 5B 重复」误读三改：
#     ① 色标范围统一为全细胞 q05/q95（同 fig5b.R）→ 两图同一把尺子
#     ② 非肿瘤细胞加深为实心深灰
#     ③ 标题写「tumour cells only … same axes as Fig 5B」
#   数据：fig5_data/umap_cells.csv.gz（官方坐标）
# ============================================================
suppressWarnings(suppressMessages({library(ggplot2); library(patchwork); library(scattermore)}))
source(file.path(file.path(Sys.getenv("LPS_FIG_ROOT", "."), "utils"), "themes.R"))
source(file.path(file.path(Sys.getenv("LPS_FIG_ROOT", "."), "utils"), "figS6_common.R"))

um <- read.csv(file.path(FIG6_ROOT, "fig5_data", "umap_cells.csv.gz"), stringsAsFactors = FALSE)
stopifnot(nrow(um) == 101890)
bg <- um[, c("umap1", "umap2")]
fg <- um[um$celltype == "Tumor cells", ]
stopifnot(nrow(fg) == 36592)
fg <- fg[order(fg$score_ams), ]
Q_ALL <- quantile(um$score_ams, c(0.05, 0.95))   # 与 fig5b.R 同口径
q <- Q_ALL
cat(sprintf("S3A: tumour cells=%d | colormap limits (all-cell q05/q95) = %.4f / %.4f\n",
            nrow(fg), q[1], q[2]))
XL <- range(um$umap1); YL <- range(um$umap2)

pA1 <- ggplot() +
  geom_scattermore(data = bg, aes(umap1, umap2), colour = "#9E9E9E",
                   pointsize = 0.62, alpha = 0.85, pixels = c(1300, 900)) +
  geom_scattermore(data = fg, aes(umap1, umap2, colour = score_ams),
                   pointsize = 0.62, alpha = 0.95, pixels = c(1300, 900)) +
  scale_colour_gradientn(colours = SCORE_GRAD, limits = c(q[1], q[2]),
                         oob = scales::squish, name = "Score A",
                         guide = guide_colourbar(barwidth = unit(5.5, "pt"),
                                                 barheight = unit(48, "pt"),
                                                 ticks.colour = "black",
                                                 frame.colour = "black",
                                                 title.position = "top")) +
  labs(x = "UMAP 1", y = "UMAP 2",
       title = "Score A in tumour cells",
       subtitle = paste0("tumour only (n = 36,592); grey = other cells\n",
                         "same axes and colour scale as Fig 5B")) +
  coord_fixed(xlim = XL, ylim = YL, expand = FALSE) +
  theme_publication(base_size = 7.5) +
  theme(legend.position = "right",
        legend.title = element_text(size = 6.4),
        legend.text = element_text(size = 6.2),
        plot.title = element_text(size = 7.5),
        plot.subtitle = element_text(size = 5.9, colour = "grey30", lineheight = 1.05),
        plot.margin = margin(3, 2, 2, 3))

save6(pA1, "figS3a")
