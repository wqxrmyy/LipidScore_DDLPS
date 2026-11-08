# ============================================================
# Fig 5B · 同坐标 UMAP 按评分 A 着色
#   与 5A 严格同坐标、同 extent（coord_fixed 同 xlim/ylim）→ 可并排对读
#   色阶：浅灰→橙→深红（非 viridis：评分高度集中于 0，viridis 会把「≈0 的大多数」压成暗色）
#   两端按 5%/95% 分位截断（防离群压平色阶；1%/99% 饱和不足，单细胞层面难辨认）
#   ★ 截断比例写入图注（不隐藏处理）
# ============================================================
suppressWarnings(suppressMessages({library(ggplot2); library(scattermore); library(viridisLite)}))
source(file.path(Sys.getenv("LPS_FIG_ROOT", "."), "utils", "colors.R"))
source(file.path(Sys.getenv("LPS_FIG_ROOT", "."), "utils", "themes.R"))
source(file.path(Sys.getenv("LPS_FIG_ROOT", "."), "utils", "fig5_common.R"))

um <- read.csv(file.path(FIG5_DATA, "umap_cells.csv.gz"), stringsAsFactors = FALSE)
stopifnot(nrow(um) == 101890)

# 画序：低分先画、高分后画（高分细胞浮在最上层，不被盖住）
um <- um[order(um$score_ams), ]
XL <- range(um$umap1); YL <- range(um$umap2)
q  <- quantile(um$score_ams, c(0.05, 0.95))   # 5%/95%：比 1%/99% 更早饱和，单细胞层面更易辨认
cat(sprintf("评分 A: 全域 [%.3f, %.3f] | 1%%-99%% 分位 [%.3f, %.3f]\n",
            min(um$score_ams), max(um$score_ams), q[1], q[2]))

p <- ggplot(um, aes(umap1, umap2, colour = score_ams)) +
  geom_scattermore(pointsize = 0.78, alpha = 0.95, pixels = c(1400, 900)) +
  # 色阶：浅灰（低）→ 橙 → 深红（高）
  #   不用 viridis：评分高度集中于 0，viridis 会把「≈0 的大多数」压成深紫黑，
  #   全场一片暗色、看不出高分区；灰-橙-红让零附近的细胞保持中性、高分细胞跳出
  scale_colour_gradientn(
    colours = c("#E5E5E5", "#FED976", "#FC8D59", "#99000D"),
    limits = c(q[1], q[2]), oob = scales::squish,
    name = "Score A",
    guide = guide_colourbar(barwidth = unit(6, "pt"), barheight = unit(62, "pt"),
                            ticks.colour = "black", frame.colour = "black",
                            title.position = "top")) +
  labs(x = "UMAP 1", y = "UMAP 2", title = "Lipid score (A)") +
  coord_fixed(xlim = XL, ylim = YL, expand = FALSE) +
  theme_publication() +
  theme(legend.position = "right",
        legend.title = element_text(size = 6.2),
        legend.text = element_text(size = 6.2),
        plot.margin = margin(3, 2, 3, 3))

export_fig(p, file.path(FIG5_OUT, "fig5b.pdf"), w = 3.2, h = 2.1)
cat(sprintf("fig5b: 色阶截断于 5%%/95%% 分位 (%.4f, %.4f)\n", q[1], q[2]))
