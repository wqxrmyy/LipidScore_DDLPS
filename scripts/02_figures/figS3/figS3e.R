# ============================================================
# Fig S3E · 跨细胞室：TAM 端样本级共变（三格散点，原 Fig 6D 上栏）
#   三格散点（n = 27 取样）x = 肿瘤细胞评分 A
#     y = 巨噬 PPARG（唯一显著）/ M2 样打分 / 巨噬 CD36（均不显著）
#   ★ 三格同轴同范围，便于直接目视比较「显著 vs 不显著」
#   数据：figS6_data/figS6d_sample_coupling.csv + figS6d_coupling_corr.csv
# ============================================================
suppressWarnings(suppressMessages({library(ggplot2); library(patchwork); library(scales)}))
source(file.path(file.path(Sys.getenv("LPS_FIG_ROOT", "."), "utils"), "themes.R"))
source(file.path(file.path(Sys.getenv("LPS_FIG_ROOT", "."), "utils"), "figS6_common.R"))

s  <- read6("figS6d_sample_coupling.csv")
cc <- read6("figS6d_coupling_corr.csv")
stopifnot(nrow(s) == 27, all(cc$n == 27))
cat(sprintf("S3E: samples=%d | patients=%d\n", nrow(s), length(unique(s$patient))))

GET <- function(y) {
  r <- cc[cc$x == "tumor_score" & cc$y == y, ]
  list(rho = r$rho[1], P = r$p[1])
}
YS <- list(list(col = "mac_PPARG",  ttl = "Macrophage PPARG", key = "mac_PPARG"),
           list(col = "mac_m2like", ttl = "M2-like score",    key = "mac_m2like"),
           list(col = "mac_CD36",   ttl = "Macrophage CD36",  key = "mac_CD36"))
long <- do.call(rbind, lapply(YS, function(y) {
  g <- GET(y$key)
  data.frame(panel = y$ttl, x = s$tumor_score, yv = s[[y$col]],
             rho = g$rho, P = g$P, stringsAsFactors = FALSE)
}))
lab <- unique(long[, c("panel", "rho", "P")])
lab$txt <- sprintf("rho = %+.3f\n%s", lab$rho, fmt_p(lab$P))
lab$panel <- factor(lab$panel, levels = sapply(YS, function(y) y$ttl))
lab$sig <- ifelse(lab$P < 0.05, "P < 0.05", "not significant")
lab$x <- Inf; lab$y <- Inf
long$panel <- factor(long$panel, levels = sapply(YS, function(y) y$ttl))

pD1 <- ggplot(long, aes(x, yv)) +
  geom_hline(yintercept = 0, linewidth = 0.3, colour = "grey80") +
  geom_point(shape = 21, size = 1.5, fill = "grey95", colour = "grey30", stroke = 0.45) +
  geom_smooth(method = "lm", se = TRUE, linewidth = 0.55, colour = SIG_COL,
              fill = SIG_COL, alpha = 0.14) +
  geom_label(data = lab, aes(x = x, y = y, label = txt, colour = sig),
             hjust = 1.02, vjust = 1.06, size = 2.2, lineheight = 0.92,
             fill = "white", alpha = 0.88, label.size = 0,
             label.padding = unit(0.6, "pt")) +
  facet_wrap(~panel, ncol = 1, scales = "free_y") +
  scale_colour_manual(values = c("P < 0.05" = SIG_COL, "not significant" = "grey35"),
                      guide = "none") +
  labs(x = "Tumour-cell score A", y = "Macrophage readout",
       title = wrap_txt("Sample coupling (n = 27)", 28)) +
  theme_publication(base_size = 7.5) +
  theme(axis.text = element_text(size = 5.8),
        axis.title = element_text(size = 6.6),
        strip.background = element_blank(),
        strip.text = element_text(size = 6.4, face = "bold"),
        plot.title = element_text(size = 7.2),
        plot.margin = margin(3, 3, 2, 3))

save6(pD1, "figS3e")

## 关键读数回显
for (i in seq_along(YS)) cat(sprintf("   %-18s rho=%+.4f P=%.4g\n",
                                     YS[[i]]$ttl, GET(YS[[i]]$key)$rho, GET(YS[[i]]$key)$P))
