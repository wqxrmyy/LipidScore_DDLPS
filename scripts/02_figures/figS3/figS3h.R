# ============================================================
# Fig S3H · 配对 WD vs DD 患者级平均伪时序连线（原 Fig 6E 下栏）
#   7 对配对患者（WD vs DD）的患者级平均伪时序连线
#   副标题给 mean change 与 Wilcoxon P；★ 须并列「1/3 谱系」与 L3 反向
#   数据：figS6_data/figS6e_pair_delta.csv + figS6e_pair_meta.csv
# ============================================================
suppressWarnings(suppressMessages({library(ggplot2); library(patchwork); library(scattermore); library(scales)}))
source(file.path(file.path(Sys.getenv("LPS_FIG_ROOT", "."), "utils"), "themes.R"))
source(file.path(file.path(Sys.getenv("LPS_FIG_ROOT", "."), "utils"), "figS6_common.R"))

d <- read6("figS6e_pair_delta.csv")
m <- read6("figS6e_pair_meta.csv")
md  <- as.numeric(m$value[m$metric == "mean_delta_DD_minus_WD"])
pw  <- as.numeric(m$value[m$metric == "wilcoxon_P_exact"])
stopifnot(nrow(d) == 7)
lg <- rbind(data.frame(patient = d$patient, tissue = "WD", pt = d$WD),
            data.frame(patient = d$patient, tissue = "DD", pt = d$DD))
lg$x <- ifelse(lg$tissue == "WD", 1, 2)
lg$col <- ifelse(lg$tissue == "WD", "WD", "DD")
cat(sprintf("S3H: n=%d | mean change %+.3f | Wilcoxon P=%.4f\n", nrow(d), md, pw))

pE2 <- ggplot(lg, aes(x, pt)) +
  geom_line(aes(group = patient), colour = "#A93226", linewidth = 0.42, alpha = 0.55) +
  geom_point(aes(fill = col), shape = 21, size = 1.9, stroke = 0.5, colour = "#A93226") +
  scale_fill_manual(values = c(WD = "white", DD = "#A93226"), guide = "none") +
  scale_x_continuous(breaks = c(1, 2), labels = c("WD", "DD"), limits = c(0.62, 2.38)) +
  labs(x = NULL, y = "Mean pseudotime",
       title = "Paired WD vs DD (n = 7)",
       subtitle = wrap_txt(sprintf("1 of 3 lineages, L3 reversed; mean delta %+.2f, P = %.3f",
                                   md, pw), 30)) +
  theme_publication(base_size = 7.5) +
  theme(axis.text.x = element_text(size = 7.0),
        axis.text.y = element_text(size = 6.2),
        plot.title = element_text(size = 7.2),
        plot.subtitle = element_text(size = 6.4, colour = "grey25"),
        plot.margin = margin(3, 3, 2, 3))

save6(pE2, "figS3h")
