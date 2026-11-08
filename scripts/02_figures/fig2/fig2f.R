# ============================================================
# Fig 2F · 结构类 vs 负荷类变量（校正纯度前后对比 分组棒棒糖）
# A 基准 d12 结论：结构类（ploidy/WGD/CancerDNA）校正后增强或转显著；
# 负荷类（CIN/mutation load/subclonal）校正后失效
# ============================================================
suppressWarnings(suppressMessages(library(ggplot2)))
source(file.path(Sys.getenv("LPS_FIG_ROOT", "."), "utils", "colors.R"))
source(file.path(Sys.getenv("LPS_FIG_ROOT", "."), "utils", "themes.R"))

load(file.path(Sys.getenv("LPS_FIG_ROOT", "."), "p5_data.RData"))
OUT <- file.path(Sys.getenv("LPS_FIG_ROOT", "."), "out")

numf <- function(x) suppressWarnings(as.numeric(as.character(x)))
y <- score_A$Score
pur <- numf(clin_real[score_A$Sample, "paper_purity"])

vars <- list(
  "CN cluster"        = "paper_CN.cluster",
  "Ploidy"            = "paper_ploidy",
  "Genome doublings"  = "paper_Genome.doublings",
  "Cancer DNA frac."  = "paper_Cancer.DNA.fraction",
  "CIN"               = "paper_CIN",
  "Mutation load"     = "paper_Mutation.Load",
  "Subclonal frac."   = "paper_Subclona.genome.fraction"
)
struct <- c("CN cluster", "Ploidy", "Genome doublings", "Cancer DNA frac.")

partial_r <- function(x, y, Z) {
  keep <- complete.cases(x, y, Z)
  if (sum(keep) < 10) return(c(NA, NA, sum(keep)))
  Zdf <- as.data.frame(Z[keep, , drop = FALSE])
  names(Zdf) <- paste0("z", seq_len(ncol(Zdf)))
  Zdf$.x <- x[keep]; Zdf$.y <- y[keep]
  rx <- residuals(lm(.x ~ . - .y, data = Zdf)); ry <- residuals(lm(.y ~ . - .x, data = Zdf))
  tt <- cor.test(rx, ry)
  c(unname(tt$estimate), tt$p.value, sum(keep))
}

rows <- NULL
for (vn in names(vars)) {
  xv <- numf(clin_real[score_A$Sample, vars[[vn]]])
  ok <- !is.na(xv)
  if (sum(ok) < 10) next
  t1 <- cor.test(xv[ok], y[ok])
  t2 <- partial_r(xv, y, data.frame(purity = pur))
  rows <- rbind(rows, data.frame(
    Variable = vn, Class = ifelse(vn %in% struct, "Structural", "Burden"),
    r_raw = unname(t1$estimate), p_raw = t1$p.value,
    r_adj = as.numeric(t2[1]), p_adj = as.numeric(t2[2])
  ))
}
# 长表
library(reshape2)
lmelt <- melt(rows[, c("Variable", "Class", "r_raw", "r_adj")], id.vars = c("Variable", "Class"))
lmelt$Adj <- ifelse(lmelt$variable == "r_raw", "Raw", "Purity-adjusted")
lmelt$Variable <- factor(lmelt$Variable, levels = rev(names(vars)))
lmelt <- lmelt[!is.na(lmelt$value), ]

p <- ggplot(lmelt, aes(x = value, y = Variable, color = Adj, shape = Class)) +
  geom_vline(xintercept = 0, linetype = "dotted", color = "grey40") +
  geom_segment(aes(x = 0, xend = value, y = Variable, yend = Variable), color = "grey70", size = 0.5) +
  geom_point(size = 3.2) +
  scale_color_manual(values = c(Raw = "#7F8C8D", "Purity-adjusted" = "#E74C3C")) +
  labs(x = "Pearson r with score (A)", y = "", color = "", shape = "",
       title = "Structural vs burden CN variables",
       # v3：原 81 字符副标题在半幅被裁断（3.25 in / 7 pt 上限约 54 字符）
       subtitle = "Structural retained, burden lost after adjustment") +
  theme_publication() +
  # v3：两个图例（颜色=校正与否，形状=结构/负荷类）横向单排共 4 键，
  #     实测在半幅 3.25 in 下末键 "Structural" 被右缘截断 → 改为上下两排。
  theme(legend.position = "bottom", legend.box = "vertical",
        legend.title = element_blank(),
        legend.text = element_text(size = 6.5),
        legend.key.size = unit(0.22, "cm"),
        legend.spacing.x = unit(0.03, "cm"),
        legend.spacing.y = unit(0.02, "cm"),
        legend.margin = margin(0, 0, 0, 0),
        ## ★★ S39：右距 4 -> 9 pt 改善右余量（原 14 px）；
        ##   下距 4 -> 6 pt 给底部两排图例留呼吸
        plot.margin = margin(4, 9, 6, 4))

export_fig(p, file.path(OUT, "fig2f.pdf"), w = 6, h = 4)
cat("fig2f done\n")
