# ============================================================
# Fig 2E · 纯度校正偏相关（原始 vs 控制纯度 vs 三重校正 棒棒糖）
# A 基准：0.5354 / 0.629 / 0.6146
# ============================================================
suppressWarnings(suppressMessages(library(ggplot2)))
source(file.path(Sys.getenv("LPS_FIG_ROOT", "."), "utils", "colors.R"))
source(file.path(Sys.getenv("LPS_FIG_ROOT", "."), "utils", "themes.R"))

load(file.path(Sys.getenv("LPS_FIG_ROOT", "."), "p5_data.RData"))
OUT <- file.path(Sys.getenv("LPS_FIG_ROOT", "."), "out")

numf <- function(x) suppressWarnings(as.numeric(as.character(x)))
cn  <- numf(clin_real[score_A$Sample, "paper_CN.cluster"])
pur <- numf(clin_real[score_A$Sample, "paper_purity"])
gra <- numf(clin_real[score_A$Sample, "paper_FNCLCC.grade"])
mit <- numf(clin_real[score_A$Sample, "paper_mitotic.rate"])
y <- score_A$Score

partial_r <- function(x, y, Z) {
  keep <- complete.cases(x, y, Z)
  X <- x[keep]; Y <- y[keep]; Zm <- as.matrix(Z[keep, , drop = FALSE])
  A <- cbind(1, Zm)
  bx <- coef(lm(X ~ A)); by <- coef(lm(Y ~ A))
  rx <- residuals(lm(X ~ A)); ry <- residuals(lm(Y ~ A))
  tt <- cor.test(rx, ry)
  c(r = unname(tt$estimate), p = tt$p.value, n = sum(keep))
}

r_raw <- cor.test(cn, y)
r_pur <- partial_r(cn, y, data.frame(purity = pur))
r_tri <- partial_r(cn, y, data.frame(purity = pur, grade = gra, mitotic = mit))

df <- data.frame(
  Model = c("Raw", "Adj. purity", "Adj. purity +\ngrade + mitotic"),
  r = c(unname(r_raw$estimate), r_pur["r"], r_tri["r"]),
  p = c(r_raw$p.value, r_pur["p"], r_tri["p"]),
  n = c(unname(r_raw$parameter) + 2, r_pur["n"], r_tri["n"])
)
df$Model <- factor(df$Model, levels = df$Model)
# v3：改为「右对齐锚在点上」，并去掉 (n = ...)。
#   原写法 hjust 默认 0.5 + 文本含 n（32 字符 ≈ 2.2 in），而半幅绘图区仅约 2.3 in，
#   点在 r≈0.62 处已占 93%，右半段被设备边缘裁断。
df$lab <- sprintf("r = %.3f, P = %.2g", df$r, df$p)

p <- ggplot(df, aes(x = Model, y = r)) +
  geom_hline(yintercept = 0, color = "grey50") +
  geom_segment(aes(x = Model, xend = Model, y = 0, yend = r), color = "#2980B9", size = 1.2) +
  geom_point(size = 4.5, color = "#E74C3C") +
  ## 0929 二次修正：vjust=-1.1 + 0.6 扩充仍被面板上缘裁断顶行（裁掉约 1/3 字高）。
  ##   改为 vjust=-0.7（贴近数据点）+ 上端扩充 0.75 类目单位
  geom_text(aes(label = lab), hjust = 1, vjust = -0.7, size = 2.9, color = "grey20") +
  scale_x_discrete(expand = expansion(add = c(0.25, 0.75))) +
  coord_flip() +
  labs(x = "", y = "Partial Pearson r (CN cluster ~ score)", title = "Purity-adjusted association") +
  ## ★★ S39：右距 4 -> 9 pt 改善右余量（原 14 px）；
  ##   上距 4 -> 12 pt 容纳 geom_text(vjust=-1.1) 在 coord_flip 后的顶部标注
  theme_publication() +
  theme(plot.margin = margin(12, 9, 6, 4))

export_fig(p, file.path(OUT, "fig2e.pdf"), w = 5.5, h = 3)
cat("fig2e done\n")
