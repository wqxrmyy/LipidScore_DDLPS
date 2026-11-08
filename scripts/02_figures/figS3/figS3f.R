# ============================================================
# Fig S3F · 效力 vs 显著性森林（原 Fig 6D 下栏，2026-09-28 重构）
#   森林点图（rho 点估计 + Fisher-z 95%CI）＋ 双阈值线：
#     ① 实线细（grey25）= P < 0.05 临界 |rho| = 0.3809（n=27, df=25）
#     ② 点线红（SIG_COL）= 80% power MDE = 0.517
#   ★ 铁律：相关系数比较用点估计 + 95%CI，禁柱状图
#   数据：figS6_data/figS6d_forest.csv
# ============================================================
suppressWarnings(suppressMessages({library(ggplot2); library(patchwork); library(scales)}))
source(file.path(file.path(Sys.getenv("LPS_FIG_ROOT", "."), "utils"), "themes.R"))
source(file.path(file.path(Sys.getenv("LPS_FIG_ROOT", "."), "utils"), "figS6_common.R"))

F6 <- read6("figS6d_forest.csv")
ORD <- c("Macrophage PPARG", "M2-like score", "Macrophage CD36", "Macrophage RICTOR",
         "A-PPARG vs CD36", "A-PPARG vs M2-like")
d <- data.frame(label = F6$label, est = F6$rho, lo = F6$lo, hi = F6$hi, cls = F6$class)
d <- d[match(ORD, d$label), ]
stopifnot(!any(is.na(d$label)))
MDE <- F6$MDE[1]
SIG_CRIT <- 0.3809
cat(sprintf("   P=0.05 critical |rho| at n=27 = %.4f | 80%% power MDE = %.3f\n", SIG_CRIT, MDE))
stopifnot(abs(d$est[d$label == "Macrophage PPARG"] - SIG_CRIT) < 0.01)

pD2 <- forest_panel(d, "Spearman rho (95% CI)", c(-0.85, 0.85),
                    base_size = 7.5, lab_size = 6.2,
                    sig_line = SIG_CRIT, mde_line = MDE,
                    sig_lab = "P=0.05", mde_lab = "MDE",
                    ## 2026-09-29：标题折行（文字不变，仅排版）——该标题≈满宽且
                    ##   plot.title.position="plot"，单行时顶到面板左缘、必然压
                    ##   面板内左上角的字母 F；折行后每行变短，居中即与字母脱开
                    title = "Effect size vs significance\nand limits")

## ★★ 2026-09-29：折行后的标题**居中**（hjust = 0.5）。
##   依据（实测，out/figS3f.png 649x540，300 dpi）：
##     - 标题第 1 行墨迹 x = 14-397、y = 12-39（左对齐顶到面板左缘）；
##     - 人工拼版字母 F 占 x = 8-32、y = 8-39（PANEL 左上角 +8/+8，10.5 pt 粗体）
##     -> 二者正面相撞，即 PI 报的「字母 F 与子图标题重叠」。
##   居中后第 1 行 x ≈ 133-516，与字母右缘（32 px）余量 ≈101 px，彻底脱开。
##   ★ 只动排版（折行 + 对齐），标题文字与字号/字重/颜色逐项保持原规则
##     （theme_publication: size = base_size = 7.5, face = "bold", colour = "black"）。
pD2 <- pD2 + theme(plot.title = element_text(hjust = 0.5, size = 7.5,
                                             face = "bold", colour = "black"))

save6(pD2, "figS3f")
