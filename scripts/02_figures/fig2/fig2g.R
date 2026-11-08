# ============================================================
# Fig 2G · 与 12q 扩增子正交性（CN 变量互相关热图 + 扩增子表达 NS）
# 上半：CN 变量间互相关矩阵；体现 ploidy~WGD=0.925、CancerDNA~purity=0.954 共线性
# 下半：CN 簇 vs 12q 扩增子基因表达（全 NS → 正交）
# ============================================================
suppressWarnings(suppressMessages(library(ggplot2)))
suppressWarnings(suppressMessages(library(reshape2)))
source(file.path(Sys.getenv("LPS_FIG_ROOT", "."), "utils", "colors.R"))
source(file.path(Sys.getenv("LPS_FIG_ROOT", "."), "utils", "themes.R"))

load(file.path(Sys.getenv("LPS_FIG_ROOT", "."), "p5_data.RData"))
OUT <- file.path(Sys.getenv("LPS_FIG_ROOT", "."), "out")

numf <- function(x) suppressWarnings(as.numeric(as.character(x)))

# --- 上：互相关矩阵（A 基准）---
vlist <- c("CN cluster" = "paper_CN.cluster", "Ploidy" = "paper_ploidy",
           "Genome dbl." = "paper_Genome.doublings", "Cancer DNA" = "paper_Cancer.DNA.fraction",
           "Purity" = "paper_purity", "CIN" = "paper_CIN", "12q CN" = "paper_avg.CN.12q13.15")
mat <- matrix(NA, length(vlist), length(vlist), dimnames = list(names(vlist), names(vlist)))
for (i in seq_along(vlist)) for (j in seq_along(vlist)) {
  a <- numf(clin_real[score_A$Sample, vlist[i]])
  b <- numf(clin_real[score_A$Sample, vlist[j]])
  ok <- !is.na(a) & !is.na(b)
  if (sum(ok) > 10 && length(unique(a[ok])) > 2) mat[i, j] <- cor(a[ok], b[ok])
}
diag(mat) <- 1
mm <- melt(mat); names(mm) <- c("V1", "V2", "r")

p1 <- ggplot(mm, aes(x = V1, y = V2, fill = r)) +
  geom_tile(color = "white") +
  geom_text(aes(label = ifelse(is.na(r), "", sprintf("%.2f", r))), size = 2.3) +
  scale_fill_gradient2(low = "#2980B9", mid = "white", high = "#E74C3C", midpoint = 0,
                       na.value = "grey95", limits = c(-1, 1)) +
  labs(x = "", y = "", fill = "r", title = "CN variable correlation matrix") +
  theme_publication() +
  ## ★★ S39：右距 4 -> 9 pt 改善右余量；下距 4 -> 6 pt；
  ##   上距 4 -> 10 pt —— 实测原 T 贴顶，标题 "CN variable correlation matrix"
  ##   首字母 C 顶部被切、colorbar 图例 1.0 刻度压顶（放大 3x 原生像素确认）
  theme(axis.text.x = element_text(angle = 45, hjust = 1),
        plot.margin = margin(10, 9, 6, 4))

# v3（2026-09-18）：改走 export_fig，尺寸由 themes.R 的 PANEL_BOX 统一接管。
# 旧写法直接 ggsave(width=5.5/4.5) 绕过版面盒 → 拼版时 1350–1650 px 宽超出半幅。
export_fig(p1, file.path(OUT, "fig2g1.pdf"), w = 5.5, h = 4.2)

# --- 下：CN 簇 vs 扩增子基因表达（MDM2/CDK4/HMGA2/FRS2，全 NS）---
log2cpm <- readRDS(file.path(Sys.getenv("LPS_DATA_ROOT", "."), "1", "data", "raw", "TCGA_LPS_real_log2cpm.rds"))
idx <- setNames(rownames(log2cpm), sub("\\.[0-9]+$", "", rownames(log2cpm)))
amp <- c(MDM2 = "ENSG00000135679", CDK4 = "ENSG00000135446", HMGA2 = "ENSG00000149948", FRS2 = "ENSG00000166225")
cn <- numf(clin_real[score_A$Sample, "paper_CN.cluster"])
rows2 <- NULL
for (g in names(amp)) {
  row <- idx[amp[g]]
  ex <- as.numeric(log2cpm[row, score_A$Sample])
  tt <- cor.test(cn, ex)
  rows2 <- rbind(rows2, data.frame(Gene = g, r = unname(tt$estimate), p = tt$p.value))
}
## ★★ S39 修订（2026-09-24，PI 指令「Fig2 也改善一下」）
##   实测根因（像素法 + 放大 2x 原生像素）：
##     ① 6 个面板右余量均仅 14 px（theme plot.margin 右距 4 pt 太紧）
##     ② g2 的 P 标注：原 geom_text(vjust = -1.2) 固定向上偏，
##        r<0 的三项（CDK4/HMGA2/MDM2）标注正好落在 y=0 参考线上被穿字；
##        r>0 的 FRS2 标注则离点过远（点 -0.15、标注 +0.25），视觉误导
##     ③ g2 竖排 y 标题 "Pearson r (CN cluster vs expression)"（36 字符）
##        在 535 px 面板高下超出 -> 顶/底两端均被切（Bmarg=0）
##     ④ g2 x 轴标签贴底（Bmarg=0）
##   修法：① 标注按 r 正负分向（正→上、负→下），并给 y 轴留安全边
##         ② y 标题折两行 ③ 下距加大 ④ g 整体右距加大
p2 <- ggplot(rows2, aes(x = Gene, y = r)) +
  geom_hline(yintercept = 0, color = "grey50") +
  geom_segment(aes(x = Gene, xend = Gene, y = 0, yend = r), color = "#7F8C8D", size = 1.2) +
  geom_point(size = 4, color = "#2980B9") +
  ## ★ 标注按 r 正负分向：正 r 标在点上方，负 r 标在点下方，
  ##   两者都背离 y=0 参考线 -> 彻底避开"标注压在零线上"（原 vjust=-1.2 固定向上）
  geom_text(aes(label = sprintf("P=%.2f", p),
                vjust = ifelse(r >= 0, -1.3, 2.3)),
            size = 3) +
  ## ★ 上下各留 0.09 安全边：原 c(-0.35, 0.35) 使负 r 项的标注无下落空间
  coord_cartesian(ylim = c(-0.42, 0.42)) +
  ## ★ y 标题折两行（36 字符竖排超 535 px 高 -> 原顶/底被切）
  labs(x = NULL,
       y = "Pearson r\n(CN cluster vs expression)",
       title = "Orthogonality to 12q amplicon expression",
       subtitle = "All NS: score tracks CN architecture, not 12q dose") +
  theme_publication() +
  ## ★ 下距 4 -> 8 pt 给 x 轴标签呼吸（原 Bmarg=0 贴底）；
  ##   右距 4 -> 10 pt 改善右余量（原 14 px）；上距 4 -> 6 pt 容纳 P 标注
  theme(plot.margin = margin(6, 10, 8, 4))

export_fig(p2, file.path(OUT, "fig2g2.pdf"), w = 4.5, h = 3.2)
cat("fig2g done\n")
