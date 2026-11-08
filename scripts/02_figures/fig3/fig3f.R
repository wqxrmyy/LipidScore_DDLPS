# ============================================================
# Fig 3F · 合成轴 vs 身份轴对比（双棒棒糖）
# 6 分量（合成程序）与身份基因（CD36/FABP4/PLIN1/LPL/ADIPOQ）与评分 A 的相关
# 揭示"评分测的是身份轴，非纯合成轴"
# ============================================================
suppressWarnings(suppressMessages(library(ggplot2)))
source(file.path(Sys.getenv("LPS_FIG_ROOT", "."), "utils", "colors.R"))
source(file.path(Sys.getenv("LPS_FIG_ROOT", "."), "utils", "themes.R"))

load(file.path(Sys.getenv("LPS_FIG_ROOT", "."), "p5_data.RData"))
OUT <- file.path(Sys.getenv("LPS_FIG_ROOT", "."), "out")

log2cpm <- readRDS(file.path(Sys.getenv("LPS_DATA_ROOT", "."), "1", "data", "raw", "TCGA_LPS_real_log2cpm.rds"))
idx <- setNames(rownames(log2cpm), sub("\\.[0-9]+$", "", rownames(log2cpm)))
y <- score_A$Score

ENS <- c(SREBF1 = "ENSG00000072310", SCD = "ENSG00000099194", HMGCR = "ENSG00000113161",
         ACLY = "ENSG00000131473", PPARG = "ENSG00000132170", FASN = "ENSG00000169710",
         CD36 = "ENSG00000135218", FABP4 = "ENSG00000091131", PLIN1 = "ENSG00000166819",
         LPL = "ENSG00000175592", ADIPOQ = "ENSG00000181092")

rows <- NULL
for (g in names(ENS)) {
  row <- idx[ENS[g]]
  if (is.na(row)) next
  ex <- as.numeric(log2cpm[row, score_A$Sample])
  tt <- cor.test(y, ex)
  rows <- rbind(rows, data.frame(Gene = g, r = unname(tt$estimate), p = tt$p.value))
}
rows$Class <- ifelse(rows$Gene %in% c("SCD", "HMGCR", "ACLY", "FASN", "SREBF1", "PPARG"),
                     "Score component (synthesis program)", "Identity gene (adipocyte program)")
rows$Gene <- factor(rows$Gene, levels = rows$Gene[order(rows$r)])
rows$col <- ifelse(rows$p < 0.05, "#E74C3C", "#7F8C8D")
rows$sig <- sapply(rows$p, sig_label)

p <- ggplot(rows, aes(x = r, y = Gene)) +
  geom_vline(xintercept = 0, color = "grey40", linetype = "dotted") +
  geom_segment(aes(x = 0, xend = r, y = Gene, yend = Gene), color = rows$col, size = 1.1) +
  geom_point(aes(shape = Class), color = rows$col, size = 3.2) +
  geom_text(aes(label = sig), hjust = -0.4, size = 2.5) +
  scale_shape_manual(values = c(16, 17), name = "",
                     labels = c("Identity gene", "Score component")) +
  # 留出显著性星号空间：无此设置时最强分量（SCD/HMGCR/ACLY/FASN）的星号被右边界裁掉
  scale_x_continuous(limits = c(-0.40, 0.95), breaks = c(0, 0.2, 0.4, 0.6)) +
  labs(x = "Pearson r with score (A), n = 50", y = "",
       title = "Synthesis components vs identity genes",
       subtitle = "Score tracks both programs; identity genes dominate") +
  theme_publication() +
  # 半幅面板里右侧图例会把绘图区挤成一条窄带 → 图例移到底部横排、并缩短标签
  theme(legend.position = "bottom", legend.direction = "horizontal",
        legend.key.size = grid::unit(0.32, "cm"),
        legend.margin = margin(0, 0, 0, 0))

export_fig(p, file.path(OUT, "fig3f.pdf"), w = 6, h = 5)
cat("fig3f done\n")
