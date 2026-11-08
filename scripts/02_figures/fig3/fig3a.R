# ============================================================
# Fig 3A · CD36 vs 评分散点（TCGA，A 基准，r=0.3574 P=0.011）
# ============================================================
suppressWarnings(suppressMessages(library(ggplot2)))
source(file.path(Sys.getenv("LPS_FIG_ROOT", "."), "utils", "colors.R"))
source(file.path(Sys.getenv("LPS_FIG_ROOT", "."), "utils", "themes.R"))

load(file.path(Sys.getenv("LPS_FIG_ROOT", "."), "p5_data.RData"))
OUT <- file.path(Sys.getenv("LPS_FIG_ROOT", "."), "out")

cd36 <- as.numeric(expr6["PPARG", ]) * 0  # 占位避免作用域问题
log2cpm <- readRDS(file.path(Sys.getenv("LPS_DATA_ROOT", "."), "1", "data", "raw", "TCGA_LPS_real_log2cpm.rds"))
idx <- setNames(rownames(log2cpm), sub("\\.[0-9]+$", "", rownames(log2cpm)))
row36 <- idx["ENSG00000135218"]
cd36 <- as.numeric(log2cpm[row36, score_A$Sample])

df <- data.frame(Score = score_A$Score, CD36 = cd36)
tt <- cor.test(df$Score, df$CD36)

p <- ggplot(df, aes(x = Score, y = CD36)) +
  geom_point(size = 1.8, alpha = 0.75, color = "#2980B9") +
  geom_smooth(method = "lm", se = TRUE, color = "#E74C3C", size = 0.8) +
  labs(x = "Lipid metabolism score (A)", y = "CD36 log2(CPM)",
       title = "CD36 association (TCGA, n = 50)",
       subtitle = sprintf("r = %.3f, P = %.3g", tt$estimate, tt$p.value)) +
  theme_publication()

export_fig(p, file.path(OUT, "fig3a.pdf"), w = 4, h = 3.5)
cat("fig3a done\n")
