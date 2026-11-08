# ============================================================
# Fig S2C · 效力分析曲线（PPARG vs CD36 的 MDIE 曲线，80% 阈值线）
# 数据：s3_power_PPARG.json / s3_power_CD36.json -> grid(nie_true, power)
# （原 Fig S4B，2026-09-28 重构统一为 Fig. S2 C）
# ============================================================
suppressWarnings(suppressMessages(library(ggplot2)))
suppressWarnings(suppressMessages(library(jsonlite)))
source(file.path(Sys.getenv("LPS_FIG_ROOT", "."), "utils", "colors.R"))
source(file.path(Sys.getenv("LPS_FIG_ROOT", "."), "utils", "themes.R"))

OUT <- file.path(Sys.getenv("LPS_FIG_ROOT", "."), "out")

pp <- fromJSON(file.path(Sys.getenv("LPS_WORK_ROOT", "."), "_work", "s3_power_PPARG.json"), simplifyVector = FALSE)
cd <- fromJSON(file.path(Sys.getenv("LPS_WORK_ROOT", "."), "_work", "s3_power_CD36.json"), simplifyVector = FALSE)

gp <- data.frame(nie_true = sapply(pp$grid, function(x) x$nie_true),
                 power = sapply(pp$grid, function(x) x$power), Mediator = "PPARG")
gc <- data.frame(nie_true = sapply(cd$grid, function(x) x$nie_true),
                 power = sapply(cd$grid, function(x) x$power), Mediator = "CD36")
df <- rbind(gp, gc)

# MDIE(80%) 由线性插值（P4 报告：PPARG 0.806, CD36 0.896）
mdie <- c(PPARG = 0.806, CD36 = 0.896)

p <- ggplot(df, aes(x = nie_true, y = power, color = Mediator)) +
  geom_hline(yintercept = 0.8, linetype = "dashed", color = "grey40") +
  geom_line(size = 1) +
  geom_point(size = 2.5) +
  geom_vline(xintercept = mdie["PPARG"], linetype = "dotted", color = "#2980B9") +
  geom_vline(xintercept = mdie["CD36"], linetype = "dotted", color = "#E67E22") +
  scale_color_manual(values = c(PPARG = "#2980B9", CD36 = "#E67E22"), name = "") +
  annotate("text", x = 0.02, y = 0.86, hjust = 0, size = 3, color = "grey30",
           label = "80% power threshold") +
  # 0929 版式重排：面板压矮后两条竖排 MDIE 标注互相重叠 → 合并为一条横排
  annotate("text", x = 1.16, y = 0.02, hjust = 1, vjust = 0, size = 3, color = "grey20",
           label = "MDIE(80%):\nPPARG 0.806, CD36 0.896") +
  labs(x = "True indirect effect (NIE)", y = "Power",
       title = "Power to detect mediation (n = 50, 23 events)",
      subtitle = "Observed NIE: PPARG 0.10, CD36 -0.63") +
  theme_publication()

export_fig(p, file.path(OUT, "figS2c.pdf"), w = 5, h = 4)
cat("figS2c done\n")
