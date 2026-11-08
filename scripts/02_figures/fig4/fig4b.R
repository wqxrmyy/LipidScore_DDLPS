# ============================================================
# Fig 4B · 患者内配对连线图（8 对 WD→DD，6/8 下降，Wilcoxon P=0.078）
# 数据：fig4b_pair_data.csv（官方 score_all53 = log2CPM + z-all-53，P4 口径）
# ============================================================
suppressWarnings(suppressMessages(library(ggplot2)))
suppressWarnings(suppressMessages(library(reshape2)))
source(file.path(Sys.getenv("LPS_FIG_ROOT", "."), "utils", "colors.R"))
source(file.path(Sys.getenv("LPS_FIG_ROOT", "."), "utils", "themes.R"))

OUT <- file.path(Sys.getenv("LPS_FIG_ROOT", "."), "out")

d <- read.csv(file.path(Sys.getenv("LPS_FIG_ROOT", "."), "fig4b_pair_data.csv"),
              stringsAsFactors = FALSE)
# 官方口径数值（硬编码，与 prep_fig4b_pair.py 校验一致）
n_down <- 6L
n_pair <- 8L
mean_delta <- -0.587
wp <- 0.078

d$patient <- factor(d$patient, levels = d$patient)
wl <- melt(d[, c("patient", "WD", "DD")], id.vars = "patient",
           variable.name = "Group", value.name = "Score")
wl$Group <- factor(wl$Group, levels = c("WD", "DD"))

p <- ggplot(wl, aes(x = Group, y = Score, group = patient)) +
  geom_line(aes(color = d$delta[match(patient, d$patient)] < 0), size = 0.8, alpha = 0.85) +
  geom_point(size = 2, aes(color = d$delta[match(patient, d$patient)] < 0)) +
  scale_color_manual(values = c(`TRUE` = "#E67E22", `FALSE` = "#2980B9"),
                     labels = c(`TRUE` = "DD < WD", `FALSE` = "DD >= WD"), name = "") +
  labs(x = "", y = "Lipid metabolism score (A)",
       title = "Within-patient paired WD vs DD",
       # v3：原句 58 字符，半幅（3.25 in / 7 pt，上限约 54 字符）会被裁断
       subtitle = sprintf("%d/%d decreased; delta = %.3f; Wilcoxon P = %.3f",
                          n_down, n_pair, mean_delta, wp)) +
  theme_publication()

export_fig(p, file.path(OUT, "fig4b.pdf"), w = 5, h = 4)
cat("fig4b done\n")
