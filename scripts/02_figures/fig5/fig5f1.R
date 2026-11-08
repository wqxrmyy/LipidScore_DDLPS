# ============================================================
# Fig 5F1 · 状态层：脂肪系谱系内平均评分 WD -> DD（7 对配对）
#   与 F2/F3 同 y 范围（F_Y_SCORE）；与 E1 同谱系同色 —— 镜像行规范 §8.4
#   度量 = 该样本「脂肪系」全部细胞的评分 A 均值（= Step 7 分解里的 m，
#          与 pseudobulk_sample_x_bucket.csv 逐值一致，见 s9a 报告 §[5]）
# ============================================================
suppressWarnings(suppressMessages({library(ggplot2)}))
source(file.path(Sys.getenv("LPS_FIG_ROOT", "."), "utils", "colors.R"))
source(file.path(Sys.getenv("LPS_FIG_ROOT", "."), "utils", "themes.R"))
source(file.path(Sys.getenv("LPS_FIG_ROOT", "."), "utils", "fig5_common.R"))

d <- fig5_read("state_sample_bucket.csv")
d <- d[d$bucket == "Adipose" & d$patient %in% FIG5_PAIRS &
       d$tissue %in% c("DDLPS-WD", "DDLPS-DD"), ]
stopifnot(nrow(d) == 14)
dl <- fig5_pair_delta(d, "mean_score_ams")
p <- fig5_pair_panel(d, "mean_score_ams", "Adipose", "Adipocytic score", F_Y_SCORE,
                     dlab = dl, mult = 1, unit = "", dgits = 4)
export_fig(p, file.path(FIG5_OUT, "fig5f1.pdf"), w = 1.62, h = 1.61)
cat(sprintf("fig5f1 done | mean change = %+.5f\n", dl))
