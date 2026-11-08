# ============================================================
# Fig 5F3 · 状态层：免疫系谱系内平均评分 WD -> DD（7 对配对）
#   与 F1/F2 同 y 范围（F_Y_SCORE）；与 E3 同谱系同色
# ============================================================
suppressWarnings(suppressMessages({library(ggplot2)}))
source(file.path(Sys.getenv("LPS_FIG_ROOT", "."), "utils", "colors.R"))
source(file.path(Sys.getenv("LPS_FIG_ROOT", "."), "utils", "themes.R"))
source(file.path(Sys.getenv("LPS_FIG_ROOT", "."), "utils", "fig5_common.R"))

d <- fig5_read("state_sample_bucket.csv")
d <- d[d$bucket == "Immune" & d$patient %in% FIG5_PAIRS &
       d$tissue %in% c("DDLPS-WD", "DDLPS-DD"), ]
stopifnot(nrow(d) == 14)
dl <- fig5_pair_delta(d, "mean_score_ams")
p <- fig5_pair_panel(d, "mean_score_ams", "Immune", "Immune score", F_Y_SCORE,
                     dlab = dl, mult = 1, unit = "", dgits = 4)
export_fig(p, file.path(FIG5_OUT, "fig5f3.pdf"), w = 1.62, h = 1.61)
cat(sprintf("fig5f3 done | mean change = %+.5f\n", dl))
