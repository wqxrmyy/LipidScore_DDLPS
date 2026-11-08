# ============================================================
# Fig 5E3 · 组成层：免疫系占比 WD -> DD（7 对配对）
#   与 E1/E2 同 y 范围（E_Y_FRAC）；与 F3 同谱系同色
# ============================================================
suppressWarnings(suppressMessages({library(ggplot2)}))
source(file.path(Sys.getenv("LPS_FIG_ROOT", "."), "utils", "colors.R"))
source(file.path(Sys.getenv("LPS_FIG_ROOT", "."), "utils", "themes.R"))
source(file.path(Sys.getenv("LPS_FIG_ROOT", "."), "utils", "fig5_common.R"))

d <- fig5_read("comp_sample_bucket.csv")
d <- d[d$bucket == "Immune" & d$patient %in% FIG5_PAIRS &
       d$tissue %in% c("DDLPS-WD", "DDLPS-DD"), ]
stopifnot(nrow(d) == 14, length(unique(d$patient)) == 7)
dl <- fig5_pair_delta(d, "frac")
p <- fig5_pair_panel(d, "frac", "Immune", "Immune fraction", E_Y_FRAC,
                     dlab = dl, mult = 100, unit = " pp")
export_fig(p, file.path(FIG5_OUT, "fig5e3.pdf"), w = 1.62, h = 1.61)
cat(sprintf("fig5e3 done | mean change = %+.2f pp\n", dl * 100))
