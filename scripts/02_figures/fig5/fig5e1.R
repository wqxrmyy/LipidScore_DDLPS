# ============================================================
# Fig 5E1 · 组成层：脂肪系占比 WD -> DD（7 对配对）
#   与 E2/E3 同 y 范围（E_Y_FRAC = 0–1，本行共享刻度常量写在 fig5_common.R）
#   与 F1 同谱系同色（LIN_COL$Adipose）—— 镜像行规范 §8.4
# ============================================================
suppressWarnings(suppressMessages({library(ggplot2)}))
source(file.path(Sys.getenv("LPS_FIG_ROOT", "."), "utils", "colors.R"))
source(file.path(Sys.getenv("LPS_FIG_ROOT", "."), "utils", "themes.R"))
source(file.path(Sys.getenv("LPS_FIG_ROOT", "."), "utils", "fig5_common.R"))

d <- fig5_read("comp_sample_bucket.csv")
d <- d[d$bucket == "Adipose" & d$patient %in% FIG5_PAIRS &
       d$tissue %in% c("DDLPS-WD", "DDLPS-DD"), ]
stopifnot(nrow(d) == 14, length(unique(d$patient)) == 7)
dl <- fig5_pair_delta(d, "frac")
p <- fig5_pair_panel(d, "frac", "Adipose", "Adipocytic fraction", E_Y_FRAC,
                     dlab = dl, mult = 100, unit = " pp")
export_fig(p, file.path(FIG5_OUT, "fig5e1.pdf"), w = 1.62, h = 1.61)
dir_up <- sum(sapply(FIG5_PAIRS, function(k)
  with(d, frac[patient == k & tissue == "DDLPS-DD"] -
         frac[patient == k & tissue == "DDLPS-WD"])) > 0)
cat(sprintf("fig5e1 done | mean change = %+.2f pp | 逐患者升 %d/7\n", dl * 100, dir_up))
