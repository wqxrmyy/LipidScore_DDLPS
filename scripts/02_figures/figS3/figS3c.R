# ============================================================
# Fig S3C · 通路层：GSVA 循环的排除（pseudobulk，去评分六基因）
#   上栏：热图 —— top 通路 x 4 细胞池（rho_clean，已剔除评分六基因）
#   下栏：PPAR 与 PI3K-AKT-mTOR 两条通路的双口径对照
#   （原 Fig 6B，2026-09-28 重构统一为 Fig. S3 C）
#   数据：figS6_data/figS6b_wide_rho_clean.csv + figS6b_wide_q_clean.csv
#         + pathway_vs_score.csv（含 rho_tum / rho_tum_noPP 双口径）
# ============================================================
suppressWarnings(suppressMessages({library(ggplot2); library(patchwork); library(scales)}))
source(file.path(file.path(Sys.getenv("LPS_FIG_ROOT", "."), "utils"), "themes.R"))
source(file.path(file.path(Sys.getenv("LPS_FIG_ROOT", "."), "utils"), "figS6_common.R"))

POOLS <- c("P1_tumor_all", "P2_tumor_prog", "P3_adipo_diff", "P4_normal_adi")
PLAB  <- c(P1_tumor_all = "all", P2_tumor_prog = "prog",
           P3_adipo_diff = "adip", P4_normal_adi = "fat")
rho <- read6("figS6b_wide_rho_clean.csv")
qq  <- read6("figS6b_wide_q_clean.csv")
stopifnot(identical(rho$pathway, qq$pathway))
stopifnot(all(POOLS %in% names(rho)))

## --- top 通路：按 4 池平均 |rho| 降序，取 20 条 ---
sc <- apply(abs(rho[, POOLS]), 1, mean)
rho$mean_abs <- sc
KEY <- c("KEGG_PPAR_SIGNALING_PATHWAY", "HALLMARK_PI3K_AKT_MTOR_SIGNALING")
stopifnot(all(KEY %in% rho$pathway))
top <- unique(c(KEY, rho$pathway[order(-sc)][1:20]))
top <- rho$pathway[order(-sc)][1:20]
sel <- rho[rho$pathway %in% top, ]
cat(sprintf("S3C up: top20 by mean|rho| (range: %.3f - %.3f)\n",
            min(sel$mean_abs), max(sel$mean_abs)))

long <- do.call(rbind, lapply(POOLS, function(p)
  data.frame(pathway = sel$pathway, pool = p, rho = sel[[p]],
             q = qq[[p]][match(sel$pathway, qq$pathway)], stringsAsFactors = FALSE)))
long$pool <- factor(long$pool, levels = POOLS)
ylev <- sel$pathway[order(sel$mean_abs)]
PW_SHORT <- c(
  HALLMARK_INTERFERON_GAMMA_RESPONSE           = "Interferon-g",
  HALLMARK_COMPLEMENT                          = "Complement",
  HALLMARK_KRAS_SIGNALING_UP                   = "Kras up",
  HALLMARK_IL6_JAK_STAT3_SIGNALING             = "IL6-JAK-STAT3",
  HALLMARK_ALLOGRAFT_REJECTION                 = "Allograft",
  HALLMARK_COAGULATION                         = "Coagulation",
  HALLMARK_INFLAMMATORY_RESPONSE               = "Inflammatory",
  HALLMARK_IL2_STAT5_SIGNALING                 = "IL2-STAT5",
  KEGG_GLYCEROLIPID_METABOLISM                 = "Glycerolipid",
  HALLMARK_INTERFERON_ALPHA_RESPONSE           = "Interferon-a",
  HALLMARK_TNFA_SIGNALING_VIA_NFKB             = "TNFA-NFKB",
  HALLMARK_PEROXISOME                          = "Peroxisome",
  HALLMARK_REACTIVE_OXYGEN_SPECIES_PATHWAY     = "ROS pathway",
  KEGG_GLYCEROPHOSPHOLIPID_METABOLISM          = "Glycerophos.",
  HALLMARK_DNA_REPAIR                          = "DNA repair",
  KEGG_BIOSYNTHESIS_OF_UNSATURATED_FATTY_ACIDS = "Unsat. FA",
  HALLMARK_MYC_TARGETS_V2                      = "MYC targets",
  HALLMARK_ANDROGEN_RESPONSE                   = "Androgen",
  HALLMARK_WNT_BETA_CATENIN_SIGNALING          = "WNT-b-catenin",
  HALLMARK_TGF_BETA_SIGNALING                  = "TGF-beta")
stopifnot(all(sel$pathway %in% names(PW_SHORT)))
sh <- unname(PW_SHORT[sel$pathway])
stopifnot(!any(is.na(sh)), max(nchar(sh)) <= 13, length(unique(sh)) == nrow(sel))
long$pw <- factor(unname(PW_SHORT[long$pathway]), levels = unname(PW_SHORT[ylev]))
long$star <- ifelse(long$q < 0.05, "*", "")
cat("  strata with q<0.05: ", sum(long$star == "*"), " / ", nrow(long), "\n", sep = "")

pB1 <- ggplot(long, aes(pool, pw)) +
  geom_tile(aes(fill = rho), colour = "white", linewidth = 0.5) +
  geom_text(aes(label = star), size = 2.6, colour = "black", vjust = 0.78, fontface = "bold") +
  scale_fill_gradientn(colours = c("#2166AC", "#F7F7F7", "#B2182B"),
                       limits = c(-0.7, 0.7), oob = scales::squish,
                       name = "rho",
                       guide = guide_colourbar(direction = "horizontal",
                                               barwidth = unit(56, "pt"),
                                               barheight = unit(5, "pt"),
                                               ticks.colour = "black",
                                               frame.colour = "black",
                                               title.position = "top")) +
  scale_x_discrete(labels = PLAB[POOLS], expand = c(0, 0)) +
  scale_y_discrete(expand = c(0, 0)) +
  labs(x = NULL, y = NULL, title = "Pathways vs score") +
  theme_publication(base_size = 7.5) +
  theme(axis.line = element_blank(), axis.ticks = element_blank(),
        axis.text.x = element_text(size = 6.4, lineheight = 0.92),
        axis.text.y = element_text(size = 5.3),
        plot.title = element_text(size = 7.5),
        legend.position = "bottom",
        legend.title = element_text(size = 6.2),
        legend.text = element_text(size = 6.0),
        legend.margin = margin(0, 0, 0, 0),
        plot.margin = margin(3, 3, 2, 3))

## ---------------------------------------------------------------- 下栏
pv <- read6("pathway_vs_score.csv")
pv <- pv[pv$pool %in% POOLS & pv$pathway %in% KEY, ]
stopifnot(nrow(pv) == 8)
pv$short <- ifelse(grepl("PPAR", pv$pathway), "PPAR", "PI3K")
dl <- rbind(
  data.frame(short = pv$short, pool = pv$pool, spec = "A + PPARG", rho = pv$rho_tum),
  data.frame(short = pv$short, pool = pv$pool, spec = "A - PPARG", rho = pv$rho_tum_noPP))
dl$pool <- factor(dl$pool, levels = POOLS)
PSHORT <- c(P1_tumor_all = "all", P2_tumor_prog = "progs",
            P3_adipo_diff = "adipo", P4_normal_adi = "fat")
dl$row  <- factor(paste(dl$short, PSHORT[dl$pool], sep = " / "),
                  levels = rev(c(paste(rep(c("PPAR", "PI3K"), each = 4),
                                       rep(PSHORT[POOLS], 2), sep = " / "))))
dl$spec <- factor(dl$spec, levels = c("A + PPARG", "A - PPARG"))

wide <- reshape(dl[, c("row", "spec", "rho")], idvar = "row", timevar = "spec", direction = "wide")
names(wide) <- c("row", "noPP", "withPP")

pB2 <- ggplot(wide) +
  geom_vline(xintercept = 0, linewidth = 0.35, colour = "grey55") +
  geom_segment(aes(x = withPP, xend = noPP, y = row, yend = row),
               colour = "grey55", linewidth = 0.5) +
  geom_point(aes(x = withPP, y = row, fill = "A + PPARG"), shape = 21,
             size = 2.1, stroke = 0.5, colour = "grey20") +
  geom_point(aes(x = noPP, y = row, fill = "A - PPARG"), shape = 21,
             size = 2.1, stroke = 0.5, colour = "grey20") +
  scale_fill_manual(values = c("A + PPARG" = "#B3B3B3", "A - PPARG" = SIG_COL),
                    name = NULL) +
  scale_x_continuous(limits = c(-0.85, 0.85), breaks = c(-0.5, 0, 0.5)) +
  scale_y_discrete(expand = expansion(add = 0.45)) +
  labs(x = "rho vs score", y = NULL,
       title = wrap_txt("PPARG: with vs without", 22)) +
  theme_publication(base_size = 7.5) +
  theme(axis.text.y = element_text(size = 6.0),
        legend.position = "bottom",
        legend.key.height = unit(7, "pt"), legend.key.width = unit(9, "pt"),
        legend.text = element_text(size = 6.0),
        legend.margin = margin(0, 0, 0, 0),
        plot.title = element_text(size = 7.5),
        plot.margin = margin(3, 3, 2, 3))

p <- figS6_stack(pB1, pB2, heights = c(1.58, 1))
save6(p, "figS3c")

## 关键读数回显
for (pw0 in c("HALLMARK_PPAR_SIGNALING", "HALLMARK_PI3K_AKT_MTOR_SIGNALING")) {
  s <- pv[pv$pathway == pw0, ]
  cat(sprintf("   %-34s P1: A=%+.3f -> A-PPARG=%+.3f\n", pw0,
              s$rho_tum[s$pool == "P1_tumor_all"], s$rho_tum_noPP[s$pool == "P1_tumor_all"]))
}
