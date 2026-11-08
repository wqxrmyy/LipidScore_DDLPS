# ============================================================
# Fig S2B · 关键证据：方向矛盾对照（TCGA 高分预后差 vs GSE221492 DD 低分）
# 左：TCGA 评分 vs OS 风险（高分=差预后，HR=1.98）
# 右：GSE221492 DD vs WD 评分（DD=低分，Δ=-0.469）
# （原 Fig S3D，2026-09-28 重构统一为 Fig. S2 B；cairo_pdf 直存）
# ============================================================
suppressWarnings(suppressMessages(library(ggplot2)))
suppressWarnings(suppressMessages(library(jsonlite)))
suppressWarnings(suppressMessages(library(gridExtra)))
source(file.path(Sys.getenv("LPS_FIG_ROOT", "."), "utils", "colors.R"))
source(file.path(Sys.getenv("LPS_FIG_ROOT", "."), "utils", "themes.R"))

OUT <- file.path(Sys.getenv("LPS_FIG_ROOT", "."), "out")

# ---- 左图：TCGA 评分箱线（按 OS 事件分组：死亡 vs 存活）----
load(file.path(Sys.getenv("LPS_FIG_ROOT", "."), "p5_data.RData"))
df_tcga <- data.frame(
  Score = score_A$Score,
  status = clin[score_A$Sample, "OS_status"]
)
df_tcga <- df_tcga[!is.na(df_tcga$status), ]
df_tcga$Event <- factor(ifelse(df_tcga$status == 1, "Died", "Censored"),
                        levels = c("Died", "Censored"))

p_left <- ggplot(df_tcga, aes(x = Event, y = Score, fill = Event)) +
  geom_boxplot(alpha = 0.6, outlier.shape = NA, color = "grey20") +
  geom_jitter(width = 0.15, size = 1.5, alpha = 0.6) +
  scale_fill_manual(values = c(Died = "#E74C3C", Censored = "#7F8C8D"), name = "") +
  labs(x = "", y = "Lipid metabolism score (A)",
       title = "TCGA-SARC: higher score = worse OS",
       subtitle = "HR = 1.98, P = 0.0013; C = 0.690") +
  theme_publication()

# ---- 右图：GSE221492 DD vs WD 评分箱线 ----
load(file.path(Sys.getenv("LPS_FIG_ROOT", "."), "p5_gse221492.RData"))
d <- gse_meta[gse_meta$Group %in% c("DD", "WD"), ]
d$Group <- factor(d$Group, levels = c("WD", "DD"))

p_right <- ggplot(d, aes(x = Group, y = Score, fill = Group)) +
  geom_boxplot(alpha = 0.6, outlier.shape = NA, color = "grey20") +
  geom_jitter(width = 0.15, size = 1.5, alpha = 0.6) +
  scale_fill_manual(values = c(WD = "#2980B9", DD = "#E67E22"), name = "") +
  labs(x = "", y = "Lipid metabolism score (A)",
       title = "GSE221492: DD component = lower score",
       subtitle = "Delta (DD - WD) = -0.469, P = 0.025 (unpaired)") +
  theme_publication()

# 合并（尺寸由 themes.R 的 PANEL_BOX 统一接管）
g <- arrangeGrob(p_left, p_right, nrow = 2)
box <- panel_box("figS2b")
cairo_pdf(file.path(OUT, "figS2b.pdf"), width = box[1], height = box[2])
grid::grid.draw(g); dev.off()
png(file.path(OUT, "figS2b.png"), width = box[1], height = box[2], units = "in",
    res = 300, bg = "white", type = "cairo")
grid::grid.draw(g); dev.off()
cat(sprintf("figS2b done | export %.2f x %.2f in (=%d x %d px)\n",
            box[1], box[2], round(box[1] * 300), round(box[2] * 300)))
