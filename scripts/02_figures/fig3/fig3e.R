# ============================================================
# Fig 3E · KM 曲线（A 评分高 vs 低，TCGA OS n=50）
# 高分组中位 OS 2.79 年 vs 低分 7.05 年
# ------------------------------------------------------------
# v3（2026-09-18）：改用 utils/themes.R 的 panel_box("fig3e") 取尺寸，
#   与其余面板同口径（按最终印刷尺寸绘制，拼版 1:1 不缩放）。
#   基准字号 11 → 8；risk.table 字号同步下调。
# ============================================================
suppressWarnings(suppressMessages({
  library(ggplot2); library(survival); library(survminer)
}))
source(file.path(Sys.getenv("LPS_FIG_ROOT", "."), "utils", "colors.R"))
source(file.path(Sys.getenv("LPS_FIG_ROOT", "."), "utils", "themes.R"))

load(file.path(Sys.getenv("LPS_FIG_ROOT", "."), "p5_data.RData"))
OUT <- file.path(Sys.getenv("LPS_FIG_ROOT", "."), "out")

df <- data.frame(
  time = clin[score_A$Sample, "OS_time"] / 365,
  status = clin[score_A$Sample, "OS_status"],
  Score = score_A$Score
)
df <- df[!is.na(df$time) & !is.na(df$status), ]
df$Risk <- ifelse(df$Score > median(df$Score), "High", "Low")
df$Risk <- factor(df$Risk, levels = c("Low", "High"))

fit <- survfit(Surv(time, status) ~ Risk, data = df)
sd <- survdiff(Surv(time, status) ~ Risk, data = df)
p_lr <- pchisq(sd$chisq, length(sd$n) - 1, lower.tail = FALSE)

box <- panel_box("fig3e")   # c(宽, 高)，单位英寸

p <- ggsurvplot(fit, data = df,
                pval = FALSE,
                risk.table = TRUE,
                risk.table.height = 0.26,
                fontsize = 2.5,          # risk.table 字号（mm 单位 ≈ 7 pt）
                palette = c("#2980B9", "#E74C3C"),
                legend.labs = c("Low score", "High score"),
                xlab = "Time (years)", ylab = "Overall survival probability",
                break.time.by = 2,
                surv.median.line = "hv",
                ggtheme = theme_classic(base_size = 8))
p$plot <- p$plot +
  annotate("text", x = 6.5, y = 0.92, size = 2.5,
           label = sprintf("Log-rank P = %.3g", p_lr)) +
  ggtitle("Overall survival by lipid metabolism score (A)") +
  theme(plot.title = element_text(size = 8, face = "bold"))

# v3：窄面板里风险表的 "Strata" 轴标题会压住组名 → 去掉轴标题
# 注意：只动 axis.title.y；axis.text.y 由 survminer 预置为向量化 element_text，
#       再叠加会触发 ggplot2「Only elements of the same class can be merged」错误。
p$table <- p$table + theme(axis.title.y = element_blank())

# newpage = FALSE：print.ggsurvplot 默认会 grid.newpage()，
# 在 pdf 设备中导致第 1 页空白、内容落到第 2 页（PNG 单页所以看似正常）
cairo_pdf(file.path(OUT, "fig3e.pdf"), width = box[1], height = box[2])
print(p, newpage = FALSE)
dev.off()
png(file.path(OUT, "fig3e.png"), width = box[1], height = box[2], units = "in", res = 300)
print(p, newpage = FALSE)
dev.off()
cat(sprintf("fig3e done | export %.2f x %.2f in (=%d x %d px)\n",
            box[1], box[2], round(box[1] * 300), round(box[2] * 300)))
