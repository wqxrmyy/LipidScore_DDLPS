# ============================================================
# Fig 1B · 队列构成（b1: 3 队列 × 病理分组堆叠柱；b2: 253 样本来源环形图）
# 数据：03_clinical_outcomes.csv
# 注：旧 b2（TCGA 亚型饼图）废弃 —— TCGA 50 例全部为 DDLPS，单分类无信息量
# ============================================================
suppressWarnings(suppressMessages(library(ggplot2)))
source(file.path(Sys.getenv("LPS_FIG_ROOT", "."), "utils", "colors.R"))
source(file.path(Sys.getenv("LPS_FIG_ROOT", "."), "utils", "themes.R"))

load(file.path(Sys.getenv("LPS_FIG_ROOT", "."), "p5_data.RData"))
OUT <- file.path(Sys.getenv("LPS_FIG_ROOT", "."), "out")

# ---- 病理大类归类（保持真实标签语义，不合并到失真）----
# DDLPS = 去分化脂肪肉瘤；WDLPS = 高分化脂肪肉瘤；
# Other sarcoma = 非脂肪肉瘤亚型（黏液/多形性/MFH/平滑肌肉瘤等）；Normal = 正常脂肪/对照
grp <- function(cohort, subtype) {
  ifelse(subtype %in% c("DDLPS", "Dedifferentiated"), "DDLPS",
  ifelse(subtype %in% c("Well_differentiated"), "WDLPS",
  ifelse(subtype %in% c("Normal adipose", "disease status: Normal control"), "Normal/control",
         "Other sarcoma")))
}
clin$PathGroup <- grp(clin$Cohort, clin$Subtype)
clin$PathGroup <- factor(clin$PathGroup,
                         levels = c("DDLPS", "WDLPS", "Other sarcoma", "Normal/control"))
clin$Cohort <- factor(clin$Cohort, levels = c("TCGA", "GSE21122", "GSE159659"))

pg_tab <- as.data.frame(table(clin$Cohort, clin$PathGroup), stringsAsFactors = FALSE)
names(pg_tab) <- c("Cohort", "PathGroup", "n")

pal_pg <- c("DDLPS" = "#E74C3C", "WDLPS" = "#E67E22",
            "Other sarcoma" = "#95A5A6", "Normal/control" = "#27AE60")

# ---- b1：堆叠柱（各队列病理构成）----
p1 <- ggplot(pg_tab, aes(x = Cohort, y = n, fill = PathGroup)) +
  geom_bar(stat = "identity", position = "stack", width = 0.62,
           color = "white", linewidth = 0.3) +
  geom_text(data = subset(pg_tab, n >= 8),
            aes(label = n), position = position_stack(vjust = 0.5),
            size = 2.8, color = "white", fontface = "bold") +
  scale_fill_manual(values = pal_pg) +
  labs(x = NULL, y = "Number of samples",
       title = "Cohort composition by pathology",
       subtitle = "253 samples across three cohorts") +
  theme_publication(base_size = 10) +
  theme(legend.position = "bottom", legend.title = element_blank(),
        legend.text = element_text(size = 8),
        plot.subtitle = element_text(size = 8.5, color = "grey30"))

# ---- b2：环形图（253 样本的队列来源），楔内标百分比 + 右侧图例标队列与 n ----
# v3 修正（2026-09-18）：原把「队列名 + n + 百分比」三行写在楔内，面板收窄到
# 半幅（3.25 in）后，相邻两块（TCGA 19.8% / GSE159659 17.8%）的标签在 12 点
# 方向互相压字。改为楔内只留百分比，队列名与 n 交给图例。
coh_df <- aggregate(n ~ Cohort, pg_tab, sum)
coh_df$Cohort <- factor(coh_df$Cohort, levels = c("TCGA", "GSE21122", "GSE159659"))
coh_df <- coh_df[order(coh_df$Cohort), ]
coh_df$pct <- coh_df$n / sum(coh_df$n) * 100
# 扇区中点（y 单位与堆叠柱一致）；文本放在环带中央 x = 1
coh_df$mid <- cumsum(coh_df$n) - coh_df$n / 2
coh_df$lab <- sprintf("%.0f%%", coh_df$pct)
pal_coh <- c("TCGA" = "#E74C3C", "GSE21122" = "#2980B9", "GSE159659" = "#27AE60")
lab_coh <- setNames(sprintf("%s (n = %d)", coh_df$Cohort, coh_df$n), as.character(coh_df$Cohort))

p2 <- ggplot(coh_df, aes(x = "", y = n, fill = Cohort)) +
  geom_bar(stat = "identity", width = 1, color = "white", linewidth = 0.6) +
  coord_polar("y", start = 0) +
  geom_text(aes(x = 1.0, y = mid, label = lab),
            size = 3.0, color = "white", fontface = "bold") +
  scale_fill_manual(values = pal_coh, labels = lab_coh, name = NULL) +
  labs(title = "Sample source (n = 253)") +
  theme_void() +
  theme(text = element_text(colour = "black"),  # theme_void 会 blank 掉所有 text（实测连显式 plot.title 也被吞），先整体恢复
        legend.position = "right",
        legend.title = element_blank(),
        legend.text = element_text(size = 6.5, colour = "black"),
        legend.key.size = unit(0.26, "cm"),
        legend.margin = margin(0, 0, 0, 0),
        plot.title = element_text(face = "bold", size = 10, hjust = 0.5,
                                  colour = "black", margin = margin(b = 6)))

export_fig(p1, file.path(OUT, "fig1b1.pdf"), w = 5.0, h = 3.9)
export_fig(p2, file.path(OUT, "fig1b2.pdf"), w = 3.4, h = 3.4)
cat("fig1b done\n")
