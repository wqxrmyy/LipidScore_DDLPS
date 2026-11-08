# ============================================================
# Fig S3B · RICTOR/CD36 六层级共定位证据矩阵（原 Fig 6A 下栏，2026-09-28 重构）
#   ★ fill 只编码「方向 x 显著性」，数值与 P 印在格内 —— 不做跨层归一化
#   数据：figS6_data/figS6a_layers.csv
# ============================================================
suppressWarnings(suppressMessages({library(ggplot2); library(patchwork)}))
source(file.path(file.path(Sys.getenv("LPS_FIG_ROOT", "."), "utils"), "themes.R"))
source(file.path(file.path(Sys.getenv("LPS_FIG_ROOT", "."), "utils"), "figS6_common.R"))

L <- read6("figS6a_layers.csv")
L$cls <- ifelse(L$sig == "yes", ifelse(L$est > 0, "sig+", "sig-"), "NS")
ord <- c("L6_bulk_adj", "L5_bulk_raw", "L4_sample", "L3_scc", "L2_fisher", "L1_cell")
# ★ y 标签保持 <= 16 字符（20 字符 ≈ 270 px，会把 panel 压到放不下数值）
NL <- c(L1_cell = "single cell", L2_fisher = "within detected",
        L3_scc = "sample x program", L4_sample = "sample (scRNA)",
        L5_bulk_raw = "TCGA-SARC raw", L6_bulk_adj = "TCGA-SARC adj")
L$row <- factor(NL[L$layer], levels = rev(unname(NL[ord])))
L$gene <- factor(L$gene, levels = c("RICTOR", "CD36"))
L$val <- ifelse(L$metric == "OR", sprintf("%.2f", L$est), sprintf("%+.3f", L$est))
L$lab1 <- paste0(ifelse(L$metric == "OR", "OR ", ""), L$val)
L$lab2 <- fmt_p(L$P)
L$cx <- as.numeric(L$gene)

pA2 <- ggplot(L, aes(cx, row)) +
  geom_tile(aes(fill = cls), colour = "white", linewidth = 0.6, width = 0.96, height = 0.94) +
  geom_text(aes(label = lab1), size = 2.35, colour = ifelse(L$cls == "NS", "grey15", "white"),
            fontface = "bold") +
  scale_fill_manual(values = c("sig+" = "#C0392B", "sig-" = "#1B6B7A", "NS" = "#EDEDED"),
                    guide = "none") +
  scale_x_continuous(breaks = c(1, 2), labels = c("RICTOR", "CD36"),
                     limits = c(0.5, 2.5), expand = c(0, 0)) +
  labs(x = NULL, y = NULL, title = "Co-localisation by level") +
  theme_publication(base_size = 7.5) +
  theme(axis.line = element_blank(), axis.ticks = element_blank(),
        axis.text.x = element_text(size = 7.0, face = "bold"),
        axis.text.y = element_text(size = 6.6),
        plot.title = element_text(size = 7.5),
        panel.background = element_blank(),
        plot.margin = margin(3, 3, 2, 3))

save6(pA2, "figS3b")

## 关键读数回显（供验收对账）
for (i in seq_len(nrow(L))) cat(sprintf("   %-22s %-7s %-9s %s\n",
  as.character(L$row[i]), as.character(L$gene[i]), L$lab1[i], L$lab2[i]))
