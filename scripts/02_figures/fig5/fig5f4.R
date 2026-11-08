# ============================================================
# Fig 5F4 · 状态层关键证据：Kitagawa 对称分解（判决的图形化）
#   内容：逐患者（细条）与合计（粗条）的组成项 C、状态项 T 对总差 dS 的贡献
#   数据：Step 7 冻结产物（s7/kitagawa_patient_L2.csv）＋ D23 主判据（T1_main_pooled.csv）
#   ⚠️ P4 是唯一反向病例（dS > 0），标签加 * 单列 —— 不得静默剔除（D23 明令）
#   ⚠️ 主判据口径 = 合计（D23 选项 C）；逐患者中位数作为异质性披露，二者并列呈现
# ============================================================
suppressWarnings(suppressMessages({library(ggplot2)}))
source(file.path(Sys.getenv("LPS_FIG_ROOT", "."), "utils", "colors.R"))
source(file.path(Sys.getenv("LPS_FIG_ROOT", "."), "utils", "themes.R"))
source(file.path(Sys.getenv("LPS_FIG_ROOT", "."), "utils", "fig5_common.R"))

kp <- read.csv(file.path(FIG5_DATA, "kitagawa_patient_F4.csv"),
               stringsAsFactors = FALSE, fileEncoding = "UTF-8-BOM")
t1 <- read.csv(file.path(W5, "s24_evidence_pack", "T1_main_pooled.csv"),
               stringsAsFactors = FALSE, fileEncoding = "UTF-8-BOM")
stopifnot(nrow(kp) == 7, nrow(t1) == 1)
stopifnot(abs(t1$ratio - 0.8663) < 5e-4)                  # D23 冻结值防漂移

kp <- kp[order(-kp$ratio), ]                              # 组成主导 -> 状态主导
med_ratio <- median(kp$ratio)
stopifnot(abs(med_ratio - 0.5370) < 5e-4)                 # D23 冻结值防漂移

LEV <- rev(kp$patient)                                    # 自下而上：P4(反向) … P1(最强组成)
stopifnot(LEV[1] == "P4")                                 # 反向病例必须在最下（D23 明令单列）
LEV <- c(LEV, "pooled")
y_of <- function(p) factor(p, levels = LEV)

# ★ 标签直接由 LEV 推导（不用手工 rev 对齐 —— 上一版正因手工 rev 而标签与数据反序）
#   变量名不得叫 labs：会遮蔽 ggplot2::labs() 使下游报错
ylab_vec <- c(paste0(LEV[1:7],
                     ifelse(kp$dS[match(LEV[1:7], kp$patient)] > 0, " *", "")), "pooled")

dp <- rbind(data.frame(y = y_of(kp$patient), val = kp[["C"]], series = "Composition"),
            data.frame(y = y_of(kp$patient), val = kp[["T"]], series = "State"))
dl <- data.frame(y = y_of("pooled"), val = c(t1$C_total, t1$T_total),
                 series = c("Composition", "State"))
sub <- sprintf("pooled %.2f | median %.2f", t1$ratio, med_ratio)

# 可证性检查：逐层级打印「标签 ↔ 实际数值」，防标签与数据错位
chk <- data.frame(level = seq_along(LEV), label = ylab_vec,
                  C = c(kp$C[match(LEV[1:7], kp$patient)], t1$C_total),
                  T = c(kp$T[match(LEV[1:7], kp$patient)], t1$T_total))
cat("--- 层级 <-> 标签 <-> 数值 ---\n"); print(chk, row.names = FALSE, digits = 5)
stopifnot(chk$label[1] == "P4 *", chk$label[8] == "pooled")

p <- ggplot() +
  geom_hline(yintercept = 7.5, colour = "grey80", linewidth = 0.3) +
  geom_vline(xintercept = 0, colour = "grey45", linewidth = 0.35) +
  geom_col(data = dp, aes(x = val, y = y, fill = series),
           position = position_dodge(width = 0.66), width = 0.58) +
  geom_col(data = dl, aes(x = val, y = y, fill = series),
           position = position_dodge(width = 1.00), width = 0.88) +
  scale_fill_manual(values = c(Composition = unname(CT_COL["C"]),
                               State = unname(CT_COL["T"])), name = NULL) +
  scale_y_discrete(labels = ylab_vec) +
  scale_x_continuous(limits = c(-0.175, 0.095), breaks = c(-0.15, -0.10, -0.05, 0)) +
  labs(x = "Contribution to dS", y = NULL,
       title = "Kitagawa decomposition", subtitle = sub) +
  theme_publication() +
  theme(axis.line.y = element_blank(), axis.ticks.y = element_blank(),
        axis.title.y = element_blank(),
        axis.text.y = element_text(size = 6.2),
        axis.text.x = element_text(size = 6.2),
        axis.title.x = element_text(size = 6.2),
        plot.subtitle = element_text(size = 6.2, colour = "grey25"),
        legend.position = "bottom", legend.direction = "horizontal",
        legend.key.size = unit(6.5, "pt"), legend.text = element_text(size = 6.2),
        legend.margin = margin(0, 0, 0, 0), legend.box.margin = margin(-3, 0, 0, 0),
        plot.margin = margin(2, 2, 2, 2))

export_fig(p, file.path(FIG5_OUT, "fig5f4.pdf"), w = 1.62, h = 1.61)
cat(sprintf("fig5f4 done | pooled C=%.5f T=%.5f ratio=%.4f | median ratio=%.4f | 自下而上: %s\n",
            t1$C_total, t1$T_total, t1$ratio, med_ratio, paste(ylab_vec, collapse = " < ")))
