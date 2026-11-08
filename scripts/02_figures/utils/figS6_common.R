# ============================================================
# FigS6 公共层（机制层：评分 <-> mTORC2 / TAM 轴）
#   依据：LAYOUT_SPEC_v3.md §8.6（v5）+ F6_caption_FINAL.md
#   ① 尺寸：一律 649 x 1120 px = 2.1656 x 3.7334 in（PANEL_BOX figS6a-f）
#      ★ 换算式：1949/300/3 = 2.16557 in（**不是** 649/300 = 2.1633 -> 那会给 648 px）
#   ② 图内文本一律 ASCII（禁 Unicode：rho / P / - / x 而非 ρ / Δ / − / ×）
#   ③ 数据只读 figS6_data/*（冻结件，哈希见 figS6_outplot_data_hash.json）
#   ④ 面板字母不在图内 —— 字母由人工拼版统一添加（同 Fig5 做法）
# ============================================================
suppressWarnings(suppressMessages({
  library(ggplot2); library(patchwork); library(scales); library(survival); library(Matrix)
}))

FIG6_ROOT <- file.path(Sys.getenv("LPS_FIG_ROOT", "."))
FIG6_DATA <- file.path(FIG6_ROOT, "figS6_data")
FIG6_OUT  <- file.path(FIG6_ROOT, "out")
W6        <- file.path(Sys.getenv("LPS_WORK_ROOT", "."), "_work")

# ---- 配色（唯一真相源；与 Fig5 谱系色不冲突）----------------------------
AXIS_COL <- c(RICTOR = "#A93226",   # 砖红：轴基因主读数（同 Fig5 Tumour 谱系）
              CD36   = "#1B6B7A",   # 深青：脂质摄取受体（同 Fig5 KEY_COL）
              MTOR   = "#7F8C8D",   # 灰：两复合物共享（被否）
              RPTOR  = "#B7950B")   # 芥黄：mTORC1 专有阴性对照
TAM_COL  <- c(TAM1 = "#F0B27A",     # 浅橙：TAM1（Myeloid 同系浅色）
              TAM3 = "#C0392B")     # 深红：TAM3 / M2 样
POOL_COL <- c(P1_tumor_all   = "#A93226",
              P2_tumor_prog  = "#CA6F1E",
              P3_adipo_diff  = "#D4A017",
              P4_normal_adi  = "#1B6B7A")
REG_COL  <- c("Up" = "#A93226", "Down" = "#1B6B7A", "NS" = "#B3B3B3")
SCORE_GRAD <- c("#E5E5E5", "#FED976", "#FC8D59", "#99000D")  # 同 Fig5B
SIG_COL    <- "#A93226"

# ---- 读冻结数据（CSV 一律 utf-8-sig）-----------------------------------
read6 <- function(f) read.csv(file.path(FIG6_DATA, f), stringsAsFactors = FALSE,
                              fileEncoding = "UTF-8-BOM", check.names = FALSE)

# ---- P 值格式化（ASCII，防 Unicode 溢出）-------------------------------
fmt_p <- function(p, digits = 3) {
  ifelse(is.na(p), "NA",
    ifelse(p < 1e-4, sprintf("P = %.1e", p), sprintf("P = %.3f", p)))
}

# ---- 通路名缩短（只缩短字面，不改含义；须在图注给对照）----------------
short_pw <- function(x, width = 20) {
  x <- sub("^(HALLMARK|KEGG|REACTOME|GOBP|WP)_", "", x)
  x <- gsub("_", " ", x)
  x <- tolower(x)
  x <- paste0(toupper(substr(x, 1, 1)), substr(x, 2, nchar(x)))
  ifelse(nchar(x) > width, paste0(substr(x, 1, width - 1), "."), x)
}

# ---- 长文本折行兜底（★ 2.166 in 面板：7.5 pt 标题上限约 28 字符）------
#   实测教训：标题/副标题比 panel 宽时，ggplot 的 gtable 会被撑宽 -> 内容溢出被裁
wrap_txt <- function(s, n = 30) {
  if (is.null(s) || !length(s)) return(s)
  if (nchar(s) <= n) return(s)
  paste(strwrap(s, width = n), collapse = "\n")
}

# ---- 上下两栏拼版（面板内复合读数；先例 figS3d/figS4d）-----------------
#   ★ 用 wrap_elements() 各自包裹：避免 patchwork 的「同列 panel 强制等宽」
#     把上栏 panel 向下栏的宽 y 轴标签看齐（实测导致上栏标题右移出界）
#   heights 比例之和决定两栏相对高度；tag 由人工拼版添加，故此处不加字母
figS6_stack <- function(up, dn, heights = c(1.25, 1)) {
  wrap_elements(up) / wrap_elements(dn) + plot_layout(heights = heights)
}

# ---- 保存（统一走 export_fig；PANEL_BOX 会覆盖为 649 x 1120）-----------
save6 <- function(p, stem) {
  export_fig(p, file.path(FIG6_OUT, paste0(stem, ".pdf")))
  cat(sprintf("  [figS6] %s -> out/%s.pdf + .png\n", stem, stem))
}

# ---- KM 曲线 -> data.frame（自绘，避免 survminer 的 print/尺寸不可控）--
km_df <- function(fit) {
  out <- do.call(rbind, lapply(names(fit$strata), function(s) {
    idx <- match(s, names(fit$strata))
    end <- cumsum(fit$strata)[idx]
    beg <- if (idx == 1) 1 else cumsum(fit$strata)[idx - 1] + 1
    data.frame(time = c(0, fit$time[beg:end]),
               surv = c(1, fit$surv[beg:end]),
               strata = sub("^[^=]*=", "", s), stringsAsFactors = FALSE)
  }))
  out
}

# 风险表：在给定时点统计 at-risk（n.risk）
km_risk <- function(d, time, group, breaks) {
  do.call(rbind, lapply(breaks, function(t0) {
    do.call(rbind, lapply(names(table(d[[group]])), function(g) {
      dd <- d[d[[group]] == g, ]
      data.frame(group = g, time = t0, n = sum(dd[[time]] >= t0), stringsAsFactors = FALSE)
    }))
  }))
}

# KM 对数秩 P（分层变量名，字符）
km_p <- function(d, group) {
  f <- as.formula(paste("Surv(OS_time, OS_status) ~", group))
  sd <- survdiff(f, data = d)
  1 - pchisq(sd$chisq, length(sd$n) - 1)
}

# ---- 自绘 KM + 风险表（单面板，无 C-index，无 Unicode）-----------------
#   col_map: named vector (strata 名 -> 颜色)
km_panel <- function(d, group, col_map, title, xmax = NULL, show_risk = TRUE,
                     base_size = 7.5, risk_h = 0.20, legend_title = NULL,
                     xbreaks = NULL) {
  f <- as.formula(paste("Surv(OS_time, OS_status) ~", group))
  fit <- survfit(f, data = d)
  km <- km_df(fit)
  lv <- names(col_map)
  km$strata <- factor(km$strata, levels = lv)
  if (is.null(xmax)) xmax <- max(d$OS_time)
  p <- ggplot(km, aes(time, surv, colour = strata)) +
    geom_step(linewidth = 0.55) +
    scale_colour_manual(values = col_map, name = legend_title, breaks = lv) +
    scale_x_continuous(limits = c(0, xmax), breaks = xbreaks,
                       expand = expansion(mult = c(0.01, 0.02))) +
    scale_y_continuous(limits = c(0, 1.02), expand = expansion(mult = c(0.01, 0.02))) +
    labs(x = "Time (days)", y = "Overall survival", title = title) +
    theme_publication(base_size = base_size) +
    theme(legend.position = c(0.72, 0.86),
          legend.background = element_rect(fill = "white", colour = NA),
          legend.key.height = unit(7, "pt"), legend.key.width = unit(9, "pt"),
          legend.text = element_text(size = base_size - 1.4),
          legend.title = element_blank(),
          plot.margin = margin(3, 3, 2, 3))
  p
}

# ---- 森林点图（点估计 + 95%CI；铁律：C-index/相关系数比较禁柱状图）----
#   ★★ 2026-09-24 修订（S36）：vline 语义歧义修复。
#   原参数 vline 只画一条线，6D 传的是 MDE(0.517)，而读者会把它当作
#   「显著性阈值」-> 把 P=0.070 的 underpowered 项误读为「无效」。
#   现改为双阈值：sig_line（P<0.05 临界 |rho|，实线细）与
#   mde_line（80% power MDE，点线）**分别可按需给**，并在图上加区分性标注；
#   只传 vline 时行为与旧版一致（向后兼容，供其他面板复用）。
forest_panel <- function(d, xlab, xlims, vline = NULL, base_size = 7.5,
                         lab_size = 6.6, title = NULL,
                         sig_line = NULL, mde_line = NULL,
                         sig_lab = NULL, mde_lab = NULL) {
  d$label <- factor(d$label, levels = rev(d$label))
  p <- ggplot(d, aes(x = est, y = label)) +
    geom_vline(xintercept = 0, linewidth = 0.35, colour = "grey55") +
    geom_errorbarh(aes(xmin = lo, xmax = hi), height = 0.16, linewidth = 0.45,
                   colour = "grey35") +
    geom_point(aes(fill = cls), shape = 21, size = 2.0, stroke = 0.5, colour = "grey20") +
    scale_fill_manual(values = c(ref = SIG_COL, neg = "#B3B3B3", NS = "#7F8C8D",
                                sens = "#1B6B7A"), guide = "none") +
    scale_x_continuous(limits = xlims) +
    labs(x = xlab, y = NULL, title = title) +
    theme_publication(base_size = base_size) +
    theme(axis.text.y = element_text(size = lab_size),
          ## ★ 2026-09-24：y 轴标签普遍较长（"A-PPARG vs M2-like" 等 57 pt），
          ##   默认 position="panel" 会把标题推向面板右缘并被裁（实测右余量 0 px）。
          ##   改 "plot" 锚定整图左缘。本面板 labs(y=NULL) 无纵排 y 标题，安全。
          plot.title.position = "plot",
          plot.margin = margin(3, 5, 3, 3))
  # ---- P < 0.05 临界线（细实线，深灰；代表「显著性阈值」）----
  if (!is.null(sig_line)) {
    p <- p + geom_vline(xintercept = c(-sig_line, sig_line), linetype = 1,
                        linewidth = 0.42, colour = "grey25", alpha = 0.90)
  }
  # ---- 80% power MDE（点线，砖红；代表「可检出的最小效应」）----
  if (!is.null(mde_line)) {
    p <- p + geom_vline(xintercept = c(-mde_line, mde_line), linetype = 3,
                        linewidth = 0.45, colour = SIG_COL, alpha = 0.85)
  }
  # ---- 仅传 vline（旧接口）：保留原行为 ----
  if (!is.null(vline) && is.null(sig_line) && is.null(mde_line)) {
    p <- p + geom_vline(xintercept = c(-vline, vline), linetype = 3,
                        linewidth = 0.38, colour = SIG_COL, alpha = 0.75)
  }
  # ---- 阈值线标注（放绘图区下沿，避免与点/CI 相撞）----
  #   ★★ 宽度实测（s37_w.R）：贴片文字须 <= 面板宽 155.9 pt，扣 y 标签占位后红线约 90 pt。
  #     故用极简式：`solid: P=0.05 (0.381)  dotted: MDE (0.517)` = 102.0 pt ✗
  #     改 `P=0.05 (0.381); MDE (0.517)` = **69.8 pt** ✓（线型由 linetype 本身区分，
  #     文字不必再重复 "solid/dotted"，避免与标题争宽）
  lab <- character(0)
  if (!is.null(sig_line) && !is.null(mde_line))
    lab <- sprintf("%s (%.3f); %s (%.3f)", sig_lab, sig_line, mde_lab, mde_line)
  else if (!is.null(sig_line))
    lab <- sprintf("%s (%.3f)", sig_lab, sig_line)
  else if (!is.null(mde_line))
    lab <- sprintf("%s (%.3f)", mde_lab, mde_line)
  if (length(lab)) {
    p <- p + labs(caption = lab) +
      theme(plot.caption = element_text(size = 5.4, colour = "grey30",
                                        hjust = 0, margin = margin(t = 2)))
  }
  p
}
