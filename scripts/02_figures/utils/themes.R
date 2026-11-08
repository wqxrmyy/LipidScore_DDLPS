# ============================================================
# P5 出图 · 全局主题（ggplot2 出版级）
# 依赖：ggplot2
# ------------------------------------------------------------
# v3 变更（2026-09-18，Cancer Science 6.5×8 in 合规整改）：
#   1) 目标期刊版面盒：165 × 203 mm = 6.5 × 8 in。
#      画布宽 1949 px @300dpi；单页高上限 2400 px。
#   2) ★ 面板一律「按最终印刷尺寸」绘制，拼版 1:1 贴图、不再缩放。
#      旧做法是画大（6.2 in / base 11pt）再缩到半幅（3.25 in）→ 字号 ×0.52
#      → 4.5–6.0 pt，跌破官方 6 pt 下限。1:1 后字号 = 作者设定值。
#   3) PANEL_BOX：按输出文件名（不含扩展名）覆盖 w/h。
#      未列入的面板沿用脚本自身的 w/h（兼容旧行为，但会给出 warning）。
#   4) theme_publication 默认 base_size 11 → 8（axis.text = 7 pt，留出下限余量）。
# ============================================================
suppressWarnings(suppressMessages(library(ggplot2)))

# ---- 版面常量 ----------------------------------------------------------
FIG_W_IN  <- 1949 / 300    # 6.4967 in = 165 mm
FIG_H_IN  <- 2400 / 300    # 8.0000 in = 203 mm（单页上限）
PAIR_W_IN <- FIG_W_IN / 2  # 3.2483 in（半幅面板，2 列行）
THIRD_W_IN <- FIG_W_IN / 3 # 2.1656 in（三列行的一列）
QUART_W_IN <- FIG_W_IN / 4 # 1.6242 in（四列行的一列）

# ---- 各面板目标 box（宽 in, 高 in）------------------------------------
# 依据：_work/p5_figures/LAYOUT_SPEC_v3.md（网格分配与 r 上限）
# ★ 列数选择有实测依据（LAYOUT_SPEC §8.2）：
#     2 列 974 px、3 列 650 px 均容得下换行后的轴标签；
#     4 列 487 px **仅限轴标签极短的图型**（配对连线 / UMAP / 热图 / 示意图）。
#   行数上限 = 4（单页 2400 px，每行扣 62 px 字母带，半幅最小实测行高 389 px）。
.half  <- function(r) c(PAIR_W_IN,  round(PAIR_W_IN  * r, 4))
.third <- function(r) c(round(THIRD_W_IN, 4), round(THIRD_W_IN * r, 4))
.quart <- function(r) c(round(QUART_W_IN, 4), round(QUART_W_IN * r, 4))
.full  <- function(r) c(FIG_W_IN,   round(FIG_W_IN   * r, 4))

PANEL_BOX <- list(
  # ---------------- Fig 1（A–F 六联图；原 1F/1G 已迁至 FigS5）------------
  #   版式：r1 = A(整行) | r2 = B(整行内 b1|b2) | r3 = C|D(半幅) | r4 = E|F(半幅)
  #   高预算：604 + 536 + 487 + 487 + 4×62(字母带) = 2362 px < 2400 ✓
  fig1a  = .full(0.31),    # 示意图，需整行宽（半幅会把 9 个文本框挤爆）
  fig1b1 = .half(0.55),    # 队列构成-堆叠柱（面板 B 左半）
  fig1b2 = .half(0.55),    # 队列构成-环形图（面板 B 右半）
  fig1c  = .half(0.50),    # 6 基因 z 分热图
  fig1d  = .half(0.50),    # 6 基因表达小提琴
  fig1e  = .half(0.50),    # 评分分布直方图
  fig1h  = .half(0.50),    # 三队列评分箱线（现为面板 F）

  # ---------------- Fig 2（A|B、C|D、E|F 半幅；G 整行）------------------
  # 0929 版式重排（对齐 FigS3 标准）：B/C 稳健性直方图上下叠、D/E 校正点线图上下叠，
  #   A/F 独柱加高至 830，小面板统一 415 —— 行内零空档；g 行不动。
  fig2a = .half(0.852), fig2b = .half(0.426), fig2c = .half(0.426),
  fig2d = .half(0.426), fig2e = .half(0.426), fig2f = .half(0.852),
  fig2g1 = .half(0.55), fig2g2 = .half(0.55),

  # ---------------- Fig 3（3 行 × 2 列）---------------------------------
  #   A–D r=0.70（守总高）；E（KM+风险表）/F（长轴标题）需更高 → r=0.85
  #   预算：2×682 + 828 + 3×62 = 2378 px < 2400 ✓
  # 0929 版式重排：C/D 为两点稀疏统计面板，压至 415 与 Fig2 小面板同规格；A/B、E/F 不动。
  fig3a = .half(0.70), fig3b = .half(0.70), fig3c = .half(0.426),
  fig3d = .half(0.426), fig3e = .half(0.85), fig3f = .half(0.85),

  # ---------------- Fig 4（A|B、C|D、E|F 半幅；G 整行）------------------
  #   v3 行高重排（2026-09-18，二轮）：按「内容最小高度」而非均分分配。
  #     实测 chrome（标题+副标题+轴标题+轴标签+图例+边距）≈163–208 px：
  #       C 有 14 个基因类别 → 绘图区需 ≥126 pt(=525 px) → 面板 ≥645 px（图例已移右）
  #       F 是 2×3 分面     → 每分面行需 ≥0.9 in → 面板 ≥585 px
  #       G 标签两行        → 4 行 × 17.4 pt → 面板 ≥453 px
  #       A/B/D 内容简单    → 390 px 足够
  #   v3b（2026-09-25）★ 行高再分配：E/F 584 → 650 px，代价由 C 支付。
  #     起因：E 是 15 行热图，584 px 下绘图区仅 352 px → 行距 23.5 px，
  #       而 6.4 pt 标签墨高 18–19 px → 行内净空只剩 3–6 px（虽然已分开，但余量很薄）。
  #     实测响应函数（fig4e，标尺 = 文字墨心的等距网格，**不是色标**）：
  #       绘图区 ≈ 面板高 − 232 px（chrome 恒定）→ 行距 = (h − 232)/15
  #       h=584 → 23.5 px（净空 4 px）；h=650 → 28.0 px（净空 9 px）
  #     预算：Fig4 原总高 2381 px（上限 2400，余 19 px），E/F 各 +66 px → 需再腾 47 px。
  #       来源 = C：701 → 653 px。依据：C 的 14 类行距实测 35.0 px（墨高 21 px、净空 14 px
  #       = 0.48 em），判据下限 31.2 px（= 21 + 0.35 em）；653 px 下绘图区 442 px
  #       → 行距 31.6 px、净空 10.6 px（0.36 em），**仍在下限之上**。
  #       不选 G：G 的余量没能测到（坐标翻转后面板的 axis.line 不在竖直方向，轴线法失效），
  #       不动没量过的余量。
  #     ★ 收益与代价都要记账（不得静默降级）：
  #       E/F：行内净空 4 → 9 px（翻倍以上）
  #       C  ：行内净空 14.0 → 10.6 px（从 0.48 em 降到 0.36 em，仍达标但余量收窄）
  #     ★ 650 px 仍差 0.5 px 未达 §7.5 的 0.35 em 舒适线（E 的墨高 19 px，需行距 ≥28.33 px）。
  #       要真正达标需 657 px，而 657 px 会把 C 压到 0.348 em（低于下限）→ **买不到**。
  #       若日后要再放宽 E，真正的杠杆不是行高，而是 x 轴那 166 px 的 45° 标签块
  #       （占 E 的 chrome 232 px 中的 166 px）→ 折行/改 0° 可回收约 124 px。
  #     ★ 另修一处静默丢组：fig4e 的 g_order 原写 "Lipoma"、数据是 "lipoma"
  #       → 图只画 6 组而副标题写 fat/lipoma，图文矛盾。已改为双向断言 stop()。
  fig4a = .half(0.40), fig4b = .half(0.40),
  # v4（2026-09-25，S65 续）：fig4d 0.40 -> 0.6704（389 -> 653 px），与同行 fig4c 同高。
  #   ★ 动机：r2 行高本来就由 fig4c（653 px）决定，fig4d 提到同高**不增加总高预算**
  #     （Fig4 仍 2399 px < 2400），等于白捡 264 px。
  #   实测病灶（974 x 389 px）：绘图区仅 176 px、4 个 y 类别 -> 类别间距 42 px，
  #     而 7 pt 文字 em 行高 29 px + 折行第二行 -> 标签互相粘连（放大裁图实证：
  #     "DD vs WD (DDLPS)" / "DD vs any WD" / "(WD+sWD+WDLPS)" 三行连成一片）；
  #     且最上数值标签 "0.249, P = 0.013" 被 geom_text(vjust=-0.8) 顶出 panel 上沿，
  #     被 ggplot 的 clip="on" 水平切掉上半（非画布裁切）。
  #   ★ 处置 = 只加高、不改内容、不改字号：653 px 下绘图区 ~440 px -> 类别间距
  #     440/4.2 = 105 px，折行标签占 ~58 px，净空 ~47 px；首类别中心距顶 63 px
  #     > 标签上探 ~13 px -> 不再出界。
  fig4c = .half(0.6704), fig4d = .half(0.6704),
  fig4e = .half(0.6676), fig4f = .half(0.6676),
  fig4g = .full(0.236),

  # ---------------- FigS1–S3（2026-09-28 重构：6 张补充图 -> 3 张纯主题）-----
  #   FigS1 = 评分构建与逐基因（5 格）；FigS2 = 配对+方向+效力（4 格）；
  #   FigS3 = 单细胞机制（9 格，R1 上下叠）
  figS1a = .half(0.70),    # 评分分布 GSE21122（原 FigS1B）
  figS1b = .half(0.70),    # 评分分布 GSE159659（原 FigS1C）
  figS1c = .half(0.72),    # 6 基因相关矩阵（原 FigS1D）
  figS1d = .half(0.72),    # 6 分量森林（原 figS2a_components）
  figS1e = .full(0.40),    # 评分 vs 分量散点矩阵（原 FigS5A）
  figS2a = .half(1.15),    # 分量基因 DD vs WD 小提琴 4x2（原 FigS3B；0929 版式重排：加高至与 B 等高填满 r1）
  figS2b = .half(1.15),    # 方向矛盾对照（原 FigS3D）
  figS2c = .half(0.65),    # 效力曲线（原 FigS4B；0929 版式重排：与 D 等高咬合填满 r2）
  figS2d = .half(0.65),    # 标志物森林（原 figS2a_markers；0929 版式重排：与 C 等高咬合填满 r2）
  # FigS3 三列行：高密度件（C/D/I：通路热图/配体-受体/KM）保持 1120 px，
  # 上下叠拆件（A/B、E/F、G/H）各 540 px
  figS3a = .third(0.8313), figS3b = .third(0.8313),
  figS3c = .third(1.724),  figS3d = .third(1.724),
  figS3e = .third(0.8313), figS3f = .third(0.8313),
  figS3g = .third(0.8313), figS3h = .third(0.8313),
  figS3i = .third(1.724),

  # ---------------- Fig 5（单细胞：组成 vs 状态判决）★设计稿 ---------------
  #   设计依据：LAYOUT_SPEC_v3.md §8.5。6 字母 / 12 图件 / 4 行（2-2-4-4 列）。
  #   ★ 两层计数：字母面板数 6（符合期刊「鼓励 ≤6」的软口径），图件数 12。
  #     r3/r4 四列共用字母 E / F，以数字下标区分 → 图件 > 字母（机制见 §8.1）。
  #   ★ r3/r4 之所以可用 4 列：这两行全是配对连线 / 简单条形，轴标签仅
  #     「WD」「DD」与 0–1 刻度，无长轴标题，落在 4 列适用区间内。
  #     **回退条件**：出图实测若裁切或压字 → 退回 3 列（图件 12→10），
  #     不得靠缩小字号解决（字号 6 pt 是期刊硬指标）。
  #   高预算：685 + 578 + 544 + 544 = 2351 px < 2400 ✓
  fig5a  = .half(0.64),    # UMAP：按细胞类型着色
  fig5b  = .half(0.64),    # UMAP：按评分 A 着色
  fig5c  = .half(0.53),    # 点图：6 分量基因 × 细胞类型（表达比例 + 平均表达）
  fig5d  = .half(0.53),    # 评分按细胞类型（箱线 / 密度）
  fig5e1 = .quart(0.99),   # 组成层 Δf：脂肪细胞系比例 WD→DD
  fig5e2 = .quart(0.99),   # 组成层 Δf：肿瘤细胞系
  fig5e3 = .quart(0.99),   # 组成层 Δf：免疫系
  fig5e4 = .quart(0.99),   # 组成层关键证据：脂肪细胞占比梯度（含 fat/lipoma 参照上限）
  fig5f1 = .quart(0.99),   # 状态层 Δm：脂肪系（与 fig5e1 同谱系 / 同色）
  fig5f2 = .quart(0.99),   # 状态层 Δm：肿瘤系（与 fig5e2 同谱系 / 同色）
  fig5f3 = .quart(0.99),   # 状态层 Δm：免疫系（与 fig5e3 同谱系 / 同色）
  fig5f4 = .quart(0.99)   # 状态层关键证据：Kitagawa 分解（组成项 vs 状态项）

  # ---------------- Fig 6（机制层）已并入 FigS3（2026-09-28 重构）----------
  #   原 figS6a–f 六面板及其上下叠拆件（figS6a1/a2 … figS6f1/f2）
  #   已统一重命名为 figS3a–i（尺寸沿用，见上方 FigS3 段）。
)

# ---- 出版级统一主题 ----------------------------------------------------
# base_size 默认 8 → axis.text 7 pt、axis.title 8 pt，均高于 6 pt 下限。
# v3 追加：plot.title 由 base_size+1 降为 base_size，并显式设定 plot.subtitle
#   = base_size-1（7 pt）。理由：面板半幅仅 3.25 in，9 pt 标题在 44 字符时宽
#   约 3.3 in > 可用 3.19 in，会被设备边缘裁断（FigS2A/FigS4C 等实测）。
theme_publication <- function(base_size = 8, base_family = "") {
  theme_classic(base_size = base_size, base_family = base_family) +
    theme(
      axis.line = element_line(size = 0.5, colour = "black"),
      axis.ticks = element_line(size = 0.5, colour = "black"),
      axis.text = element_text(colour = "black", size = base_size - 1),
      axis.title = element_text(colour = "black", size = base_size),
      legend.text = element_text(colour = "black", size = base_size - 1),
      legend.title = element_text(colour = "black", size = base_size),
      plot.title = element_text(colour = "black", size = base_size, face = "bold"),
      plot.subtitle = element_text(colour = "grey30", size = base_size - 1,
                                   lineheight = 1.0),
      # ★ v3 备查：曾试 `plot.title.position = "plot"`（让标题对齐整图而非绘图区，
      # 以解决「y 轴标签一长就把标题挤到面板右半、被右缘裁断」）——**已撤回**。
      # 原因：标题移到整图左缘后，会与溢出的纵排 y 轴标题叠字
      # （fig2a/fig2e/fig2g2、fig4a/fig4b 实测），代价大于收益。
      # 正确做法：保留默认 "panel"，改为**压缩标题字符数**或**把长 y 轴标签折成两行**
      # （后者同时给标题腾出宽度，见 fig4d.R）。
      panel.grid = element_blank(),
      plot.margin = margin(4, 4, 4, 4)
    )
}

# ---- 取面板 box（供不走 export_fig 的面板，如 ggsurvplot / grid 直存）----
panel_box <- function(stem, default = c(4, 3.5)) {
  b <- PANEL_BOX[[stem]]
  if (is.null(b)) {
    warning(sprintf("PANEL_BOX 未覆盖 %s ---- 沿用默认 %.2f x %.2f in", stem, default[1], default[2]))
    return(default)
  }
  b
}

# ---- 导出函数 ----------------------------------------------------------
# 统一输出 PDF + PNG（300 dpi）。
# ★ v3：若输出文件名命中 PANEL_BOX，则用查表尺寸覆盖脚本自带的 w/h，
#   保证面板按最终印刷尺寸绘制（拼版 1:1，不再缩放）。
# PNG 固定白底 + cairo 设备（Windows 默认设备会乱码 Δ/−/ρ 等 Unicode）
export_fig <- function(p, outfile, w = 4, h = 3.5) {
  stem <- sub("\\.(pdf|png)$", "", basename(outfile))
  if (!is.null(PANEL_BOX[[stem]])) {
    w <- PANEL_BOX[[stem]][1]
    h <- PANEL_BOX[[stem]][2]
  } else {
    warning(sprintf("PANEL_BOX 未覆盖 %s ---- 沿用脚本自带 w=%s h=%s", stem, w, h))
  }
  dir.create(dirname(outfile), recursive = TRUE, showWarnings = FALSE)
  ggsave(sub("\\.(pdf|png)$", ".pdf", outfile), plot = p, width = w, height = h, device = cairo_pdf)
  ggsave(sub("\\.(png|pdf)$", ".png", outfile), plot = p, width = w, height = h,
         dpi = 300, bg = "white", type = "cairo")
  cat(sprintf("  [export] %-8s %.2f x %.2f in  (=%d x %d px @300dpi)\n",
              stem, w, h, round(w * 300), round(h * 300)))
  invisible(outfile)
}
