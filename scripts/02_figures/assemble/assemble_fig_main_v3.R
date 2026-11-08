# ============================================================
# 拼版正文主图 Fig 1–4（v3 · Cancer Science 6.5×8 in 合规版）
# ------------------------------------------------------------
# ★★ 已停用（2026-09-21，PI 裁定）：拼版不再走脚本，改由人工完成。
#     交付物见 重构手稿/v2/figures/manual_layout/（坐标清单 + 面板件 + 坐标预览）。
#     保留本脚本仅作几何算法参考；若日后恢复自动拼版，须以该坐标清单为验收基准。
#     面板渲染管线（themes.R 的 PANEL_BOX + 各 figXx.R）不受影响，仍然有效。
# ------------------------------------------------------------
# v3 与旧版（assemble_fig_main.R）的区别：
#   1) 画布宽 1949 px（165 mm），旧版 2126 px = 180 mm 已超宽；
#   2) ★ 面板 1:1 贴图 —— 面板已在 themes.R 的 PANEL_BOX 中按最终印刷尺寸
#      绘制，缩放会破坏字号（旧版把 11 pt 压到 4.5–6 pt，跌破官方 6 pt 下限）；
#   3) 行内字母带 62 px（旧版 90 px）；
#   4) 超宽/超页直接 stop，不再只警告（见 utils/assemble_lib.R）。
# 版式（各图高度预算均 < 2400 px）：
#   Fig1  r1 A(整行)|r2 B1|B2(整行内)|r3 C|D|r4 E|F          2362
#   Fig2  r1–r3 (A|B)(C|D)(E|F)|r4 G=g1|g2(整行)             2392
#   Fig3  r1–r3 (A|B)(C|D)(E|F)                              2378
#   Fig4  同 Fig2（r4 G 整行）                                2372
# 依据：_work/p5_figures/LAYOUT_SPEC_v3.md
# 输出：out/Fig1.png + Fig1.pdf ... Fig4.png + Fig4.pdf
# ============================================================
source(file.path(Sys.getenv("LPS_FIG_ROOT", "."), "utils", "assemble_lib.R"))

WHICH <- c("Fig1", "Fig2", "Fig4", "Fig3")

# —— Fig 1：A–F 六联图（原 1F/1G 已迁至 FigS5）————————————
if ("Fig1" %in% WHICH) {
  cat("Fig1 ...\n")
  r1 <- make_row("fig1a.png",          "A")
  r2 <- make_row(c("fig1b1.png", "fig1b2.png"), c("B1", "B2"))
  r3 <- make_row(c("fig1c.png",  "fig1d.png"),  c("C", "D"))
  r4 <- make_row(c("fig1e.png",  "fig1h.png"),  c("E", "F"))
  compose(list(r1, r2, r3, r4), "Fig1")
}

# —— Fig 2：A|B / C|D / E|F / G(g1|g2)整行 ——————————————
if ("Fig2" %in% WHICH) {
  cat("Fig2 ...\n")
  r1 <- make_row(c("fig2a.png", "fig2b.png"), c("A", "B"))
  r2 <- make_row(c("fig2c.png", "fig2d.png"), c("C", "D"))
  r3 <- make_row(c("fig2e.png", "fig2f.png"), c("E", "F"))
  r4 <- make_row(c("fig2g1.png", "fig2g2.png"), c("G", ""))
  compose(list(r1, r2, r3, r4), "Fig2")
}

# —— Fig 4：A|B / C|D / E|F / G整行 ——————————————
if ("Fig3" %in% WHICH) {
  cat("Fig3 ...\n")
  r1 <- make_row(c("fig3a.png", "fig3b.png"), c("A", "B"))
  r2 <- make_row(c("fig3c.png", "fig3d.png"), c("C", "D"))
  r3 <- make_row(c("fig3e.png", "fig3f.png"), c("E", "F"))
  compose(list(r1, r2, r3), "Fig3")
}

# —— Fig 3：A|B / C|D / E|F ——————————————
if ("Fig4" %in% WHICH) {
  cat("Fig4 ...\n")
  r1 <- make_row(c("fig4a.png", "fig4b.png"), c("A", "B"))
  r2 <- make_row(c("fig4c.png", "fig4d.png"), c("C", "D"))
  r3 <- make_row(c("fig4e.png", "fig4f.png"), c("E", "F"))
  r4 <- make_row("fig4g.png", "G")
  compose(list(r1, r2, r3, r4), "Fig4")
}

cat("ALL MAIN FIGURES DONE\n")
