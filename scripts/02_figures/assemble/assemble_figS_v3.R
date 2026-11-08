# ============================================================
# 拼版补充图 FigS1–S5（v3 · Cancer Science 6.5×8 in 合规版）
# ------------------------------------------------------------
# ★★ 已停用（2026-09-21，PI 裁定）：拼版不再走脚本，改由人工完成。
#     交付物见 重构手稿/v2/figures/manual_layout/（坐标清单 + 面板件 + 坐标预览）。
#     保留本脚本仅作几何算法参考。
# ------------------------------------------------------------
# v3 与旧版（assemble_figS.R）的区别：
#   1) 画布宽 1500→1949 以上统一为 1949 px（165 mm）；旧版输出 3000 px = 254 mm 超宽；
#   2) 面板 1:1 贴图（面板已按 PANEL_BOX 最终尺寸绘制），不再等比缩放到统一宽；
#   3) 超宽/超页直接 stop。
# 版式：
#   FigS1–S4  2×2（半幅面板，各 974 px 宽）
#   FigS5     r1 A 整行（散点矩阵）| r2 B 整行（棒棒糖）
# 输出：out/FigS1.png ~ FigS5.png（+ PDF）
# ============================================================
source(file.path(Sys.getenv("LPS_FIG_ROOT", "."), "utils", "assemble_lib.R"))

WHICH <- c("FigS1", "FigS2", "FigS3", "FigS4", "FigS5")

grid22 <- function(stems, outname) {          # stems: 4 个面板文件名
  cat(outname, "...\n")
  r1 <- make_row(stems[1:2], c("A", "B"))
  r2 <- make_row(stems[3:4], c("C", "D"))
  compose(list(r1, r2), outname)
}

if ("FigS1" %in% WHICH) grid22(paste0("figS1", c("a", "b", "c", "d"), ".png"), "FigS1")
if ("FigS2" %in% WHICH) grid22(paste0("figS2", c("a", "b", "c", "d"), ".png"), "FigS2")
if ("FigS3" %in% WHICH) grid22(paste0("figS3", c("a", "b", "c", "d"), ".png"), "FigS3")
if ("FigS4" %in% WHICH) grid22(paste0("figS4", c("a", "b", "c", "d"), ".png"), "FigS4")

# —— FigS5：承接原 Fig1F/1G ——————————————
if ("FigS5" %in% WHICH) {
  cat("FigS5 ...\n")
  r1 <- make_row("figS5a.png", "A")
  r2 <- make_row("figS5b.png", "B")
  compose(list(r1, r2), "FigS5")
}

cat("ALL SUPPLEMENTARY FIGURES DONE\n")
