# ============================================================
# 拼版公共库（v3 · Cancer Science 6.5×8 in 合规版）
# ------------------------------------------------------------
# 约定（依据 _work/p5_figures/LAYOUT_SPEC_v3.md）：
#   · 画布宽 W = 1949 px（165 mm @300 dpi），单页高上限 2400 px（203 mm）
#   · 面板「按最终印刷尺寸」绘制 → 此处 1:1 贴图，禁止二次缩放
#   · 行内字母带 62 px，字母 44 px ≈ 10.5 pt
#   · 行宽不足时以 round(i*W/n) 定位，允许 ≤1 px 重叠
# ============================================================
suppressWarnings(suppressMessages(library(magick)))

FIG_ROOT <- file.path(Sys.getenv("LPS_FIG_ROOT", "."))
OUT      <- file.path(FIG_ROOT, "out")
W        <- 1949   # 165 mm @300 dpi
HMAX     <- 2400   # 203 mm @300 dpi
PAD      <- 62     # 字母带高度
LSIZE    <- 44     # 字母字号 px ≈ 10.5 pt

rd <- function(f) image_read(file.path(OUT, f))

# ---- 一行拼合 ----------------------------------------------------------
# panels : 图块文件名（PNG）列表，应已按最终尺寸绘制
# letters: 与 panels 等长的字母；空串 "" 表示该位不标字母
make_row <- function(panels, letters = LETTERS[seq_along(panels)]) {
  n <- length(panels)
  imgs <- lapply(panels, rd)
  hs <- vapply(imgs, function(p) image_info(p)$height, numeric(1))
  row_h <- max(hs) + PAD
  canvas <- image_blank(W, row_h, color = "white")
  for (i in seq_len(n)) {
    x <- round((i - 1) * W / n)
    wi <- image_info(imgs[[i]])$width
    if (x + wi > W + 1) stop(sprintf(
      "面板 %s 宽 %d px，在 %d 列布局下超出画布（起点 %d，右缘 %d > %d）",
      panels[i], wi, n, x, x + wi, W))
    canvas <- image_composite(canvas, imgs[[i]], operator = "over",
                              offset = paste0("+", x, "+", PAD))
    if (nzchar(letters[i])) {
      canvas <- image_annotate(canvas, letters[i], size = LSIZE, font = "Arial",
                               weight = 700, gravity = "northwest",
                               location = paste0("+", x + 8, "+8"), color = "black")
    }
  }
  canvas
}

# ---- 整图输出 ----------------------------------------------------------
compose <- function(rows, outname, quiet = FALSE) {
  full <- image_append(image_join(rows), stack = TRUE)
  info <- image_info(full)
  if (info$width > W)  stop(sprintf("%s 宽 %d px 超出画布上限 %d", outname, info$width, W))
  if (info$height > HMAX) stop(sprintf(
    "%s 高 %d px 超出单页上限 %d（差 %d px）", outname, info$height, HMAX, info$height - HMAX))
  image_write(full, file.path(OUT, paste0(outname, ".png")),
              format = "png", density = "300x300")
  image_write(full, file.path(OUT, paste0(outname, ".pdf")),
              format = "pdf", density = "300x300")
  if (!quiet) cat(sprintf("  [%-5s] %5d x %5d px = %.2f x %.2f in @300dpi  (剩余高度 %d px)\n",
                          outname, info$width, info$height,
                          info$width / 300, info$height / 300, HMAX - info$height))
  invisible(info)
}

# ---- 供面板脚本反算字号下限 -------------------------------------------
# 说明：面板按最终尺寸绘制，base_size 即实际磅值；axis.text = base_size - 1。
pt_min <- function(base_size = 8) base_size - 1
