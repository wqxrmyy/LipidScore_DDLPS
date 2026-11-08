# ============================================================
# Fig 4E · 脂肪细胞身份基因热图（7 组 × 基因，组均值 z；官方 log2CPM）
# 数据：fig4_genes.csv（prep_fig4_gse_official.py，13 基因 log2CPM，CD36 均值=S5）
# ============================================================
suppressWarnings(suppressMessages(library(ggplot2)))
suppressWarnings(suppressMessages(library(reshape2)))
source(file.path(Sys.getenv("LPS_FIG_ROOT", "."), "utils", "colors.R"))
source(file.path(Sys.getenv("LPS_FIG_ROOT", "."), "utils", "themes.R"))

OUT <- file.path(Sys.getenv("LPS_FIG_ROOT", "."), "out")

d <- read.csv(file.path(Sys.getenv("LPS_FIG_ROOT", "."), "fig4_genes.csv"),
              stringsAsFactors = FALSE)
# 宽表：基因 × 样本
wl <- reshape(d[, c("sample", "gene", "expr")], idvar = "sample", timevar = "gene",
              direction = "wide")
names(wl) <- sub("^expr\\.", "", names(wl))
genes <- names(wl)[names(wl) != "sample"]
m <- as.matrix(wl[, genes])
rownames(m) <- wl$sample
# 每基因跨 53 样本 z 分
z <- t(scale(t(m)))

# 组均值
grp <- d[!duplicated(d$sample), c("sample", "group")]
grp <- grp[match(rownames(z), grp$sample), ]
zg <- data.frame(Group = grp$group, z, check.names = FALSE)
gm <- aggregate(. ~ Group, data = zg, FUN = mean)
rownames(gm) <- gm$Group
gm$Group <- NULL
gm <- as.matrix(t(gm))  # 基因 × 组

# 组排序：Fat > lipoma > Hibernoma > WDLPS > sWD > WD > DD（身份递减）
# ★ 2026-09-25 修复静默丢组：原写 "Lipoma"（大写），而 fig4_genes.csv 里是小写
#   "lipoma"（n=75，5 样本），被下一行 `g_order[g_order %in% colnames(gm)]`
#   **静默丢弃** —— 图上只画了 6 组，而副标题却写着 "declines from fat/lipoma to DD"，
#   图文自相矛盾。
#   处置：不用 `%in%` 过滤（它会吃掉错误），改成**双向断言** —— 少一个组、多一个组都
#   立即 stop()。宁可出图中断，也不要出一张"少了一组却看起来正常"的图。
g_order <- c("Fat", "lipoma", "Hibernoma", "WDLPS", "sWD", "WD", "DD")
.miss <- setdiff(g_order, colnames(gm))
.extra <- setdiff(colnames(gm), g_order)
if (length(.miss))
  stop("fig4e: g_order 里有数据中不存在的组（检查大小写/拼写）：", paste(.miss, collapse = ", "))
if (length(.extra))
  stop("fig4e: 数据里有未列入 g_order 的组（会被静默丢掉）：", paste(.extra, collapse = ", "))
gm <- gm[, g_order]
cat(sprintf("fig4e groups (%d): %s\n", ncol(gm), paste(colnames(gm), collapse = " | ")))

mm <- melt(gm); names(mm) <- c("Gene", "Group", "Z")

p <- ggplot(mm, aes(x = Group, y = Gene, fill = Z)) +
  geom_tile(color = "white") +
  # v4：2.6 mm = 7.37 pt（墨高约 22 px）> 原行距 24 px - 余量 → 数字上下相撞。
  #     降到 2.2 mm = 6.24 pt（墨高约 18.7 px），仍在 6 pt 下限之上。
  geom_text(aes(label = sprintf("%.1f", Z)), size = 2.2) +
  scale_fill_gradient2(low = "#2980B9", mid = "white", high = "#E74C3C", midpoint = 0) +
  labs(x = "", y = "", fill = "z", title = "Adipocyte identity genes across groups",
       # v3：原句 62 字符，半幅会被裁断；队列与量纲已由图注承载
       subtitle = "Identity declines from fat/lipoma to DD") +
  theme_publication() +
  # ============================================================
  # v5（2026-09-25，PI 裁决后）· 行距再放宽 + 补回静默丢弃的 lipoma 组
  # ------------------------------------------------------------
  # 一、行距：584 → 650 px（PANEL_BOX fig4e .half(0.6676)）
  #   实测响应函数（本主题下 chrome 恒定 232 px）：
  #     绘图区 ≈ 面板高 − 232 px  →  行距 = (h − 232) / 15
  #     h=584 → 行距 23.5 px、净空 4 px
  #     h=650 → 行距 28.0 px、净空 9 px   ← 现值（y 轴与数字两处均 9 px）
  #   ★ 预算来源：Fig4 原总高 2381 px（上限 2400，只余 19 px），E/F 各 +66 px
  #     → 另需 47 px，由 fig4c 支付（701 → 653 px）。fig4c 的 14 类行距
  #     35.0 → 31.5 px、净空 14 → 11 px（0.48 → 0.38 em，仍在下限 0.35 em 之上）。
  #     合计 Fig4 = 2399 px ≤ 2400 ✓（只余 1 px，再动任何一行都必须重算）
  #   ※ 650 px 仍差 0.5 px 未达 0.35 em 舒适线（墨高 19 px → 需行距 28.33 px）；
  #     真正达标需 657 px，但那会把 fig4c 压到 0.348 em（低于下限）→ 买不到。
  #     若日后还要放宽 E，杠杆不是行高，而是本题下方那 166 px 的 45° x 标签块
  #     （占 chrome 232 px 中的 166 px）→ 折行/改 0° 可回收约 124 px。
  # ------------------------------------------------------------
  # 二、★ 补回被静默丢弃的组（本文件顶部 g_order 处，已改成双向断言）
  #   原 g_order 写 "Lipoma"（大写），数据里是小写 "lipoma"（n=75 / 5 样本），
  #   被 `g_order[g_order %in% colnames(gm)]` 静默吃掉 → 只画了 6 组，
  #   而副标题写着 "declines from fat/lipoma to DD" —— **图文自相矛盾**。
  #   现为 7 组。几何影响已核验：行数不变（行=基因=15），列数 6→7 →
  #   格宽 101 → 87 px；45° 标签矩形不相交判据为 across = 格宽·sin45 ≥ 墨高，
  #   87×0.707 = 61.6 px ≥ 19 px ✓ 仍安全（判据见技能 §7.5，别用"占位宽度"当判据）。
  #   副作用（须记账）：lipoma 的 CD8A = −2.0 使 z 值域下限扩到 −2.0，
  #   色标刻度由 {−1, 0, 1} 变为 {−2, −1, 0, 1}。图注若引用了刻度需同步。
  # ------------------------------------------------------------
  # ★ 实测验收（口径：行墨心等距网格，**不是色标**）
  #   改前 565 px / 7 pt / 数字 2.6 mm：y 轴标签粘成 287 px 一整段、
  #                                     格内数字粘成 327 px 一整段 → FAIL
  #   只加高 584 px 其余不动      ：粘成 304 / 346 px → 仍 FAIL（证明只加高不够）
  #   现值 650 px / 6.4 pt / 2.2 mm：15 段等距，行距 28.0 px，净空 9 px → PASS
  # ※ 标尺口径：**不要用右侧色标量绘图区高度**。色标与绘图区尺寸解耦
  #   （实测偏移 20 px、高差 7 px），按它裁剪会把最上一行切掉，伪造出"14/15 粘连"。
  #   校验脚本：~/.workbuddy/skills/ggplot-publication-figure-qa/scripts/check_heatmap_overlap.py
  # ------------------------------------------------------------
  # ※ x 轴标签保持 45°：本题是「旋转是被支配策略」的**反例** —— 组名最长 Hibernoma
  #   在 6.4 pt 下墨宽约 120 px，0° 时格宽只有 87 px → 必然横向相撞。
  #   （判据：0° 时 across = 0 < 墨高 且 along = 格宽 87 < 标签宽 120 → 相交。）
  # ============================================================
  theme(
    axis.text.x = element_text(colour = "black", size = 6.4, angle = 45, hjust = 1),
    axis.text.y = element_text(colour = "black", size = 6.4),
    # chrome 压到极限：上下 margin 1 pt、ticks 1.5 pt、标题/副标题间距近乎为 0。
    # 这些是纯留白，不承载信息；每省 1 px 都直接变成行距。
    axis.ticks.length = unit(1.5, "pt"),
    plot.title = element_text(colour = "black", size = 7, face = "bold",
                              margin = margin(b = 0.5)),
    plot.subtitle = element_text(colour = "grey30", size = 6.4, lineheight = 1.0,
                                 margin = margin(t = 1, b = 0)),
    plot.margin = margin(1, 1, 1, 1)
  )

export_fig(p, file.path(OUT, "fig4e.pdf"), w = 6, h = 5)
cat("fig4e done\n")
