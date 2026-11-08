# ============================================================
# Fig5 公共层：配色 + 数据读取 + 配对面板辅助函数
#   依据：LAYOUT_SPEC_v3.md §8.4（镜像行配色规范）
#     ① 同一谱系在 r3(E) 与 r4(F) 必须同色 → LIN_COL 只此一处定义
#     ② E1–E3 / F1–F3 的 WD/DD 用「空心/实心」区分（x 轴已有 WD/DD 标签，故不需图例）
#     ③ 第 4 列（E4/F4 关键证据）单独配色，不与谱系色冲突 → KEY_COL
#   数据一律来自 fig5_data/*.csv（由 s39_s9a_prep.py / s39_s9b_dotplot.R 生成）
# ============================================================
suppressWarnings(suppressMessages({library(ggplot2); library(scales)}))

FIG5_ROOT <- file.path(Sys.getenv("LPS_FIG_ROOT", "."))
FIG5_DATA <- file.path(FIG5_ROOT, "fig5_data")
FIG5_OUT  <- file.path(FIG5_ROOT, "out")
W5        <- file.path(Sys.getenv("LPS_WORK_ROOT", "."), "_work")

# ---- 谱系色（r3/r4 逐列同色；全图唯一真相源）----------------------------
LIN_COL <- c(Adipose = "#D4A017",   # 脂肪系（金）
             Tumor   = "#A93226",   # 肿瘤系（砖红）
             Immune  = "#7D3C98")   # 免疫系（紫）
LIN_LAB <- c(Adipose = "Adipocytic", Tumor = "Tumour", Immune = "Immune")

# ---- 第 4 列关键证据：中性深青，不与上面的谱系色冲突 --------------------
KEY_COL <- "#1B6B7A"

# ---- 组成项 / 状态项（仅 F4 用；两者同现，不存在与 WD/DD 混读）---------
CT_COL <- c(C = "#1B7A8C", T = "#C0392B")

# ---- L1 细胞类型色（面板 A / D 共用：同一细胞类型同色）------------------
L1_COL <- c("Normal adipocytes" = "#D4A017",
            "Endothelial cells" = "#1ABC9C",
            "Pericytes"         = "#2E86C1",
            "Tumor cells"       = "#A93226",
            "RBC"               = "#78281F",
            "Myeloid cells"     = "#E67E22",
            "Mast cells"        = "#8E6E53",
            "Cycling cells"     = "#D35DB3",
            "Lymphoid cells"    = "#7D3C98",
            "Mesothelial cells" = "#7F8C8D")

# ---- 读入并固定 L1 顺序（按平均评分降序，图例/点图/箱线共用同一顺序）----
fig5_l1_order <- function() {
  um <- read.csv(file.path(FIG5_DATA, "umap_cells.csv.gz"), stringsAsFactors = FALSE)
  names(um)[names(um) == "celltype"] <- "L1"
  ag <- aggregate(score_ams ~ L1, data = um, FUN = mean)
  ag$L1[order(-ag$score_ams)]
}

fig5_read <- function(f) read.csv(file.path(FIG5_DATA, f), stringsAsFactors = FALSE,
                                  fileEncoding = "UTF-8-BOM")

# ---- 桶色（面板 D 用；比 LIN_COL 多一个 Stroma，r3/r4 未出现该谱系列）----
BUCKET_COL <- c(Adipose = "#D4A017", Tumor = "#A93226", Immune = "#7D3C98",
                Stroma = "#2E86C1")

# ---- 面板 C/D 的 L2 类清单（★ 唯一规则，两面板共用，防两处漂移）--------
#   规则：占全体细胞 ≥3% 的类（10 个，合计 80.4%）
#         ∪ 类平均评分最高的 5 个类（补入 Normal adipocytes / Macrophages / Pericytes）
#   理由：单纯按细胞数取会漏掉评分最高的 Normal adipocytes（仅占 0.97%）→ 图会误导
#   ⚠️ 未纳入的其余 11 类（合计约 15%）须在图注中披露，不得静默省略
fig5_l2_classes <- function(min_pct = 3, n_top = 5) {
  cs <- read.csv(file.path(W5, "s6", "cell_scores.csv.gz"), stringsAsFactors = FALSE)
  ag <- aggregate(score_ams ~ bucket + L2, data = cs,
                  FUN = function(x) c(n = length(x), m = mean(x)))
  ag <- data.frame(bucket = ag$bucket, L2 = ag$L2,
                   n = ag$score_ams[, "n"], mean_score = ag$score_ams[, "m"])
  ag$pct <- 100 * ag$n / nrow(cs)
  keep <- ag$pct >= min_pct | rank(-ag$mean_score) <= n_top
  ag <- ag[keep, ]
  list(tab = ag, coverage = sum(ag$pct),
       omitted = setdiff(unique(cs$L2), ag$L2))
}

# 显示用短标签（仅缩短字面，不改变类别含义；与 atlas 注释一一对应）
L2_SHORT <- c("Normal adipocytes" = "Normal adipocytes",
              "Endothelial cells" = "Endothelial",
              "ECM remodeling"   = "ECM remod.",
              "Adipo diff"       = "Adipo diff",
              "Macrophages"      = "Macrophages",
              "Pericytes"        = "Pericytes",
              "Invasion"         = "Invasion",
              "Hypoxia"          = "Hypoxia",
              "Monocytes"        = "Monocytes",
              "Stemness"         = "Stemness",
              "MDSC"             = "MDSC",
              "T cells"          = "T cells",
              "NK cells"         = "NK")

# ---- ★ 面板 C/D 的「定义轴」分组（PI 裁定 2026-09-22；唯一真相源）--------
#   规范（PI，2026-09-22）：「同一 y 轴或 x 轴的参数须同一类型；细胞是细胞，
#     与抽象参数不能混用」→ 13 个 L2 类不能并列在一根轴上，须按定义轴分子轴。
#   分轴依据 = **表型/来源口径**（与 BUCKET_COL 的色标对齐，避免轴分组与色标打架）：
#     - Adipocytic(2)：Adipo diff, Normal adipocytes —— 脂肪细胞样表型
#     - Immune & stromal(7)：7 个非脂肪谱系身份（BUCKET_COL 里为 Immune + Stroma 两桶）
#     - Tumour programs(4)：4 个肿瘤细胞转录程序（源注释母级均为 Tumor cells）
#   ★ Adipo diff 归 Adipocytic 侧（PI 方案 X，表型口径）——
#     注意与源注释 L1（Tumor cells）不同，属再归属，见 LIM-P3-015，
#     须在 caption / Methods 声明。这也使面板 D 的轴分组与桶色（金）自洽。
AXIS_OF_L2 <- c("Adipo diff"        = "Adipocytic",
                "Normal adipocytes" = "Adipocytic",
                "T cells"           = "Immune & stromal",
                "NK cells"          = "Immune & stromal",
                "Monocytes"         = "Immune & stromal",
                "Macrophages"       = "Immune & stromal",
                "MDSC"              = "Immune & stromal",
                "Endothelial cells" = "Immune & stromal",
                "Pericytes"         = "Immune & stromal",
                "Hypoxia"           = "Tumour programs",
                "Invasion"          = "Tumour programs",
                "Stemness"          = "Tumour programs",
                "ECM remodeling"    = "Tumour programs")
AXIS_LEV <- c("Adipocytic", "Immune & stromal", "Tumour programs")

# ---- ★ r3/r4 的共享刻度（§8.4「同一行内 y 轴范围必须一致」）-------------
#   写成常量而非在各脚本里重复 → 结构性保证不会漂移
E_Y_FRAC  <- c(0, 1)                  # 组成层：谱系占比（0–1）
F_Y_SCORE <- c(-0.10, 0.175)          # 状态层：谱系内平均评分（据实测值域 [−0.081, +0.149] 定）

# ---- 配对面板（E1–E3 / F1–F3 通用）-------------------------------------
#   编码：连线与点边 = 谱系色；空心 = WD、实心 = DD（x 轴已标 WD/DD，故不需图例）
#   数值全部来自数据；平均变化量实时计算，不硬编码
#   ★ 平均变化量放**副标题**而非面板内 annotate —— 面板内注释在部分谱系会与数据点
#     相撞（E3 免疫系实测：DD 最高点 0.93 与注释重叠）。副标题在面板外，结构性免疫。
#   ★ 图内文本一律 ASCII（项目铁律：图内禁 Unicode，避免 PDF/PNG 渲染乱码）
fig5_pair_panel <- function(d, ycol, lin, title, ylims, dlab, mult = 1, unit = "", dgits = 1) {
  d <- d[order(d$sample), ]
  d$x <- ifelse(d$tissue == "DDLPS-WD", 1, 2)
  col <- unname(LIN_COL[[lin]])
  ann <- sprintf("mean change %+.*f%s", dgits, dlab * mult, unit)
  ggplot(d, aes(x = x, y = .data[[ycol]])) +
    geom_line(aes(group = patient), colour = col, linewidth = 0.4, alpha = 0.55) +
    geom_point(aes(fill = tissue), shape = 21, size = 1.9, stroke = 0.5, colour = col) +
    scale_fill_manual(values = c("DDLPS-WD" = "white", "DDLPS-DD" = col)) +
    scale_x_continuous(breaks = c(1, 2), labels = c("WD", "DD"),
                       limits = c(0.62, 2.38)) +
    coord_cartesian(ylim = ylims) +
    labs(x = NULL, y = NULL, title = title, subtitle = ann) +
    theme_publication() +
    theme(legend.position = "none",
          plot.subtitle = element_text(size = 6.4, colour = "grey25"),
          axis.text.x = element_text(size = 6.6),
          axis.text.y = element_text(size = 6.2),
          plot.margin = margin(2, 2, 2, 2))
}

# 配对平均变化量（WD → DD），由数据实时计算
fig5_pair_delta <- function(d, ycol) {
  d <- d[d$patient %in% c("P1", "P4", "P7", "P10", "P14", "P15", "P16"), ]
  w <- d[d$tissue == "DDLPS-WD", c("patient", ycol)]
  dd <- d[d$tissue == "DDLPS-DD", c("patient", ycol)]
  m <- merge(w, dd, by = "patient", suffixes = c(".w", ".d"))
  mean(m[[paste0(ycol, ".d")]] - m[[paste0(ycol, ".w")]])
}

# 7 对患者（严格成对，顺序固定）
FIG5_PAIRS <- c("P1", "P4", "P7", "P10", "P14", "P15", "P16")

