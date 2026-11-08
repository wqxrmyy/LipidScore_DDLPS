suppressPackageStartupMessages({
  library(Seurat); library(Matrix); library(harmony); library(future); library(FNN)
})
# ★ Step 3 实测：Windows 上 plan("multisession") 会因 future globals 超限崩溃，
#   且 2 GB 级对象反复序列化反而慢 20 倍。全程单进程。
plan("sequential")

W   <- file.path(Sys.getenv("LPS_WORK_ROOT", "."), "_work")
LOG <- file.path(W, "s39_s4_report.txt")
sink(LOG, split = TRUE)
say <- function(...) cat(sprintf(...), "\n")
t_start <- Sys.time()

say("== Step 4  批次校正（batch = 患者，非样本）==")
say("start: %s", format(t_start))
say("")
say("设计依据：sample = patient x tissue。harmony 按 patient 校正整体偏移；")
say("若按 sample 校正，会把同一患者 WD/DD 之间的真实生物学差异一并抹掉。")
say("本步定位：稳健性验证 + 可视化底图。核心判决（Q2 Kitagawa 分解）不依赖嵌入。")
say("工具：FNN::get.knn（KD-tree），避免展开 n x n 距离矩阵。")

so <- readRDS(file.path(W, "s39_seurat.s3.rds"))
say("\nloaded: %d genes x %d cells", nrow(so), ncol(so))

# ================= [0] metadata 核对 =================
say("\n-- [0] metadata 核对 --")
say("cols:       %s", paste(colnames(so@meta.data), collapse = ", "))
say("assays:     %s", paste(names(so@assays), collapse = ", "))
say("reductions: %s", paste(names(so@reductions), collapse = ", "))
say("本步使用：  patient / tissue / L1 / L2 / L3")

gd <- so@meta.data[["grade"]]
if (!is.null(gd)) {
  uv <- sort(unique(as.character(gd)))
  say("\n※ 记账：metadata 的 grade 列取值 %d 个，例：%s",
      length(uv), paste(utils::head(uv, 10), collapse = ", "))
  say("  这不是 1/2/3 病理分级（取值与 L3 细胞类型数一一对应，如 11777 = Adipo diff_1）；")
  say("  疑为注释文件里的类型编码被误命名为 grade。本步不使用该列；")
  say("  原始 barcodes 列定义已另行核对（见 s39_s4_barcode_cols.txt）。")
}

# ================= [1] 指标函数 =================
DIMS <- 1:30
KK   <- 30   # kNN 邻居数

mix_stats <- function(emb, batch, l1, l3, tissue, k = KK) {
  nn <- FNN::get.knn(emb, k = k)$nn.index
  n  <- nrow(emb)
  f  <- function(v) {
    v <- as.character(v)
    mean(rowMeans(matrix(v[nn], nrow = n) == v))
  }
  c(same_patient = f(batch), same_L1 = f(l1), same_L3 = f(l3), same_tissue = f(tissue))
}
rand_base <- function(x) { p <- table(x) / length(x); sum(p^2) }

# 二分类 kNN 判别，返回 acc / 多数类基线 / 增益
classify_knn <- function(emb, lab, k = 15) {
  lab <- as.character(lab)
  if (length(unique(lab)) < 2) return(NULL)
  nn   <- FNN::get.knn(emb, k = k)$nn.index
  pred <- apply(matrix(lab[nn], nrow = nrow(emb)), 1,
                function(x) { tb <- table(x); names(tb)[which.max(tb)] })
  acc  <- mean(pred == lab)
  base <- max(prop.table(table(lab)))
  c(acc = acc, base = base, gain = acc - base,
    ngain = (acc - base) / max(1 - base, 1e-9), n = length(lab))
}

say("\n-- [1] 随机基线（按当前细胞组成）--")
say("  同患者 %.4f | 同 L1 %.4f | 同 L3 %.4f | 同组织 %.4f",
    rand_base(so$patient), rand_base(so$L1), rand_base(so$L3), rand_base(so$tissue))

# ================= [2] 校正前 =================
pca_emb <- Embeddings(so, "pca")[, DIMS]
say("\n-- [2] 校正前（PCA 1:30）kNN(%d) --", KK)
t0 <- Sys.time()
m0 <- mix_stats(pca_emb, so$patient, so$L1, so$L3, so$tissue)
say("  同患者 %.4f | 同 L1 %.4f | 同 L3 %.4f | 同组织 %.4f   [%.1fs]",
    m0["same_patient"], m0["same_L1"], m0["same_L3"], m0["same_tissue"],
    as.numeric(difftime(Sys.time(), t0, units = "secs")))

# ================= [3] harmony =================
say("\n-- [3] RunHarmony (group.by.vars = patient) --")
say("  harmony %s；RunHarmony 为 S3 泛型，参数名用候选链实测：",
    as.character(packageVersion("harmony")))
cands <- list(
  list(group.by.vars = "patient", reduction = "pca", dims.use = DIMS,
       reduction.save = "harmony", ncores = 8, verbose = FALSE),
  list(group.by.vars = "patient", reduction = "pca", dims.use = DIMS,
       reduction.save = "harmony", verbose = FALSE),
  list(group.by.vars = "patient", reduction = "pca", dims.use = DIMS,
       reduction.name = "harmony", verbose = FALSE),
  list(group.by.vars = "patient", reduction = "pca", dims.use = DIMS,
       reduction.name = "harmony"),
  list(group.by.vars = "patient", reduction = "pca", dims.use = DIMS),
  list(group.by.vars = "patient", reduction = "pca", dims = DIMS,
       reduction.name = "harmony"),
  list(group.by.vars = "patient", reduction = "pca")
)
t0 <- Sys.time(); so_new <- NULL; hit <- NA
for (i in seq_along(cands)) {
  r <- tryCatch(do.call(RunHarmony, c(list(object = so), cands[[i]])),
                error = function(e) {
                  say("  [try %d] 失败: %s", i, substr(conditionMessage(e), 1, 150)); NULL
                })
  if (!is.null(r)) { so_new <- r; hit <- i; break }
}
if (is.null(so_new)) stop("RunHarmony 全部候选失败")
so <- so_new
say("  [OK] 候选 %d 成功；reductions: %s", hit, paste(names(so@reductions), collapse = ", "))
HRED <- if ("harmony" %in% names(so@reductions)) "harmony" else
        setdiff(names(so@reductions), c("pca", "umap"))[1]
say("  使用 reduction: %s | 耗时 %.1fs", HRED,
    as.numeric(difftime(Sys.time(), t0, units = "secs")))

# ================= [4] 校正后 =================
har_emb <- Embeddings(so, HRED)[, DIMS]
say("\n-- [4] 校正后（harmony 1:30）kNN(%d) --", KK)
t0 <- Sys.time()
m1 <- mix_stats(har_emb, so$patient, so$L1, so$L3, so$tissue)
say("  同患者 %.4f | 同 L1 %.4f | 同 L3 %.4f | 同组织 %.4f   [%.1fs]",
    m1["same_patient"], m1["same_L1"], m1["same_L3"], m1["same_tissue"],
    as.numeric(difftime(Sys.time(), t0, units = "secs")))

say("\n-- [5] 对照表 --")
say("  指标                 校正前    校正后   随机基线")
say("  同患者邻居比例 (低好) %.4f    %.4f    %.4f", m0["same_patient"], m1["same_patient"], rand_base(so$patient))
say("  同 L1 邻居比例 (高好)  %.4f    %.4f    %.4f", m0["same_L1"], m1["same_L1"], rand_base(so$L1))
say("  同 L3 邻居比例 (高好)  %.4f    %.4f    %.4f", m0["same_L3"], m1["same_L3"], rand_base(so$L3))
say("  同组织邻居比例 (高好)  %.4f    %.4f    %.4f", m0["same_tissue"], m1["same_tissue"], rand_base(so$tissue))

# ================= [6] 过度校正检验 =================
PAIRS <- c("P1", "P4", "P7", "P10", "P14", "P15", "P16")
say("\n-- [6] ★ 过度校正检验：同患者内 WD vs DD 可分性 --")
say("  判读：acc 应显著高于多数类基线 base；若 acc 回落到 base 附近 => 真实差异被抹平。")
say("  （不可用 0.5 作基线：配对患者 WD/DD 细胞数不均衡）")
say("  %-5s %7s %7s | %9s %9s %9s | %9s %9s %9s", "pat", "nWD", "nDD",
    "acc_pre", "base", "gain_pre", "acc_post", "gain_post", "ngain_post")
PRE <- POST <- NULL
for (p in PAIRS) {
  sel <- which(so$patient == p & so$tissue %in% c("DDLPS-WD", "DDLPS-DD"))
  lab <- ifelse(so$tissue[sel] == "DDLPS-WD", "WD", "DD")
  a <- classify_knn(pca_emb[sel, , drop = FALSE], lab)
  b <- classify_knn(har_emb[sel, , drop = FALSE], lab)
  say("  %-5s %7d %7d | %9.3f %9.3f %9.3f | %9.3f %9.3f %9.3f",
      p, sum(lab == "WD"), sum(lab == "DD"),
      a["acc"], a["base"], a["gain"], b["acc"], b["gain"], b["ngain"])
  PRE  <- rbind(PRE,  a); POST <- rbind(POST, b)
}
say("  平均：acc %.3f -> %.3f | base %.3f | 增益 gain %.3f -> %.3f | 归一化增益 %.3f -> %.3f",
    mean(PRE[, "acc"]), mean(POST[, "acc"]), mean(PRE[, "base"]),
    mean(PRE[, "gain"]), mean(POST[, "gain"]),
    mean(PRE[, "ngain"]), mean(POST[, "ngain"]))

# ================= [7] Adipo diff 内部（描述性） =================
say("\n-- [7] ★ 描述性观察：Adipo diff 内部 WD vs DD 可分性（不作判决）--")
say("  说明：若 WD/DD 差异纯由「组成」驱动，则 Adipo diff 类内部应难以区分；")
say("  本项仅供 Step 6/7 对照，判决以事先登记的分解线为准。")
adipo <- which(so$L2 == "Adipo diff" & so$patient %in% PAIRS &
               so$tissue %in% c("DDLPS-WD", "DDLPS-DD"))
say("  Adipo diff（7 配对患者内）n = %d", length(adipo))
say("  %-5s %7s %7s | %9s %9s %9s | %9s %9s", "pat", "nWD", "nDD",
    "acc_pre", "acc_post", "base", "gain_pre", "gain_post")
gp <- gq <- numeric(0)
for (p in PAIRS) {
  sel <- adipo[so$patient[adipo] == p]
  if (length(sel) < 100) { say("  %-5s  n=%d <100，跳过", p, length(sel)); next }
  lab <- ifelse(so$tissue[sel] == "DDLPS-WD", "WD", "DD")
  if (length(unique(lab)) < 2) { say("  %-5s  单一组织，跳过", p); next }
  a <- classify_knn(pca_emb[sel, , drop = FALSE], lab)
  b <- classify_knn(har_emb[sel, , drop = FALSE], lab)
  say("  %-5s %7d %7d | %9.3f %9.3f %9.3f | %9.3f %9.3f",
      p, sum(lab == "WD"), sum(lab == "DD"),
      a["acc"], b["acc"], a["base"], a["gain"], b["gain"])
  gp <- c(gp, a["gain"]); gq <- c(gq, b["gain"])
}
if (length(gp)) say("  平均增益：校正前 %+.3f -> 校正后 %+.3f", mean(gp), mean(gq))

say("\n-- [7b] 跨患者可迁移性（留一患者：其余 6 例作参考，预测目标患者）--")
transfer <- function(emb, pool, k = 15) {
  res <- NULL
  for (p in PAIRS) {
    tr <- pool[so$patient[pool] != p]; te <- pool[so$patient[pool] == p]
    if (length(te) < 100 || length(tr) < 100) next
    ltr <- ifelse(so$tissue[tr] == "DDLPS-WD", "WD", "DD")
    lte <- ifelse(so$tissue[te] == "DDLPS-WD", "WD", "DD")
    if (length(unique(ltr)) < 2 || length(unique(lte)) < 2) next
    nn <- FNN::get.knnx(emb[tr, , drop = FALSE], emb[te, , drop = FALSE], k = k)$nn.index
    pred <- apply(matrix(ltr[nn], nrow = length(te)), 1,
                  function(x) { tb <- table(x); names(tb)[which.max(tb)] })
    res <- rbind(res, c(acc = mean(pred == lte),
                        base = max(prop.table(table(lte))), n = length(te)))
  }
  if (is.null(res)) return(NULL)
  c(acc = mean(res[, "acc"]), base = mean(res[, "base"]),
    gain = mean(res[, "acc"] - res[, "base"]), n_test = sum(res[, "n"]))
}
pool_all <- which(so$patient %in% PAIRS & so$tissue %in% c("DDLPS-WD", "DDLPS-DD"))
tp0 <- transfer(pca_emb, pool_all); tp1 <- transfer(har_emb, pool_all)
say("  全部 WD/DD 细胞：校正前 acc %.3f (base %.3f, gain %+.3f) -> 校正后 acc %.3f (base %.3f, gain %+.3f)",
    tp0["acc"], tp0["base"], tp0["gain"], tp1["acc"], tp1["base"], tp1["gain"])
ta0 <- transfer(pca_emb, adipo); ta1 <- transfer(har_emb, adipo)
if (!is.null(ta0) && !is.null(ta1)) {
  say("  仅 Adipo diff   ：校正前 acc %.3f (base %.3f, gain %+.3f) -> 校正后 acc %.3f (base %.3f, gain %+.3f)",
      ta0["acc"], ta0["base"], ta0["gain"], ta1["acc"], ta1["base"], ta1["gain"])
}

# ================= [8] harmony-UMAP =================
say("\n-- [8] RunUMAP (reduction = %s, dims = 1:30) --", HRED)
t0 <- Sys.time()
so <- RunUMAP(so, reduction = HRED, dims = DIMS,
              reduction.name = "umap.harmony", verbose = FALSE)
say("  done in %.1fs", as.numeric(difftime(Sys.time(), t0, units = "secs")))

purity <- function(emb, lab, n = 3000, k = 10, seed = 42) {
  set.seed(seed)
  ii <- sample(seq_len(nrow(emb)), min(n, nrow(emb)))
  e <- emb[ii, , drop = FALSE]; v <- as.character(lab[ii])
  nn <- FNN::get.knn(e, k = k)$nn.index
  mean(rowMeans(matrix(v[nn], nrow = length(ii)) == v))
}
uh <- Embeddings(so, "umap.harmony")
base_l3 <- rand_base(so$L3)
say("\n  L3 kNN 纯度（抽样3000, k=10）：随机基线 %.4f | harmony-UMAP %.4f (%.1fx)",
    base_l3, purity(uh, so$L3), purity(uh, so$L3) / base_l3)
say("  L1 kNN 纯度                  ：随机基线 %.4f | harmony-UMAP %.4f (%.1fx)",
    rand_base(so$L1), purity(uh, so$L1), purity(uh, so$L1) / rand_base(so$L1))

# ================= [9] 保存 =================
out <- file.path(W, "s39_seurat.s4.rds")
say("\nsaving -> %s", out)
t0 <- Sys.time()
saveRDS(so, out, compress = "gzip")
say("saved: %.1f MB in %.1fs", file.size(out) / 1048576,
    as.numeric(difftime(Sys.time(), t0, units = "secs")))
say("\ntotal elapsed: %.1f min", as.numeric(difftime(Sys.time(), t_start, units = "mins")))
say("\nDONE_S4")
sink()
