suppressPackageStartupMessages({library(Seurat); library(Matrix); library(future)})
plan("sequential")                      # ★ 关键：multisession 会因 globals 超限而崩
options(future.globals.maxSize = 4 * 1024^3)

W <- file.path(Sys.getenv("LPS_WORK_ROOT", "."), "_work")
LOG <- file.path(W, "s39_s3b1_report.txt")
sink(LOG, split = TRUE)
say <- function(...) cat(sprintf(...), "\n")
tt <- function(t0) as.numeric(difftime(Sys.time(), t0, units = "secs"))

say("== Step 3b-1  过滤 + 归一化 + HVG + Scale + PCA ==")
say("start: %s   [plan = sequential]", format(Sys.time()))

so <- readRDS(file.path(W, "s39_seurat.qc.rds"))
say("loaded: %d genes x %d cells", nrow(so), ncol(so))

# ---------- 1. QC 过滤 ----------
NF_MIN <- 200; MT_MAX <- 20
md0 <- so@meta.data
n0 <- ncol(so)
keep <- md0$nFeature_RNA >= NF_MIN & md0$percent.mt <= MT_MAX
so <- subset(so, cells = colnames(so)[keep])
rm(md0); gc(verbose = FALSE)
say("\ncriterion: nFeature_RNA >= %d AND percent.mt <= %d", NF_MIN, MT_MAX)
say("kept %d / %d (%.2f%%)", ncol(so), n0, 100 * ncol(so) / n0)

# ---------- 2. 归一化 ----------
say("\n-- NormalizeData (LogNormalize, 1e4) --")
t0 <- Sys.time()
so <- NormalizeData(so, normalization.method = "LogNormalize", scale.factor = 1e4, verbose = FALSE)
say("done in %.1fs", tt(t0))

# ---------- 3. HVG（多取后剔除技术性/谱系私有基因）----------
say("\n-- FindVariableFeatures (vst, n=3000) --")
t0 <- Sys.time()
so <- FindVariableFeatures(so, selection.method = "vst", nfeatures = 3000, verbose = FALSE)
say("done in %.1fs", tt(t0))
hv_all <- VariableFeatures(so)
say("原始 top 30: %s", paste(head(hv_all, 30), collapse = ", "))

EXCL <- "^(IG[HKL]|MT-|HB[ABDEGMQZ]|RP[SL][0-9]|MTRNR|RNA5S|RNA18S|RNA28S)"
bad <- grepl(EXCL, hv_all)
say("被剔除的技术性/谱系私有基因 %d 个（%s...）", sum(bad),
    paste(head(hv_all[bad], 12), collapse = ", "))
hv <- hv_all[!bad][1:2000]
VariableFeatures(so) <- hv
say("最终 HVG %d 个，top 30: %s", length(hv), paste(head(hv, 30), collapse = ", "))

score6 <- c("SREBF1", "SCD", "HMGCR", "ACLY", "PPARG", "FASN")
say("评分 6 基因入选情况: %s", paste(sprintf("%s=%s", score6, score6 %in% hv), collapse = ", "))

# ---------- 4. Scale + PCA ----------
say("\n-- ScaleData --")
t0 <- Sys.time()
so <- ScaleData(so, features = hv, verbose = FALSE)
say("done in %.1fs", tt(t0))

say("\n-- RunPCA (npcs = 50) --")
t0 <- Sys.time()
so <- RunPCA(so, features = hv, npcs = 50, verbose = FALSE)
say("done in %.1fs", tt(t0))
ve <- so[["pca"]]@stdev^2
pct <- 100 * ve / sum(ve)
say("PC1-10 方差解释(%%): %s", paste(sprintf("%.1f", pct[1:10]), collapse = ", "))
say("累计 PC1-20 = %.1f%% ; PC1-30 = %.1f%% ; PC1-50 = %.1f%%",
    sum(pct[1:20]), sum(pct[1:30]), sum(pct[1:50]))
say("PC1 前 15 基因: %s", paste(head(rownames(Loadings(so, "pca")), 15), collapse = ", "))

# ---------- 5. 落盘 ----------
out <- file.path(W, "s39_seurat.s3.rds")
say("\nsaving -> %s", out)
t0 <- Sys.time()
saveRDS(so, out, compress = "gzip")
say("saved: %.1f MB in %.1fs", file.size(out) / 1048576, tt(t0))
say("\nDONE_S3B1")
sink()
