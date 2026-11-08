## =============================================================================
##  S34 · 6C：CellChat 细胞通讯（肿瘤程序 <-> TAM <-> 瘤旁脂肪）
##  ---------------------------------------------------------------------------
##  定位（依据 s34_6c_probe.txt 实测）：
##    CellChatDB.human 不含 CD36/PPARG/FABP4（脂质介导非蛋白互作）→ 6C **不做脂质轴**；
##    命中的是 SPP1-CD44 / MIF-CD74 / TGFB1 / IL10 / CCL2-CCR2 / CXCL12 / VEGF / ANGPTL4
##    → 6C = 「TAM <-> 肿瘤 的免疫抑制性通讯图谱」，用于补强 N3（M2 样命题）
##  粒度：L2 冻结口径；仅把 Macrophages 按 L3 拆为 TAM1 / TAM3（N3 需要两者对比）
##  下采样：每类按患者分层，每患者上限 = floor(CAP/患者数)，小类全取（保留 TAM1/TAM3 真实比例）
##  并行：future 保持 sequential（Windows 下 multisession 会崩）
##  输出：s34_6c_report.txt / s34_cellchat.rds / fig6_data/cellchat_*.csv
## =============================================================================
options(stringsAsFactors = FALSE)
suppressPackageStartupMessages({ library(Seurat); library(Matrix); library(CellChat); library(future) })
options(future.globals.maxSize = 6 * 1024^3)
future::plan("sequential")   ## ★ Windows 下 multisession 会崩；一律 sequential

W <- file.path(Sys.getenv("LPS_WORK_ROOT", "."), "_work")
FIG <- file.path(W, "p5_figures/fig6_data"); dir.create(FIG, showWarnings = FALSE)
LOG <- file.path(W, "s34_6c_report.txt"); if (file.exists(LOG)) file.remove(LOG)
A <- function(...) { s <- paste0(...); cat(s, "\n", sep = ""); cat(s, "\n", sep = "", file = LOG, append = TRUE) }
P <- function(x, ...) { capture.output(print(x, ...)) -> o; for (l in o) A(l) }

set.seed(20260924)
A("===============================================================================")
A("S34 6C CellChat (tumor programs <-> TAM <-> peritumoral fat)   ", format(Sys.time(), "%Y-%m-%d %H:%M:%S"))
A("===============================================================================")

obj <- readRDS(file.path(W, "s39_seurat.s6.rds")); md <- obj[[]]
A("loaded: ", paste(dim(obj), collapse = " x "))

## ---- labels：L2 冻结粒度，Macrophages 拆 TAM1/TAM3 -------------------------
lab <- as.character(md$L2)
lab[md$L2 == "Macrophages" & md$L3 == "TAM1"] <- "TAM1"
lab[md$L2 == "Macrophages" & md$L3 == "TAM3"] <- "TAM3"
md$lab <- lab

KEEP <- c("Adipo diff", "Invasion", "ECM remodeling", "Stemness", "Hypoxia",
          "TAM1", "TAM3", "Monocytes", "MDSC", "TAN", "Normal adipocytes", "Endothelial cells")
CAP <- 800
A("")
A(sprintf("-- 0. subsampling (CAP=%d per class, stratified by patient) --", CAP))
sel <- integer(0)
for (g in KEEP) {
  cells <- which(md$lab == g); pats <- md$sample[cells]; up <- unique(pats)
  per <- max(1, floor(CAP / length(up)))
  pick <- unlist(lapply(up, function(p) { v <- cells[pats == p]; if (length(v) > per) sample(v, per) else v }))
  A(sprintf("   %-20s total=%6d  patients=%2d  per_cap=%3d  kept=%5d", g, length(cells), length(up), per, length(pick)))
  sel <- c(sel, pick)
}
sel <- sort(sel)
A(sprintf("   -> selected %d cells (%.1f%% of library)", length(sel), 100 * length(sel) / ncol(obj)))

data.input <- GetAssayData(obj, assay = "RNA", layer = "data")[, sel]
meta <- data.frame(labels = factor(md$lab[sel], levels = KEEP), row.names = colnames(obj)[sel])
A(sprintf("-- input: %d genes x %d cells ; labels = %d", nrow(data.input), ncol(data.input), nlevels(meta$labels)))
P(table(meta$labels))

## ---- CellChat 主流程 -------------------------------------------------------
A("")
A("-- 1. createCellChat / DB --")
cc <- createCellChat(object = data.input, meta = meta, group.by = "labels")
cc@DB <- CellChatDB.human
A(sprintf("   DB = CellChatDB.human ; %d interactions ; %d pathways",
          nrow(cc@DB$interaction), length(unique(cc@DB$interaction$pathway_name))))
cc <- subsetData(cc)
A(sprintf("   after subsetData: %d signaling genes", nrow(cc@data.signaling)))

A("-- 2. identifyOverExpressedGenes / Interactions --")
cc <- identifyOverExpressedGenes(cc, thresh.p = 0.05)
cc <- identifyOverExpressedInteractions(cc)
A(sprintf("   over-expressed LR: %d", nrow(cc@LR$LRsig)))

A("-- 3. computeCommunProb (triMean, nboot=100) --")
t0 <- Sys.time()
cc <- computeCommunProb(cc, type = "triMean", nboot = 100)
A(sprintf("   done in %.1f min", as.numeric(difftime(Sys.time(), t0, units = "mins"))))
cc <- filterCommunication(cc, min.cells = 10)
cc <- computeCommunProbPathway(cc)
cc <- aggregateNet(cc)
saveRDS(cc, file.path(W, "s34_cellchat.rds"))
A("   saved: s34_cellchat.rds")

g <- levels(cc@idents)
A("")
A("-- A. group sizes --")
A(paste(sprintf("%s=%d", g, as.numeric(table(cc@idents))), collapse = " | "))

net <- cc@net$weight; dimnames(net) <- list(g, g)
cnt <- cc@net$count;  dimnames(cnt) <- list(g, g)
A("")
A("-- B. net$weight : sender(rows) x receiver(cols) --")
P(round(net, 3))
write.csv(round(net, 5), file.path(FIG, "cellchat_net_weight.csv"))
write.csv(round(cnt, 2), file.path(FIG, "cellchat_net_count.csv"))

tg <- intersect(c("Adipo diff", "Invasion", "ECM remodeling", "Stemness", "Hypoxia"), g)
mg <- intersect(c("TAM1", "TAM3"), g)
A("")
A("-- C. tumor programs -> TAM --")
P(round(net[tg, mg, drop = FALSE], 3))
A("-- C2. TAM -> tumor programs --")
P(round(net[mg, tg, drop = FALSE], 3))
A(sprintf("   rowSum(tumor->TAM) = %.3f ; colSum(TAM) = %.3f",
          sum(net[tg, mg]), sum(net[, mg])))

## ---- 通路层：source/target 强度 -------------------------------------------
prob <- cc@netP$prob   ## [source, target, pathway]
pws <- dimnames(prob)[[3]]
A("")
A(sprintf("-- D. pathways with signal (%d) --", length(pws)))
A(paste(pws, collapse = ", "))

## 每个通路的总强度（全群对求和），按降序
tot <- apply(prob, 3, function(m) sum(as.numeric(m), na.rm = TRUE))
o <- order(-tot)
A("")
A("-- D2. pathway total strength (top 30) --")
P(data.frame(pathway = pws[o][1:min(30, length(o))], strength = round(tot[o][1:min(30, length(o))], 3)),
  row.names = FALSE)
write.csv(data.frame(pathway = pws, strength = as.numeric(tot)),
          file.path(FIG, "cellchat_pathway_strength.csv"), row.names = FALSE)

## 重点通路：肿瘤<->TAM 双向
FOCUS <- intersect(c("SPP1", "MIF", "TGFb", "IL10", "CCL", "CXCL", "VEGF", "ANGPTL",
                     "ADIPONECTIN", "VISFATIN", "INSULIN", "GALECTIN", "MK", "THBS", "FN1",
                     "COLLAGEN", "LAMININ", "CD39", "IL1", "TNF", "PDGF", "EGF", "IGF", "ApoE"), pws)
A("")
A(sprintf("-- E. tumor <-> TAM : focus pathways (%d of focus set present) --", length(FOCUS)))
A(sprintf("%-14s %-12s %12s %12s %12s", "pathway", "dir", "AdipoDiff", "Invasion", "ECMremod"))
for (p in FOCUS) {
  for (dir in c("t2m", "m2t")) {
    strs <- character(0)
    for (a in tg) {
      v <- if (dir == "t2m") prob[a, mg, p] else prob[mg, a, p]
      strs <- c(strs, sprintf("%12.4f", sum(as.numeric(v), na.rm = TRUE)))
    }
    A(sprintf("%-14s %-12s %s", if (dir == "t2m") p else "", dir, paste(strs, collapse = " ")))
  }
}

## ---- LR 对层：TAM3 -> tumor 的 top LR -------------------------------------
lrprob <- cc@net$prob; dimnames(lrprob) <- list(g, g, dimnames(cc@net$prob)[[3]])
lrpw <- dimnames(cc@net$prob)[[3]]
A("")
A("-- F. top LR pairs: TAM3 -> tumor programs --")
agg <- apply(lrprob["TAM3", tg, , drop = FALSE], 3, function(x) sum(as.numeric(x), na.rm = TRUE))
oo <- order(-agg)
for (i in oo[1:min(25, length(oo))]) A(sprintf("   %-40s %.4f", lrpw[i], agg[i]))
A("")
A("-- G. top LR pairs: tumor -> TAM3 --")
agg2 <- apply(lrprob[tg, "TAM3", , drop = FALSE], 3, function(x) sum(as.numeric(x), na.rm = TRUE))
oo2 <- order(-agg2)
for (i in oo2[1:min(25, length(oo2))]) A(sprintf("   %-40s %.4f", lrpw[i], agg2[i]))

## ---- TAM1 vs TAM3：作为发送者的总强度 -------------------------------------
A("")
A("-- H. TAM1 vs TAM3 as SENDER (outgoing) / RECEIVER (incoming) --")
A(sprintf("   TAM1 out=%.3f  in=%.3f   TAM3 out=%.3f  in=%.3f",
          sum(net["TAM1", ]), sum(net[, "TAM1"]), sum(net["TAM3", ]), sum(net[, "TAM3"])))
A(sprintf("   per-cell out: TAM1=%.5f  TAM3=%.5f",
          sum(net["TAM1", ]) / as.numeric(table(cc@idents)["TAM1"]),
          sum(net["TAM3", ]) / as.numeric(table(cc@idents)["TAM3"])))
A("")
A("-- I. per-pathway TAM3 vs TAM1 outgoing (top by |diff|) --")
d <- sapply(pws, function(p) sum(as.numeric(prob["TAM3", , p]), na.rm = TRUE) -
                                sum(as.numeric(prob["TAM1", , p]), na.rm = TRUE))
od <- order(-abs(d))
for (i in od[1:min(25, length(od))])
  A(sprintf("   %-30s TAM3-TAM1 = %+8.4f  (TAM3=%.4f TAM1=%.4f)", pws[i], d[i],
            sum(as.numeric(prob["TAM3", , pws[i]]), na.rm = TRUE),
            sum(as.numeric(prob["TAM1", , pws[i]]), na.rm = TRUE)))

A("")
A("DONE")
