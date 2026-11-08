# s211 — D3(a) 用 TCGA-SARC 官方 gene-level 拷贝数**独立重导**聚类，替代"未重导"的 limitation。
#
# 输入
# ----
#   _work/s199_cnv/<participant>.tsv   GDC ASCAT2 gene-level copy number
#    列：gene_id / gene_name / chromosome / start / end / copy_number / min_copy_number / max_copy_number
#    copy_number 为**整数总拷贝数**
#
# v4 修订（2026-09-28）
# -------------------
#   [1] 基因覆盖体检（关键，决定了样本取舍）
#       GDC 这批 ASCAT2 文件并非都完整。实测：
#         TCGA-DX-A3U6  26,147 基因，只覆盖 chr1-chr8          → **截断**（末行残行）
#         TCGA-DX-AB36  24,848 基因，只覆盖 chr1-chr8          → **截断**（末行 = ENSG000001，chromosome 缺失）
#         另 2 例        约 60,200 基因，覆盖 chr1-chr22 但缺 chrX
#         余 46 例       覆盖 chr1-chr22 + chrX（含 chrY）
#       ⇒ **不可**用「全 50 例交集」当主分析：交集被截断样本压到 24,234 基因且全部落在 chr1-chr8，
#          整段丢掉 chr12q（DDLPS 的 12q13-15 扩增、MDM2/CDK4 所在）——口径不可用。
#       ⇒ 主分析 46 例（染色体覆盖完整）× 交集；敏感性 48 例（再纳入仅缺 chrX 者）。
#   [2] 标签重编号 bug 修正
#       `order(tapply(...))` 丢 names → 原写法 `names(ordc)[ordc]` 恒为 NULL → 重导标签全 NA。
#   [3] **全 k 矩阵报告（防选择性报告）**
#       两种特征口径 × k = 2..6 全部评估：规模 / r(score) / KW / Cramér's V vs 已发表 / ARI vs 已发表。
#       特征口径：(a) 基因级 top-5000 高变异基因；(b) 染色体级（23 条染色体中位 CN）——贴近 DDLPS
#       的臂级扩增生物学，也是 SARC 类 CN 分型的常规尺度。
#   [4] 非有限值防护：CN = 0 时 log2 出 -Inf，整行剔除并报数。
#
# 用法：Rscript s211_D3_cnv_consensus_clustering.R
suppressMessages({
  library(ConsensusClusterPlus)
})

OUT <- "s211_out"
dir.create(OUT, showWarnings = FALSE)
CN_DIR <- "s199_cnv"
CLIN <- file.path(Sys.getenv("LPS_DATA_ROOT", "."), "1", "data", "raw", "TCGA_clinical_real.csv")
SCORE <- file.path(Sys.getenv("LPS_DATA_ROOT", "."), "1", "data", "raw", "TCGA_Score.csv")
SEED <- 20260928L
TOPN <- 5000L
KSET <- 2:6

say <- function(...) cat(..., "\n", sep = "")
rule <- function(t) cat("\n", strrep("=", 78), "\n== ", t, "\n", strrep("=", 78), "\n", sep = "")

ari <- function(a, b) {
  tb <- table(a, b); n <- sum(tb)
  sij <- sum(choose(tb, 2))
  sa <- sum(choose(rowSums(tb), 2)); sb <- sum(choose(colSums(tb), 2))
  exp <- sa * sb / choose(n, 2); mx <- (sa + sb) / 2
  (sij - exp) / (mx - exp)
}
pacf <- function(M) { m <- M[lower.tri(M)]; mean(m > 0.1 & m < 0.9) }
norm_chr <- function(x) ifelse(grepl("^chr", x), x, paste0("chr", x))

## ---------------- 0. 读取 + 覆盖体检 ----------------
rule("0. 读取拷贝数矩阵与基因覆盖体检")
files <- list.files(CN_DIR, pattern = "[.]tsv$", full.names = TRUE)
parts <- sub("[.]tsv$", "", basename(files))
say("文件数 = ", length(files))

lst <- vector("list", length(files)); names(lst) <- parts
rows_v <- integer(length(files)); names(rows_v) <- parts
nchr_v <- integer(length(files)); names(nchr_v) <- parts
covers_auto <- logical(length(files)); names(covers_auto) <- parts
covers_x <- logical(length(files)); names(covers_x) <- parts
chrom_map <- NULL
for (i in seq_along(files)) {
  d <- read.delim(files[i], stringsAsFactors = FALSE, check.names = FALSE)
  rows_v[i] <- nrow(d)
  v <- suppressWarnings(as.numeric(d$copy_number))
  ok <- !is.na(v)
  if (any(duplicated(d$gene_id[ok]))) ok <- ok & !duplicated(d$gene_id)
  x <- setNames(v[ok], d$gene_id[ok])
  lst[[i]] <- x
  ch <- unique(d$chromosome[ok])
  nchr_v[i] <- length(ch); covers_auto[i] <- "chr22" %in% ch; covers_x[i] <- "chrX" %in% ch
  if (is.null(chrom_map) || length(x) > length(chrom_map)) {
    chrom_map <- setNames(d$chromosome[ok], d$gene_id[ok])
  }
  say(sprintf("  %-16s 行=%-6d 有效基因=%-6d 染色体=%-3d 到=%-6s %s",
              parts[i], nrow(d), length(x), nchr_v[i],
              if (nchr_v[i] == 0) "-" else ch[length(ch)],
              if (!covers_auto[i]) "**截断（未到 chr22）**" else if (!covers_x[i]) "缺 chrX" else ""))
}

trunc <- parts[!covers_auto]
part_pri <- parts[covers_auto & covers_x]
part_s48 <- parts[covers_auto]
say("\n截断样本 = ", length(trunc), "：", paste(trunc, collapse = ", "))
say("主分析 = ", length(part_pri), " 例（染色体完整）；敏感性 = ", length(part_s48), " 例（再纳入仅缺 chrX 者）")
say("※ 全 50 例口径被弃用：交集被截断样本压到 chr1-chr8，会整段丢失 chr12q 扩增")

mk <- function(ids) {
  g <- Reduce(intersect, lapply(ids, function(p) names(lst[[p]])))
  m <- sapply(ids, function(p) lst[[p]][g]); rownames(m) <- g; m
}
cn <- mk(part_pri); chrom <- chrom_map[rownames(cn)]
cn48 <- mk(part_s48); chrom48 <- chrom_map[rownames(cn48)]
say("\n【主】", length(part_pri), " 例 × 交集 = ", nrow(cn), " 基因；覆盖 ",
    length(unique(norm_chr(chrom))), " 条染色体")
say("【敏感性】", length(part_s48), " 例 × 交集 = ", nrow(cn48), " 基因；覆盖 ",
    length(unique(norm_chr(chrom48))), " 条染色体")

## ---------------- 1. 临床与评分 ----------------
rule("1. 临床标签与评分")
clin <- read.csv(CLIN, fileEncoding = "UTF-8-BOM", check.names = FALSE)
sc <- read.csv(SCORE, stringsAsFactors = FALSE)
say("临床 n = ", nrow(clin), " | 评分 n = ", nrow(sc))
part_of <- function(bc) paste(strsplit(bc, "-")[[1]][1:3], collapse = "-")
clin$participant <- vapply(clin$barcode, part_of, "")
sc$participant <- vapply(sc$Sample, part_of, "")
pub_map <- setNames(as.numeric(clin$paper_CN.cluster), clin$participant)
sc_map <- setNames(sc$Score, sc$participant)
say("已发表 CN cluster 分布（全 50 例）：",
    paste(names(table(clin$paper_CN.cluster)), as.integer(table(clin$paper_CN.cluster)),
          sep = "=", collapse = "  "))
say("被剔除样本的标签：", paste(sprintf("%s=%s", c(trunc, setdiff(parts, c(part_pri, trunc))),
                                     pub_map[c(trunc, setdiff(parts, c(part_pri, trunc)))]), collapse = "  "))
say("主分析 n = ", length(part_pri), "：", paste(names(table(pub_map[part_pri])),
                                               as.integer(table(pub_map[part_pri])), sep = "=", collapse = "  "))
say("※ SARC 原始 CN cluster 记号为 {2, 3, 6}，Pearson r 按其数值编码（沿用现稿口径）")

## ---------------- 2. 复现闸门 ----------------
rule("2. 复现闸门（已发表标签 vs 评分，全 50 例）")
g_pub <- pub_map[names(sc_map)]; g_scv <- sc_map[names(sc_map)]
r_pub <- cor(g_scv, g_pub); p_pub <- cor.test(g_scv, g_pub)$p.value
rho_pub <- cor(g_scv, g_pub, method = "spearman")
say(sprintf("  r(score, published cluster) = %.4f  P = %.3g  (现稿报 0.535 / 6.2e-05)", r_pub, p_pub))
say(sprintf("  Spearman rho = %.4f  (现稿报 0.453)", rho_pub))
ok_gate <- abs(r_pub - 0.535) < 0.001 && abs(rho_pub - 0.453) < 0.001
say("  复现闸门：", ifelse(ok_gate, "OK", "FAIL"))
if (!ok_gate) stop("复现闸门未通过：数据对齐有问题，先查再跑。")

## ---------------- 3. 特征构造 ----------------
rule("3. 特征构造")
base <- apply(cn, 2, median, na.rm = TRUE)
say("样本内中位 CN 范围 = ", sprintf("%.1f - %.1f", min(base), max(base)), "（反映全基因组倍性）")
L <- log2(sweep(cn, 2, base, "/"))
nf <- !is.finite(L)
if (any(nf)) { keep <- rowSums(nf) == 0; say("非有限值（CN = 0 等）剔除基因 = ", sum(!keep)); L <- L[keep, , drop = FALSE]; chrom <- chrom[keep] }
L <- L - rowMeans(L)
v <- apply(L, 1, var); o <- order(v, decreasing = TRUE); sel <- o[seq_len(min(TOPN, length(o)))]
Xgene <- L[sel, , drop = FALSE]
say("特征 A（基因级）：top-", nrow(Xgene), " 高变异基因（方差下界 = ", sprintf("%.4f", v[sel[length(sel)]]), "）")
chrk <- norm_chr(chrom)
ord_chr <- c(paste0("chr", 1:22), "chrX", "chrY")
chrk <- factor(chrk, levels = ord_chr[ord_chr %in% unique(chrk)])
Xchr <- do.call(rbind, lapply(levels(chrk), function(a) colMeans(L[chrk == a, , drop = FALSE])))
rownames(Xchr) <- levels(chrk)
say("特征 B（染色体级）：", nrow(Xchr), " 条染色体的基因 CN 中位（", paste(levels(chrk), collapse = ","), "）")

## ---------------- 4. 全 k 矩阵 ----------------
rule("4. 共识聚类（两种特征口径 × k = 2..6）")
evalk <- function(res, k, cn) {
  cl <- res[[k]]$consensusClass
  sn <- names(cl); if (is.null(sn)) sn <- colnames(cn)
  y <- sc_map[sn]; p <- pub_map[sn]
  tv <- tapply(y, cl, mean, na.rm = TRUE); o2 <- order(tv)
  lk <- setNames(seq_along(o2), names(tv)[o2])
  clo <- setNames(as.integer(lk[as.character(cl)]), sn)
  tab <- table(derived = clo, published = p)
  tab <- tab[rowSums(tab) > 0, colSums(tab) > 0, drop = FALSE]
  ch <- tryCatch(chisq.test(tab, simulate.p.value = TRUE, B = 10000), error = function(e) NULL)
  cramv <- if (is.null(ch)) NA_real_ else sqrt(as.numeric(ch$statistic) / (sum(tab) * (min(dim(tab)) - 1)))
  chp <- if (is.null(ch)) NA_real_ else ch$p.value
  kw <- kruskal.test(y ~ factor(clo))
  list(k = k, cluster = clo, samples = sn, tab = tab, sizes = as.integer(table(clo)),
       r = cor(y, clo), p = cor.test(y, clo)$p.value, KW_H = as.numeric(kw$statistic), KW_P = kw$p.value,
       cramv = cramv, chi_p = chp, ari = ari(clo, p))
}
run_set <- function(X, cn, tag) {
  say("\n-- ", tag, "：", nrow(X), " 特征 x ", ncol(X), " 样本")
  res <- ConsensusClusterPlus(as.matrix(X), maxK = max(KSET), reps = 1000, pItem = 0.8, pFeature = 0.8,
                              clusterAlg = "hc", distance = "euclidean",
                              innerLinkage = "ward.D2", finalLinkage = "ward.D2",
                              seed = SEED, plot = NULL, verbose = FALSE)
  pacv <- sapply(KSET, function(k) pacf(res[[k]]$consensusMatrix)); names(pacv) <- paste0("k", KSET)
  say("   PAC：", paste(sprintf("%s=%.4f", names(pacv), pacv), collapse = "  "))
  ev <- lapply(KSET, function(k) evalk(res, k, cn)); names(ev) <- paste0("k", KSET)
  for (k in KSET) {
    e <- ev[[paste0("k", k)]]
    say(sprintf("   k=%d | 规模 %s | r(score)=%.4f (P=%.3g) | KW P=%.3g | Cramér's V=%.4f (P=%.3g) | ARI=%.4f",
                k, paste(e$sizes, collapse = "/"), e$r, e$p, e$KW_P, e$cramv, e$chi_p, e$ari))
  }
  list(res = res, pac = pacv, ev = ev, tag = tag, nX = nrow(X), cn = cn, samples = colnames(X))
}

Rg <- run_set(Xgene, cn, "特征 A：基因级 top-5000")
Rc <- run_set(Xchr, cn, "特征 B：染色体级 23 特征")

## 敏感性：48 例（基因级）
base48 <- apply(cn48, 2, median, na.rm = TRUE)
L48 <- log2(sweep(cn48, 2, base48, "/")); nf48 <- !is.finite(L48)
if (any(nf48)) L48 <- L48[rowSums(nf48) == 0, , drop = FALSE]
L48 <- L48 - rowMeans(L48)
v48 <- apply(L48, 1, var); o48 <- order(v48, decreasing = TRUE)
Xg48 <- L48[o48[seq_len(min(TOPN, length(o48)))], , drop = FALSE]
Rs <- run_set(Xg48, cn48, "敏感性：48 例 × 基因级 top-5000")

## ---------------- 5. 落盘 ----------------
rule("5. 产物")
row_of <- function(RS, setlab, ns) {
  do.call(rbind, lapply(KSET, function(k) {
    e <- RS$ev[[paste0("k", k)]]
    data.frame(set = setlab, k = k, n_samples = ns, n_features = RS$nX,
               sizes = paste(e$sizes, collapse = "/"),
               r_score = e$r, p_score = e$p, KW_P = e$KW_P,
               cramersV_vs_pub = e$cramv, chi_p = e$chi_p, ARI_vs_pub = e$ari,
               pac = RS$pac[[paste0("k", k)]], stringsAsFactors = FALSE)
  }))
}
summ <- rbind(row_of(Rg, "gene_top5000_primary46", length(part_pri)),
              row_of(Rc, "chromosome23_primary46", length(part_pri)),
              row_of(Rs, "gene_top5000_sens48", length(part_s48)))
write.csv(summ, file.path(OUT, "D3_cnv_summary.csv"), row.names = FALSE)

mkout <- function(RS) data.frame(participant = RS$samples, score = sc_map[RS$samples],
                                 published_cluster = pub_map[RS$samples],
                                 derived_gene_k3 = RS$ev[["k3"]]$cluster,
                                 derived_gene_kPACmin = RS$ev[[paste0("k", as.integer(sub("k", "", names(which.min(RS$pac)))) )]]$cluster,
                                 stringsAsFactors = FALSE)
write.csv(mkout(Rg), file.path(OUT, "D3_cnv_derived_clusters_gene.csv"), row.names = FALSE)
o4 <- data.frame(participant = Rc$samples, score = sc_map[Rc$samples],
                 published_cluster = pub_map[Rc$samples],
                 derived_chrom_k2 = Rc$ev[["k2"]]$cluster, derived_chrom_k3 = Rc$ev[["k3"]]$cluster,
                 derived_chrom_k4 = Rc$ev[["k4"]]$cluster, stringsAsFactors = FALSE)
write.csv(o4, file.path(OUT, "D3_cnv_derived_clusters_chromosome.csv"), row.names = FALSE)
write.csv(as.data.frame.matrix(Rg$ev[["k3"]]$tab), file.path(OUT, "D3_cnv_contingency_gene_k3.csv"))
write.csv(as.data.frame.matrix(Rc$ev[["k3"]]$tab), file.path(OUT, "D3_cnv_contingency_chrom_k3.csv"))

cov <- data.frame(participant = parts, rows_raw = as.integer(rows_v),
                  n_genes_effective = vapply(lst, length, 0L), n_chromosomes = as.integer(nchr_v),
                  covers_autosomes = covers_auto, covers_chrX = covers_x,
                  in_primary46 = parts %in% part_pri, in_sens48 = parts %in% part_s48,
                  stringsAsFactors = FALSE)
write.csv(cov, file.path(OUT, "D3_cnv_gene_coverage.csv"), row.names = FALSE)
saveRDS(list(gene = Rg, chrom = Rc, sens48 = Rs), file.path(OUT, "D3_cnv_ccp.rds"))

say(" 已写：", paste(list.files(OUT), collapse = ", "))
say("\n完成。")
