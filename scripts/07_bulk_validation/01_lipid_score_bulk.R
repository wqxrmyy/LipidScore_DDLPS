# ============================================================
# scripts/07_bulk_validation/01_lipid_score_bulk.R
# 目的：GSE221492 bulk RNA-seq 六基因脂质代谢评分 + 配对比较
#
# 本脚本是 **两套并存实现之一**（另一套：scripts/07_bulk_validation/01_lipid_score_gse221492.py）。
# 两者必须给出**逐位一致**的结果；由 scripts/09_qc/assert_manuscript_numbers.py 做交叉校验。
# 本文报告的每个数字都由这两套实现共同复现（容差见 QC 脚本）。
#
# ★ 2026-09-24 修正两处口径错误（旧版会把本文报告的第 3 位小数算错）：
#   [FIX-1] library size 必须在**合并重复 hgnc_symbol 之前**计算。
#           合并丢弃 12,873 行（28.3%）的计数 -> 分母偏小 0.47% -> CPM 系统性偏高
#           -> 评分 Δ 由 −0.587 变为 −0.586（报告值被改错）。
#   [FIX-2] DDLPS_sWD 不得并入 DDLPS_WD。sWD 是独立标签；把它当 WD 会把配对集
#           从「严格 8 对」放大到 10 对，得到本文从未报告过的 n=10 结果。
#           严格 8 对为**主分析**，含 sWD 的 10 对降为敏感性分析。
# ============================================================

# ---------- 0. 项目根（可用环境变量 LPS_ROOT 覆盖；默认由脚本位置推断） ----------
Root <- Sys.getenv("LPS_ROOT", unset = NA)
if (is.na(Root)) {
  args <- commandArgs(trailingOnly = FALSE)
  fa <- grep("^--file=", args, value = TRUE)
  if (length(fa) == 1) {
    Root <- normalizePath(file.path(dirname(sub("^--file=", "", fa)), "..", ".."))
  } else {
    stop("无法推断项目根目录：请设置环境变量 LPS_ROOT")
  }
}
cat("项目根:", Root, "\n\n")

suppressPackageStartupMessages({
  library(data.table); library(ggplot2); library(dplyr); library(tidyr)
})

# ---------- 1. 评分基因（读口径锁，不硬编码） ----------
score_lock <- fread(file.path(Root, "governance/metadata_lock/score_genes_v1.tsv"),
                    data.table = FALSE)
stopifnot("include_in_scoring" %in% names(score_lock),
          is.logical(score_lock$include_in_scoring))
score_genes <- score_lock$gene_symbol[score_lock$include_in_scoring]
cat("评分基因:", paste(score_genes, collapse = ", "), "\n\n")

# ---------- 2. 读取 counts ----------
raw_file <- file.path(Root, "data/GSE221492_bulk/raw",
                      "GSE221492_bulk_LPS_atlas_raw_counts_matrix.tsv.gz")
counts <- fread(raw_file, data.table = FALSE)
counts_mat <- as.matrix(counts[, -c(1, 2)])
rownames(counts_mat) <- counts$hgnc_symbol

# ---- [FIX-1] library size 必须在合并重复符号之前算 ----
# 合并重复行会丢弃 12,873 行的计数（占全部 45,473 行的 28.3%），
# 使列和不再等于该样本的总计数。library size 应始终取自**完整矩阵**。
lib_size <- colSums(counts_mat, na.rm = TRUE)
stopifnot(all(lib_size > 0))
cat(sprintf("完整矩阵 library size: 首个样本 %.0f（合并后为 %.0f，差 %.2f%%）\n\n",
            lib_size[1], sum(counts_mat[!duplicated(rownames(counts_mat)), 1], na.rm = TRUE),
            (sum(counts_mat[!duplicated(rownames(counts_mat)), 1], na.rm = TRUE) / lib_size[1] - 1) * 100))

# ---------- 3. 数据清洗 ----------
# 3.1 重复基因名：保留总表达最高的一行（仅用于**基因取值**，不影响上面的 library size）
n_before <- nrow(counts_mat)
if (any(duplicated(rownames(counts_mat)))) {
  row_sums <- rowSums(counts_mat, na.rm = TRUE)
  ord <- order(rownames(counts_mat), -row_sums)
  counts_mat <- counts_mat[ord, ]
  counts_mat <- counts_mat[!duplicated(rownames(counts_mat)), ]
  cat(sprintf("✓ 重复基因名已处理: %d -> %d 行（丢弃 %d 行）\n",
              n_before, nrow(counts_mat), n_before - nrow(counts_mat)))
}

# 3.2 NA 检查与处理
na_before <- sum(is.na(counts_mat))
if (na_before > 0) {
  na_samples <- colnames(counts_mat)[colSums(is.na(counts_mat)) > 0]
  na_genes_in_score <- intersect(score_genes,
                                 rownames(counts_mat)[rowSums(is.na(counts_mat)) > 0])
  cat("检测到 NA:", na_before, "个\n")
  cat("涉及样本:", length(na_samples), "个\n")
  cat("评分基因中是否含 NA:", if (length(na_genes_in_score) > 0)
    paste(na_genes_in_score, collapse = ",") else "无", "\n")
  if (length(na_genes_in_score) > 0) {
    stop("评分基因含 NA，需要进一步处理，不能直接替换为 0")
  }
  counts_mat[is.na(counts_mat)] <- 0
  cat("✓ NA 已替换为 0（评分基因无 NA，安全）\n\n")
}
cat("矩阵:", nrow(counts_mat), "x", ncol(counts_mat), "\n\n")

# ---------- 4. 计算评分（与 Python 实现同式） ----------
cpm <- sweep(counts_mat, 2, lib_size, "/") * 1e6          # [FIX-1] 用完整矩阵的 library size
log_cpm <- log2(cpm + 1)
z_mat <- t(scale(t(log_cpm[score_genes, , drop = FALSE])))
score <- colMeans(z_mat, na.rm = FALSE)

cat("评分为 NaN 的样本数:", sum(is.nan(score)), "\n")
cat("评分为 NA 的样本数:", sum(is.na(score)), "\n\n")
stopifnot(!anyNA(score))

# ---------- 5. 标签解析（sWD 独立，不得并入 WD） ----------
classify <- function(x) {
  if (grepl("_DDLPS_DD$", x))              return("DDLPS_DD")
  if (grepl("_DDLPS_WD$", x))              return("DDLPS_WD")     # [FIX-2] 只认严格 WD
  if (grepl("_DDLPS_sWD$", x))             return("DDLPS_sWD")    # [FIX-2] sWD 独立
  if (grepl("_WDLPS$", x))                 return("WDLPS")
  return("Control")
}
score_df <- data.frame(
  sample  = names(score),
  patient = sub("^(P[0-9]+)_.*", "\\1", names(score)),
  group   = sapply(names(score), classify),
  score   = as.numeric(score),
  stringsAsFactors = FALSE
)
print(table(score_df$group))

# ---------- 6. 配对分析 ----------
# 主分析：严格配对 = 同一患者同时拥有真正的 DDLPS_WD 与 DDLPS_DD 样本
wide_all <- score_df %>%
  filter(group %in% c("DDLPS_DD", "DDLPS_WD")) %>%
  select(patient, group, score) %>%
  pivot_wider(names_from = group, values_from = score) %>%
  filter(!is.na(DDLPS_WD), !is.na(DDLPS_DD),
         !is.nan(DDLPS_WD), !is.nan(DDLPS_DD)) %>%
  mutate(delta = DDLPS_DD - DDLPS_WD) %>%
  arrange(patient)

# 敏感性：把 sWD 也当作 WD（旧版误当作主分析）。
# 同一患者若同时有真正的 WD 与 sWD 样本，优先取严格 WD（rank=1）。
dd_part <- score_df %>% filter(group == "DDLPS_DD") %>%
  transmute(patient, DDLPS_DD = score)
wd_part <- score_df %>% filter(group %in% c("DDLPS_WD", "DDLPS_sWD")) %>%
  mutate(rank = ifelse(group == "DDLPS_WD", 1L, 2L)) %>%
  arrange(patient, rank) %>%
  group_by(patient) %>%
  summarise(WD_any = first(score), .groups = "drop")
sens <- inner_join(dd_part, wd_part, by = "patient") %>%
  mutate(delta = DDLPS_DD - WD_any) %>%
  arrange(patient)

cat("\n=== 主分析（严格配对，sWD 除外）===\n")
cat("有效配对数:", nrow(wide_all), " 患者:", paste(wide_all$patient, collapse = ", "), "\n")
print(as.data.frame(wide_all))
cat("\n=== 敏感性分析（WD 含 sWD）===\n")
cat("有效配对数:", nrow(sens), "\n")

run_test <- function(w, label, wd_col = "DDLPS_WD") {
  wt <- wilcox.test(w$DDLPS_DD, w[[wd_col]], paired = TRUE)
  res <- list(label = label, n_pairs = nrow(w),
              V = unname(wt$statistic), p = wt$p.value,
              mean_delta = mean(w$delta), median_delta = median(w$delta),
              n_dd_lt_wd = sum(w$delta < 0))
  cat(sprintf("\n%s: n=%d  V=%.1f  P=%.6f\n", label, res$n_pairs, res$V, res$p))
  cat(sprintf("  mean Δ   = %.6f  <- 本文报告值（正文写 Δ = -0.587）\n", res$mean_delta))
  cat(sprintf("  median Δ = %.6f  <- 不要把中位数当成本文报告的 Δ\n", res$median_delta))
  cat(sprintf("  DD<WD    = %d/%d\n", res$n_dd_lt_wd, res$n_pairs))
  res
}
main_res <- run_test(wide_all, "主分析")
sens_res <- run_test(sens, "敏感性（WD 含 sWD）", wd_col = "WD_any")

# ---------- 7. 主分析断言：钉死本文报告的数字 ----------
# 若数据/口径变动导致任一报告位改变，必须失败并人工复核，不得静默通过。
stopifnot(
  main_res$n_pairs      == 8,
  main_res$n_dd_lt_wd   == 6,
  abs(main_res$p        - 0.078125) < 1e-9,
  round(main_res$mean_delta, 3) == -0.587
)
cat("\n✓ 主分析断言通过：n=8 / 6-of-8 / P=0.078125 / mean Δ=-0.587\n")

# ---------- 8. 可视化（仅主分析） ----------
score_paired_clean <- score_df %>%
  filter(group %in% c("DDLPS_DD", "DDLPS_WD"), patient %in% wide_all$patient)
p1 <- ggplot(score_paired_clean, aes(group, score, fill = group)) +
  geom_boxplot(alpha = 0.6, width = 0.5) +
  geom_point(size = 2) +
  geom_line(aes(group = patient), color = "grey50", alpha = 0.6) +
  scale_fill_manual(values = c("DDLPS_WD" = "#4fb04c", "DDLPS_DD" = "#d8311f")) +
  labs(title = "Six-gene lipid metabolism score",
       subtitle = paste0("Paired DDLPS (n=", nrow(wide_all), ")"),
       x = "", y = "Score") +
  theme_classic(base_size = 12) + theme(legend.position = "none")

dir.create(file.path(Root, "results/figures"), recursive = TRUE, showWarnings = FALSE)
ggsave(file.path(Root, "results/figures/bulk_lipid_score_paired.pdf"),
       p1, width = 4, height = 5)
ggsave(file.path(Root, "results/figures/bulk_lipid_score_paired.png"),
       p1, width = 4, height = 5, dpi = 300)

# ---------- 9. 保存结果 ----------
dir.create(file.path(Root, "results/tables"), recursive = TRUE, showWarnings = FALSE)
write.csv(score_df, file.path(Root, "results/tables/bulk_lipid_score_all_samples.csv"),
          row.names = FALSE)
write.csv(as.data.frame(wide_all),
          file.path(Root, "results/tables/bulk_lipid_score_paired_delta.csv"),
          row.names = FALSE)
write.csv(data.frame(
  analysis = c("Main (strict pairs)", "Sensitivity (WD includes sWD)"),
  n_pairs = c(main_res$n_pairs, sens_res$n_pairs),
  V_statistic = c(main_res$V, sens_res$V),
  p_value = c(main_res$p, sens_res$p),
  mean_delta = c(main_res$mean_delta, sens_res$mean_delta),
  median_delta = c(main_res$median_delta, sens_res$median_delta),
  n_dd_lt_wd = c(main_res$n_dd_lt_wd, sens_res$n_dd_lt_wd),
  reported_in_manuscript = c("YES (n=8, Delta=-0.587, P=0.078, 6/8)", "no"),
  stringsAsFactors = FALSE),
  file.path(Root, "results/tables/bulk_lipid_score_summary.csv"), row.names = FALSE)

# ---------- 10. 记录限制（追加，不覆盖历史） ----------
lim <- c(
  paste0("issue_id: LIM-", format(Sys.time(), "%Y%m%d%H%M")),
  "analysis_module: 01_lipid_score_bulk (R)",
  paste0("issue: 原始矩阵含 ", na_before, " 个 NA，分布于 ", length(na_samples), " 个样本"),
  paste0("na_samples: ", paste(na_samples, collapse = ", ")),
  "root_cause: 622 个基因在 9 个样本中缺失，呈批次特异性模式，疑似不同建库批次基因注释差异",
  "resolution: NA 替换为 0（评分基因无 NA，不影响评分计算）",
  paste0("dup_gene_rows_dropped: ", n_before - nrow(counts_mat)),
  "library_size_source: full_matrix_before_symbol_collapse  # [FIX-1]",
  "sWD_handling: kept_as_separate_label_not_merged_into_WD  # [FIX-2]",
  paste0("effective_n_pairs_main: ", main_res$n_pairs),
  paste0("effective_n_pairs_sensitivity: ", sens_res$n_pairs),
  "status: active"
)
writeLines(lim, file.path(Root, "governance/limitations/LIM_lipid_score_bulk.yaml"))

cat("\n✓ 完成\n")
