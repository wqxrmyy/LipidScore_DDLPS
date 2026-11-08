# s148 · 概率杠杆实测 —— TCGA-SARC DDLPS 多终点可用性探查
# 目的：手稿现仅用 OS(n=50, ev=23)。本脚本核查同一批 50 例是否还有 DSS / RFS 可用，
#       若有，则"换/加终点"为零新数据成本的杠杆，需实测其效应量后再判断值不值。
# 口径：与 S4 一致 —— Harrell C 用风险约定（高分=高风险），survival::concordance(Surv ~ score)
# 幂等：全部输出为覆盖写。
# 日期：2026-09-25

suppressMessages(library(survival))

base <- file.path(Sys.getenv("LPS_DATA_ROOT", "."), "1", "data")
outdir <- file.path(Sys.getenv("LPS_WORK_ROOT", "."), "_work", "s148_out")
dir.create(outdir, showWarnings = FALSE, recursive = TRUE)

clin <- read.csv(file.path(base, "raw/TCGA_clinical_real.csv"),
                 fileEncoding = "UTF-8-BOM", check.names = FALSE)
sc <- read.csv(file.path(base, "raw/TCGA_Score.csv"),
               fileEncoding = "UTF-8-BOM", check.names = FALSE)

cat("== 载入 ==\n")
cat("clin rows =", nrow(clin), " score rows =", nrow(sc), "\n\n")

cat("== 候选终点列可用性（非空计数）==\n")
cand <- c("barcode", "OS_time", "OS_status", "days_to_death", "days_to_last_follow_up",
          "paper_OS.days", "paper_OS.status", "paper_DSS.status",
          "paper_RFS.status", "paper_RFS.days",
          "paper_new.tumor.event.dx.indicator", "paper_days.to.new.tumor.event",
          "paper_local.recurrence", "paper_local.recurrence.days",
          "paper_distant.recurrence", "paper_distant.recurrence.days",
          "paper_tumor.status", "paper_treatment.outcome.at.tcga.followup")
for (c in cand) {
  if (!c %in% names(clin)) { cat(sprintf("%-44s MISSING COL\n", c)); next }
  v <- clin[[c]]
  ne <- sum(!is.na(v) & as.character(v) != "")
  u <- paste(head(unique(as.character(v[!is.na(v) & as.character(v) != ""])), 5), collapse = "|")
  cat(sprintf("%-44s nonNA=%3d  nuniq=%3d  eg=[%s]\n", c, ne, length(unique(v)), u))
}

cat("\n== 评分与临床的合并键检查 ==\n")
cat("score Sample eg  :", head(sc$Sample, 2), "\n")
cat("clin barcode eg  :", head(clin$barcode, 2), "\n")
key_s <- substr(sc$Sample, 1, 12)
key_c <- substr(clin$barcode, 1, 12)
cat("score 唯一患者 =", length(unique(key_s)), " clin 唯一患者 =", length(unique(key_c)), "\n")
cat("交集患者 =", length(intersect(key_s, key_c)), "\n")
cat("clin 独有 =", length(setdiff(key_c, key_s)), "\n")

cat("\n== OS 基线核对（应与手稿 ev=23, n=50 一致）==\n")
clin$key <- key_c
sc$key <- key_s
sc_agg <- aggregate(Score ~ key, data = sc, FUN = mean)
m <- merge(clin, sc_agg, by = "key", all = FALSE)
cat("merged n =", nrow(m), "\n")
ok <- !is.na(m$OS_time) & !is.na(m$OS_status)
cat("OS 可用 n =", sum(ok), " events =", sum(m$OS_status[ok] == 1), "\n")

sink(file.path(outdir, "s148_probe.txt"))
cat("clin rows =", nrow(clin), "score rows =", nrow(sc), "\n")
cat("merged n =", nrow(m), " OS n =", sum(ok), " ev =", sum(m$OS_status[ok] == 1), "\n")
sink()
cat("\n[OK] probe written\n")
