# s149 · 概率杠杆实测 —— TCGA-SARC DDLPS 多终点效应量
# 上游：s148（探查确认 RFS n=50 完整、DSS n=47）
# 目的：量化"加终点"这一零新数据成本杠杆的收益。手稿现仅报 OS(n=50, ev=23)。
# 口径：与 S4 严格一致 ——
#   Harrell C : survival::concordance(Surv(time,status) ~ score)，风险约定（高分=高风险，C>0.5 表示与叙事一致）
#   HR        : coxph(Surv(time,status) ~ score_z)，score 按 SD 标准化 → "per-SD HR"
#   全部为单变量、未校正（与手稿 §3.4 口径一致）
# 幂等：全部输出覆盖写。
# 日期：2026-09-25

suppressMessages(library(survival))

base   <- file.path(Sys.getenv("LPS_DATA_ROOT", "."), "1", "data")
outdir <- file.path(Sys.getenv("LPS_WORK_ROOT", "."), "_work", "s148_out")
dir.create(outdir, showWarnings = FALSE, recursive = TRUE)

clin <- read.csv(file.path(base, "raw/TCGA_clinical_real.csv"),
                 fileEncoding = "UTF-8-BOM", check.names = FALSE)
sc   <- read.csv(file.path(base, "raw/TCGA_Score.csv"),
                 fileEncoding = "UTF-8-BOM", check.names = FALSE)

clin$key <- substr(clin$barcode, 1, 12)
sc$key   <- substr(sc$Sample, 1, 12)
sc_agg   <- aggregate(Score ~ key, data = sc, FUN = mean)
m        <- merge(clin, sc_agg, by = "key", all = FALSE)

# ---- 端点构造 -------------------------------------------------------------
# OS : 原生 OS_time / OS_status(0,1)
# DSS: 疾病特异死亡；TCGA 未单列 DSS 时间列，采用 OS 时间 + DSS 状态（标准做法）
# RFS: 原生 paper_RFS.days / paper_RFS.status(Relapse, Non-Relapse)
# NTE: 新肿瘤事件 paper_days.to.new.tumor.event / paper_new.tumor.event.dx.indicator
#      → 无事件者的删失时间需借用随访时间，定义较弱，仅作参考不计入主判读

m$dss_status <- suppressWarnings(as.numeric(as.character(m$paper_DSS.status)))
m$rfs_status <- ifelse(m$paper_RFS.status == "Relapse", 1,
                ifelse(m$paper_RFS.status == "Non-Relapse", 0, NA))
m$nte_status <- ifelse(m$paper_new.tumor.event.dx.indicator == "YES", 1,
                ifelse(m$paper_new.tumor.event.dx.indicator == "NO", 0, NA))
# NTE 删失时间：以 OS 随访时间兜底
m$nte_time   <- ifelse(is.na(m$paper_days.to.new.tumor.event), m$OS_time,
                       m$paper_days.to.new.tumor.event)

cat("== DSS 状态值域 ==\n"); print(table(m$paper_DSS.status, useNA = "ifany"))
cat("\n== RFS 状态值域 ==\n"); print(table(m$paper_RFS.status, useNA = "ifany"))

run_ep <- function(label, time, status) {
  ok <- !is.na(time) & !is.na(status) & time > 0
  n  <- sum(ok); ev <- sum(status[ok] == 1)
  S  <- Surv(time[ok], status[ok])
  s  <- m$Score[ok]
  sz <- (s - mean(s)) / sd(s)

  cx  <- coxph(S ~ sz)
  sm  <- summary(cx)
  # reverse = TRUE → 风险约定（高分 = 高风险）；与 S4 文档 §一 的基准测试口径一致
  cc  <- concordance(S ~ s, reverse = TRUE)
  HR  <- sm$conf.int[1, 1]; lo <- sm$conf.int[1, 3]; hi <- sm$conf.int[1, 4]
  pv  <- sm$coef[1, 5]
  C   <- cc$concordance
  Cse <- sqrt(cc$var)

  # log-rank（中位分割）作为方向旁证
  grp <- ifelse(s > median(s), "High", "Low")
  lr  <- survdiff(S ~ grp)
  lrp <- 1 - pchisq(lr$chisq, df = 1)

  data.frame(
    endpoint = label, n = n, events = ev,
    HR_perSD = round(HR, 3),
    CI_low = round(lo, 3), CI_high = round(hi, 3),
    P_cox = signif(pv, 3),
    C_index = round(C, 3),
    C_CI_low = round(C - 1.96 * Cse, 3), C_CI_high = round(C + 1.96 * Cse, 3),
    dir_consistent = ifelse(C > 0.5, "YES(higher score = worse)", "NO"),
    P_logrank_mediansplit = signif(lrp, 3),
    stringsAsFactors = FALSE
  )
}

res <- rbind(
  run_ep("OS",  m$OS_time,        m$OS_status),
  run_ep("DSS", m$OS_time,        m$dss_status),
  run_ep("RFS", m$paper_RFS.days, m$rfs_status),
  run_ep("NTE", m$nte_time,       m$nte_status)
)

cat("\n== 结果（风险约定：C>0.5 = 高分→预后差，与手稿叙事一致）==\n")
print(res, row.names = FALSE)

write.csv(res, file.path(outdir, "s149_tcga_endpoints.csv"),
          row.names = FALSE, fileEncoding = "UTF-8")

# 中位生存（仅报告，用于叙事）
cat("\n== 中位生存（月，按评分中位分割）==\n")
for (ep in list(list("OS", m$OS_time, m$OS_status),
                list("RFS", m$paper_RFS.days, m$rfs_status))) {
  lab <- ep[[1]]; tt <- ep[[2]]; ss <- ep[[3]]
  ok <- !is.na(tt) & !is.na(ss) & tt > 0
  grp <- ifelse(m$Score[ok] > median(m$Score[ok]), "High", "Low")
  f <- survfit(Surv(tt[ok], ss[ok]) ~ grp)
  print(f)
}

cat("\n[OK] s149 done\n")
