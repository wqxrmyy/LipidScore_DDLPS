## S34 6E 患者级重算：把「细胞当独立样本」的伪重复改为患者级
## 输入：s34_6e_slingshot.rds（V1/V2 的 pt + cd）＋ cell_meta_tumor.csv（patient 映射）
## 输出：s34_6e_patientlevel.txt + fig6_data/slingshot_patientlevel.csv
suppressPackageStartupMessages({ library(Matrix); library(Seurat) })
options(stringsAsFactors = FALSE)

`%+%` <- function(a, b) paste0(a, b)

W   <- file.path(Sys.getenv("LPS_WORK_ROOT", "."), "_work")
FIG <- file.path(W, "p5_figures", "fig6_data")
REP <- file.path(W, "s34_6e_patientlevel.txt")
con <- file(REP, "w", encoding = "UTF-8")
A <- function(...) { s <- paste0(...); cat(s, "\n", sep = ""); cat(s, "\n", sep = "", file = con) }

A("===============================================================================")
A("S34 6E patient-level re-analysis (fix pseudoreplication)   ", format(Sys.time()))
A("===============================================================================")

rds <- file.path(W, "s34_6e_slingshot.rds")
if (!file.exists(rds)) { A("ERRS: rds not found -> ", rds); close(con); quit(status = 1) }
obj <- readRDS(rds)
A("loaded: ", paste(names(obj), collapse = ", "))

meta <- read.csv(file.path(FIG, "cell_meta_tumor.csv"), fileEncoding = "UTF-8-BOM",
                 colClasses = c(barcode = "character", patient = "character"))
pt_map <- setNames(meta$patient, meta$barcode)
A("meta: ", nrow(meta), " cells ; patients = ", length(unique(meta$patient)))

PAIRED <- c("P1", "P4", "P7", "P10", "P14", "P15", "P16")
rows <- list()

for (nm in names(obj)) {
  o <- obj[[nm]]; if (is.null(o)) next
  pt <- o$pt; cd <- o$cd
  pat <- pt_map[rownames(cd)]
  score <- as.numeric(cd$score); tissue <- as.character(cd$tissue)
  A("")
  A("---- ", nm, " : ", nrow(cd), " cells / ", ncol(pt), " lineages ----")
  A("   patient mapped: ", sum(!is.na(pat)), " / ", length(pat), " (unmapped dropped)")
  keep <- !is.na(pat) & !is.na(score)
  for (l in seq_len(ncol(pt))) {
    x <- pt[, l]; k <- keep & !is.na(x)
    if (sum(k) < 50) { A("   L", l, " : too few cells, skip"); next }
    ## --- cell-level (pseudoreplicated, for contrast only) ---
    ct <- suppressWarnings(cor.test(x[k], score[k], method = "spearman"))
    ## --- patient-level (correct) ---
    dp <- data.frame(pat = pat[k], x = x[k], y = score[k])
    ag <- aggregate(cbind(x, y) ~ pat, data = dp, FUN = mean)
    ag <- ag[order(ag$x), ]
    rp <- suppressWarnings(cor.test(ag$x, ag$y, method = "spearman"))
    npat <- nrow(ag)
    z95 <- qnorm(0.975); z80 <- qnorm(0.80)
    mde <- tanh((z95 + z80) / sqrt(max(npat - 3, 1)))
    qq <- quantile(ag$x, c(0.2, 0.8)); lo <- ag$y[ag$x <= qq[1]]; hi <- ag$y[ag$x >= qq[2]]
    A(sprintf("   L%d | cell-level rho=%+.4f (P=%.2e, n=%d) [PSEUDOREPLICATED]",
              l, ct$estimate, ct$p.value, sum(k)))
    A(sprintf("        patient-level rho=%+.4f (P=%.4f, n=%d patients) | MDE|rho|(80%% power)=%.3f | d_top20_minus_bot20=%+.4f",
              rp$estimate, rp$p.value, npat, mde, mean(hi) - mean(lo)))
    rows[[length(rows) + 1]] <- data.frame(ver = nm, lineage = l,
        cell_rho = round(ct$estimate, 4), cell_P = signif(ct$p.value, 3), cell_n = sum(k),
        pat_rho = round(rp$estimate, 4), pat_P = round(rp$p.value, 4), pat_n = npat,
        MDE = round(mde, 3), d_delta = round(mean(hi) - mean(lo), 4))
    ## --- tissue-level mean (this lineage) ---
    tis <- tissue[k]
    dd_m <- mean(x[k][tis == "DDLPS-DD"], na.rm = TRUE)
    wd_m <- mean(x[k][tis == "DDLPS-WD"], na.rm = TRUE)
    A(sprintf("        tissue mean pt: DD=%.3f WD=%.3f diff(DD-WD)=%+.3f", dd_m, wd_m, dd_m - wd_m))
    ## --- paired 7-pair WD vs DD at patient level ---
    dpt <- data.frame(pat = pat[k], tis = tis, x = x[k])
    dpt <- dpt[dpt$tis %in% c("DDLPS-DD", "DDLPS-WD"), ]
    if (nrow(dpt) > 0) {
      pm <- aggregate(x ~ pat + tis, data = dpt, FUN = mean)
      wide <- reshape(pm, idvar = "pat", timevar = "tis", direction = "wide")
      cDD <- "x.DDLPS-DD"; cWD <- "x.DDLPS-WD"
      if (all(c(cDD, cWD) %in% names(wide))) {
        sub <- wide[!is.na(wide[[cDD]]) & !is.na(wide[[cWD]]), ]
        sub <- sub[sub$pat %in% PAIRED, ]
        if (nrow(sub) >= 4) {
          wt <- suppressWarnings(wilcox.test(sub[[cDD]], sub[[cWD]], paired = TRUE))
          A(sprintf("        paired 7-pair WD vs DD: n=%d  mean diff(DD-WD)=%+.3f  Wilcoxon P=%.4f",
                    nrow(sub), mean(sub[[cDD]] - sub[[cWD]]), wt$p.value))
          A("           pairs: ", paste(sub$pat, collapse = ","))
        }
      }
    }
  }
}

if (length(rows)) {
  out <- do.call(rbind, rows)
  write.csv(out, file.path(FIG, "slingshot_patientlevel.csv"), row.names = FALSE, fileEncoding = "UTF-8")
  A(""); A("wrote: fig6_data/slingshot_patientlevel.csv (", nrow(out), " rows)")
}
A(""); A("DONE")
close(con)
cat("EXIT_OK\n")
