suppressPackageStartupMessages({
  library(Seurat); library(harmony); library(Matrix); library(future)
})
options(future.globals.maxSize = 16 * 1024^3)
W <- file.path(Sys.getenv("LPS_WORK_ROOT", "."), "_work")
sink(file.path(W, "s39_s4_precheck.txt"), split = TRUE)
say <- function(...) cat(sprintf(...), "\n")

say("== S4 precheck ==")
say("R      %s", paste(R.version$major, R.version$minor, sep = "."))
say("Seurat  %s", as.character(packageVersion("Seurat")))
say("harmony %s", as.character(packageVersion("harmony")))
say("RunHarmony formals: %s", paste(names(formals(harmony::RunHarmony)), collapse = ", "))
for (p in c("FNN", "RANN", "BiocNeighbors", "uwot", "AUCell")) {
  say("pkg %-15s %s", p, requireNamespace(p, quietly = TRUE))
}

so <- readRDS(file.path(W, "s39_seurat.s3.rds"))
say("\ndim: %d genes x %d cells", nrow(so), ncol(so))
say("assays:    %s", paste(names(so@assays), collapse = ", "))
say("reductions: %s", paste(names(so@reductions), collapse = ", "))
for (r in names(so@reductions)) say("   %-14s ncol=%d", r, ncol(Embeddings(so, r)))
md <- colnames(so@meta.data)
say("meta (%d cols): %s", length(md), paste(md, collapse = ", "))

for (v in c("patient", "tissue", "L1", "L2", "L3", "grade", "sample")) {
  if (v %in% md) {
    say("\n-- %s --", v)
    print(table(so@meta.data[[v]], useNA = "ifany"))
  }
}

say("\n-- paired WD/DD cell counts (majority-class baseline) --")
PAIRS <- c("P1", "P4", "P7", "P10", "P14", "P15", "P16")
for (p in PAIRS) {
  sel <- which(so$patient == p & so$tissue %in% c("DDLPS-WD", "DDLPS-DD"))
  tb <- table(droplevels(factor(so$tissue[sel])))
  say("%-5s n=%6d  WD=%6d DD=%6d  majority_base=%.3f",
      p, length(sel), sum(tb["DDLPS-WD"]), sum(tb["DDLPS-DD"]),
      max(prop.table(tb)))
}

say("\n-- tissue x patient (top tissues) --")
tb <- table(so$tissue, so$patient)
print(tb)
say("\nPRE_DONE")
sink()
