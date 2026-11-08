suppressPackageStartupMessages({ library(Seurat) })
W <- file.path(Sys.getenv("LPS_WORK_ROOT", "."), "_work")
LOG <- file.path(W, "s39_s4_fixcol_report.txt")
sink(LOG, split = TRUE)
say <- function(...) cat(sprintf(...), "\n")

say("== 修正 metadata 列名：grade -> subtype_code ==")
say("依据：原始注释文件 10 列实为 barcode/patient/tissue/site/cluster/L1/L2/L3/umap1/umap2，")
say("      无 grade 列；col5 是 43 簇亚型编号（与 L1/L2/L3 交叉纯度均为 1.000000）。")
say("时间: %s", format(Sys.time()))

f <- file.path(W, "s39_seurat.s4.rds")
so <- readRDS(f)
say("\nloaded: %d genes x %d cells", nrow(so), ncol(so))
say("改前 cols: %s", paste(colnames(so@meta.data), collapse = ", "))

stopifnot("grade" %in% colnames(so@meta.data))
old_vals <- as.character(so@meta.data[["grade"]])
say("\n待改列 grade：uniq = %d", length(unique(old_vals)))

so[["subtype_code"]] <- factor(old_vals)
so[["grade"]] <- NULL

say("改后 cols: %s", paste(colnames(so@meta.data), collapse = ", "))

# 回读验证
new_vals <- as.character(so@meta.data[["subtype_code"]])
say("\n校验：")
say("  值完全一致        : %s", identical(sort(unique(old_vals)), sort(unique(new_vals))))
say("  逐位一致          : %s", all(old_vals == new_vals))
say("  subtype_code uniq : %d", length(unique(new_vals)))
say("  grade 已移除      : %s", !("grade" %in% colnames(so@meta.data)))
say("  subtype_code 与 L3 的列联（前 5 组示例）:")
tb <- table(so@meta.data[["subtype_code"]], so@meta.data[["L3"]])
top <- order(-rowSums(tb))[1:5]
for (i in top) {
  r <- tb[i, ]; r <- r[r > 0]
  say("    %-10s -> %s", rownames(tb)[i],
      paste(sprintf("%s(%d)", names(r), as.integer(r)), collapse = "; "))
}

saveRDS(so, f, compress = "gzip")
say("\nresaved: %.1f MB", file.size(f) / 1048576)

# 回读磁盘确认
chk <- readRDS(f)
say("回读磁盘：%d cells | cols = %s", ncol(chk),
    paste(colnames(chk@meta.data), collapse = ", "))
say("  grade 不再存在: %s", !("grade" %in% colnames(chk@meta.data)))
say("  subtype_code 存在: %s", "subtype_code" %in% colnames(chk@meta.data))
say("\nDONE_FIXCOL")
sink()
