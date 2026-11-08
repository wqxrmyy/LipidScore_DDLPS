suppressPackageStartupMessages({ library(Seurat); library(harmony) })
cat("methods(RunHarmony):\n"); print(methods(RunHarmony))
cat("\n--- formals of each method ---\n")
for (m in methods(RunHarmony)) {
  f <- tryCatch(formals(m), error = function(e) NULL)
  cat("\n[", m, "]\n", sep = "")
  if (!is.null(f)) print(names(f)) else cat("  <no formals>\n")
}
cat("\n--- harmony 版本 & 顶层函数 ---\n")
print(packageVersion("harmony"))
cat("exports:", paste(sort(getNamespaceExports("harmony")), collapse = ", "), "\n")
