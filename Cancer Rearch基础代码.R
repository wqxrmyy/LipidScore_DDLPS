# ==================== 设置工作路径 ====================
setwd("C:/Users/Tao25/Documents/Liposarcoma/2")

# ==================== 第1步：创建缺失的必要文件 ====================

# ---- 1.1 创建 .gitignore ----
writeLines(
  c(
    "# R临时文件",
    ".Rproj.user/",
    ".Rhistory",
    ".RData",
    ".Ruserdata",
    "*.Rproj",
    "",
    "# 大型数据文件（不提交到GitHub）",
    "*.rds",
    "*.RDS",
    "raw/TCGA_LPS_real_counts.rds",
    "raw/TCGA_LPS_real_counts_filtered.rds",
    "raw/TCGA_LPS_real_log2cpm.rds",
    "raw/TCGA_6genes_expr.rds",
    "raw/TCGA_Lipid_Score.rds",
    "raw/GSE21122_6genes_expr.rds",
    "raw/GSE21122_Lipid_Score.rds",
    "raw/GSE21122_log2expr.rds",
    "raw/GSE21122_sample_names.rds",
    "raw/GSE159659_6genes_expr.rds",
    "raw/GSE159659_Lipid_Score.rds",
    "",
    "# 原始未处理文件",
    "raw/不用的原始文件/",
    "raw/GSE21653/",
    "raw/GSE30929/",
    "",
    "# 安装程序",
    "raw/Git-*.exe",
    "",
    "# 临时文件",
    "临时.txt",
    "MANIFEST.txt",
    "",
    "# 系统文件",
    ".DS_Store",
    "Thumbs.db"
  ),
  ".gitignore"
)

# ---- 1.2 创建 README.md ----
writeLines(
  c(
    "# A Six-Gene Lipid Metabolism Score Stratifies Prognosis in Dedifferentiated Liposarcoma",
    "",
    "## Overview",
    "This repository contains all analysis code and processed data for the development and validation of a lipid metabolism-based prognostic score in dedifferentiated liposarcoma (DDLPS).",
    "",
    "## Data Sources",
    "",
    "| Dataset | n | Source | Accession |",
    "|---------|---|--------|-----------|",
    "| TCGA-SARC | 50 | The Cancer Genome Atlas | https://portal.gdc.cancer.gov/ |",
    "| GSE21122 | 158 | Gene Expression Omnibus | https://www.ncbi.nlm.nih.gov/geo/query/acc.cgi?acc=GSE21122 |",
    "| GSE159659 | 45 | Gene Expression Omnibus | https://www.ncbi.nlm.nih.gov/geo/query/acc.cgi?acc=GSE159659 |",
    "",
    "## Core Genes",
    "The six-gene lipid metabolism score comprises: **FASN**, **ACLY**, **SCD**, **SREBF1**, **HMGCR**, and **PPARG**.",
    "",
    "## Repository Structure",
    "",
    "```",
    "├── scripts/              # All analysis scripts (Figures 1-6 + Supplementary)",
    "├── results/",
    "│   ├── figures/          # Main and supplementary figures (PDF/PNG/TIFF)",
    "│   └── tables/           # Supplementary tables (S1-S4) + Key Resources Table",
    "├── raw/                  # Processed expression matrices and metadata",
    "├── data/                 # Final gene lists and pathway results",
    "├── sessionInfo.txt       # R session information",
    "├── .gitignore",
    "└── README.md",
    "```",
    "",
    "## Requirements",
    "",
    "- **R** >= 4.2.0",
    "- R packages are documented in `sessionInfo.txt`",
    "",
    "## Quick Start",
    "",
    "```r",
    "# Set working directory",
    "setwd(\"path/to/Liposarcoma/2\")",
    "",
    "# Run scripts in order",
    "source(\"scripts/Figure1_script.R\")   # Study overview and cohort description",
    "source(\"scripts/Figure2_script.R\")   # Six-gene score development",
    "source(\"scripts/Figure3_script.R\")   # Survival and ROC analyses",
    "source(\"scripts/Figure4_script.R\")   # Pathway enrichment and immune analysis",
    "source(\"scripts/Figure5_script.R\")   # Nomogram and calibration",
    "source(\"scripts/Figure6_script.R\")   # External validation",
    "```",
    "",
    "## Supplementary Figures",
    "",
    "```r",
    "source(\"scripts/Supplementary_Figure_S1.R\")",
    "source(\"scripts/Supplementary_Figure_S2.R\")",
    "source(\"scripts/Supplementary_Figure_S3.R\")",
    "source(\"scripts/Supplementary_Figure_S4-S7.R\")",
    "```",
    "",
    "## License",
    "",
    "MIT License",
    "",
    "## Citation",
    "",
    "If you use this code or data, please cite: [Manuscript citation - to be added upon publication]"
  ),
  "README.md"
)

# ---- 1.3 更新 sessionInfo.txt（如果已有则跳过）----
if (!file.exists("sessionInfo.txt") || file.info("sessionInfo.txt")$size < 100) {
  sink("sessionInfo.txt")
  cat("R Session Information for DDLPS Lipid Score Analysis\n")
  cat(paste("Date:", Sys.Date(), "\n\n"))
  sessionInfo()
  sink()
}

# ==================== 第2步：清理不应存放的文件 ====================

# 删除Git安装程序（198 MB）
if (file.exists("raw/Git-2.55.0.3-64-bit.exe")) {
  file.remove("raw/Git-2.55.0.3-64-bit.exe")
  cat("已删除 Git 安装程序\n")
}

# 删除临时文件
for (f in c("临时.txt", "MANIFEST.txt")) {
  if (file.exists(f)) {
    file.remove(f)
    cat("已删除", f, "\n")
  }
}

# ==================== 第3步：移动散落的图表文件至正确位置 ====================

# 将 raw/ 中的图表文件移至 results/figures/
chart_files <- c(
  "raw/Calibration_1yr.png",
  "raw/Calibration_3yr.png",
  "raw/Calibration_5yr.png",
  "raw/Figure_circular_tree_heatmap.pdf",
  "raw/Figure_circular_tree_heatmap.png",
  "raw/Figure_circular_tree_heatmap.svg",
  "raw/Figure_KM_Curve.png",
  "raw/pheatmap.pdf"
)

for (f in chart_files) {
  if (file.exists(f)) {
    dest <- file.path("results/figures", basename(f))
    file.rename(f, dest)
    cat("已移动", f, "→", dest, "\n")
  }
}

# ==================== 第4步：整理根目录散落文件 ====================

# 将 df.rds, results.rds 移入 data/ 或 raw/（按内容决定）
for (f in c("df.rds", "results.rds")) {
  if (file.exists(f)) {
    file.rename(f, file.path("raw", f))
    cat("已移动", f, "→ raw/", f, "\n")
  }
}

# ==================== 第5步：创建最终目录结构确认 ====================

cat("\n========== 整理后的文件结构 ==========\n")
key_files <- list.files(".", recursive = TRUE, pattern = "\\.(R|md|gitignore|txt|csv|pdf)$",
                        full.names = FALSE)
# 排除不用的原始文件
key_files <- key_files[!grepl("不用的原始文件", key_files)]
cat(paste("  ", key_files), sep = "\n")

cat("\n========== 补齐完成！ ==========\n")
cat("下一步：\n")
cat("1. 检查 README.md 内容是否需要修改\n")
cat("2. 运行 Git 初始化和推送（见下方命令）\n")