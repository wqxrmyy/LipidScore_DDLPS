# ============================================================
# DDLPS LipidScore Extended Pipeline
# 功能：一键执行全部分析（多数据库版本）
# ============================================================
setwd("C:/Users/Tao25/Documents/Liposarcoma/2")

cat("\n========================================\n")
cat("  Stage 1 Extended: Multi-database download\n")
cat("========================================\n")
source("scripts/stage1_extended_data.R")

cat("\n========================================\n")
cat("  Stage 2 Advanced: Multi-cohort validation\n")
cat("========================================\n")
source("scripts/stage2_advanced_prognostic.R")

cat("\n========================================\n")
cat("  Stage 3 Multiomics: Genomics + scRNA + Immune\n")
cat("========================================\n")
source("scripts/stage3_multiomics_mechanism.R")

cat("\n========================================\n")
cat("  Stage 4 Translation: Drug sensitivity + clinical\n")
cat("========================================\n")
source("scripts/stage2_prognostic_model/05_GSE30929_validation.R")
source("scripts/stage2_prognostic_model/05_GSE30929_validation.R")
source("scripts/stage4_translation_extended.R")

cat("\n========================================\n")
cat("  Pipeline completed.\n")
cat("========================================\n")
