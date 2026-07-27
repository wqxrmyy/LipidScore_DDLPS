# ======================================================
# Stage1 扩展版：多数据库批量下载、整合、批次校正
# 新增队列：GSE30929、ICGC-SARC、GTEx正常脂肪对照
# 原有队列：TCGA-SARC、GSE21122、GSE159659
# ======================================================
# 0. 加载依赖包
if (!require("BiocManager", quietly = TRUE)) install.packages("BiocManager")
pkg_list = c("GEOquery","TCGAbiolinks","ICGCportal","sva","ComBat","tidyverse","limma","rtracklayer","readr")
for(pkg in pkg_list){
  if(!require(pkg,character.only = T)) BiocManager::install(pkg,update=F)
}
library(GEOquery)
library(TCGAbiolinks)
library(sva)
library(sva)
library(limma)
library(tidyverse)

# 路径配置（匹配原仓库目录结构）
root_path = getwd()
raw_exp_path = file.path(root_path,"raw","bulk")
raw_clin_path = file.path(root_path,"raw")
dir.create(raw_exp_path,recursive = T,showWarnings = F)
dir.create(raw_clin_path,recursive = T,showWarnings = F)

# 1. 批量下载GEO队列函数
download_geo_cohort <- function(geo_id){
  gse <- getGEO(geo_id, destdir = raw_exp_path, getGPL = F)
  expr_raw <- exprs(gse[[1]])
  clin_raw <- Biobase::pData(gse[[1]])
  saveRDS(expr_raw,file.path(raw_exp_path,paste0(geo_id,"_expr.rds")))
  saveRDS(clin_raw,file.path(raw_exp_path,paste0(geo_id,"_clin.rds")))
  message(paste0(geo_id," 下载完成"))
  return(list(expr=expr_raw,clin=clin_raw))
}

# 批量拉取新增+原有GEO数据集
geo_ids <- c("GSE21122","GSE159659","GSE30929")
geo_list <- lapply(geo_ids,download_geo_cohort)
names(geo_list) <- geo_ids

# 2. TCGA-SARC 下载与标准化
# 匹配原项目TCGA软组织肉瘤队列
query_tcga <- GDCquery(project = "TCGA-SARC",
                       data.category = "Transcriptome Profiling",
                       data.type = "Gene Expression Quantification",
                       workflow.type = "STAR - Counts")
GDCdownload(query_tcga,directory = raw_exp_path)
tcga_data <- GDCprepare(query_tcga,directory = raw_exp_path)
tcga_expr <- assay(tcga_data,"unstranded")
# TPM归一化处理
tcga_tpm <- convertCountstoTPM(tcga_expr)
tcga_clin <- colData(tcga_data) %>% as.data.frame()
saveRDS(tcga_tpm,file.path(raw_exp_path,"TCGA-SARC_expr.rds"))
saveRDS(tcga_clin,file.path(raw_exp_path,"TCGA-SARC_clin.rds"))

# 3. ICGC-SARC 肉瘤队列读取
# ICGC官网提前下载count矩阵后读取
icgc_raw <- readRDS(file.path(raw_exp_path,"ICGC-SARC_raw.rds"))
icgc_tpm <- convertCountstoTPM(icgc_raw$count)
icgc_clin <- icgc_raw$clinical
saveRDS(icgc_tpm,file.path(raw_exp_path,"ICGC-SARC_expr.rds"))
saveRDS(icgc_clin,file.path(raw_exp_path,"ICGC-SARC_clin.rds"))

# 4. GTEx正常脂肪组织对照（DDLPS正常对照基线）
gtex_adipose <- readRDS(file.path(raw_exp_path,"GTEx_Adipose_expr.rds"))
gtex_clin <- readRDS(file.path(raw_exp_path,"GTEx_Adipose_clin.rds"))

# 5. 基因ID统一转换（所有队列统一Gene Symbol）
# 加载基因注释文件（原仓库自带gene_anno.csv）
gene_anno <- read_csv(file.path(root_path,"data","gene_anno.csv"))
unify_gene_symbol <- function(mat,anno){
  mat$gene <- rownames(mat)
  mat <- inner_join(mat,anno,by="gene")
  mat <- mat %>% group_by(Symbol) %>% summarise_all(mean,na.rm=T)
  mat <- column_to_rownames(mat,"Symbol")
  return(as.matrix(mat))
}

# 批量统一所有队列基因名
all_expr_list <- list(
  TCGA=tcga_tpm,
  GSE21122=geo_list$GSE21122$expr,
  GSE159659=geo_list$GSE159659$expr,
  GSE30929=geo_list$GSE30929$expr,
  ICGC=icgc_tpm,
  GTEx=gtex_adipose
)
all_expr_unified <- lapply(all_expr_list,unify_gene_symbol,anno=gene_anno)

# 6. 取交集基因，合并全部表达矩阵
common_genes <- Reduce(intersect,lapply(all_expr_unified,rownames))
expr_merge_list <- lapply(all_expr_unified,function(x) x[common_genes,])
# 标记队列来源batch标签
for(b in names(expr_merge_list)){
  colnames(expr_merge_list[[b]]) <- paste0(b,"_",colnames(expr_merge_list[[b]]))
}
expr_all_raw <- do.call(cbind,expr_merge_list)
batch_label <- rep(names(expr_merge_list),times=sapply(expr_merge_list,ncol))

# 7. ComBat跨平台批次校正（消除芯片/测序平台差异）
expr_corrected <- ComBat(dat=expr_all_raw,batch=batch_label,mod=NULL)
saveRDS(expr_corrected,file.path(raw_exp_path,"exp_all_cohorts_corrected.rds"))
message("批次校正完成，整合表达矩阵已保存")

# 8. 整合全部临床信息，统一临床字段
clean_clin <- function(clin_df,cohort_name){
  clin_df$sample_id <- paste0(cohort_name,"_",rownames(clin_df))
  clin_df$cohort <- cohort_name
  # 统一关键临床字段：OS、OS.time、status、age、gender、tumor_size、MDM2_amp
  standard_col <- c("sample_id","cohort","age","gender","OS.time","OS","tumor_size","MDM2_amp")
  for(col in standard_col){
    if(!col %in% colnames(clin_df)) clin_df[[col]] <- NA
  }
  clin_df <- clin_df[,standard_col]
  return(clin_df)
}

# 批量清洗各队列临床数据
clin_list <- list(
  TCGA=tcga_clin,
  GSE21122=geo_list$GSE21122$clin,
  GSE159659=geo_list$GSE159659$clin,
  GSE30929=geo_list$GSE30929$clin,
  ICGC=icgc_clin
)
clin_clean <- map2_dfr(clin_list,names(clin_list),clean_clin)
# 保存合并临床表型
write_csv(clin_clean,file.path(raw_clin_path,"clinical_merged_all_cohorts.csv"))
saveRDS(clin_clean,file.path(raw_clin_path,"clinical_merged_all_cohorts.rds"))
message("全部队列临床信息整合完毕")

# 9. 输出队列统计表格（写入data目录用于论文附表）
cohort_stat <- clin_clean %>% 
  group_by(cohort) %>% 
  summarise(
    total_sample=n(),
    death_case=sum(OS==1,na.rm=T),
    median_age=median(age,na.rm=T),
    MDM2_amp_rate=mean(MDM2_amp=="yes",na.rm=T)
  )
write_csv(cohort_stat,file.path(root_path,"data","cohort_basic_stat.csv"))

message("==== Stage1 全部数据预处理完成 ====")
message("输出文件：")
message("1. raw/bulk/exp_all_cohorts_corrected.rds 校正后合并表达矩阵")
message("2. raw/clinical_merged_all_cohorts.csv 统一临床信息表")
message("3. data/cohort_basic_stat.csv 各队列样本统计附表")