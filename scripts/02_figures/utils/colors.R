# ============================================================
# P5 出图 · 全局配色方案（对齐旧 PDF v1.0 的 colors.R 规范）
# ============================================================
COLORS <- list(
  # 组别色
  group = c(
    Ctrl      = "#7F8C8D",  # 灰色
    High      = "#E74C3C",  # 红色（高评分）
    Low       = "#2980B9",  # 蓝色（低评分）
    WD        = "#2980B9",  # 蓝色（well-differentiated 分化好）
    DD        = "#E67E22",  # 橙色（dedifferentiated 去分化）
    Inhibitor = "#E67E22",  # 橙色
    Rescue    = "#27AE60",  # 绿色
    MDM2_amp  = "#C0392B"   # 深红
  ),
  # 细胞类型色（保留以备扩展）
  celltype = c(
    Tumor      = "#E74C3C",
    Fibroblast = "#3498DB",
    Endothelial= "#1ABC9C",
    T_cell     = "#9B59B6",
    NK_cell    = "#F39C12",
    B_cell     = "#2ECC71",
    Myeloid    = "#E67E22",
    Mast       = "#D35400",
    Adipocyte  = "#F1C40F"
  ),
  # 连续色阶
  heatmap   = c("#2980B9", "#FFFFFF", "#E74C3C"),  # 蓝-白-红
  umap_grad = c("#5DADE2", "#F7DC6F", "#EC7063")   # 蓝-黄-红
)

# 统计标注辅助
sig_label <- function(p) {
  if (is.na(p)) return("ns")
  if (p < 0.0001) return("****")
  if (p < 0.001)  return("***")
  if (p < 0.01)   return("**")
  if (p < 0.05)   return("*")
  "ns"
}
