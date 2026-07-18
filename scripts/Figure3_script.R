# ============================================================
# Figure 3: Subgroup Survival Analysis (Tower-style)
# 出版级规范版本 - 字号全面上调，保证可读性
# 子图: A-E（每行一个亚组，含标签+HR+KM+ROC）
# 输出: PDF（矢量）, TIFF（600 dpi LZW）, PNG（300 dpi）
# ============================================================

rm(list = ls())
setwd("C:/Users/Tao25/Documents/Liposarcoma/1")

library(ggplot2)
library(survival)
library(patchwork)
library(pROC)
library(grid)

# ============================================================
# 1. 版本兼容性检测（消除 size/linewidth 警告）
# ============================================================

ggplot2_version <- packageVersion("ggplot2")
if (ggplot2_version >= "3.4.0") {
  LINE_SIZE <- "linewidth"
} else {
  LINE_SIZE <- "size"
}

# ============================================================
# 2. 颜色系统
# ============================================================

dark_blue <- "#0F4C81"
bright_red <- "#C93312"
light_green <- "#8DAA7D"
warm_orange <- "#E69F00"
indigo <- "#6053a0"

dir.create("figures", recursive = TRUE, showWarnings = FALSE)

# ============================================================
# 3. 数据准备
# ============================================================

scores <- read.csv("data/results/02_feature_scores.csv", row.names = 1)
clinical <- read.csv("data/results/03_clinical_outcomes.csv", row.names = 1)

common_rows <- intersect(rownames(scores), rownames(clinical))
scores <- scores[common_rows, ]
clinical <- clinical[common_rows, ]
valid_idx <- !is.na(scores$Lipid_Metabolism_Score)
scores_valid <- scores[valid_idx, ]
clinical_valid <- clinical[valid_idx, ]

tcga_idx <- clinical_valid$Cohort == "TCGA"
tcga_data <- clinical_valid[tcga_idx, ]
tcga_scores <- scores_valid[tcga_idx, ]

tcga_ok <- !is.na(tcga_data$OS_time) & !is.na(tcga_data$OS_status) & tcga_data$OS_time > 0
tcga_data <- tcga_data[tcga_ok, ]
tcga_scores <- tcga_scores[tcga_ok, ]

tcga_data$AgeGroup <- ifelse(tcga_data$Age <= 60, "Age ≤ 60", "Age > 60")
tcga_data$GenderGroup <- tcga_data$Gender

# ============================================================
# 4. 定义亚组
# ============================================================

floors <- list()
floors[[1]] <- list(name = "All Patients", n = nrow(tcga_data), data = tcga_data, scores = tcga_scores$Lipid_Metabolism_Score)

idx2 <- tcga_data$AgeGroup == "Age ≤ 60"
if (sum(idx2) >= 10) {
  floors[[length(floors)+1]] <- list(name = "Age ≤ 60", n = sum(idx2), data = tcga_data[idx2, ], scores = tcga_scores$Lipid_Metabolism_Score[idx2])
}

idx3 <- tcga_data$AgeGroup == "Age > 60"
if (sum(idx3) >= 10) {
  floors[[length(floors)+1]] <- list(name = "Age > 60", n = sum(idx3), data = tcga_data[idx3, ], scores = tcga_scores$Lipid_Metabolism_Score[idx3])
}

idx4 <- tcga_data$GenderGroup == "female"
if (sum(idx4) >= 10) {
  floors[[length(floors)+1]] <- list(name = "Female", n = sum(idx4), data = tcga_data[idx4, ], scores = tcga_scores$Lipid_Metabolism_Score[idx4])
}

idx5 <- tcga_data$GenderGroup == "male"
if (sum(idx5) >= 10) {
  floors[[length(floors)+1]] <- list(name = "Male", n = sum(idx5), data = tcga_data[idx5, ], scores = tcga_scores$Lipid_Metabolism_Score[idx5])
}

# ============================================================
# 5. 绘图函数（字号全面上调）
# ============================================================

# 5.1 标签面板（无字母，字母由外部统一添加）
plot_label <- function(name, n) {
  p <- ggplot() +
    annotate("text", x = 0.5, y = 0.6, label = name,
             size = 5.5, fontface = "plain", family = "serif", color = dark_blue) +
    annotate("text", x = 0.5, y = 0.3, label = paste0("n = ", n),
             size = 4.5, fontface = "plain", family = "serif", color = dark_blue) +
    xlim(0, 1) + ylim(0, 1) +
    theme_void() +
    theme(panel.border = element_rect(color = "gray70", fill = NA, linewidth = 0.5))
  return(p)
}

# 5.2 HR面板
plot_hr <- function(data, scores) {
  if (nrow(data) < 5) {
    return(ggplot() + annotate("text", x = 0.5, y = 0.5, label = "NA", size = 4.5, family = "serif") + theme_void())
  }
  
  med <- median(scores)
  df <- data.frame(time = data$OS_time, status = data$OS_status, group = ifelse(scores > med, "High", "Low"))
  df$group <- factor(df$group, levels = c("Low", "High"))
  
  if (length(unique(df$group)) < 2) {
    return(ggplot() + annotate("text", x = 0.5, y = 0.5, label = "NA", size = 4.5, family = "serif") + theme_void())
  }
  
  fit <- tryCatch(coxph(Surv(time, status) ~ group, data = df), error = function(e) NULL)
  if (is.null(fit)) {
    return(ggplot() + annotate("text", x = 0.5, y = 0.5, label = "Error", size = 4.5, family = "serif") + theme_void())
  }
  
  hr <- round(exp(coef(fit)["groupHigh"]), 2)
  ci <- round(exp(confint(fit)["groupHigh", ]), 2)
  p_val <- summary(fit)$coefficients["groupHigh", "Pr(>|z|)"]
  
  hr_text <- sprintf("HR = %.2f\n(%.2f-%.2f)", hr, ci[1], ci[2])
  p_text <- if (p_val < 0.001) "P < 0.001" else if (p_val < 0.01) "P < 0.01" else if (p_val < 0.05) "P < 0.05" else sprintf("P = %.2f", p_val)
  color <- if (hr > 1) bright_red else dark_blue
  
  ggplot() +
    annotate("rect", xmin = 0, xmax = 1, ymin = 0, ymax = 1, fill = color, alpha = 0.12) +
    annotate("text", x = 0.5, y = 0.65, label = hr_text,
             size = 4.5, fontface = "bold", family = "serif", color = dark_blue) +
    annotate("text", x = 0.5, y = 0.30, label = p_text,
             size = 4.0, fontface = "plain", family = "serif", color = dark_blue) +
    xlim(0, 1) + ylim(0, 1) +
    theme_void() +
    theme(panel.border = element_rect(color = "gray70", fill = NA, linewidth = 0.5))
}

# 5.3 KM曲线
plot_km <- function(data, scores) {
  if (nrow(data) < 5) {
    return(ggplot() + annotate("text", x = 0.5, y = 0.5, label = "Insufficient data", size = 4, family = "serif") + theme_void())
  }
  
  med <- median(scores)
  df <- data.frame(time = data$OS_time, status = data$OS_status, group = ifelse(scores > med, "High", "Low"))
  df$group <- factor(df$group, levels = c("Low", "High"))
  
  fit <- survfit(Surv(time, status) ~ group, data = df)
  sdiff <- survdiff(Surv(time, status) ~ group, data = df)
  p_val <- 1 - pchisq(sdiff$chisq, 1)
  p_text <- if (p_val < 0.001) "P < 0.001" else if (p_val < 0.01) "P < 0.01" else if (p_val < 0.05) "P < 0.05" else sprintf("P = %.2f", p_val)
  
  km_data <- data.frame(
    time = fit$time,
    surv = fit$surv,
    group = rep(c("Low", "High"), times = fit$strata)
  )
  
  p <- ggplot(km_data, aes(x = time / 365, y = surv, color = group)) +
    geom_step(linewidth = 0.9) +
    scale_color_manual(values = c(light_green, bright_red)) +
    annotate("text", x = max(km_data$time / 365) * 0.55, y = 0.20,
             label = p_text, size = 4.0, fontface = "plain", family = "serif", color = dark_blue, hjust = 0) +
    labs(x = "", y = "") +
    coord_cartesian(ylim = c(0, 1)) +
    theme_minimal(base_size = 9) +
    theme(
      plot.title = element_blank(),
      legend.position = "none",
      panel.grid.minor = element_blank(),
      panel.grid.major = element_line(linewidth = 0.2, color = "#E8E8E8"),
      axis.text = element_text(size = 8, color = "#333333", family = "serif"),
      panel.border = element_rect(color = "gray70", fill = NA, linewidth = 0.5)
    )
  return(p)
}

# 5.4 ROC曲线
plot_roc <- function(data, scores) {
  if (nrow(data) < 5) {
    return(ggplot() + annotate("text", x = 0.5, y = 0.5, label = "Insufficient data", size = 4, family = "serif") + theme_void())
  }
  
  df <- data.frame(time = data$OS_time, status = data$OS_status)
  times <- c(365, 1095, 1825)
  roc_data <- data.frame()
  
  for (t in times) {
    status_t <- ifelse(df$time <= t & df$status == 1, 1, 0)
    if (length(unique(status_t)) > 1 && sum(status_t) > 0 && sum(1 - status_t) > 0) {
      roc_obj <- roc(status_t, scores, quiet = TRUE)
      auc_val <- round(auc(roc_obj), 2)
      roc_df <- data.frame(
        FP = 1 - roc_obj$specificities,
        TP = roc_obj$sensitivities,
        Time = paste0(t / 365, "y (AUC=", auc_val, ")")
      )
      roc_data <- rbind(roc_data, roc_df)
    }
  }
  
  if (nrow(roc_data) == 0) {
    return(ggplot() + annotate("text", x = 0.5, y = 0.5, label = "ROC N/A", size = 4, family = "serif") + theme_void())
  }
  
  p <- ggplot(roc_data, aes(x = FP, y = TP, color = Time)) +
    geom_line(linewidth = 0.8) +
    geom_abline(intercept = 0, slope = 1, linetype = "dashed", color = "gray50", linewidth = 0.4) +
    scale_color_manual(values = c(light_green, warm_orange, bright_red)) +
    labs(x = "", y = "") +
    coord_cartesian(xlim = c(0, 1), ylim = c(0, 1)) +
    theme_minimal(base_size = 8) +
    theme(
      plot.title = element_blank(),
      legend.position = "none",
      panel.grid.minor = element_blank(),
      panel.grid.major = element_line(linewidth = 0.2, color = "#E8E8E8"),
      axis.text = element_text(size = 7, color = "#333333", family = "serif"),
      panel.border = element_rect(color = "gray70", fill = NA, linewidth = 0.5)
    )
  return(p)
}

# ============================================================
# 6. 字母模板（独立左列，字号增大）
# ============================================================

letter_template <- function(label) {
  ggplot() +
    annotate("text", x = 0.5, y = 0.5, label = label,
             size = 8.0, fontface = "plain", family = "serif", color = dark_blue) +
    xlim(0, 1) + ylim(0, 1) +
    theme_void()
}

# ============================================================
# 7. 构建楼层（每个楼层：字母 | 标签 | HR | KM | ROC）
# ============================================================

row_labels <- c("A", "B", "C", "D", "E")
row_indices <- 1:5

floor_plots <- list()
for (i in seq_along(row_indices)) {
  idx <- row_indices[i]
  f <- floors[[idx]]
  
  letter_panel <- letter_template(row_labels[i])
  label_panel <- plot_label(f$name, f$n)
  hr_panel <- plot_hr(f$data, f$scores)
  km_panel <- plot_km(f$data, f$scores)
  roc_panel <- plot_roc(f$data, f$scores)
  
  floor_row <- (letter_panel | label_panel | hr_panel | km_panel | roc_panel) +
    plot_layout(widths = c(0.3, 0.5, 0.9, 1.2, 1.2))
  
  floor_plots[[i]] <- floor_row
}

# ============================================================
# 8. 组合
# ============================================================

figure3 <- wrap_plots(floor_plots, ncol = 1) +
  plot_annotation(
    title = "Figure 3: Subgroup Analysis of Lipid Metabolism Score (TCGA Cohort)",
    subtitle = "Tower-style visualization of prognostic performance",
    caption = paste(
      "Note: High-risk group (red); Low-risk group (green). HR > 1 indicates higher risk.",
      "Subgroups: A: All Patients, B: Age ≤ 60, C: Age > 60, D: Female, E: Male.",
      "KM: Kaplan-Meier curves with log-rank P values. ROC: Time-dependent ROC curves.",
      sep = "\n"
    ),
    theme = theme(
      plot.title = element_text(hjust = 0.5, face = "bold", size = 14, color = dark_blue,
                                family = "serif", margin = margin(b = 6)),
      plot.subtitle = element_text(hjust = 0.5, face = "plain", size = 10, color = indigo,
                                   family = "serif", margin = margin(b = 4)),
      plot.caption = element_text(hjust = 0, size = 7, color = "gray50",
                                  family = "serif", margin = margin(t = 6), lineheight = 1.3)
    )
  )

# ============================================================
# 9. 出版级输出（多格式 + 高分辨率）
# ============================================================

total_height <- 0.8 + length(floors) * 2.3

ggsave("figures/Figure3.pdf", figure3, width = 14, height = total_height, units = "in", dpi = 600, limitsize = FALSE, device = cairo_pdf)
ggsave("figures/Figure3.tiff", figure3, width = 14, height = total_height, units = "in", dpi = 600, compression = "lzw")
ggsave("figures/Figure3.png", figure3, width = 14, height = total_height, units = "in", dpi = 300)

cat("\n============================================================\n")
cat("✅ Figure 3 已保存（出版级，字号已上调）\n")
cat("   PDF:  figures/Figure3.pdf\n")
cat("   TIFF: figures/Figure3.tiff\n")
cat("   PNG:  figures/Figure3.png\n")
cat(sprintf("总高度: %.1f 英寸\n", total_height))
cat(sprintf("楼层数: %d\n", length(floors)))
cat("\n字号设置:\n")
cat("  字母编号: 8.0pt\n")
cat("  亚组名称: 5.5pt\n")
cat("  样本量: 4.5pt\n")
cat("  HR值: 4.5pt\n")
cat("  p值: 4.0pt\n")
cat("  KM轴标签: 8pt\n")
cat("  ROC轴标签: 7pt\n")
cat("============================================================\n")