# ============================================================
# Fig S3I · 临床加固：评分 A 的 KM 分层（主图 = 2 层；两版 4 层进补充）
#   （原 Fig 6F，2026-09-28 重构统一为 Fig. S3 I）
#   上栏：评分 A 2 层 KM（中位切分）+ 风险表 —— 唯一显著生存读数（P = 0.017）
#   下栏：两版 4 层 KM 对照（Score x RICTOR / Score x CD36）—— 均不显著
#   ★ 曲线自绘（不依赖 survminer），确保像素与字号可控
#   数据：Liposarcoma/1/data/raw/TCGA_Score.csv + TCGA_clinical_real.csv
#         + TCGA_LPS_real_log2cpm.rds（RICTOR = ENSG00000164327 权威 ID）
# ============================================================
suppressWarnings(suppressMessages({library(ggplot2); library(patchwork); library(survival); library(scales)}))
source(file.path(file.path(Sys.getenv("LPS_FIG_ROOT", "."), "utils"), "themes.R"))
source(file.path(file.path(Sys.getenv("LPS_FIG_ROOT", "."), "utils"), "figS6_common.R"))

Root <- file.path(Sys.getenv("LPS_DATA_ROOT", "."), "1")
sc <- read.csv(file.path(Root, "data/raw/TCGA_Score.csv"), stringsAsFactors = FALSE)
cl <- read.csv(file.path(Root, "data/raw/TCGA_clinical_real.csv"), check.names = FALSE)
x  <- readRDS(file.path(Root, "data/raw/TCGA_LPS_real_log2cpm.rds"))
rn <- sub("\\..*$", "", rownames(x))
ENS <- c(RICTOR = "ENSG00000164327", CD36 = "ENSG00000135218")   # D27 权威 ID
miss <- setdiff(ENS, rn); stopifnot(!length(miss))
E <- as.data.frame(t(x[match(ENS, rn), , drop = FALSE])); colnames(E) <- names(ENS)
E$Sample <- rownames(E)
df <- merge(sc, E, by = "Sample")
df <- merge(df, cl[, c("barcode", "OS_time", "OS_status")], by.x = "Sample", by.y = "barcode")
df <- df[!is.na(df$OS_time) & df$OS_time > 0 & !is.na(df$OS_status), ]
cat(sprintf("S3I: n=%d | events=%d\n", nrow(df), sum(df$OS_status == 1)))
stopifnot(nrow(df) == 50, sum(df$OS_status == 1) == 23)

med <- function(v) median(v)
df$S <- factor(ifelse(df$Score  > med(df$Score),  "Score high", "Score low"),  c("Score low", "Score high"))
df$R <- factor(ifelse(df$RICTOR > med(df$RICTOR), "RICTOR high", "RICTOR low"), c("RICTOR low", "RICTOR high"))
df$C <- factor(ifelse(df$CD36   > med(df$CD36),   "CD36 high",  "CD36 low"),   c("CD36 low", "CD36 high"))
mk4 <- function(a, b, an, bn) {
  ka <- ifelse(endsWith(as.character(a), "high"), "hi", "lo")
  kb <- ifelse(endsWith(as.character(b), "high"), "hi", "lo")
  mp <- c("lo/lo" = paste0(an, "-lo/", bn, "-lo"), "lo/hi" = paste0(an, "-lo/", bn, "-hi"),
          "hi/lo" = paste0(an, "-hi/", bn, "-lo"), "hi/hi" = paste0(an, "-hi/", bn, "-hi"))
  factor(unname(mp[paste(ka, kb, sep = "/")]),
         levels = unname(mp[c("lo/lo", "lo/hi", "hi/lo", "hi/hi")]))
}
df$SR <- mk4(df$S, df$R, "S", "R")
df$SC <- mk4(df$S, df$C, "S", "C")
stopifnot(all(nchar(levels(df$SR)) == 9), all(nchar(levels(df$SC)) == 9))
stopifnot(all(table(df$SR) > 0), all(table(df$SC) > 0), nlevels(df$SR) == 4)

XMAX <- 3150
BREAKS <- c(0, 1000, 2000, 3000)

## ---------------------------------------------------------------- 上栏（KM + 风险表）
KM_LM <- margin(3, 8, 2, 0)
COL2 <- c("Score low" = "#2980B9", "Score high" = "#C0392B")
pS2 <- km_panel(df, "S", COL2, "Score A, median split", xmax = XMAX, base_size = 7.5,
                xbreaks = BREAKS)
P2 <- km_p(df, "S")
cat(sprintf("   2-strata log-rank P = %.4f\n", P2))
pS2 <- pS2 + labs(subtitle = wrap_txt(
    sprintf("25 vs 25 patients, 8 vs 15 events; log-rank P = %.3f", P2), 30)) +
  theme(plot.subtitle = element_text(size = 6.4, colour = "grey25"),
        legend.position = c(0.80, 0.84),
        plot.margin = KM_LM)

rk <- km_risk(df, "OS_time", "S", BREAKS)
rk$group <- factor(rk$group, levels = levels(df$S))
pRk <- ggplot(rk, aes(time, group, label = n)) +
  geom_text(size = 2.15, colour = "grey15") +
  scale_x_continuous(limits = c(0, XMAX), breaks = BREAKS, expand = c(0, 0)) +
  scale_y_discrete(expand = c(0.35, 0.35)) +
  labs(x = NULL, y = NULL, title = "Number at risk") +
  theme_publication(base_size = 7.5) +
  theme(axis.line = element_blank(), axis.ticks = element_blank(),
        axis.text.x = element_blank(),
        axis.text.y = element_text(size = 6.0),
        plot.title = element_text(size = 6.4, face = "plain"),
        plot.margin = margin(1, 8, 2, 0))

## ---------------------------------------------------------------- 下栏（两版 4 层）
COL4 <- c("#A9CCE3", "#2980B9", "#F5B7B1", "#C0392B")
SHORT_LAB <- c("S-lo/R-lo" = "lo/lo", "S-lo/R-hi" = "lo/hi",
               "S-hi/R-lo" = "hi/lo", "S-hi/R-hi" = "hi/hi",
               "S-lo/C-lo" = "lo/lo", "S-lo/C-hi" = "lo/hi",
               "S-hi/C-lo" = "hi/lo", "S-hi/C-hi" = "hi/hi")

m4 <- function(g, ttl) {
  lv <- levels(df[[g]])
  cm <- setNames(COL4[seq_along(lv)], lv)
  p <- km_panel(df, g, cm, ttl, xmax = XMAX, base_size = 7.0,
                xbreaks = c(0, 1500, 3000))
  pp <- km_p(df, g)
  p <- p + scale_x_continuous(limits = c(0, XMAX), breaks = c(0, 1500, 3000),
                              expand = expansion(mult = c(0, 0.02)))
  p + labs(subtitle = wrap_txt(sprintf("P = %.3f", pp), 28),
          y = NULL, x = "Time (days)") +
    scale_colour_manual(values = cm, name = NULL, breaks = lv,
                        labels = unname(SHORT_LAB[lv])) +
    theme(plot.title = element_text(size = 6.6),
          plot.title.position = "plot",
          plot.subtitle = element_text(size = 6.0, colour = "grey25"),
          legend.position = c(0.045, 0.045),
          legend.justification = c(0, 0),
          legend.direction = "vertical",
          legend.background = element_rect(fill = "white", colour = NA),
          legend.text = element_text(size = 5.2),
          legend.key.height = unit(5.4, "pt"), legend.key.width = unit(7, "pt"),
          legend.spacing.y = unit(0.6, "pt"),
          axis.text = element_text(size = 5.6),
          axis.title = element_text(size = 6.2),
          plot.margin = margin(3, 8, 1, 0))
}
pSR <- m4("SR", "Score x RICTOR")
pSC <- m4("SC", "Score x CD36")
cat(sprintf("   4-strata P: SR=%.4f  SC=%.4f\n", km_p(df, "SR"), km_p(df, "SC")))

p <- pS2 + pRk + pSR + pSC +
  plot_layout(design = "
AAAA
BBBB
CCDD
", heights = c(3.15, 1, 2.55))
save6(p, "figS3i")

## 中位 OS（供图注核对）
for (g in c("S", "SR", "SC")) {
  fit <- survfit(as.formula(paste("Surv(OS_time, OS_status) ~", g)), data = df)
  cat(sprintf("   %-3s median OS: %s\n", g,
              paste(sprintf("%.0f", summary(fit)$table[, "median"]), collapse = " / ")))
}
