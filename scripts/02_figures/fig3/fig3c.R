# ============================================================
# Fig 3C · 免疫层阴性披露（正确 CD8A 校正前后 NS，棒棒糖）
# S18 修正后：n=50 r=-0.2311 P=0.106；控纯度 r=-0.0328 P=0.835
# ============================================================
suppressWarnings(suppressMessages(library(ggplot2)))
source(file.path(Sys.getenv("LPS_FIG_ROOT", "."), "utils", "colors.R"))
source(file.path(Sys.getenv("LPS_FIG_ROOT", "."), "utils", "themes.R"))

load(file.path(Sys.getenv("LPS_FIG_ROOT", "."), "p5_data.RData"))
OUT <- file.path(Sys.getenv("LPS_FIG_ROOT", "."), "out")

numf <- function(x) suppressWarnings(as.numeric(as.character(x)))
log2cpm <- readRDS(file.path(Sys.getenv("LPS_DATA_ROOT", "."), "1", "data", "raw", "TCGA_LPS_real_log2cpm.rds"))
idx <- setNames(rownames(log2cpm), sub("\\.[0-9]+$", "", rownames(log2cpm)))
cd8a <- as.numeric(log2cpm[idx["ENSG00000153563"], score_A$Sample])
pur <- numf(clin_real[score_A$Sample, "paper_purity"])
y <- score_A$Score

t1 <- cor.test(cd8a, y)
keep <- !is.na(pur)
rx <- residuals(lm(cd8a[keep] ~ pur[keep])); ry <- residuals(lm(y[keep] ~ pur[keep]))
t2 <- cor.test(rx, ry)

df <- data.frame(
  Model = c("Raw (n = 50)", "Purity-adjusted (n = 44)"),
  r = c(unname(t1$estimate), unname(t2$estimate)),
  p = c(t1$p.value, t2$p.value)
)
df$Model <- factor(df$Model, levels = df$Model)
df$lab <- sprintf("r = %.3f, P = %.2f (ns)", df$r, df$p)

p <- ggplot(df, aes(x = Model, y = r)) +
  geom_hline(yintercept = 0, color = "grey50") +
  geom_segment(aes(x = Model, xend = Model, y = 0, yend = r), color = "#7F8C8D", size = 1.2) +
  geom_point(size = 4.5, color = "#7F8C8D") +
  geom_text(aes(label = lab), vjust = 1.8, size = 3) +
  coord_cartesian(ylim = c(-0.5, 0.15)) +
  labs(x = "", y = "Pearson r",
       title = "Immune layer: negative disclosure",
       subtitle = "Corrected CD8A probe: no robust association") +
  theme_publication()

export_fig(p, file.path(OUT, "fig3c.pdf"), w = 4.8, h = 3.5)
cat("fig3c done\n")
