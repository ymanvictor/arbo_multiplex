# --- Libraries ---------------------------------------------------------------
library(tidyverse)
library(cowplot)
library(pROC)
library(ggrepel)
library(ggsci)

# --- Palette ----------------------------------------------------------------
cols <- pal_jama()(6)

# --- Data -------------------------------------------------------------------
# Posterior stats table from GMMM fitted to CHIKV VLP data

post_stats <- read_csv("output/gmm_posterior_stats_CHIKV_VLP.csv") %>%
  rename(param = 1)

# Named vector of posterior medians for direct access like post_med["mu[1]"]

post_med <- setNames(post_stats$post_med, post_stats$param)

# Normalised MFI data for control samples

ctrls <- read_csv("raw_data/CHIKV_VLP.csv") %>% 
  filter(sampletype != "U")

# Estimated empirical and analytical AUROCS for all antibody responses
aucs <- read_csv("output/empirical_vs_analytical_aucs.csv")

# --- Helpers ----------------------------------------------------------------

# CDFs for negative (D2) and positive (D1) classes under posterior medians

# (The lowest and middle components represent sero-negatives for CHIKV VLP)

CDF_D2_med <- function(x) {
  w1 <- post_med["theta[1]"]; w2 <- post_med["theta[2]"]
  m1 <- post_med["mu[1]"];    s1 <- post_med["sigma[1]"]
  m2 <- post_med["mu[2]"];    s2 <- post_med["sigma[2]"]
  1 - ((w1/(w1 + w2)) * pnorm(x, m1, s1) + (w2/(w1 + w2)) * pnorm(x, m2, s2))
}
CDF_D1_med <- function(x) {
  m3 <- post_med["mu[3]"]; s3 <- post_med["sigma[3]"]
  1 - pnorm(x, m3, s3)
}

# Trapezoidal AUC; sorts by FPR, pads with (0,0) and (1,1)

auc_trap <- function(tpr, fpr) {
  stopifnot(length(tpr) == length(fpr))
  ok <- is.finite(tpr) & is.finite(fpr)
  tpr <- pmin(pmax(tpr[ok], 0), 1)
  fpr <- pmin(pmax(fpr[ok], 0), 1)
  
  o <- order(fpr, tpr)        # ensure ascending FPR (tie-break by TPR)
  x <- c(0, fpr[o], 1)        # FPR axis
  y <- c(0, tpr[o], 1)        # TPR axis
  
  sum( (y[-1] + y[-length(y)]) * diff(x) / 2 )
}

# Invert CDF_D2_med to find threshold at desired specificity (1 - FPR)

find_quantile <- function(prob, lower = 0, upper = 20) {
  uniroot(function(x) CDF_D2_med(x) - prob, lower = lower, upper = upper)$root
}

# --- Analytical ROC (posterior medians) -------------------------------------

thresholds <- seq(0, 6, length.out = 2000)
analytical_roc_med <- tibble(
  threshold = thresholds,
  tpr = sapply(thresholds, CDF_D1_med),
  fpr = sapply(thresholds, CDF_D2_med)
) %>% arrange(fpr)

# Calculate AUROC

analyt_auc_CHIKV_VLP <- auc_trap(analytical_roc_med$tpr, analytical_roc_med$fpr)

# 97.7% specificity point (threshold)
serothr_977 <- find_quantile(1 - 0.97725)
roc_serothr_977 <- tibble(
  tpr = CDF_D1_med(serothr_977),
  fpr = CDF_D2_med(serothr_977)
)

# --- Plot: Analytical ROC ----------------------------------------------------

plot_analytical_roc <-
  ggplot(analytical_roc_med, aes(x = fpr, y = tpr)) +
  geom_path(linewidth = 1.5, color = cols[3]) +
  geom_point(data = roc_serothr_977, size = 3, color = "darkgreen") +
  geom_label_repel(
    data = roc_serothr_977,
    aes(label = "Threshold of 97.7 % specificity"),
    color = "darkgreen", nudge_y = -0.05
  ) +
  annotate("label", x = 0.5, y = 0.5,
           label = paste0("AUC: ", round(analyt_auc_CHIKV_VLP, 3)), size = 6) +
  labs(
    title = "CHIKV_VLP: Analytical ROC Curve",
    x = "False positivity rate",
    y = "True positivity rate"
  ) +
  theme_bw() +
  theme(
    plot.title = element_text(hjust = 0.5, size = 22),
    axis.text  = element_text(size = 20),
    axis.title = element_text(size = 22),
    aspect.ratio = 1
  )

# --- Empirical ROC from controls --------------------------------------------

roc_obj <- pROC::roc(sampletype ~ log_mfi, data = ctrls, ci = TRUE)

roc_empirical_df <- tibble(
  fpr = 1 - roc_obj$specificities,
  tpr = roc_obj$sensitivities,
  threshold = roc_obj$thresholds
) %>% arrange(fpr)

empir_auc_CHIKV_VLP <- roc_obj$ci

plot_empirical_roc <-
  ggplot(roc_empirical_df, aes(x = fpr, y = tpr)) +
  geom_path(linewidth = 1.5, color = cols[2]) +
  annotate("label", x = 0.5, y = 0.5,
           label = paste0("AUC: ", round(empir_auc_CHIKV_VLP[2], 3)), size = 6) +
  labs(
    title = "CHIKV_VLP: Empirical ROC curve",
    x = "False positivity rate",
    y = "True positivity rate"
  ) +
  theme_bw() +
  theme(
    plot.title = element_text(hjust = 0.5, size = 22),
    axis.text  = element_text(size = 20),
    axis.title = element_text(size = 22),
    aspect.ratio = 1
  )

roc_plot <- plot_grid(plot_empirical_roc, plot_analytical_roc, ncol = 1)

# --- Analytical AUROC scatter ------------------------------------------------

scatterplot_analyt_auroc <-
  ggplot(aucs, aes(x = med_an, y = reorder(antigen, med_an), color = group)) +
  geom_errorbarh(aes(xmin = lb_an, xmax = ub_an), height = 0) +
  geom_point(size = 5) +
  facet_grid(rows = vars(group), scales = "free", space = "free") +
  scale_colour_brewer(palette = "Accent") +
  labs(
    title = "Analytical AUROC based on Gaussian FMM",
    x = "AUROC", y = NULL
  ) +
  theme_bw() +
  theme(
    strip.background = element_blank(),
    strip.text = element_text(size = 22),
    legend.position = "none",
    plot.title = element_text(hjust = 0.5, size = 22),
    axis.text  = element_text(size = 20),
    axis.title = element_text(size = 22)
  )

# --- Final Figure ------------------------------------------------------------

fig_2 <- plot_grid(roc_plot, scatterplot_analyt_auroc)

# --- Save to file ------------------------------------------------------------

ggsave(fig_2, file = "figures/Fig_2.png", width = 15, height = 10, scale = 1.2, bg = "white")
