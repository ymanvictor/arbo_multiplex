###############################################################################
# CHIKV_E2 Analytical ROC and AUROC (Gaussian Mixture Model, k = 3)
#
# This script computes analytical ROC curves, empirical ROC curves, and
# AUROC summaries for CHIKV_E2 antibody responses using a Gaussian 
# finite mixture model (GMM) with three components.
#
# ***Model structure***
#   - Component 1 (lowest mean)  → seronegative
#   - Component 2 (middle mean)  → seronegative
#   - Component 3 (highest mean) → seropositive
#
# Thus the NEGATIVE distribution is a mixture of components 1 and 2:
#
#     f_neg(x) = w1 * N(mu1, sigma1^2) + w2 * N(mu2, sigma2^2)
#     where w1 = theta1 / (theta1 + theta2)
#
# The POSITIVE distribution is the upper component:
#
#     f_pos(x) = N(mu3, sigma3^2)
#
# ***Closed-form AUROC***
# AUROC = P(X_pos > X_neg) under the above model. Because neg is a mixture:
#
#     AUROC = w1 * Φ((mu3 - mu1) / sqrt(sigma3^2 + sigma1^2)) +
#             w2 * Φ((mu3 - mu2) / sqrt(sigma3^2 + sigma2^2))
#
# This closed-form AUROC is used instead of trapezoidal numerical integration.
#
# Thresholds for fixed specificity (e.g. 97.7%) are obtained by numerically
# inverting the negative CDF:  P(X_neg > thr) = tail_prob
#
###############################################################################

# --- Libraries ---------------------------------------------------------------
library(tidyverse)
library(cowplot)
library(pROC)
library(ggrepel)
library(ggsci)

# --- Palette -----------------------------------------------------------------
cols <- pal_jama()(6)

# --- Data --------------------------------------------------------------------
post_stats <- read_csv("data/source_data_fig_2a.csv") %>%
  rename(param = 1)

post_med <- setNames(post_stats$post_med, post_stats$param)

ctrls <- read_csv("data/source_data_fig_1.csv") %>% 
  filter(sampletype != "U")

aucs <- read_csv("data/source_data_fig_2bc.csv")

# --- Helpers -----------------------------------------------------------------

CDF_neg_med <- function(x) {
  w1 <- post_med["theta[1]"]; w2 <- post_med["theta[2]"]
  m1 <- post_med["mu[1]"];    s1 <- post_med["sigma[1]"]
  m2 <- post_med["mu[2]"];    s2 <- post_med["sigma[2]"]
  w1t <- w1 / (w1 + w2)
  w2t <- 1 - w1t
  1 - (w1t * pnorm(x, m1, s1) + w2t * pnorm(x, m2, s2))
}

CDF_pos_med <- function(x) {
  1 - pnorm(x, post_med["mu[3]"], post_med["sigma[3]"])
}

# Closed-form AUROC for mid-negative scenario
auc_closed_form_midneg <- function(post_med) {
  w1 <- post_med["theta[1]"]; w2 <- post_med["theta[2]"]
  m1 <- post_med["mu[1]"];    m2 <- post_med["mu[2]"];    m3 <- post_med["mu[3]"]
  s1 <- post_med["sigma[1]"]; s2 <- post_med["sigma[2]"]; s3 <- post_med["sigma[3]"]
  
  w1t <- w1 / (w1 + w2)
  w2t <- 1 - w1t
  
  z31 <- (m3 - m1) / sqrt(s3^2 + s1^2)
  z32 <- (m3 - m2) / sqrt(s3^2 + s2^2)
  
  w1t * pnorm(z31) + w2t * pnorm(z32)
}

find_quantile <- function(tail_prob, lower = 0, upper = 20) {
  uniroot(
    function(x) CDF_neg_med(x) - tail_prob,
    interval = c(lower, upper)
  )$root
}
# --- Analytical ROC ----------------------------------------------------------

thresholds <- seq(0, 10, length.out = 2000)

analytical_roc_med <- tibble(
  threshold = thresholds,
  tpr = sapply(thresholds, CDF_pos_med),
  fpr = sapply(thresholds, CDF_neg_med)
) %>% arrange(fpr)

analyt_auc_CHIKV_E2 <- auc_closed_form_midneg(post_med)

# Threshold with 97.7% specificity
tail_977 <- 1 - 0.97725
serothr_977 <- find_quantile(tail_977)

roc_serothr_977 <- tibble(
  tpr = CDF_pos_med(serothr_977),
  fpr = CDF_neg_med(serothr_977)
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
  annotate(
    "label", x = 0.5, y = 0.5,
    label = paste0("AUC: ", round(analyt_auc_CHIKV_E2, 3)),
    size = 6
  ) +
  labs(
    title = "CHIKV_E2: Analytical ROC Curve",
    x = "False positivity rate",
    y = "True positivity rate"
  ) +
  theme_bw() +
  theme(
    plot.title   = element_text(hjust = 0.5, size = 22),
    axis.text    = element_text(size = 20),
    axis.title   = element_text(size = 22),
    aspect.ratio = 1
  )

# --- Empirical ROC -----------------------------------------------------------

roc_obj <- pROC::roc(sampletype ~ log_mfi, data = ctrls, ci = TRUE)

roc_empirical_df <- tibble(
  fpr = 1 - roc_obj$specificities,
  tpr = roc_obj$sensitivities,
  threshold = roc_obj$thresholds
) %>% arrange(fpr)

plot_empirical_roc <-
  ggplot(roc_empirical_df, aes(x = fpr, y = tpr)) +
  geom_path(linewidth = 1.5, color = cols[2]) +
  annotate(
    "label", x = 0.5, y = 0.5,
    label = paste0("AUC: ", round(roc_obj$ci[2], 3)),
    size = 6
  ) +
  labs(
    title = "CHIKV_E2: Empirical ROC curve",
    x = "False positivity rate",
    y = "True positivity rate"
  ) +
  theme_bw() +
  theme(
    plot.title   = element_text(hjust = 0.5, size = 22),
    axis.text    = element_text(size = 20),
    axis.title   = element_text(size = 22),
    aspect.ratio = 1
  )

roc_plot <- plot_grid(plot_empirical_roc, plot_analytical_roc, ncol = 1)

# --- Analytical AUROC scatter (from full dataset) ----------------------------

scatterplot_analyt_auroc <-
  ggplot(aucs, aes(x = med_an, y = reorder(antigen, med_an), color = group)) +
  geom_errorbarh(aes(xmin = lb_an, xmax = ub_an), height = 0) +
  geom_point(size = 5) +
  facet_grid(rows = vars(group), scales = "free", space = "free") +
  scale_colour_brewer(palette = "Accent") +
  labs(
    title = "Analytical AUROC based on Gaussian FMM",
    x = "AUROC",
    y = NULL
  ) +
  theme_bw() +
  theme(
    strip.background = element_blank(),
    strip.text       = element_text(size = 22),
    legend.position  = "none",
    plot.title       = element_text(hjust = 0.5, size = 22),
    axis.text        = element_text(size = 20),
    axis.title       = element_text(size = 22)
  )

# --- Final Figure ------------------------------------------------------------

fig_2 <- plot_grid(roc_plot, scatterplot_analyt_auroc)

ggsave(fig_2, file = "figures/Fig_2.png", width = 15, height = 10, scale = 1.2, bg = "white")
