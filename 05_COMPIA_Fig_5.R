# --- Libraries ----------------------------------------------------------------

library(drc)         # load first so tidyverse masks MASS::select if you reload it
library(tidyverse)
library(ggsci)
library(viridis)
library(cowplot)

# --- Data ---------------------------------------------------------------------

# Load competitive immunoassay results
chikv_compia <- readr::read_csv("raw_data/chikv_compia.csv")
mayv_compia  <- readr::read_csv("raw_data/mayv_compia.csv")

# set factor level etc.
prep_labels <- function(df) {
  df %>%
    dplyr::mutate(
      MAYV_CHIKV = factor(
        MAYV_CHIKV,
        levels = c("0 1", "1 0", "1 1"),
        labels = c("Suspected CHIKV", "Suspected MAYV", "Suspected MAYV & CHIKV")
      ),
      facet_label = paste0(MAYV_CHIKV, ": ", anon_id)
    )
}

chikv_compia <- prep_labels(chikv_compia)
mayv_compia  <- prep_labels(mayv_compia)

# --- Summary stats ------------------------------------------------------------

chikv_aggr <- chikv_compia %>%
  dplyr::filter(anon_id != "S0025") %>%
  dplyr::group_by(MAYV_CHIKV, comp_antigen, comp_conc) %>%
  dplyr::summarise(
    median_ratio = median(mfi_ratio),
    min_ratio = quantile(mfi_ratio, probs = 0.05),
    max_ratio = quantile(mfi_ratio, probs = 0.95),
    .groups = "drop"
  )

mayv_aggr <- mayv_compia %>%
  dplyr::filter(anon_id != "S0025") %>%
  dplyr::group_by(MAYV_CHIKV, comp_antigen, comp_conc) %>%
  dplyr::summarise(
    median_ratio = median(mfi_ratio),
    min_ratio = quantile(mfi_ratio, probs = 0.05),
    max_ratio = quantile(mfi_ratio, probs = 0.95),
    .groups = "drop"
  )

#--------------------------------------
# Fit 4-PL model to comp results
#--------------------------------------

# drop outlier 

chikv_fit <- chikv_compia
chikv_fit$mfi_ratio[
  which(chikv_fit$anon_id == "S0025" & chikv_fit$comp_antigen == "MAYV_E2")
] <- NA

mayv_fit <- mayv_compia
mayv_fit$mfi_ratio[
  which(mayv_fit$anon_id == "S0025" & mayv_fit$comp_antigen == "MAYV_E2")
] <- NA

# Prediction grid
newdata_grid <- tibble::tibble(
  comp_conc = seq(0.001, 10, length.out = 10000)
)

# Safe prediction function (unchanged behavior)

safe_predict_drm <- purrr::possibly(
  function(model) {
    preds <- stats::predict(model, newdata = data.frame(comp_conc = newdata_grid$comp_conc))
    tibble::tibble(
      comp_conc = newdata_grid$comp_conc,
      pred_mfi_ratio = as.vector(preds)
    )
  },
  otherwise = tibble::tibble(comp_conc = newdata_grid$comp_conc, pred_mfi_ratio = NA_real_)
)

# Fit per group (same as your original, just namespaced + explicit skip for BSA)

chikv_fit <- chikv_fit %>%
  dplyr::group_by(comp_antigen, MAYV_CHIKV) %>%
  tidyr::nest() %>%
  dplyr::filter(comp_antigen != "BSA") %>%  # ADDED HACK (same as original)
  dplyr::mutate(
    model = purrr::map(
      data,
      ~ drc::drm(mfi_ratio ~ (comp_conc), data = .x, fct = drc::LL.4())
    ),
    predictions = purrr::map(model, safe_predict_drm)
  ) %>%
  tidyr::unnest(predictions)

mayv_fit <- mayv_fit %>%
  dplyr::group_by(comp_antigen, MAYV_CHIKV) %>%
  tidyr::nest() %>%
  dplyr::filter(comp_antigen != "BSA") %>%  # ADDED HACK (same as original)
  dplyr::mutate(
    model = purrr::map(
      data,
      ~ drc::drm(mfi_ratio ~ (comp_conc), data = .x, fct = drc::LL.4())
    ),
    predictions = purrr::map(model, safe_predict_drm)
  ) %>%
  tidyr::unnest(predictions)

# --- Plots --------------------------------------------------------------------

plot_chikv_compia <- ggplot() +
  geom_line(data = chikv_fit, aes(x = comp_conc, y = pred_mfi_ratio, colour = comp_antigen), linewidth = 1.2) +
  geom_linerange(data = chikv_aggr, aes(x = comp_conc, ymin = min_ratio, ymax = max_ratio), colour = "black", alpha = 0.4) +
  geom_point(data = chikv_aggr, aes(x = comp_conc, y = median_ratio, fill = comp_antigen), size = 5, shape = 23, show.legend = FALSE) +
  geom_jitter(data = chikv_compia, aes(x = comp_conc, y = mfi_ratio, colour = comp_antigen), alpha = 0.5, width = 0.2) +
  facet_wrap(~ MAYV_CHIKV) +
  scale_x_log10(
    breaks = c(0.001, 0.01, 0.1, 1, 10),
    labels = c("0", "0.01", "0.1", "1", "10"),
    limits = c(0.001, 10)
  ) +
  scale_y_continuous(breaks = c(0, 0.25, 0.5, 0.75, 1.0, 1.25))+#limits = c(0, 1.8)) +
  labs(x = "Competitor Concentration", y = "Predicted MFI Ratio") +
  theme_bw() +
  ggtitle("Competitive Immunoassay: IgG reactivity towards CHIKV E2") +
  scale_colour_nejm() +
  scale_fill_nejm() +
  guides(colour = guide_legend(override.aes = list(size = 6))) +
  theme(
    plot.title = element_text(hjust = 0.5, size = 16),
    strip.background = element_blank(),
    strip.text = element_text(size = 13),
    legend.position = "bottom",
    legend.text = element_text(size = 14),
    legend.title = element_text(size = 15),
    axis.title = element_text(size = 13),
    axis.text = element_text(size = 14)
  ) +
  labs(
    colour = "Competing antigen",
    x = "Concentration of competing antigen (ug/mL)",
    y = "CHIKV E2: Proportion of max MFI"
  )

plot_mayv_compia <- ggplot() +
  geom_line(data = mayv_fit, aes(x = comp_conc, y = pred_mfi_ratio, colour = comp_antigen), linewidth = 1.2) +
  geom_linerange(data = mayv_aggr, aes(x = comp_conc, ymin = min_ratio, ymax = max_ratio), colour = "black", alpha = 0.4) +
  geom_point(data = mayv_aggr, aes(x = comp_conc, y = median_ratio, fill = comp_antigen), size = 5, shape = 23, show.legend = FALSE) +
  geom_jitter(data = mayv_compia, aes(x = comp_conc, y = mfi_ratio, colour = comp_antigen), alpha = 0.5, width = 0.2) +
  facet_wrap(~ MAYV_CHIKV) +
  scale_x_log10(
    breaks = c(0.001, 0.01, 0.1, 1, 10),
    labels = c("0", "0.01", "0.1", "1", "10"),
    limits = c(0.001, 10)
  ) +
  scale_y_continuous(breaks = c(0, 0.25, 0.5, 0.75, 1.0, 1.25))+#limits = c(0, 1.8)) +
  labs(x = "Competitor Concentration", y = "Predicted MFI Ratio") +
  theme_bw() +
  ggtitle("Competitive Immunoassay: IgG reactivity towards mayv E2") +
  scale_colour_nejm() +
  scale_fill_nejm() +
  guides(colour = guide_legend(override.aes = list(size = 6))) +
  theme(
    plot.title = element_text(hjust = 0.5, size = 16),
    strip.background = element_blank(),
    strip.text = element_text(size = 13),
    legend.position = "bottom",
    legend.text = element_text(size = 14),
    legend.title = element_text(size = 15),
    axis.title = element_text(size = 13),
    axis.text = element_text(size = 14)
  ) +
  labs(
    colour = "Competing antigen",
    x = "Concentration of competing antigen (ug/mL)",
    y = "mayv E2: Proportion of max MFI"
  )

# --- Final combined plot ------------------------------------------------------

plot_compia <- cowplot::plot_grid(plot_mayv_compia, plot_chikv_compia, nrow = 2)

# --- Save ---------------------------------------------------------------------
ggsave(
  filename = "figures/Fig_5.png",
  plot = plot_compia,
  bg = "white",
  width = 6, height = 6, scale = 1.5
)
