# --- Libraries ---------------------------------------------------------------
library(tidyverse)
library(cmdstanr)
library(posterior)
library(bayesplot)
library(cowplot)

# --- Compile Stan model ------------------------------------------------------
mdl <- cmdstan_model(stan_file = "src/gmm_hybrid.stan")

# --- Data --------------------------------------------------------------------
data <- read_csv("raw_data/CHIKV_VLP.csv") %>%
  mutate(MFI = log_mfi)# - min(log_mfi) + 0.1)

antigen_name <- data$antigen %>% unique() %>% first()

set.seed(12)
MFI_train <- data %>% filter(sampletype == "U") %>% pull(MFI) %>% sample(1000)
MFI_neg   <- data %>% filter(sampletype == "N") %>% pull(MFI)
MFI_pos   <- data %>% filter(sampletype == "P") %>% pull(MFI)
MFI_un    <- data %>% filter(sampletype == "U") %>% pull(MFI)

MFI_max   <- max(data$MFI)
MFI_min   <- min(data$MFI)

stan_data <- list(
  N      = length(MFI_train),
  y      = MFI_train,
  N_neg  = length(MFI_neg),
  y_neg  = MFI_neg,
  y_max  = MFI_max,
  y_min  = MFI_min,
  K      = 3
)

# --- Fit (4 chains, 4000 warmup + 4000 sampling) ----------------------------

N_warmup <- 4000
N_iter   <- 4000

fit_k3 <- mdl$sample(
  data            = stan_data,
  iter_warmup     = N_warmup,
  iter_sampling   = N_iter,
  chains          = 4,
  parallel_chains = 4,
  seed            = 12
)

# --- Draws & posterior stats -------------------------------------------------

vars <- c("mu", "sigma", "theta")
draws_k3 <- fit_k3$draws(vars) |> as_draws_df()

params_k3_tbl <- draws_k3 %>%
  as_tibble() %>%
  select(all_of(grep(paste0("^(", paste(vars, collapse="|"), ")\\["), names(.), value = TRUE))) %>%
  pivot_longer(everything(), names_to = "param", values_to = "value") %>%
  group_by(param) %>%
  summarise(
    post_med = median(value),
    post_lb  = quantile(value, 0.025),
    post_ub  = quantile(value, 0.975),
    .groups  = "drop"
  )

post_med <- setNames(params_k3_tbl$post_med, params_k3_tbl$param)

# --- Convergence diagnostics -------------------------------------------------

summary_k3 <- fit_k3$summary(vars)
conv_diagn_k3 <- summary_k3 %>%
  select(variable, ess_bulk, ess_tail, rhat)

# --- Diagnostics plots -------------------------------------------------------

base_theme_small <- theme_bw() +
  theme(strip.background = element_blank(), axis.text.x = element_text(size = 5))

trace_plot <- mcmc_trace(fit_k3$draws(c(vars, "total_log_likelihood"))) +
  base_theme_small +
  ggtitle(paste0(antigen_name, ": Gaussian (k=3)"))

dens_overlay_plot <- mcmc_dens_overlay(fit_k3$draws(vars)) +
  base_theme_small +
  ggtitle(paste0(antigen_name, ": Gaussian (k=3)"))

# --- Mixture helpers ---------------------------------------------------------
mu_vec    <- unname(c(post_med["mu[1]"],    post_med["mu[2]"],    post_med["mu[3]"]))
sigma_vec <- unname(c(post_med["sigma[1]"], post_med["sigma[2]"], post_med["sigma[3]"]))
theta_vec <- unname(c(post_med["theta[1]"], post_med["theta[2]"], post_med["theta[3]"]))

# component i density; weighted=TRUE multiplies by theta[i], weighted=FALSE integrates to 1
comp_fun <- function(i, weighted = TRUE) {
  function(x) {
    w <- if (weighted) theta_vec[i] else 1
    dnorm(x, mean = mu_vec[i], sd = sigma_vec[i]) * w
  }
}

# full mixture density (weighted)
mix_fun <- function(x) {
  theta_vec[1]*dnorm(x, mu_vec[1], sigma_vec[1]) +
    theta_vec[2]*dnorm(x, mu_vec[2], sigma_vec[2]) +
    theta_vec[3]*dnorm(x, mu_vec[3], sigma_vec[3])
}

# closures for stat_function
comp1_unw_fun <- comp_fun(1, weighted = FALSE)  # unit-integral component 1
comp1_w_fun   <- comp_fun(1, weighted = TRUE)
comp2_w_fun   <- comp_fun(2, weighted = TRUE)
comp3_w_fun   <- comp_fun(3, weighted = TRUE)

# --- Plot panels -------------------------------------------------------------
p_theme <- theme_bw() + theme(plot.title = element_text(size = 8))

# Use small tibbles for ggplot inputs to keep aes() simple
df_neg <- tibble(x = MFI_neg)
df_pos <- tibble(x = MFI_pos)
df_un  <- tibble(x = MFI_un)

gmm_neg <- ggplot(df_neg, aes(x = x)) +
  geom_histogram(aes(y = after_stat(density)), bins = 30, fill = "blue", alpha = 0.5) +
  stat_function(fun = comp1_unw_fun, color = "blue", linewidth = 1) +   # <-- unit-integral
  p_theme +
  ggtitle(paste0(antigen_name, ": Gaussian (k=3) - Unexposed controls")) +
  xlab("log(MFI)") +
  scale_x_continuous(limits = c(2.5, MFI_max + 0.5))

gmm_pos <- ggplot(df_pos, aes(x = x)) +
  geom_histogram(aes(y = after_stat(density)), bins = 30, fill = "red", alpha = 0.5) +
  p_theme +
  ggtitle(paste0(antigen_name, ": Gaussian (k=3) - Exposed controls")) +
  xlab("log(MFI)") +
  scale_x_continuous(limits = c(2.5, MFI_max + 0.5))

gmm_un <- ggplot(df_un, aes(x = x)) +
  geom_histogram(aes(y = after_stat(density)), bins = 30, alpha = 0.5) +
  stat_function(fun = mix_fun,     color = "black",  linewidth = 1) +   # mixture (weighted)
  stat_function(fun = comp1_w_fun, color = "blue",   linewidth = 1) +
  stat_function(fun = comp2_w_fun, color = "purple", linewidth = 1) +
  stat_function(fun = comp3_w_fun, color = "red",    linewidth = 1) +
  p_theme +
  ggtitle(paste0(antigen_name, ": Gaussian (k=3) - Unknown samples")) +
  xlab("log(MFI)") +
  scale_x_continuous(limits = c(2.5, MFI_max + 0.5))

plot_gmm_k3 <- plot_grid(gmm_un, gmm_neg, gmm_pos, ncol = 1)

# --- Save figures ------------------------------------------------------------
ggsave(plot_gmm_k3, filename = "figures/Fig_1.png", width = 5, height = 7)

# Optionally save diagnostics
# ggsave(trace_plot,        filename = "figures/trace_k3.png", width = 12, height = 8)
# ggsave(dens_overlay_plot, filename = "figures/dens_overlay_k3.png", width = 12, height = 8)

# --- Outputs -----------------------------------------------------------------
write_csv(params_k3_tbl, "output/gmm_posterior_stats_CHIKV_VLP.csv")
write_csv(conv_diagn_k3, "output/gmm_ess_rhat_CHIKV_VLP.csv")
