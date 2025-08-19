library(tidyverse)
library(cmdstanr)
library(posterior)
library(bayesplot)
library(cowplot)

# source (and compile) stan models

mdl <- cmdstan_model(stan_file = "src/gmm_hybrid.stan")

# load data
  
data <- read_csv(file = "raw_data/CHIKV_VLP.csv")

antigen_name <- unique(data$antigen)

print(antigen_name)

MFI <- data$log_mfi

MFI_train <- sample(MFI[which(data$sampletype == "U")], 1000)

MFI_neg <- MFI[which(data$sampletype=="N")]

MFI_max <- max(MFI)

MFI_min <- min(MFI)

N_train <- length(MFI_train) # no. of observations

N_neg <- length(MFI_neg)


##############################################################################

# fit models (4 chains - sequential, 10000 iterations of which 2000 warmup)

N_warmup <-  4000
N_iter <- 4000


fit_k3 <- mdl$sample(data = list(N = N_train, ## num observations
                                 y = MFI_train,
                                 N_neg = N_neg,
                                 y_neg = MFI_neg,
                                 y_max = MFI_max,
                                 y_min = MFI_min,
                                 K = 3),
                     iter_warmup = N_warmup,
                     iter_sampling = N_iter,
                     chains = 4,
                     parallel_chains = 4,
                     seed = 12
                     
)


draws_k3 <- as_draws_df(fit_k3$draws(c("mu", "sigma", "theta"))) %>% #dplyr::filter(.chain!=3 & .chain != 4) |>
  data.frame()

# posterior stats

median_params_k3 <- sapply(draws_k3, median)[1:9]

lb_params_k3 <- sapply(draws_k3, quantile, 0.025)[1:9]

ub_params_k3 <- sapply(draws_k3, quantile, 0.975)[1:9]

params_k3 <- cbind(median_params_k3, lb_params_k3, ub_params_k3)

# model convergence diagnostics

summary_k3 <- fit_k3$summary(c("mu", "sigma", "theta"))

ess_k3 <- summary_k3[, c("variable", "ess_bulk", "ess_tail")]
rhat_k3 <- summary_k3[, c("variable", "rhat")]

conv_diagn_k3 <- ess_k3 %>% left_join(rhat_k3)

##########

# model diagnostic plots

# trace plots

mcmc_trace(fit_k3$draws(c("mu", "sigma", "theta", "total_log_likelihood")))+ # k = 3
  
  theme_bw()+
  
  theme(
    strip.background = element_blank(),
    axis.text.x = element_text(size = 5)
  )+
  
  ggtitle(paste0(antigen_name, ": Gaussian (k=3)"))

# posterior draws plot

mcmc_dens_overlay(fit_k3$draws(c("mu", "sigma", "theta")))+
  
  theme_bw()+
  
  theme(
    
    strip.background = element_blank(),
    axis.text.x = element_text(size = 5)
  )+
  
  ggtitle(paste0(antigen_name, ": Gaussian (k=3)"))



########################################

# Fig. 01


# gaussian k = 3

gmm_neg <- ggplot()+geom_histogram(aes(y = after_stat(density), x = MFI[which(data$sampletype == "N")]), fill = "blue", bins = 30, alpha = 0.5)+ # for hybrid model with negatives only fit components to MFI_train
  
  stat_function(fun = function(x) dnorm(x, mean = median_params_k3["mu.1."],
                                        sd = median_params_k3["sigma.1."])
                , color = "blue", , linewidth = 1)+
  
  theme_bw()+
  
  ggtitle(paste0(antigen_name, ": Gaussian (k=3) - Unexposed controls"))+
  
  xlab("log(MFI)")+
  
  theme(plot.title = element_text(size = 8))+
  
  scale_x_continuous(limits= c(0, MFI_max+0.5))


####################

# positive controls

gmm_pos <- ggplot()+geom_histogram(aes(y = after_stat(density), x = MFI[which(data$sampletype == "P")]),
                                           fill = "red", bins = 30, alpha = 0.5)+ # for hybrid model with negatives only fit components to MFI_train

  theme_bw()+
  
  ggtitle(paste0(antigen_name, ": Gaussian (k=3) - Exposed controls"))+
  
  xlab("log(MFI)")+
  
  theme(plot.title = element_text(size = 8))+
  
  scale_x_continuous(limits= c(0, MFI_max+0.5))

####################

# Unknown samples

gmm_un <- ggplot()+geom_histogram(aes(y = after_stat(density), x = MFI[which(data$sampletype == "U")]), bins = 30, alpha = 0.5)+ # for hybrid model with negatives only fit components to MFI_train
  
  stat_function(fun = function(x) dnorm(x, mean = median_params_k3["mu.1."],
                                        sd = median_params_k3["sigma.1."])*median_params_k3["theta.1."]+
                  dnorm(x, mean = median_params_k3["mu.2."],
                        sd = median_params_k3["sigma.2."])*median_params_k3["theta.2."]+
                  dnorm(x, mean = median_params_k3["mu.3."],
                        sd = median_params_k3["sigma.3."])*median_params_k3["theta.3."],
                color = "black", linewidth = 1)+
  
  stat_function(fun = function(x) dnorm(x, mean = median_params_k3["mu.1."],
                                        sd = median_params_k3["sigma.1."])*median_params_k3["theta.1."]
                , color = "blue", , linewidth = 1)+
  
  stat_function(fun = function(x) dnorm(x, mean = median_params_k3["mu.2."],
                                        sd = median_params_k3["sigma.2."])*median_params_k3["theta.2."]
                , color = "purple", , linewidth = 1)+
  
  stat_function(fun = function(x) dnorm(x, mean = median_params_k3["mu.3."],
                                        sd = median_params_k3["sigma.3."])*median_params_k3["theta.3."]
                , color = "red", linewidth = 1)+
  
  theme_bw()+
  
  ggtitle(paste0(antigen_name, ": Gaussian (k=3) - Unknown samples"))+
  
  xlab("log(MFI)")+
  
  theme(plot.title = element_text(size = 8))+
  
  scale_x_continuous(limits= c(0, MFI_max+0.5))

plot_gmm_k3 <- plot_grid(gmm_un, gmm_neg, gmm_pos,  ncol = 1)


########################################

ggsave(plot_gmm_k3, filename = "figures/Fig_1.png", width = 12.5, height = 7)


##############################################################################


# output

########################

# posterior statistics

write.csv(params_k3, file = "output/gmm_posterior_stats.csv")

# model diagnostics

write_csv(conv_diagn_k3, file = "output/gmm_ess_rhat.csv")


  




 
