# ============================================================
# 
# Rebuilds GMM plots from:
#   - output/histogram<antigen>histo_{neg,pos,un}.csv
#   - output/Table_S1.csv (GMM parameter estimates)
# ============================================================

library(tidyverse)
library(patchwork)
library(ggrepel)
library(viridis)

# ========================== Supplementary Fig S1. =============================

# -------------------------------
# Load parameter estimates
# -------------------------------

gmm_pars <- read_csv("data/source_data_fig_S1.csv")

# Expect columns: antigen, param, median, lb, ub
pars_nest <- gmm_pars %>%
  group_by(antigen) %>%
  nest(pars = c(param, median, lb, ub))

# -------------------------------
# Load histogram data
# -------------------------------

hist_files <- list.files(
  "data/histograms_fig_s1",
  pattern = "histo_(neg|pos|un)\\.csv$",
  full.names = TRUE
)

hist_tbl <- hist_files %>%
  set_names(basename(.) %>% str_remove("\\.csv$")) %>%
  map(read_csv)

# Names look like "DENV1 VLPhisto_neg", etc.

# -------------------------------
# Theme
# -------------------------------

p_theme <- theme_bw() +
  theme(
    plot.title  = element_text(size = 16),
    axis.title  = element_text(size = 14),
    axis.text   = element_text(size = 12),
    strip.text  = element_text(size = 14)
  )

# -------------------------------
# GMM function builder (2- or 3-comp)
# -------------------------------

make_gmm_funs <- function(pars, antigen_name = NA_character_) {
  
  get1 <- function(p) {
    v <- pars %>% filter(param == p) %>% pull(median)
    if (length(v) != 1L || is.na(v)) {
      stop(
        paste0("Missing parameter '", p,
               "' for antigen '", antigen_name, "'"),
        call. = FALSE
      )
    }
    v
  }
  
  k <- sum(grepl("^mu\\.", pars$param))
  
  if (k == 2) {
    mu1 <- get1("mu.1.")
    mu2 <- get1("mu.2.")
    
    s1  <- get1("sigma.1.")
    s2  <- get1("sigma.2.")
    
    th1 <- get1("theta.1.")
    th2 <- get1("theta.2.")
    
    w1 <- th1 / (th1 + th2)
    w2 <- th2 / (th1 + th2)
    
    comp1_unw_fun <- function(x) dnorm(x, mu1, s1)
    
    comp1_w_fun <- function(x) w1 * dnorm(x, mu1, s1)
    comp2_w_fun <- function(x) w2 * dnorm(x, mu2, s2)
    
    mix_fun <- function(x) comp1_w_fun(x) + comp2_w_fun(x)
    
    return(list(
      k = 2,
      comp1_unw_fun = comp1_unw_fun,
      comp1_w_fun   = comp1_w_fun,
      comp2_w_fun   = comp2_w_fun,
      mix_fun       = mix_fun
    ))
  }
  
  if (k == 3) {
    mu1 <- get1("mu.1.")
    mu2 <- get1("mu.2.")
    mu3 <- get1("mu.3.")
    
    s1  <- get1("sigma.1.")
    s2  <- get1("sigma.2.")
    s3  <- get1("sigma.3.")
    
    th1 <- get1("theta.1.")
    th2 <- get1("theta.2.")
    th3 <- get1("theta.3.")
    
    w1 <- th1 / (th1 + th2 + th3)
    w2 <- th2 / (th1 + th2 + th3)
    w3 <- th3 / (th1 + th2 + th3)
    
    comp1_unw_fun <- function(x) dnorm(x, mu1, s1)
    
    comp1_w_fun <- function(x) w1 * dnorm(x, mu1, s1)
    comp2_w_fun <- function(x) w2 * dnorm(x, mu2, s2)
    comp3_w_fun <- function(x) w3 * dnorm(x, mu3, s3)
    
    mix_fun <- function(x) comp1_w_fun(x) + comp2_w_fun(x) + comp3_w_fun(x)
    
    return(list(
      k = 3,
      comp1_unw_fun = comp1_unw_fun,
      comp1_w_fun   = comp1_w_fun,
      comp2_w_fun   = comp2_w_fun,
      comp3_w_fun   = comp3_w_fun,
      mix_fun       = mix_fun
    ))
  }
  
  stop("Unsupported k for antigen: ", antigen_name, call. = FALSE)
}

# -------------------------------
# Single-antigen plot constructor
# -------------------------------

make_gmm_plot <- function(antigen_name) {
  
  # names in hist_tbl are like "<antigen_name>histo_neg"
  key_neg <- paste0(antigen_name, "histo_neg")
  key_pos <- paste0(antigen_name, "histo_pos")
  key_un  <- paste0(antigen_name, "histo_un")
  
  h_neg <- hist_tbl[[key_neg]]
  h_pos <- hist_tbl[[key_pos]]
  h_un  <- hist_tbl[[key_un]]
  
  pars_row <- pars_nest %>% filter(antigen == antigen_name)
  if (nrow(pars_row) != 1L) {
    stop("No unique parameter row found for antigen: ", antigen_name)
  }
  funs <- make_gmm_funs(pars_row$pars[[1]], antigen_name)
  
  # Define x-limits from available hist data
  x_all <- c(h_neg$xmin, h_neg$xmax, h_pos$xmin, h_pos$xmax, h_un$xmin, h_un$xmax)
  x_all <- x_all[is.finite(x_all)]
  if (length(x_all) == 0L) {
    stop("No finite histogram x-range for antigen: ", antigen_name)
  }
  x_lim <- range(x_all)
  
  # Unknown samples: histogram + mixture + components
  gmm_un <- ggplot() +
    geom_rect(
      data = h_un,
      aes(xmin = xmin, xmax = xmax, ymin = 0, ymax = density),
      alpha = 0.5
    ) +
    stat_function(fun = funs$mix_fun,       color = "black", linewidth = 1) +
    stat_function(fun = funs$comp1_w_fun,   color = "blue",  linewidth = 1)
  
  # Special coloring for WNV EDIII component 2 (to mimic your original code)
  if (funs$k == 2 && antigen_name == "WNV EDIII") {
    gmm_un <- gmm_un +
      stat_function(fun = funs$comp2_w_fun, color = "red", linewidth = 1)
  } else {
    gmm_un <- gmm_un +
      stat_function(fun = funs$comp2_w_fun, color = "purple", linewidth = 1)
  }
  
  if (funs$k == 3) {
    gmm_un <- gmm_un +
      stat_function(fun = funs$comp3_w_fun, color = "red", linewidth = 1)
  }
  
  gmm_un <- gmm_un +
    p_theme +
    ggtitle(paste0(antigen_name, " - Unknown samples")) +
    xlab("log(MFI)") +
    ylab("Density") +
    scale_x_continuous(limits = x_lim)
  
  # Unexposed controls: histogram + negative component (unweighted)
  gmm_neg <- ggplot() +
    geom_rect(
      data = h_neg,
      aes(xmin = xmin, xmax = xmax, ymin = 0, ymax = density),
      fill = "blue", alpha = 0.5
    ) +
    stat_function(fun = funs$comp1_unw_fun, color = "blue", linewidth = 1) +
    p_theme +
    ggtitle(paste0(antigen_name, " - Unexposed controls")) +
    xlab("log(MFI)") +
    ylab("Density") +
    scale_x_continuous(limits = x_lim)
  
  # Exposed controls: histogram only
  gmm_pos <- ggplot() +
    geom_rect(
      data = h_pos,
      aes(xmin = xmin, xmax = xmax, ymin = 0, ymax = density),
      fill = "red", alpha = 0.5
    ) +
    p_theme +
    ggtitle(paste0(antigen_name, " - Exposed controls")) +
    xlab("log(MFI)") +
    ylab("Density") +
    scale_x_continuous(limits = x_lim)
  
  # Combine vertically as in your original
  gmm_full <- wrap_elements(
    (gmm_un / gmm_neg / gmm_pos) +
      plot_annotation(
        title = antigen_name,
        theme = theme(plot.title = element_text(
          hjust = 0.5, size = 20
        ))
      )
  )
  
  gmm_full
}

# -------------------------------
# Build full figure layout
# (you can customise the order/blocks to match p1–p4)
# -------------------------------

# Derive antigen list from parameter table or from histogram names
antigen_levels <- gmm_pars %>%
  distinct(antigen) %>%
  pull(antigen)

# If you want a specific order, define it here instead:
antigen_order <- factor(
  antigen_levels,
  levels = c(
    "DENV1 VLP", "DENV2 VLP", "DENV3 VLP", "DENV4 VLP",
    "DENV1 NS1", "DENV2 NS1", "DENV3 NS1", "DENV4 NS1",
    "DENV1 EDIII","DENV2 EDIII","DENV3 EDIII","DENV4 EDIII",
    "ZIKV VLP", "ZIKVSU NS1","ZIKVU NS1","ZIKVAS EDIII",
    "JEV E","JEV NS1","USUV NS1","WNV NS1",
    "YFV E","YFV NS1","WNV EDIII",
    "CHIKV VLP","CHIKV E2","ONNV VLP","MAYV E2","MAYV NSP1"
  )
)

antigen_order <- as.character(antigen_order[!is.na(antigen_order)])

plot_list <- lapply(antigen_order, make_gmm_plot)
names(plot_list) <- antigen_order

# Example layout mimicking your p1–p4
p2 <- wrap_plots(
  plot_list[c(
    "DENV1 VLP","DENV2 VLP","DENV3 VLP","DENV4 VLP",
    "DENV1 NS1","DENV2 NS1","DENV3 NS1","DENV4 NS1"
  )],
  ncol = 4
)

p3 <- wrap_plots(
  plot_list[c(
    "DENV1 EDIII","DENV2 EDIII","DENV3 EDIII","DENV4 EDIII",
    "ZIKV VLP","ZIKVSU NS1","ZIKVU NS1","ZIKVAS EDIII"
  )],
  ncol = 4
)

p4 <- wrap_plots(
  list(
    plot_list[["JEV E"]], plot_list[["JEV NS1"]],
    plot_list[["USUV NS1"]], plot_list[["WNV NS1"]],
    plot_list[["YFV E"]], plot_list[["YFV NS1"]],
    patchwork::plot_spacer(),
    plot_list[["WNV EDIII"]]
  ),
  ncol = 4
)

p1 <- wrap_plots(
  list(
    plot_list[["CHIKV VLP"]], plot_list[["CHIKV E2"]],
    plot_list[["ONNV VLP"]],  plot_list[["MAYV E2"]],
    patchwork::plot_spacer(), patchwork::plot_spacer(),
    patchwork::plot_spacer(),
    plot_list[["MAYV NSP1"]]
  ),
  ncol = 4
)

pdf("figures/Fig_S1.pdf",
    width = 21, height = 25, onefile = TRUE)
print(p1)
print(p2)
print(p3)
print(p4)
dev.off()

cat("✅ Figure reproduced from histogram data and GMM parameters\n")


# ========================== Supplementary Fig S2. =============================

aucs <- read_csv(file = "data/source_data_fig_s2.csv")

alpha <- c("CHIKV_VLP", "CHIKV_E2", "ONNV_VLP", "MAYV_E2", "MAYV_NSP1")

# gen virus var (alpha flavi)

aucs <- aucs %>% 
  mutate(
    virus = ifelse(antigen %in% alpha, "Alphavirus", "Flavivirus"),
  )

cor <- cor(aucs$med_em, aucs$med_an, method = "pearson")

# plot

fig_s2 <- aucs %>% 
  ggplot(aes(x = med_em, y = med_an, colour = virus))+
  geom_linerange(aes(xmin = lb_em, xmax = ub_em, colour = virus), alpha = 0.9)+
  geom_linerange(aes(ymin = lb_an, ymax = ub_an, colour = virus), alpha = 0.9)+
  annotate(geom = "label", label = paste0("Pearson r = ", round(cor, 3)), x = 0.7, y = 0.95, size = 5)+
  geom_point(size = 4)+
  theme_bw()+
  theme(
    legend.position = "bottom", 
    aspect.ratio = 1, 
    legend.title = element_blank(), 
    legend.text = element_text(size = 12)
  )+
  geom_label_repel(aes(label = antigen), max.overlaps = 22, size = 1.4, colour = "black")+
  xlab("Empirical AUROC with 95% confidence interval")+
  ylab("Analytical AUROC with 95% credible interval")+
  scale_colour_brewer(palette = "Accent")+
  scale_x_continuous(limits = c(0.5925,1.02))+
  scale_y_continuous(limits = c(0.6,1.02))

ggsave(fig_s2, file = "figures/Fig_S2.png", scale = 1.2)

# ========================== Supplementary Fig S3. =============================

# load ab data

hist_file <- "data/source_data_fig_s3_1.csv"
thr_file  <- "data/source_data_fig_s3_2.csv"

hist_df <- read_csv(hist_file, show_col_types = FALSE)
thr_df  <- read_csv(thr_file, show_col_types = FALSE)

# ---- compute stacking geometry: ymin/ymax for each site within each bin ----
# Choose a consistent stacking order (matches viridis legend order if desired)
site_levels <- sort(unique(hist_df$site))

hist_stacked <- hist_df %>%
  mutate(site = factor(site, levels = c("Peru", "Senegal", "French Guiana", "New Caledonia"))) %>%
  arrange(antigen, bin_id, desc(site)) %>%
  group_by(antigen, bin_id) %>%
  mutate(
    ymin = lag(cumsum(count), default = 0),
    ymax = cumsum(count)
  ) %>%
  ungroup()

x_rng <- range(c(hist_stacked$xmin, hist_stacked$xmax), na.rm = TRUE)

fig_s3 <- ggplot(hist_stacked, aes(fill = site)) +
  geom_rect(
    aes(xmin = xmin, xmax = xmax, ymin = ymin, ymax = ymax),
    alpha = 0.6,
    colour = NA
  ) +
  geom_vline(
    data = thr_df,
    aes(xintercept = seropos_threshold),
    linetype = "dashed",
    inherit.aes = FALSE
  ) +
  scale_fill_viridis_d() +
  facet_wrap(~ antigen, ncol = 7) +
  scale_x_continuous(limits = x_rng) +
  labs(x = "mfi", y = "count") +
  theme(
    legend.position = "bottom",
    legend.title = element_blank()
  )


ggsave(fig_s3, file = "figures/Fig_S3.png", scale = 2, width = 6)

# ========================== Supplementary Fig S4B. =============================

# load FMM estimated thresholds

thresholds <- read_csv(file = "data/source_data_fig_s4.csv") %>%
  pivot_wider(names_from = model, values_from = seropos_threshold)
  

# vector of alphavirus antigen names

alpha <- c("CHIKV VLP", "CHIKV E2", "ONNV VLP", "MAYV E2", "MAYV NSP1")

# gen virus (alpha vs flavi) var 

thresholds <- thresholds%>% 
  mutate(
    virus = ifelse(antigen %in% alpha, "Alphavirus", "Flavivirus")
  )

# estimate correlation

cor(thresholds$GMM, thresholds$`Weibull FMM`)

# plot (Supplementary Figure S4B)

fig_s4 <- thresholds %>%  
  ggplot(aes(x = GMM, y = `Weibull FMM`, color = virus)) + 
  geom_point(size = 6)+
  theme_bw()+
  theme(
    aspect.ratio = 1,
    legend.position = "bottom", 
    legend.title = element_blank(),
    legend.text = element_text(size = 12)
  )+
  coord_fixed()+
  geom_label_repel(aes(label = antigen), max.overlaps = 20, size = 2, colour = "black")+
  xlab("Threshold of seropositivity obtained from Guassian model")+
  ylab("Threshold of seropositivity obtained from Weibull model")+
  annotate(geom = "label", label = "Pearson r = 0.94", x = 5.5, y = 8, size = 4)+
  scale_colour_brewer(palette = "Accent")+
  scale_x_continuous(limits = c(4.5, 8.7))+
  scale_y_continuous(limits = c(4.5, 8.6))

ggsave(fig_s4, file = "figures/Fig_S4B.png", scale = 1.5)


# ========================== Supplementary Fig S4B. =============================

ab_box_data_sex <- read_csv(file = "data/source_data_fig_s5.csv") %>% 
  mutate(
    age_cat = factor(age_cat, levels = c("<5",  "5-9", "10-17", "18-39", "40-59", ">59"))
    )

fig_s5 <- 
  ab_box_data_sex %>% 
  ggplot(
  aes(
    x = factor(age_cat),
    ymin = ymin, lower = lower, middle = middle, upper = upper, ymax = ymax,
    fill = male
  )
) +
  geom_boxplot(stat = "identity", position = position_dodge(width = 0.8)) +
  labs(x = "Age group", y = "MFI", fill = "Sex") +
  theme_bw()+
  theme(axis.text.x = element_text(angle = 45, hjust = 1, vjust = 1))+
  facet_wrap(~site+antigen, ncol = 14)

ggsave(fig_s5, file = "figures/Fig_S5.png", scale = 4)


# ========================== Supplementary Fig S6. =============================

ab_box_data<- read_csv(file = "data/source_data_fig_s6.csv") %>% 
  mutate(
    age_cat = factor(age_cat, levels = c("<5",  "5-9", "10-17", "18-39", "40-59", ">59"))
  )

fig_s6 <- 
  ab_box_data %>% 
  ggplot(
    aes(
      x = factor(age_cat),
      ymin = ymin, lower = lower, middle = middle, upper = upper, ymax = ymax,
      fill = site
    )
  ) +
  geom_boxplot(stat = "identity", position = position_dodge(width = 0.8)) +
  labs(x = "Site", y = "Value", fill = "Age group") +
  theme_bw()+
  theme(axis.text.x = element_text(angle = 45, hjust = 1, vjust = 1))+
  facet_wrap(~antigen+site, ncol = 16)+
  scale_fill_viridis_d()

ggsave(fig_s6, file = "figures/Fig_S6.png", scale = 4)


# ========================== Supplementary Fig S7. =============================

seroprev_data <- read_csv(file = "data/source_data_fig_s7.csv") %>% 
  mutate(
    age_cat = factor(age_cat, levels = c("<5",  "5-9", "10-17", "18-39", "40-59", ">59"))
  )

fig_s7 <- seroprev_data %>% 
  ggplot(aes(y = adjusted_prev, x = factor(age_cat), colour = site_year, group = site_year))+
  geom_line(linewidth = 0.5, colour = "black")+
  geom_point(size = 3)+
  facet_wrap(~antigen, ncol = 7)+
  scale_color_viridis_d(alpha = 0.95)+
  scale_y_continuous(limits = c(0,1))+
  theme_bw()+
  ylab("Proportion seropositive")+
  xlab("Age-group (years)")+
  guides(colour = guide_legend(override.aes = list(size=10)))+
  theme(
    axis.text.x = element_text(angle = 45, hjust = 1, size = 8, margin = margin(t = 0, b = 0)),
    axis.text.y = element_text(angle = 0, size = 10),
    axis.title = element_text(margin = margin(b = 0, t = 0, r = 0, l = 0), size = 8),
    legend.position = "bottom", 
    legend.title = element_blank(),
    legend.text = element_text(size = 15),
    strip.background = element_blank(),
    strip.text = element_text(size = 10, margin = margin()),
    plot.title = element_text(hjust = 0.5,  size = 12, margin=margin(b = 0, t = 0)), 
    panel.spacing = unit(0.1, "lines")
  )


ggsave(fig_s7, file = "figures/Fig_S7.png", scale = 2)
