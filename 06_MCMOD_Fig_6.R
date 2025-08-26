# ---- packages ----
library(dplyr)
library(tidyr)
library(readr)
library(purrr)
library(ggplot2)
library(stringr)
library(scales)

# =========================
# 1) Load posterior CSV
# =========================
# CSV: param, stats ∈ {mean, lower95, upper95}, value
path_to_csv <- "output/MC_mod_params.csv"
post <- read_csv(path_to_csv, show_col_types = FALSE)

# =========================
# 2) Extract α (FOI) summaries using MEAN
# =========================
site_map <- tibble(
  code = c("Iquitos","Maroni","Dielmo","Ndiop","NC"),
  site = c("Iquitos (Peru)","Maroni (French Guiana)",
           "Dielmo (Senegal)","Ndiop (Senegal)","New Caledonia")
)

alphas <- post %>%
  filter(str_detect(param, "^alpha[MC]_")) %>%
  mutate(virus = if_else(str_starts(param, "alphaM_"), "MAYV", "CHIKV"),
         code  = str_replace(param, "^alpha[MC]_", "")) %>%
  left_join(site_map, by = "code") %>%
  select(site, virus, stats, value) %>%
  pivot_wider(names_from = stats, values_from = value) %>%
  transmute(site, virus,
            alpha_mid = mean,      # central tendency = MEAN
            alpha_lo  = lower95,   # 95% CrI lower
            alpha_hi  = upper95)   # 95% CrI upper

# =========================
# 3) FOI assumptions (edit if needed)
# =========================
# CHIKV: outbreak in Iquitos / Maroni / New Caledonia (T=2014, tau=0.5);
#        constant in Dielmo & Ndiop (Senegal).
# MAYV : constant in Iquitos & Maroni; absent in Senegal & New Caledonia.
foi_model <- tribble(
  ~site,                    ~virus,  ~model,
  "Iquitos (Peru)",         "CHIKV", "outbreak",
  "Maroni (French Guiana)", "CHIKV", "outbreak",
  "Dielmo (Senegal)",       "CHIKV", "constant",
  "Ndiop (Senegal)",        "CHIKV", "constant",
  "New Caledonia",          "CHIKV", "outbreak",
  "Iquitos (Peru)",         "MAYV",  "constant",
  "Maroni (French Guiana)", "MAYV",  "constant",
  "Dielmo (Senegal)",       "MAYV",  "absent",
  "Ndiop (Senegal)",        "MAYV",  "absent",
  "New Caledonia",          "MAYV",  "absent"
)

T_peak_tbl <- tribble(
  ~site,                    ~virus,  ~T_peak,
  "Iquitos (Peru)",         "CHIKV", 2014,
  "Maroni (French Guiana)", "CHIKV", 2014,
  "New Caledonia",          "CHIKV", 2014
)
tau_outbreak <- 0.5

# Survey years (two for Maroni)
survey_years <- tribble(
  ~site,                    ~survey_year,
  "Iquitos (Peru)",         2016,
  "Maroni (French Guiana)", 2015,
  "Maroni (French Guiana)", 2019,
  "Dielmo (Senegal)",       2015,
  "Ndiop (Senegal)",        2015,
  "New Caledonia",          2016
)

# =========================
# 4) Helpers
# =========================
gauss_w <- function(y, T, tau) exp(-((y - T)^2) / (tau^2))

# Share of total α accrued by years [S-age, S-1] for an outbreak with peak T and width tau
outbreak_share <- function(T, tau, age, S, window = 30L) {
  if (age <= 0) return(0)
  yrs_full <- (T - window):(T + window)
  w        <- gauss_w(yrs_full, T, tau)
  denom    <- sum(w)
  yrs_sub  <- (S - age):(S - 1)
  num      <- sum(w[match(yrs_sub, yrs_full, nomatch = 0)])
  num / denom
}

# =========================
# 5) Build grid & compute seroprevalence (mean and 95% CrI)
# =========================
age_max <- 80

grid <- survey_years %>%
  tidyr::expand_grid(virus = c("MAYV","CHIKV")) %>%
  left_join(foi_model, by = c("site","virus")) %>%
  left_join(alphas,    by = c("site","virus")) %>%
  mutate(
    alpha_mid = if_else(model == "absent", 0, alpha_mid),
    alpha_lo  = if_else(model == "absent", 0, alpha_lo),
    alpha_hi  = if_else(model == "absent", 0, alpha_hi)
  ) %>%
  left_join(T_peak_tbl, by = c("site","virus")) %>%
  tidyr::expand_grid(age = 0:age_max) %>%
  rowwise() %>%
  mutate(
    f_share = if (model == "outbreak")
      outbreak_share(T_peak, tau_outbreak, age, survey_year)
    else NA_real_,
    # cumulative hazard Λ(age)
    Lambda_mid = dplyr::case_when(
      model == "absent"   ~ 0,
      model == "constant" ~ alpha_mid * age,
      model == "outbreak" ~ alpha_mid * f_share
    ),
    Lambda_lo = dplyr::case_when(
      model == "absent"   ~ 0,
      model == "constant" ~ alpha_lo * age,
      model == "outbreak" ~ alpha_lo * f_share
    ),
    Lambda_hi = dplyr::case_when(
      model == "absent"   ~ 0,
      model == "constant" ~ alpha_hi * age,
      model == "outbreak" ~ alpha_hi * f_share
    ),
    seroprev_mid = 1 - exp(-Lambda_mid),
    seroprev_lo  = 1 - exp(-Lambda_lo),
    seroprev_hi  = 1 - exp(-Lambda_hi)
  ) %>%
  ungroup() %>%
  mutate(
    site = factor(site,
                  levels = c("Iquitos (Peru)",
                             "Dielmo (Senegal)",
                             "Ndiop (Senegal)",
                             "Maroni (French Guiana)",
                             "New Caledonia"))
  ) %>%
  mutate(
    grp = interaction(site, virus, survey_year, drop = TRUE),
    # short dashes ONLY for CHIKV in Maroni (2015); others solid
    lt  = if_else(site == "Maroni (French Guiana)" &
                    virus == "CHIKV" &
                    survey_year == 2015,
                  "22",  # short dashes pattern
                  "solid")
  )

plot_df <- grid %>% filter(model != "absent")

# =========================
# 6) Plot
# =========================
cols <- c("MAYV" = "purple", "CHIKV" = "orange")

MC_mod_plot <- ggplot() +
  # 95% CrI ribbons
  geom_ribbon(
    data = plot_df,
    aes(x = age, ymin = seroprev_lo, ymax = seroprev_hi,
        fill = virus, group = grp),
    alpha = 0.18, color = NA
  ) +
  # mean lines
  geom_path(
    data = plot_df,
    aes(x = age, y = seroprev_mid,
        color = virus, linetype = lt, group = grp),
    linewidth = 1
  ) +
  facet_wrap(~ site, ncol = 5) +
  scale_color_manual(values = cols) +
  scale_fill_manual(values  = cols) +
  scale_linetype_identity() +                 # use provided linetypes
  scale_y_continuous(labels = percent_format(accuracy = 1), limits = c(0, 1)) +
  labs(x = "Age (years)",
       y = "Seroprevalence",
       color = "Virus",
       title = "Age-specific seroprevalence by site (mean) with 95% CrI")+
  guides(fill = "none") +
  theme_bw(base_size = 12) +
  theme(panel.grid.minor = element_blank(),
        strip.background = element_blank(),
        plot.title = element_text(face = "bold"),
        legend.position = "bottom", 
        aspect.ratio = 1)

# --- Save to file ------------------------------------------------------------

ggsave(MC_mod_plot, file = "figures/Fig_6.png", bg = "white", width = 10, height = 5, scale = 1.2)
