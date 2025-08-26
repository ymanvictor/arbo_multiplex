# --- Libraries ---------------------------------------------------------------

library(tidyverse)
library(cowplot)

# --- Data --------------------------------------------------------------------

# Summary statistics to recreate box-plot
# (Site- and age-specific antibody levels for selected antigens:
# CHIKV VLP, CHIKV E2, ONNV VLP, MAYV E2, 
# DENV1 VLP, DENV NS1, ZIKVU NS1, YFV NS1)

box_stats <- read_csv(file = "output/site_age_Ab_levels.csv")

# order factor

box_stats <- box_stats %>%
  mutate(
    age_cat = factor(
      age_cat,
      levels = c("<5", "5-9", "10-17", "18-39", "40-59", ">59")
       
    )
  )

# titles

titles <- levels(factor(box_stats$title))

#--- Helpers ----------------------------------------------------------------

# 1) rebuild a single plot from stats (same styling, no hlines)
rebuild_one <- function(stats_boxes, antigen_title = "") {
  
  antigen_title <- unique(stats_boxes$title)
  if (length(antigen_title) > 1) {
    antigen_title <- antigen_title[1] # safeguard if duplicates
  }
  
  ggplot(stats_boxes, aes(x = age_cat, fill = study)) +
    geom_boxplot(
      aes(ymin = ymin, lower = lower, middle = middle, upper = upper, ymax = ymax,
          group = interaction(study, age_cat, drop = TRUE)),
      stat = "identity"
    ) +
    scale_fill_viridis_d(alpha = 0.7) +
    facet_wrap(~ study, ncol = 7) +
    theme_bw() +
    scale_y_log10(limits = c(10, 30000), breaks = c(100, 1000, 10000)) +
    labs(
      x = "Age-group (years)",
      y = "Normalised-MFI",
      title = antigen_title
    ) +
    theme(
      axis.text.x = element_text(angle = 45, hjust = 1, size = 8, margin = margin(t = 0, b = 0)),
      axis.text.y = element_text(angle = 90, size = 8, hjust = 0.5),
      axis.title  = element_text(margin = margin(b = 0, t = 0, r = 0, l = 0)),
      legend.position = "none",
      strip.background = element_blank(),
      strip.text = element_text(size = 10, margin = margin()),
      plot.title = element_text(hjust = 0.5, size = 12, margin = margin(b = 0, t = 0)),
      panel.spacing = unit(0.1, "lines"),
      aspect.ratio = 1
    )
}

# --- Generate each figure element ------

# 2) rebuild all and combine like before

rebuilt_list <- map(
  seq_along(1:8),
  ~{
    id <- paste0("plot_", .x)
    rebuild_one(
      stats_boxes   = filter(box_stats, plot_id == id),
      antigen_title = titles[.x]
    )
  }
)

# --- Final Figure ------------------------------------------------------------

fig_3 <- plot_grid(plotlist = rebuilt_list, ncol = 2)

# --- Save to file ------------------------------------------------------------

ggsave(fig_3, file = "figures/Fig_3.png", bg = "white", 
       width = 10, height = 8)
