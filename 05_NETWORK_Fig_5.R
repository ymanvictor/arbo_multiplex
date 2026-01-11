# ---- packages ----
library(dplyr)
library(tidyr)
library(ggplot2)
library(scales)
library(igraph)
library(ggraph)
library(ggrepel)
library(graphlayouts)
library(patchwork)

# ===================== 1) LOAD DERIVED INPUTS ========================

# 1) residual correlation matrix + antigen names
net_input <- readRDS("data/ns1_resid_cor.rds")
C_resid   <- net_input$C_resid
ant_cols  <- net_input$ant_cols
n_resid <- net_input$n_resid

# 2) PCA scores (PC1-3 + phenotype only)
pca_scores   <- readRDS("data/ns1_pca_scores_pc1_3.rds")

# 3) PCA loadings (long format: antigen, virus_group, PC, loading)
pca_loadings_long <- readRDS("data/ns1_pca_loadings_pc1_3.rds")

# 4) Phenotype proportions by site × age-group (saved with base::save)
phenotype_props_age_site <- readRDS("data/ns1_phenotype_props.rds")  # creates phenotype_props_age_site


# ===================== 2) NETWORK PLOT (PANEL A) =====================
virus_palette <- c(
  DENV="#7B3294", ZIKV="#1B9E77", YFV="#E69F00", WNV="#56B4E9",
  JEV ="#D55E00", USUV="#66A61E", CHIKV="#E7298A", MAYV="#4E79A7",
  ONNV="#A6D854", Other="grey60"
)

virus_order <- c("CHIKV","MAYV","ONNV","DENV","JEV","USUV","WNV","ZIKV","YFV")

gamma_ebic   <- 1
plot_thr     <- 0.1   # visual threshold: |rho| >= 0.1

# ---- EBICglasso on residual correlation matrix ----
P <- qgraph::EBICglasso(
  S     = C_resid,
  n     = n_resid,
  gamma = gamma_ebic
)

# Clean but DO NOT threshold for layout
P[!is.finite(P)] <- 0
diag(P)          <- 0
P[P < 0]         <- 0   # drop negative edges if you like

# ---- build igraph object from full EBICglasso network ----
g <- graph_from_adjacency_matrix(
  P,
  mode     = "undirected",
  weighted = TRUE,
  diag     = FALSE
)
E(g)$w <- E(g)$weight

# virus tags from antigen names
vnames <- V(g)$name
virus <- ifelse(grepl("^DENV", vnames),  "DENV",
                ifelse(grepl("^ZI[KV]", vnames), "ZIKV",
                       ifelse(grepl("^YFV", vnames),    "YFV",
                              ifelse(grepl("^WNV", vnames),    "WNV",
                                     ifelse(grepl("^JEV", vnames),    "JEV",
                                            ifelse(grepl("^USUV", vnames),   "USUV",
                                                   ifelse(grepl("^CHIKV", vnames),  "CHIKV",
                                                          ifelse(grepl("^MAYV", vnames),   "MAYV",
                                                                 ifelse(grepl("^ONNV", vnames),   "ONNV", "Other")))))))))
V(g)$virus <- factor(virus, levels = virus_order)

# node size by degree
rescale_safe <- function(x, to = c(7, 18)) {
  if (!length(x) || all(!is.finite(x))) return(rep(mean(to), length(x)))
  rng <- range(x[is.finite(x)])
  if (!is.finite(rng[1]) || rng[1] == rng[2]) return(rep(mean(to), length(x)))
  scales::rescale(x, to = to)
}
deg <- degree(g)
V(g)$size <- rescale_safe(deg, to = c(7, 18))

# ---- layout on FULL network (no 0.1 threshold here) ----
set.seed(12)

w_pos <- if (ecount(g)) pmax(1e-6, abs(E(g)$w)) else numeric(0)
xy    <- if (vcount(g)) layout_with_fr(g, weights = w_pos, niter = 4000) else matrix(0, 0, 2)

if (length(xy)) {
  V(g)$x <- xy[,1]
  V(g)$y <- xy[,2]
} else {
  V(g)$x <- numeric(0)
  V(g)$y <- numeric(0)
}

# ---- visual threshold: only draw edges with |w| >= 0.1 ----
E(g)$w_plot <- ifelse(abs(E(g)$w) >= plot_thr, E(g)$w, NA_real_)

p_network_A <-
  ggraph(g, layout = "manual", x = V(g)$x, y = V(g)$y) +
  geom_edge_link0(
    aes(edge_width = abs(w_plot)),
    colour       = "black",
    alpha        = 0.42,
    show.legend  = TRUE,
    na.rm        = TRUE   # drop edges with NA width from plotting
  ) +
  scale_edge_width(
    range  = c(0.25, 5),
    limits = c(0, 1),
    breaks = c(0.1, 0.2, 0.4, 0.6, 0.8, 1),
    name   = "Partial correlation (|ρ| ≥ 0.1)"
  ) +
  geom_node_point(
    aes(colour = virus, size = size),
    stroke      = 0.5,
    show.legend = c(colour = TRUE, size = FALSE)
  ) +
  scale_colour_manual(
    values = virus_palette,
    limits = virus_order,
    breaks = virus_order,
    name   = "Virus"
  ) +
  scale_size_identity() +
  geom_label_repel(
    aes(x = x, y = y, label = name),
    size          = 2.5,
    max.overlaps  = 200,
    point.padding = 0.75,
    box.padding   = 0.5,
    segment.alpha = 0.25
  ) +
  labs(
    title = expression(bold("A"))
  ) +
  theme_graph(base_family = "Helvetica") +
  theme(
    legend.position      = "bottom",
    legend.direction     = "horizontal",
    legend.box           = "vertical",
    legend.justification = "left",
    plot.title           = element_text(size = 20)
  )#+coord_flip()

# ===================== 3) PCA SCATTER (PANEL B) ======================

scores_full <- pca_scores

pca_cols <- c(
  "Unexposed-like (all low)"        = "#3300CC",
  "Other"                           = "grey80",
  "High DENV-only"                  = "#7B3294",
  "High ZIKV-only"                  = "#1B9E77",
  "High Non-DENV/Non-ZIKV-Flavi"    = "#E6AB02",
  "High DENV+ZIKV-only"             = "#3E7FBF",
  "Pan-Flavivirus-high"             = "#D450B4"
)

p_pc12 <- ggplot() +
  geom_point(
    data = subset(scores_full, pheno == "Other"),
    aes(PC1, PC2),
    color = "grey80", alpha = 0.8, size = 1
  ) +
  geom_point(
    data = subset(scores_full, pheno != "Other"),
    aes(PC1, PC2, color = pheno),
    alpha = 0.8, size = 1.7
  ) +
  scale_colour_manual(values = pca_cols, 
                      guide  = guide_legend(nrow = 1, byrow = TRUE))+
  theme_bw()+
  guides(colour = guide_legend(ncol=1))+
  theme(legend.position = "bottom", 
        panel.grid = element_blank(), 
        plot.title = element_text(size = 20))+
  labs(colour = "Sero-phenotype", 
       title = expression(bold("B")))
  
p_pc13 <- ggplot() +
  geom_point(
    data = subset(scores_full, pheno == "Other"),
    aes(PC1, PC3),
    color = "grey80", alpha = 0.8, size = 1
  ) +
  geom_point(
    data = subset(scores_full, pheno != "Other"),
    aes(PC1, PC3, color = pheno),
    alpha = 0.8, size = 1.7
  ) +
  scale_colour_manual(values = pca_cols, guide = "none") +
  guides(colour = guide_legend(ncol = 1))+
  theme_bw() +
  theme(
    legend.position = "none",
    panel.grid      = element_blank(),
    plot.title      = element_text(size = 14)
  )

p_pc23 <- ggplot() +
  geom_point(
    data = subset(scores_full, pheno == "Other"),
    aes(PC2, PC3),
    color = "grey80", alpha = 0.8, size = 1
  ) +
  geom_point(
    data = subset(scores_full, pheno != "Other"),
    aes(PC2, PC3, color = pheno),
    alpha = 0.8, size = 1.7
  ) +
  scale_colour_manual(values = pca_cols, guide = "none") +
  theme_bw() +
  theme(
    legend.position = "none",
    panel.grid      = element_blank(),
    plot.title      = element_text(size = 14)
  )

design_B <- "
AB
CD
"

p_scores_B <-
  p_pc12 + theme(legend.title.position = "top") +
  guide_area() +
  p_pc13 +
  p_pc23 +
  plot_layout(design = design_B, guides = "collect")


# ===================== 4) LOADINGS LOLLIPOP (PANEL C) =================



virus_cols <- c(
  DENV              = "#7B3294",
  ZIKV              = "#1B9E77",
  `Other Flavivirus` = "#E6AB02",
  Other             = "grey60"
)

p_loadings_C <- ggplot(pca_loadings_long,
                       aes(x = antigen, y = loading, colour = virus_group)) +
  geom_hline(yintercept = 0, linetype = 3, colour = "grey70") +
  geom_segment(aes(xend = antigen, y = 0, yend = loading), linewidth = 0.6) +
  geom_point(size = 2) +
  coord_flip() +
  facet_wrap(~ PC, ncol = 1, scales = "free_y") +
  scale_colour_manual(values = virus_cols) +
  theme_bw() +
  theme(
    strip.background     = element_blank(),
    legend.position      = "none",
    legend.title.position = "top",
    plot.title           = element_text(size = 20)
  ) +
  labs(
    x      = "Antigen",
    y      = "Loading",
    colour = "Virus group",
    title  = expression(bold("C"))
  )


# ===================== 5) PHENOTYPE PROPS (PANEL D) ===================

fill_vals <- c(
  "Other"                        = "grey90",
  "Unexposed-like (all low)"     = "#3300CC",
  "High DENV-only"               = "#7B3294",
  "High ZIKV-only"               = "#1B9E77",
  "High Non-DENV/Non-ZIKV-Flavi" = "#E6AB02",
  "High DENV+ZIKV-only"          = "#3E7FBF",
  "Pan-Flavivirus-high"          = "#D450B4"
)

p_phenotype_D <-
  phenotype_props_age_site %>%
  filter(!site %in% c("Alpha", "Flavi"), !is.na(age_cat)) %>%
  ggplot(aes(x = age_cat, y = prop, fill = pheno_bar)) +
  geom_bar(stat = "identity", position = "fill", alpha = 0.85) +
  theme_bw() +
  ylab("Proportion") +
  scale_fill_manual(values = fill_vals,
                    guide  = guide_legend(nrow = 1, byrow = TRUE)) +
  facet_wrap(~site, nrow = 1) +
  theme(
    axis.text.x      = element_text(angle = 90, vjust = 0.5, hjust = 1),
    legend.position  = "bottom",
    strip.background = element_blank(),
    plot.title       = element_text(size = 20)
  ) +
  labs(
    title = expression(bold("D")),
    x     = "Age-group (yrs)",
    fill  = "Sero-phenotype"
  )


# ===================== 6) COMBINE A–D INTO p_net_pca ==================

design2 <- "
A
B
"

p_net_pca <- ((((p_network_A+#labs(title = "A")+
                   # coord_flip()+
                   theme(legend.position= "bottom"))/guide_area())+plot_layout(heights = c(4,1), design = design2, guides="collect")) |((((p_scores_B | p_loadings_C)+plot_layout(widths = c(4,1)))  / 
                                                                                                                                           p_phenotype_D)+plot_layout(heights = c(4,1))))+
  plot_layout(widths = c(1, 2))

# Save
ggsave(
  p_net_pca,
  file   = "figures/p_net_pca.png",
  dpi    = 300,
  bg     = "white",
  units  = "cm",
  width  = 28,
  height = 21,
  scale  = 1.5
)

