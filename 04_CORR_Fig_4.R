# --- Libraries ---------------------------------------------------------------
library(tidyverse)
library(ggcorrplot)
library(cowplot)
library(viridis)

# --- Data --------------------------------------------------------------------

# load correlation matrices

cor_flavi <- read_csv("output/corr_matrix_flavi.csv")
cor_alpha <- read_csv("output/corr_matrix_alpha.csv")


# flavi

cor_flavi_mtrx <- as.matrix(cor_flavi[, 2:ncol(cor_flavi)])
rownames(cor_flavi_mtrx) <- cor_flavi$term

# alpha

cor_alpha_mtrx <- as.matrix(cor_alpha[, 2:ncol(cor_alpha)])
rownames(cor_alpha_mtrx) <- cor_alpha$term


#-------------------------------------------------------------------------------


# Figure 4

# Panel A

plot_flavi_corr <- ggcorrplot(cor_flavi_mtrx,
                                type = "upper", 
                                outline.color = "white")+
  
  theme(panel.grid.major = element_blank(), 
        legend.position = "right", 
        axis.text.x = element_text(angle = 45, vjust = 0, hjust = 0, size = 15) , 
        axis.text.y = element_text(size = 15)
  )+
  
  labs(fill= "Spearman rho")+
  scale_fill_viridis(option = "B", limits = c(-0.1,1))+
  scale_x_discrete(position = "top")




# panel B

plot_alpha_corr <- ggcorrplot(cor_alpha_mtrx,
                              type = "upper", 
                              outline.color = "white")+
  
  theme(panel.grid.major = element_blank(), 
        legend.position = "right", 
        axis.text.x = element_text(angle = 45, vjust = 0, hjust = 0, size = 15) , 
        axis.text.y = element_text(size = 15)
  )+
  
  labs(fill= "Spearman rho")+
  scale_fill_viridis(option = "B", limits = c(-0.1,1))+
  scale_x_discrete(position = "top")



# combine panels

fig_4 <- plot_grid(plot_alpha_corr, plot_flavi_corr, labels = c("A", "B"), label_size = 36)


ggsave(fig_4, file = "figures/Fig_4.png", width = 15, height = 8, bg = "white", scale = 1.2)
