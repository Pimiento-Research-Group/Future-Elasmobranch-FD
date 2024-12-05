###############################################################################################################
# 20. AquaMaps shift maps
## This R code plots all shift maps by metric
## it produces Figures 2, S6 and S7
###############################################################################################################

# Import packages
library(tidyverse)
library(reshape2)
library(ggsci)
library(tibble)
library(scales)
library(cowplot)
library(RColorBrewer)
library(maps)
library(viridis)
library(ggthemes)
library(wesanderson)

# Form map template
world <- map_data("world")
gmap <- ggplot() +
  geom_map(data = world, map = world,
           aes(long, lat, map_id = region),
           color = "black", fill = "lightgray", size = 0.2) +
  xlab("Longitude") +
  ylab("Latitude") +
  ggtitle("") +
  theme(plot.title = element_text(hjust = 0.5)) +
  theme_minimal() +
  theme_bw() +
  theme(panel.grid = element_blank(), 
        panel.border = element_rect(fill = "transparent"))

# Load all shift maps
## Absolute shifts
load(file = "C:/Users/2022207/Dropbox/Jack's PhD/Chapter 3. Future shark FD/Analyses/Final paper R codes/data/AquaMaps/RCP26_Spp_Rich_absolute.RData")
load(file = "C:/Users/2022207/Dropbox/Jack's PhD/Chapter 3. Future shark FD/Analyses/Final paper R codes/data/AquaMaps/RCP26_FRic_absolute.RData")
load(file = "C:/Users/2022207/Dropbox/Jack's PhD/Chapter 3. Future shark FD/Analyses/Final paper R codes/data/AquaMaps/RCP26_FUn_absolute.RData")
load(file = "C:/Users/2022207/Dropbox/Jack's PhD/Chapter 3. Future shark FD/Analyses/Final paper R codes/data/AquaMaps/RCP26_FSp_absolute.RData")
load(file = "C:/Users/2022207/Dropbox/Jack's PhD/Chapter 3. Future shark FD/Analyses/Final paper R codes/data/AquaMaps/RCP45_Spp_Rich_absolute.RData")
load(file = "C:/Users/2022207/Dropbox/Jack's PhD/Chapter 3. Future shark FD/Analyses/Final paper R codes/data/AquaMaps/RCP45_FRic_absolute.RData")
load(file = "C:/Users/2022207/Dropbox/Jack's PhD/Chapter 3. Future shark FD/Analyses/Final paper R codes/data/AquaMaps/RCP45_FUn_absolute.RData")
load(file = "C:/Users/2022207/Dropbox/Jack's PhD/Chapter 3. Future shark FD/Analyses/Final paper R codes/data/AquaMaps/RCP45_FSp_absolute.RData")
load(file = "C:/Users/2022207/Dropbox/Jack's PhD/Chapter 3. Future shark FD/Analyses/Final paper R codes/data/AquaMaps/RCP85_Spp_Rich_absolute.RData")
load(file = "C:/Users/2022207/Dropbox/Jack's PhD/Chapter 3. Future shark FD/Analyses/Final paper R codes/data/AquaMaps/RCP85_FRic_absolute.RData")
load(file = "C:/Users/2022207/Dropbox/Jack's PhD/Chapter 3. Future shark FD/Analyses/Final paper R codes/data/AquaMaps/RCP85_FUn_absolute.RData")
load(file = "C:/Users/2022207/Dropbox/Jack's PhD/Chapter 3. Future shark FD/Analyses/Final paper R codes/data/AquaMaps/RCP85_FSp_absolute.RData")
load(file = "C:/Users/2022207/Dropbox/Jack's PhD/Chapter 3. Future shark FD/Analyses/Final paper R codes/data/AquaMaps/RCP26ext_Spp_Rich_absolute.RData")
load(file = "C:/Users/2022207/Dropbox/Jack's PhD/Chapter 3. Future shark FD/Analyses/Final paper R codes/data/AquaMaps/RCP26ext_FRic_absolute.RData")
load(file = "C:/Users/2022207/Dropbox/Jack's PhD/Chapter 3. Future shark FD/Analyses/Final paper R codes/data/AquaMaps/RCP26ext_FUn_absolute.RData")
load(file = "C:/Users/2022207/Dropbox/Jack's PhD/Chapter 3. Future shark FD/Analyses/Final paper R codes/data/AquaMaps/RCP26ext_FSp_absolute.RData")
load(file = "C:/Users/2022207/Dropbox/Jack's PhD/Chapter 3. Future shark FD/Analyses/Final paper R codes/data/AquaMaps/RCP45ext_Spp_Rich_absolute.RData")
load(file = "C:/Users/2022207/Dropbox/Jack's PhD/Chapter 3. Future shark FD/Analyses/Final paper R codes/data/AquaMaps/RCP45ext_FRic_absolute.RData")
load(file = "C:/Users/2022207/Dropbox/Jack's PhD/Chapter 3. Future shark FD/Analyses/Final paper R codes/data/AquaMaps/RCP45ext_FUn_absolute.RData")
load(file = "C:/Users/2022207/Dropbox/Jack's PhD/Chapter 3. Future shark FD/Analyses/Final paper R codes/data/AquaMaps/RCP45ext_FSp_absolute.RData")
load(file = "C:/Users/2022207/Dropbox/Jack's PhD/Chapter 3. Future shark FD/Analyses/Final paper R codes/data/AquaMaps/RCP85ext_Spp_Rich_absolute.RData")
load(file = "C:/Users/2022207/Dropbox/Jack's PhD/Chapter 3. Future shark FD/Analyses/Final paper R codes/data/AquaMaps/RCP85ext_FRic_absolute.RData")
load(file = "C:/Users/2022207/Dropbox/Jack's PhD/Chapter 3. Future shark FD/Analyses/Final paper R codes/data/AquaMaps/RCP85ext_FUn_absolute.RData")
load(file = "C:/Users/2022207/Dropbox/Jack's PhD/Chapter 3. Future shark FD/Analyses/Final paper R codes/data/AquaMaps/RCP85ext_FSp_absolute.RData")

# Plot Figures (absolute changes)
## Figure 2 - RCP 4.5
Fig_2_map_RCP4.5 <- plot_grid(Map_RCP45ext_Ab_sprich,Map_RCP45ext_Ab_fric,
                       Map_RCP45ext_Ab_fun,Map_RCP45ext_Ab_fsp,
                       labels = c("(a)", "(b)", "(c)", "(d)"), 
                       label_size = 12, align = "hv", label_fontface = "bold", nrow = 4)

## Supplementary figure - climate change only
CC_map <- plot_grid(Map_RCP26_Ab_sprich, 
                       Map_RCP45_Ab_sprich,
                       Map_RCP85_Ab_sprich,
                       labels = c("(a)", "(b)", "(c)"), 
                       label_size = 12, align = "hv", label_fontface = "bold", nrow = 3)

## Supplementary figure - extinction species and FD maps (RCP 2.6, 8.5)
Fig_S7_map <- plot_grid(Map_RCP26ext_Ab_sprich,Map_RCP85ext_Ab_sprich,
                        Map_RCP26ext_Ab_fric,Map_RCP85ext_Ab_fric,
                        Map_RCP26ext_Ab_fun,Map_RCP85ext_Ab_fun,
                        Map_RCP26ext_Ab_fsp,Map_RCP85ext_Ab_fsp,
                        labels = c("(a)", "(b)", "(c)", "(d)", "(e)", "(f)", "(g)", "(h)"), 
                        label_size = 12, align = "hv", label_fontface = "bold", nrow = 4)
