###############################################################################################################
# 13. AquaMaps current maps
## This R code plots the global distribution of present-day functional diversity
## it produces Figure S5
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

# Current distribution - species richness and FD
load(file = "C:/Users/2022207/Dropbox/Jack's PhD/Chapter 3. Future shark FD/Analyses/Final paper R codes/data/AquaMaps/FD_AqMap_map_df.RData")

# Load current species richness map
load(file = "C:/Users/2022207/Dropbox/Jack's PhD/Chapter 3. Future shark FD/Analyses/Final paper R codes/data/AquaMaps/FD_AqMap_Species_richness_map.RData")

final_map_df_AqMap <- final_map_df_AqMap %>%
  drop_na()

# Species richness map
world <- map_data("world")

Map1 <- ggplot() +
  geom_point(
    data = AqMap_Occ_map,
    aes(CenterLong, CenterLat, color = Spp_Richn),
    size = 0.01,
    alpha = 0.3,
    inherit.aes = FALSE
  ) +
  geom_map(data = world, map = world,
           aes(long, lat, map_id = region),
           color = "black", fill = "lightgray", size = 0.2) +
  scale_color_viridis(discrete = FALSE, name = "SR") +
  labs(title = "", x = "Longitude", y = "Latitude") +
  scale_x_continuous(breaks = seq(-180, 180, by = 45)) +
  scale_y_continuous(breaks = seq(-90, 90, by = 45)) +  
  theme_minimal() +
  theme_bw() +
  theme(
    panel.grid = element_blank(), 
    panel.border = element_rect(fill = "transparent"),
    legend.position = "right",
    plot.title = element_text(hjust = 0.5)
  )

# FD Maps
Map2 <- ggplot() +
  geom_point(
    data = final_map_df_AqMap,
    aes(CenterLong, CenterLat, color = fric),
    size = 0.01,
    alpha = 0.3,
    inherit.aes = FALSE
  ) +
  geom_map(data = world, map = world,
           aes(long, lat, map_id = region),
           color = "black", fill = "lightgray", size = 0.2) +
  scale_color_viridis(discrete = FALSE, name = "FRic") +
  labs(title = "", x = "Longitude", y = "Latitude") +
  scale_x_continuous(breaks = seq(-180, 180, by = 45)) + 
  scale_y_continuous(breaks = seq(-90, 90, by = 45)) +  
  theme_minimal() +
  theme_bw() +
  theme(
    panel.grid = element_blank(), 
    panel.border = element_rect(fill= "transparent"),
    legend.position = "right",
    plot.title = element_text(hjust = 0.5)
  )

Map3 <- ggplot() +
  geom_point(
    data = final_map_df_AqMap,
    aes(CenterLong, CenterLat, color = fun),
    size = 0.01,
    alpha = 0.3,
    inherit.aes = FALSE
  ) +
  geom_map(data = world, map = world,
           aes(long, lat, map_id = region),
           color = "black", fill = "lightgray", size = 0.2) +
  scale_color_viridis(discrete = FALSE, name = "FUn") +
  labs(title = "", x = "Longitude", y = "Latitude") +
  scale_x_continuous(breaks = seq(-180, 180, by = 45)) + 
  scale_y_continuous(breaks = seq(-90, 90, by = 45)) +  
  theme_minimal() +
  theme_bw() +
  theme(
    panel.grid = element_blank(), 
    panel.border = element_rect(fill= "transparent"),
    legend.position = "right",
    plot.title = element_text(hjust = 0.5)
  )

Map4 <- ggplot() +
  geom_point(
    data = final_map_df_AqMap,
    aes(CenterLong, CenterLat, color = fsp),
    size = 0.01,
    alpha = 0.3,
    inherit.aes = FALSE
  ) +
  geom_map(data = world, map = world,
           aes(long, lat, map_id = region),
           color = "black", fill = "lightgray", size = 0.2) +
  scale_color_viridis(discrete = FALSE, name = "FSp") +
  labs(title = "", x = "Longitude", y = "Latitude") +
  scale_x_continuous(breaks = seq(-180, 180, by = 45)) +  
  scale_y_continuous(breaks = seq(-90, 90, by = 45)) +    
  theme_minimal() +
  theme_bw() +
  theme(
    panel.grid = element_blank(), 
    panel.border = element_rect(fill= "transparent"),
    legend.position = "right",
    plot.title = element_text(hjust = 0.5)
  )

# Plot 4 maps together - produces Figure S8
Aquamap <- plot_grid(Map1, Map2, Map3, Map4,
                     labels = c("(a)", "(b)", "(c)", "(d)"), 
                     label_size = 12, align = "hv", label_fontface = "bold", nrow = 2, ncol = 2)