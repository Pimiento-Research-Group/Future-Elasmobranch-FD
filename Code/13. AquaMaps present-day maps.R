###############################################################################################################
# 13. AquaMaps current maps
## This R code plots the global distribution of present-day functional diversity
## it produces Figure S8
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

# Current distribution - species richness and FD
load(file = "C:/Users/2022207/Dropbox/Jack's PhD/Chapter 3. Future shark FD/Data/R data/FD_AqMap_map_df.RData")

# Load current species richness map
load(file = "C:/Users/2022207/Dropbox/Jack's PhD/Chapter 3. Future shark FD/Data/R data/FD_AqMap_Species_richness_map.RData")

# Species richness map
Map1 <- gmap +
  geom_point(
    data = AqMap_Occ_map,
    aes(CenterLong, CenterLat, color = Spp_Richn),
    size = 0.01,
    alpha = 0.2,
    inherit.aes = FALSE
  ) +
  scale_color_viridis(discrete = FALSE) +
  labs(title = "Current species richness distribution", x = "Longitude", y = "Latitude") +
  theme_minimal() +
  theme_bw() +
  theme(
    panel.grid = element_blank(), 
    panel.border = element_rect(fill= "transparent"),
    legend.position = "right",
    plot.title = element_text(hjust = 0.5)
  )

# Rename columns for legend
final_map_df_AqMap <- final_map_df_AqMap %>%
  rename(FRic = fric, FSp = fspe)

# FD maps
Map2 <- gmap +
  geom_point(
    data = final_map_df_AqMap,
    aes(CenterLong, CenterLat, color = FRic),
    size = 0.01,
    alpha = 0.2,
    inherit.aes = FALSE
  ) +
  scale_color_viridis(discrete = FALSE) +
  labs(title = "Current functional richness distribution", x = "Longitude", y = "Latitude") +
  theme_minimal() +
  theme_bw() +
  theme(
    panel.grid = element_blank(), 
    panel.border = element_rect(fill= "transparent"),
    legend.position = "right",
    plot.title = element_text(hjust = 0.5)
  )

Map3 <- gmap +
  geom_point(
    data = final_map_df_AqMap,
    aes(CenterLong, CenterLat, color = FUn),
    size = 0.01,
    alpha = 0.2,
    inherit.aes = FALSE
  ) +
  scale_color_viridis(discrete = FALSE) +
  labs(title = "Current functional uniqueness distribution", x = "Longitude", y = "Latitude") +
  theme_minimal() +
  theme_bw() +
  theme(
    panel.grid = element_blank(), 
    panel.border = element_rect(fill= "transparent"),
    legend.position = "right",
    plot.title = element_text(hjust = 0.5)
  )

Map4 <- gmap +
  geom_point(
    data = final_map_df_AqMap,
    aes(CenterLong, CenterLat, color = FSp),
    size = 0.01,
    alpha = 0.2,
    inherit.aes = FALSE
  ) +
  scale_color_viridis(discrete = FALSE) +
  labs(title = "Current functional specialisation distribution", x = "Longitude", y = "Latitude") +
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
