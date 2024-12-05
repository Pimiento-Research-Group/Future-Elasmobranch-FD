###############################################################################################################
# 16. RCP 4.5 maps
## This R code plots the functional diversity shift under RCP 4.5 climate change
## It produces data necessary for all map plotting
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
library(wesanderson)
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

## Load present-day to calculate shifts
load(file = "C:/Users/2022207/Dropbox/Jack's PhD/Chapter 3. Future shark FD/Analyses/Final paper R codes/data/AquaMaps/FD_AqMap_Species_richness_map.RData")
load(file = "C:/Users/2022207/Dropbox/Jack's PhD/Chapter 3. Future shark FD/Analyses/Final paper R codes/data/AquaMaps/FD_AqMap_map_df.RData")

# Load RCP 4.5
load(file = "C:/Users/2022207/Dropbox/Jack's PhD/Chapter 3. Future shark FD/Analyses/Final paper R codes/data/AquaMaps/FD_RCP45_Species_richness_map.RData")
load(file = "C:/Users/2022207/Dropbox/Jack's PhD/Chapter 3. Future shark FD/Analyses/Final paper R codes/data/AquaMaps/FD_RCP45_map_df.RData")

# RCP 4.5 maps
MapRCP45.1 <- gmap +
  geom_point(
    data = RCP45_Occ_map,
    aes(CenterLong, CenterLat, color = Spp_Richn),
    size = 0.01,
    alpha = 0.2,
    inherit.aes = FALSE
  ) +
  scale_color_viridis(discrete = FALSE) +
  labs(title = "", x = "Longitude", y = "Latitude") +
  theme_minimal() +
  theme_bw() +
  theme(
    panel.grid = element_blank(), 
    panel.border = element_rect(fill= "transparent"),
    legend.position = "right",
    plot.title = element_text(hjust = 0.5)
  )

MapRCP45.2 <- gmap +
  geom_point(
    data = final_map_df_RCP45,
    aes(CenterLong, CenterLat, color = fric),
    size = 0.01,
    alpha = 0.2,
    inherit.aes = FALSE
  ) +
  scale_color_viridis(discrete = FALSE) +
  labs(title = "", x = "Longitude", y = "Latitude") +
  theme_minimal() +
  theme_bw() +
  theme(
    panel.grid = element_blank(), 
    panel.border = element_rect(fill= "transparent"),
    legend.position = "right",
    plot.title = element_text(hjust = 0.5)
  )

MapRCP45.3 <- gmap +
  geom_point(
    data = final_map_df_RCP45,
    aes(CenterLong, CenterLat, color = fun),
    size = 0.01,
    alpha = 0.2,
    inherit.aes = FALSE
  ) +
  scale_color_viridis(discrete = FALSE) +
  labs(title = "", x = "Longitude", y = "Latitude") +
  theme_minimal() +
  theme_bw() +
  theme(
    panel.grid = element_blank(), 
    panel.border = element_rect(fill= "transparent"),
    legend.position = "right",
    plot.title = element_text(hjust = 0.5)
  )

MapRCP45.4 <- gmap +
  geom_point(
    data = final_map_df_RCP45,
    aes(CenterLong, CenterLat, color = fsp),
    size = 0.01,
    alpha = 0.2,
    inherit.aes = FALSE
  ) +
  scale_color_viridis(discrete = FALSE) +
  labs(title = "", x = "Longitude", y = "Latitude") +
  theme_minimal() +
  theme_bw() +
  theme(
    panel.grid = element_blank(), 
    panel.border = element_rect(fill= "transparent"),
    legend.position = "right",
    plot.title = element_text(hjust = 0.5)
  )

# Plot 4 maps together
RCP45_map <- plot_grid(MapRCP45.1, MapRCP45.2, MapRCP45.3, MapRCP45.4,
                       labels = c("(a)", "(b)", "(c)", "(d)"), 
                       label_size = 12, align = "hv", label_fontface = "bold", nrow = 2, ncol = 2)

# Join the dataframes based on common columns & calculate shifts
merged_data_RCP45_FD <- left_join(final_map_df_AqMap, final_map_df_RCP45, by = c("Grid"))
merged_data_RCP45_FD <- merged_data_RCP45_FD %>%
  mutate(FRic_Diff = fric.y - fric.x) %>% 
  mutate(FUn_Diff = fun.y - fun.x) %>% 
  mutate(FSp_Diff = fsp.y - fsp.x)

merged_data_RCP45_FD <- merged_data_RCP45_FD %>%
  mutate(FRic_Percent_Change = (FRic_Diff/fric.x)*100) %>% 
  mutate(FUn_Percent_Change = (FUn_Diff/fun.x)*100) %>% 
  mutate(FSp_Percent_Change = (FSp_Diff/fsp.x)*100)

merged_data_RCP45_FD <- merged_data_RCP45_FD %>%
  mutate(
    FRic_Shift = case_when(
      FRic_Diff > 0 ~ "Increase",
      FRic_Diff == 0 ~ "No change",
      FRic_Diff < 0 ~ "Decrease",
      TRUE ~ NA_character_
    )) %>% 
  mutate(
    FUn_Shift = case_when(
      FUn_Diff > 0 ~ "Increase",
      FUn_Diff == 0 ~ "No change",
      FUn_Diff == 0.000000e+00 ~ "No change",
      FUn_Diff < 0 ~ "Decrease",
      TRUE ~ NA_character_
    )) %>%  
  mutate(
    FSp_Shift = case_when(
      FSp_Diff > 0 ~ "Increase",
      FSp_Diff == 0 ~ "No change",
      FSp_Diff < 0 ~ "Decrease",
      TRUE ~ NA_character_
    )
  ) %>% 
  na.omit()
merged_data_RCP45_FD$FRic_Shift <- ordered(merged_data_RCP45_FD$FRic_Shift,levels=c("Decrease", "Increase"))
merged_data_RCP45_FD$FUn_Shift <- ordered(merged_data_RCP45_FD$FUn_Shift,levels=c("Decrease", "Increase"))
merged_data_RCP45_FD$FSp_Shift <- ordered(merged_data_RCP45_FD$FSp_Shift,levels=c("Decrease", "Increase"))

# Calculate absolute shifts
merged_data_RCP45_SR <- left_join(AqMap_Occ_map, RCP45_Occ_map, by = c("Grid"))
merged_data_RCP45_SR <- merged_data_RCP45_SR %>%
  mutate(
    CenterLat.y = coalesce(CenterLat.y, CenterLat.x),
    CenterLong.y = coalesce(CenterLong.y, CenterLong.x),
    Spp_Richn.y = coalesce(Spp_Richn.y, Spp_Richn.x)
  )

merged_data_RCP45_SR <- merged_data_RCP45_SR %>%
  mutate(Spp_Richn_Diff = Spp_Richn.y - Spp_Richn.x)

merged_data_RCP45_SR <- merged_data_RCP45_SR %>%
  mutate(SR_Percent_Change = (Spp_Richn_Diff/Spp_Richn.x)*100)

merged_data_RCP45_SR <- merged_data_RCP45_SR %>%
  mutate(
    Spp_Richn_Shift = case_when(
      Spp_Richn_Diff > 0 ~ "Increase",
      Spp_Richn_Diff == 0 ~ "No change",
      Spp_Richn_Diff < 0 ~ "Decrease",
      TRUE ~ NA_character_
    )
  )
merged_data_RCP45_SR$Spp_Richn_Shift <- ordered(merged_data_RCP45_SR$Spp_Richn_Shift,levels=c("Decrease","No change","Increase"))

# Absolute shift maps
RCP45_colors <- c("#e31a1c", "#ff7f00", "#ffffff", "#6a3d9a", "#1f78b4")
SR_breaks <- c(-149,-75,0,1,37,74)
FRic_breaks <- c(-0.8,-0.4,0,0.0000001,0.4,0.8)
FUn_breaks <- c(-0.27, -0.13, 0, 0.0000000001, 0.16, 0.32)
FSp_breaks <- c(-0.28, -0.14, 0, 0.00000001, 0.12, 0.23)

Map_RCP45_Ab_sprich <- ggplot() +
  geom_point(
    data = merged_data_RCP45_SR,
    aes(CenterLong.x, CenterLat.x, color = Spp_Richn_Diff),
    size = 0.01,
    alpha = 0.5,
    inherit.aes = FALSE
  ) +
  geom_map(data = world, map = world,
           aes(long, lat, map_id = region),
           color = "black", fill = "lightgray", size = 0.2) +
  scale_color_gradientn(
    colors = RCP45_colors,
    values = scales::rescale(SR_breaks, to = c(0, 1)),
    limits = c(-149, 74),
    breaks = c(-149, -75, 0, 37, 74),
    name = expression(Delta ~ SR)
  ) +
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
save(Map_RCP45_Ab_sprich, file = "C:/Users/2022207/Dropbox/Jack's PhD/Chapter 3. Future shark FD/Analyses/Final paper R codes/data/AquaMaps/RCP45_Spp_Rich_absolute.RData")

Map_RCP45_Ab_fric <- ggplot() +
  geom_point(
    data = merged_data_RCP45_FD,
    aes(CenterLong.x, CenterLat.x, color = FRic_Diff),
    size = 0.01,
    alpha = 0.5,
    inherit.aes = FALSE
  ) +
  geom_map(data = world, map = world,
           aes(long, lat, map_id = region),
           color = "black", fill = "lightgray", size = 0.2) +
  scale_color_gradientn(
    colors = RCP45_colors,
    values = scales::rescale(FRic_breaks, to = c(0, 1)),
    limits = c(-0.8, 0.8),
    breaks = c(-0.8, -0.4, 0, 0.4, 0.8),
    labels = scales::label_number(accuracy = 0.01),
    name = expression(Delta ~ FRic)
  ) +
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
save(Map_RCP45_Ab_fric, file = "C:/Users/2022207/Dropbox/Jack's PhD/Chapter 3. Future shark FD/Analyses/Final paper R codes/data/AquaMaps/RCP45_FRic_absolute.RData")

Map_RCP45_Ab_fun <- ggplot() +
  geom_point(
    data = merged_data_RCP45_FD,
    aes(CenterLong.x, CenterLat.x, color = FUn_Diff),
    size = 0.01,
    alpha = 0.5,
    inherit.aes = FALSE
  ) +
  geom_map(data = world, map = world,
           aes(long, lat, map_id = region),
           color = "black", fill = "lightgray", size = 0.2) +
  scale_color_gradientn(
    colors = RCP45_colors,
    values = scales::rescale(FUn_breaks, to = c(0, 1)),
    limits = c(-0.27, 0.32),
    breaks = c(-0.27, -0.13, 0, 0.16, 0.32),
    labels = scales::label_number(accuracy = 0.01),
    name = expression(Delta ~ FUn)
  ) +
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
save(Map_RCP45_Ab_fun, file = "C:/Users/2022207/Dropbox/Jack's PhD/Chapter 3. Future shark FD/Analyses/Final paper R codes/data/AquaMaps/RCP45_FUn_absolute.RData")

Map_RCP45_Ab_fsp <- ggplot() +
  geom_point(
    data = merged_data_RCP45_FD,
    aes(CenterLong.x, CenterLat.x, color = FSp_Diff),
    size = 0.01,
    alpha = 0.5,
    inherit.aes = FALSE
  ) +
  geom_map(data = world, map = world,
           aes(long, lat, map_id = region),
           color = "black", fill = "lightgray", size = 0.2) +
  scale_color_gradientn(
    colors = RCP45_colors,
    values = scales::rescale(FSp_breaks, to = c(0, 1)),
    limits = c(-0.28, 0.23),
    breaks = c(-0.28, -0.14, 0, 0.00000001, 0.12, 0.23),
    labels = scales::label_number(accuracy = 0.01),
    name = expression(Delta ~ FSp)
  ) +
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
save(Map_RCP45_Ab_fsp, file = "C:/Users/2022207/Dropbox/Jack's PhD/Chapter 3. Future shark FD/Analyses/Final paper R codes/data/AquaMaps/RCP45_FSp_absolute.RData")

# Plot 4 maps together
RCP45_map_changes <- plot_grid(Map_RCP45_Ab_sprich, Map_RCP45_Ab_fric, Map_RCP45_Ab_fun, Map_RCP45_Ab_fsp,
                               labels = c("(a)", "(b)", "(c)", "(d)"), 
                               label_size = 12, align = "hv", label_fontface = "bold", nrow = 2, ncol = 2)
