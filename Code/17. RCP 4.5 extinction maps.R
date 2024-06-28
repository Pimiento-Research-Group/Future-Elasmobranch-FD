###############################################################################################################
# 17. RCP 4.5 extinction maps
## This R code plots the functional diversity shift under RCP 4.5 climate change
## while accounting for mode extinction times recorded by IUCN_sim
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

## Load present-day to calculate shifts
load(file = "C:/Users/2022207/Dropbox/Jack's PhD/Chapter 3. Future shark FD/Data/R data/FD_AqMap_map_df.RData")
load(file = "C:/Users/2022207/Dropbox/Jack's PhD/Chapter 3. Future shark FD/Data/R data/FD_AqMap_Species_richness_map.RData")

## Load RCP45 with extinctions
load(file = "C:/Users/2022207/Dropbox/Jack's PhD/Chapter 3. Future shark FD/Data/R data/FD_RCP45_extinctions_map_df.RData")
load(file = "C:/Users/2022207/Dropbox/Jack's PhD/Chapter 3. Future shark FD/Data/R data/FD_RCP45_extinctions_Species_richness_map.RData")

## Distribution maps
Map1 <- gmap +
  geom_point(
    data = RCP45_ext_Occ_map,
    aes(CenterLong, CenterLat, color = Spp_Richn),
    size = 0.01,
    alpha = 0.2,
    inherit.aes = FALSE
  ) +
  scale_color_viridis(discrete = FALSE) +
  labs(title = "RCP 4.5 extinctions species richness", x = "Longitude", y = "Latitude") +
  theme_minimal() +
  theme_bw() +
  theme(
    panel.grid = element_blank(), 
    panel.border = element_rect(fill= "transparent"),
    legend.position = "right",
    plot.title = element_text(hjust = 0.5)
  )

Map2 <- gmap +
  geom_point(
    data = final_map_df_RCP45_extinctions,
    aes(CenterLong, CenterLat, color = fric),
    size = 0.01,
    alpha = 0.2,
    inherit.aes = FALSE
  ) +
  scale_color_viridis(discrete = FALSE) +
  labs(title = "RCP 4.5 extinctions functional richness", x = "Longitude", y = "Latitude") +
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
    data = final_map_df_RCP45_extinctions,
    aes(CenterLong, CenterLat, color = FUn),
    size = 0.01,
    alpha = 0.2,
    inherit.aes = FALSE
  ) +
  scale_color_viridis(discrete = FALSE) +
  labs(title = "RCP 4.5 extinctions functional uniqueness", x = "Longitude", y = "Latitude") +
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
    data = final_map_df_RCP45_extinctions,
    aes(CenterLong, CenterLat, color = fspe),
    size = 0.01,
    alpha = 0.2,
    inherit.aes = FALSE
  ) +
  scale_color_viridis(discrete = FALSE) +
  labs(title = "RCP 4.5 extinctions functional specialisation", x = "Longitude", y = "Latitude") +
  theme_minimal() +
  theme_bw() +
  theme(
    panel.grid = element_blank(), 
    panel.border = element_rect(fill= "transparent"),
    legend.position = "right",
    plot.title = element_text(hjust = 0.5)
  )

# Plot 4 maps together
RCP45_extinctions_map <- plot_grid(Map1, Map2, Map3, Map4,
                                   labels = c("(a)", "(b)", "(c)", "(d)"), 
                                   label_size = 12, align = "hv", label_fontface = "bold", nrow = 2, ncol = 2)

# Join the dataframes based on common columns & calculate metric shifts
merged_data_RCP45ext_FD <- left_join(final_map_df_AqMap, final_map_df_RCP45_extinctions, by = c("Grid"))
merged_data_RCP45ext_FD <- merged_data_RCP45ext_FD %>%
  mutate(FRic_Diff = fric.y - fric.x) %>% 
  mutate(FUn_Diff = FUn.y - FUn.x) %>% 
  mutate(FSp_Diff = fspe.y - fspe.x)
merged_data_RCP45ext_FD <- merged_data_RCP45ext_FD %>%
  mutate(FRic_Percent_Change = (FRic_Diff/fric.x)*100) %>%
  mutate(FUn_Percent_Change = (FUn_Diff/FUn.x)*100) %>%
  mutate(FSp_Percent_Change = (FSp_Diff/fspe.x)*100)

merged_data_RCP45ext_SR <- left_join(AqMap_Occ_map, RCP45_ext_Occ_map, by = c("Grid"))
merged_data_RCP45ext_SR <- merged_data_RCP45ext_SR %>%
  mutate(
    CenterLat.y = coalesce(CenterLat.y, CenterLat.x),
    CenterLong.y = coalesce(CenterLong.y, CenterLong.x),
    Spp_Richn.y = coalesce(Spp_Richn.y, Spp_Richn.x)
  ) 

merged_data_RCP45ext_SR <- merged_data_RCP45ext_SR %>%
  mutate(Spp_Richn_Diff = Spp_Richn.y - Spp_Richn.x)
merged_data_RCP45ext_SR <- merged_data_RCP45ext_SR %>%
  mutate(SR_Percent_Change = (Spp_Richn_Diff/Spp_Richn.x)*100)
merged_data_RCP45ext_SR <- merged_data_RCP45ext_SR %>%
  mutate(
    Spp_Richn_Shift = case_when(
      Spp_Richn_Diff > 0 ~ "Increase",
      Spp_Richn_Diff == 0 ~ "No change",
      Spp_Richn_Diff < 0 ~ "Decrease",
      TRUE ~ NA_character_
    )
  )
merged_data_RCP45ext_SR$Spp_Richn_Shift <- ordered(merged_data_RCP45ext_SR$Spp_Richn_Shift,levels=c("Decrease","No change","Increase"))

# Generalise shifts
merged_data_RCP45ext_FD <- merged_data_RCP45ext_FD %>%
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
merged_data_RCP45ext_FD$FRic_Shift <- ordered(merged_data_RCP45ext_FD$FRic_Shift,levels=c("Decrease","Increase"))
merged_data_RCP45ext_FD$FUn_Shift <- ordered(merged_data_RCP45ext_FD$FUn_Shift,levels=c("Decrease","Increase"))
merged_data_RCP45ext_FD$FSp_Shift <- ordered(merged_data_RCP45ext_FD$FSp_Shift,levels=c("Decrease","Increase"))

## Absolute shift
RCP45_colors <- c("#e31a1c", "#ff7f00", "#ffffff", "#6a3d9a", "#1f78b4")
SR_breaks <- c(-69,-34,0,1,25,51)
FRic_breaks <- c(-0.46,-0.24,0,0.0000001,0.15,0.31)
FUn_breaks <- c(-0.0010766740, -0.000538337, 0, 0.00000001, 0.00077391, 0.001547820)
FSp_breaks <- c(-0.15922954, -0.07961477, 0, 0.00000001, 0.1189179, 0.2378358)

Map_RCP45ext_Ab_sprich <- gmap +
  geom_point(
    data = merged_data_RCP45ext_SR,
    aes(CenterLong.x, CenterLat.x, color = Spp_Richn_Diff),
    size = 0.01,
    alpha = 0.3,
    inherit.aes = FALSE
  ) +
  scale_color_gradientn(
    colors = RCP45_colors,
    values = scales::rescale(SR_breaks, to = c(0, 1)),
    limits = c(-69, 51),
    breaks = c(-69,-34,0,25,51),
    name = "Shift"
  ) +
  labs(title = "RCP 4.5 extinctions: Shift in Species Richness", x = "Longitude", y = "Latitude") +
  theme_minimal() +
  theme_bw() +
  theme(
    panel.grid = element_blank(),
    panel.border = element_rect(fill = "transparent"),
    legend.position = "right",
    plot.title = element_text(hjust = 0.5)
  )
save(Map_RCP45ext_Ab_sprich, file = "C:/Users/2022207/Dropbox/Jack's PhD/Chapter 3. Future shark FD/Data/R data/RCP45ext_Spp_Rich_absolute.RData")

Map_RCP45ext_Ab_fric <- gmap +
  geom_point(
    data = merged_data_RCP45ext_FD,
    aes(CenterLong.x, CenterLat.y, color = FRic_Diff),
    size = 0.01,
    alpha = 0.3,
    inherit.aes = FALSE
  ) +
  scale_color_gradientn(
    colors = RCP45_colors,
    values = scales::rescale(FRic_breaks, to = c(0, 1)),
    limits = c(-0.46, 0.31),
    breaks = c(-0.46,-0.24,0,0.15,0.31),
    labels = scales::label_number(accuracy = 0.01),
    name = "Shift"
  ) +
  labs(title = "RCP 4.5 extinctions: Shift in Functional Richness", x = "Longitude", y = "Latitude") +
  theme_minimal() +
  theme_bw() +
  theme(
    panel.grid = element_blank(),
    panel.border = element_rect(fill = "transparent"),
    legend.position = "right",
    plot.title = element_text(hjust = 0.5)
  )
save(Map_RCP45ext_Ab_fric, file = "C:/Users/2022207/Dropbox/Jack's PhD/Chapter 3. Future shark FD/Data/R data/RCP45ext_FRic_absolute.RData")

Map_RCP45ext_Ab_fun <- gmap +
  geom_point(
    data = merged_data_RCP45ext_FD,
    aes(CenterLong.x, CenterLat.x, color = FUn_Diff),
    size = 0.01,
    alpha = 0.3,
    inherit.aes = FALSE
  ) +
  scale_color_gradientn(
    colors = RCP45_colors,
    values = scales::rescale(FUn_breaks, to = c(0, 1)),
    limits = c(-0.0010766740, 0.001547820),
    breaks = c(-0.0010766740, -0.000538337, 0, 0.00077391, 0.001547820),
    labels = scales::label_scientific(),
    name = "Shift"
  ) +
  labs(title = "RCP 4.5 extinctions: Shift in Functional Uniqueness", x = "Longitude", y = "Latitude") +
  theme_minimal() +
  theme_bw() +
  theme(
    panel.grid = element_blank(),
    panel.border = element_rect(fill = "transparent"),
    legend.position = "right",
    plot.title = element_text(hjust = 0.5)
  )
save(Map_RCP45ext_Ab_fun, file = "C:/Users/2022207/Dropbox/Jack's PhD/Chapter 3. Future shark FD/Data/R data/RCP45ext_FUn_absolute.RData")

Map_RCP45ext_Ab_fsp <- gmap +
  geom_point(
    data = merged_data_RCP45ext_FD,
    aes(CenterLong.y, CenterLat.y, color = FSp_Diff),
    size = 0.01,
    alpha = 0.3,
    inherit.aes = FALSE
  ) +
  scale_color_gradientn(
    colors = RCP45_colors,
    values = scales::rescale(FSp_breaks, to = c(0, 1)),
    limits = c(-0.15922954, 0.2378358),
    breaks = c(-0.15922954, -0.07961477, 0, 0.1189179, 0.2378358),
    labels = scales::label_number(accuracy = 0.01),
    name = "Shift"
  ) +
  labs(title = "RCP 4.5 extinctions: Shift in Functional Specialisation", x = "Longitude", y = "Latitude") +
  theme_minimal() +
  theme_bw() +
  theme(
    panel.grid = element_blank(),
    panel.border = element_rect(fill = "transparent"),
    legend.position = "right",
    plot.title = element_text(hjust = 0.5)
  )
save(Map_RCP45ext_Ab_fsp, file = "C:/Users/2022207/Dropbox/Jack's PhD/Chapter 3. Future shark FD/Data/R data/RCP45ext_FSp_absolute.RData")

# Plot 4 maps together
RCP45_extinctions_map_changes <- plot_grid(Map_RCP45ext_Ab_sprich, Map_RCP45ext_Ab_fric, Map_RCP45ext_Ab_fun, Map_RCP45ext_Ab_fsp,
                                           labels = c("(a)", "(b)", "(c)", "(d)"), 
                                           label_size = 12, align = "hv", label_fontface = "bold", nrow = 2, ncol = 2)


# General shift maps
shift_colors <- c("Decrease" = "#F21A00", "No change" = "#EBCC2A", "Increase" = "#3B9AB2")
FD_shift_colors <- c("Decrease" = "#F21A00", "Increase" = "#3B9AB2")
Map_RCP45ext_sprich <- gmap +
  geom_point(
    data = merged_data_RCP45ext_SR,
    aes(CenterLong.x, CenterLat.x, color = Spp_Richn_Shift),
    size = 0.01,
    alpha = 0.2,
    inherit.aes = FALSE
  ) +
  scale_color_manual(values = shift_colors,
                     name = "Shift",
                     labels = c("Decrease", "No change", "Increase")) +
  labs(title = "RCP 4.5 extinctions: shift in species richness", x = "Longitude", y = "Latitude") +
  theme_minimal() +
  theme_bw() +
  theme(
    panel.grid = element_blank(), 
    panel.border = element_rect(fill= "transparent"),
    legend.position = "bottom",
    plot.title = element_text(hjust = 0.5),
    legend.key.size = unit(3, "mm")) +
  guides(color = guide_legend(
    override.aes = list(size = 5, alpha = 1, fill = shift_colors, color = shift_colors, 
                        shape = 21, stroke = 0.5))
  )
save(Map_RCP45ext_sprich, file = "C:/Users/2022207/Dropbox/Jack's PhD/Chapter 3. Future shark FD/Data/R data/RCP45ext_Spp_Rich_shift.RData")

Map_RCP45ext_fric <- gmap +
  geom_point(
    data = merged_data_RCP45ext_FD,
    aes(CenterLong.y, CenterLat.y, color = FRic_Shift),
    size = 0.01,
    alpha = 0.2,
    inherit.aes = FALSE
  ) +
  scale_color_manual(values = FD_shift_colors,
                     name = "Shift",
                     labels = c("Decrease", "Increase")) +
  labs(title = "RCP 4.5 extinctions: shift in functional richness", x = "Longitude", y = "Latitude") +
  theme_minimal() +
  theme_bw() +
  theme(
    panel.grid = element_blank(), 
    panel.border = element_rect(fill= "transparent"),
    legend.position = "bottom",
    plot.title = element_text(hjust = 0.5),
    legend.key.size = unit(3, "mm")) +
  guides(color = guide_legend(
    override.aes = list(size = 5, alpha = 1, fill = FD_shift_colors, color = FD_shift_colors, 
                        shape = 21, stroke = 0.5))
  )
save(Map_RCP45ext_fric, file = "C:/Users/2022207/Dropbox/Jack's PhD/Chapter 3. Future shark FD/Data/R data/RCP45ext_FRic_shift.RData")

Map_RCP45ext_fun <- gmap +
  geom_point(
    data = merged_data_RCP45ext_FD,
    aes(CenterLong.y, CenterLat.y, color = FUn_Shift),
    size = 0.01,
    alpha = 0.2,
    inherit.aes = FALSE
  ) +
  scale_color_manual(values = FD_shift_colors,
                     name = "Shift",
                     labels = c("Decrease", "Increase")) +
  labs(title = "RCP 4.5 extinctions: shift in functional uniqueness", x = "Longitude", y = "Latitude") +
  theme_minimal() +
  theme_bw() +
  theme(
    panel.grid = element_blank(), 
    panel.border = element_rect(fill= "transparent"),
    legend.position = "bottom",
    plot.title = element_text(hjust = 0.5),
    legend.key.size = unit(3, "mm")) +
  guides(color = guide_legend(
    override.aes = list(size = 5, alpha = 1, fill = FD_shift_colors, color = FD_shift_colors, 
                        shape = 21, stroke = 0.5))
  )
save(Map_RCP45ext_fun, file = "C:/Users/2022207/Dropbox/Jack's PhD/Chapter 3. Future shark FD/Data/R data/RCP45ext_FUn_shift.RData")

Map_RCP45ext_fsp <- gmap +
  geom_point(
    data = merged_data_RCP45ext_FD,
    aes(CenterLong.y, CenterLat.y, color = FSp_Shift),
    size = 0.01,
    alpha = 0.2,
    inherit.aes = FALSE
  ) +
  scale_color_manual(values = FD_shift_colors,
                     name = "Shift",
                     labels = c("Decrease", "Increase")) +
  labs(title = "RCP 4.5 extinctions: shift in functional specialisation", x = "Longitude", y = "Latitude") +
  theme_minimal() +
  theme_bw() +
  theme(
    panel.grid = element_blank(), 
    panel.border = element_rect(fill= "transparent"),
    legend.position = "bottom",
    plot.title = element_text(hjust = 0.5),
    legend.key.size = unit(3, "mm")) +
  guides(color = guide_legend(
    override.aes = list(size = 5, alpha = 1, fill = FD_shift_colors, color = FD_shift_colors, 
                        shape = 21, stroke = 0.5))
  )
save(Map_RCP45ext_fsp, file = "C:/Users/2022207/Dropbox/Jack's PhD/Chapter 3. Future shark FD/Data/R data/RCP45ext_FSp_shift.RData")

# Plot 4 maps together
RCP45_extinctions_map_shift <- plot_grid(Map_RCP45ext_sprich, Map_RCP45ext_fric, Map_RCP45ext_fun, Map_RCP45ext_fsp,
                                         labels = c("(a)", "(b)", "(c)", "(d)"), 
                                         label_size = 12, align = "hv", label_fontface = "bold", nrow = 2, ncol = 2)

## % changes maps
RCP4.5_colors <- c("#e31a1c", "#ff7f00","#ffa07a", "#ffffff","#b19cd9", "#6a3d9a", "#1f78b4")
SR_breaks <- c(-100, -67, -33, 0, 1, 100, 500, 900)
FRic_breaks <- c(-100, -67, -33, 0, 0.01, 100, 1000, 100000)
FUn_breaks <- c(-100, -67, -33, 0, 0.0000000001, 100, 1000, 60000)
FSp_breaks <- c(-40, -25, -10, 0, 0.0000000001, 33, 67, 100)

Map_RCP45ext_percent_sprich <- gmap +
  geom_point(
    data = merged_data_RCP45ext_SR,
    aes(CenterLong.x, CenterLat.x, color = SR_Percent_Change),
    size = 0.01,  
    alpha = 0.3,  
    inherit.aes = FALSE
  ) +
  scale_color_gradientn(
    colors = RCP4.5_colors, 
    values = scales::rescale(SR_breaks, to = c(0, 1)), 
    limits = c(-100, 700), 
    name = "Shift"  
  ) +
  labs(title = "RCP 4.5 extinctions: Shift in Species Richness", x = "Longitude", y = "Latitude") +
  theme_minimal() +
  theme_bw() +
  theme(
    panel.grid = element_blank(),
    panel.border = element_rect(fill = "transparent"),
    legend.position = "right",
    plot.title = element_text(hjust = 0.5)
  )
save(Map_RCP45ext_percent_sprich, file = "C:/Users/2022207/Dropbox/Jack's PhD/Chapter 3. Future shark FD/Data/R data/RCP45ext_Spp_Rich_percent.RData")

Map_RCP45ext_percent_fric <- gmap +
  geom_point(
    data = merged_data_RCP45ext_FD,
    aes(CenterLong.x, CenterLat.x, color = FRic_Percent_Change),
    size = 0.01, 
    alpha = 0.3,  
    inherit.aes = FALSE
  ) +
  scale_color_gradientn(
    colors = RCP4.5_colors,  
    values = scales::rescale(FRic_breaks, to = c(0, 1)),
    limits = c(-100, 100000), 
    breaks = c(-100, -67, -33, 0, 0.01, 100, 1000, 100000),
    labels = scales::label_number(accuracy = 0.01),
    name = "Shift"  
  ) +
  labs(title = "RCP 4.5 extinctions: Shift in Functional Richness", x = "Longitude", y = "Latitude") +
  theme_minimal() +
  theme_bw() +
  theme(
    panel.grid = element_blank(),
    panel.border = element_rect(fill = "transparent"),
    legend.position = "right",
    plot.title = element_text(hjust = 0.5)
  )
save(Map_RCP45ext_percent_fric, file = "C:/Users/2022207/Dropbox/Jack's PhD/Chapter 3. Future shark FD/Data/R data/RCP45ext_FRic_percent.RData")

Map_RCP45ext_percent_fun <- gmap +
  geom_point(
    data = merged_data_RCP45ext_FD,
    aes(CenterLong.x, CenterLat.x, color = FUn_Percent_Change),
    size = 0.01,  
    alpha = 0.3,  
    inherit.aes = FALSE
  ) +
  scale_color_gradientn(
    colors = RCP4.5_colors,  
    values = scales::rescale(FUn_breaks, to = c(0, 1)),  
    limits = c(-100, 60000), 
    breaks = c(-100, -67, -33, 0, 0.0000000001, 100, 1000, 60000),
    labels = scales::label_number(accuracy = 0.01),
    name = "Shift"  
  ) +
  labs(title = "RCP 4.5 extinctions: Shift in Functional Uniqueness", x = "Longitude", y = "Latitude") +
  theme_minimal() +
  theme_bw() +
  theme(
    panel.grid = element_blank(),
    panel.border = element_rect(fill = "transparent"),
    legend.position = "right",
    plot.title = element_text(hjust = 0.5)
  )
save(Map_RCP45ext_percent_fun, file = "C:/Users/2022207/Dropbox/Jack's PhD/Chapter 3. Future shark FD/Data/R data/RCP45ext_FUn_percent.RData")

Map_RCP45ext_percent_fsp <- gmap +
  geom_point(
    data = merged_data_RCP45ext_FD,
    aes(CenterLong.x, CenterLat.x, color = FSp_Percent_Change),
    size = 0.01,  
    alpha = 0.3,  
    inherit.aes = FALSE
  ) +
  scale_color_gradientn(
    colors = RCP4.5_colors,  
    values = scales::rescale(FSp_breaks, to = c(0, 1)),  
    limits = c(-40, 100), 
    breaks = c(-40, -25, -10, 0, 0.0000000001, 33, 67, 100),
    labels = scales::label_number(accuracy = 0.01),
    name = "Shift"  
  ) +
  labs(title = "RCP 4.5 extinctions: Shift in Functional Specialisation", x = "Longitude", y = "Latitude") +
  theme_minimal() +
  theme_bw() +
  theme(
    panel.grid = element_blank(),
    panel.border = element_rect(fill = "transparent"),
    legend.position = "right",
    plot.title = element_text(hjust = 0.5)
  )
save(Map_RCP45ext_percent_fsp, file = "C:/Users/2022207/Dropbox/Jack's PhD/Chapter 3. Future shark FD/Data/R data/RCP45ext_FSp_percent.RData")

# Plot 4 maps together
RCP45ext_percent_changes <- plot_grid(Map_RCP45ext_percent_sprich, Map_RCP45ext_percent_fric, 
                                   Map_RCP45ext_percent_fun, Map_RCP45ext_percent_fsp,
                                   labels = c("(a)", "(b)", "(c)", "(d)"), 
                                   label_size = 12, align = "hv", label_fontface = "bold", nrow = 2, ncol = 2)