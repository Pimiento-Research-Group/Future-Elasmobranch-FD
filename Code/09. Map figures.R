###############################################################################################################
# 9. Map figures
## This R code plots all maps of species and functional richness
## Also includes latitudinal gradients 
## it produces the map components of Figure 1, 2, S4
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

# Load world
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

# Load map dataframes
load(file="C:/Users/2022207/Dropbox/Jack's PhD/Chapter 3. Future shark FD/Data/R data/Maps/FRic_map_present.RData")
load(file="C:/Users/2022207/Dropbox/Jack's PhD/Chapter 3. Future shark FD/Data/R data/Maps/FRic_map_RCP45.RData")
load(file="C:/Users/2022207/Dropbox/Jack's PhD/Chapter 3. Future shark FD/Data/R data/Maps/FRic_map_RCP45ext.RData")
load(file="C:/Users/2022207/Dropbox/Jack's PhD/Chapter 3. Future shark FD/Data/R data/Maps/SR_map_present.RData")
load(file="C:/Users/2022207/Dropbox/Jack's PhD/Chapter 3. Future shark FD/Data/R data/Maps/SR_map_RCP45.RData")
load(file="C:/Users/2022207/Dropbox/Jack's PhD/Chapter 3. Future shark FD/Data/R data/Maps/SR_map_RCP45ext.RData")

# Plot raw value maps (SR and FRic)
## Present day
Map_SR <- ggplot() +
  geom_point(
    data = SR_AqMap_map,
    aes(x = CenterLong, y = CenterLat, color = species_richness),
    size = 0.01,
    alpha = 0.3,
    inherit.aes = FALSE
  ) +
  geom_map(
    data = world, map = world,
    aes(x = long, y = lat, map_id = region),
    color = "black", fill = "lightgray", size = 0.2
  ) +
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

Map_FRic <- ggplot() +
  geom_point(
    data = FD_AqMap,
    aes(x = CenterLong, y = CenterLat, color = fric),
    size = 0.01,
    alpha = 0.3,
    inherit.aes = FALSE
  ) +
  geom_map(
    data = world, map = world,
    aes(x = long, y = lat, map_id = region),
    color = "black", fill = "lightgray", size = 0.2
  ) +
  scale_color_viridis(discrete = FALSE, name = "FRic") +
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

## Raw maps - Figures S4
Fig_S4 <- plot_grid(Map_SR, Map_FRic,
                    labels = c("A", "B"), 
                    label_size = 12, align = "hv", label_fontface = "bold", nrow = 2)

## RCP 4.5 (climate change)
Map_SR_RCP45 <- ggplot() +
  geom_point(
    data = SR_RCP45_map,
    aes(x = CenterLong, y = CenterLat, color = species_richness),
    size = 0.01,
    alpha = 0.3,
    inherit.aes = FALSE
  ) +
  geom_map(
    data = world, map = world,
    aes(x = long, y = lat, map_id = region),
    color = "black", fill = "lightgray", size = 0.2
  ) +
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

Map_FRic_RCP45 <- ggplot() +
  geom_point(
    data = FD_RCP45,
    aes(x = CenterLong, y = CenterLat, color = fric),
    size = 0.01,
    alpha = 0.3,
    inherit.aes = FALSE
  ) +
  geom_map(
    data = world, map = world,
    aes(x = long, y = lat, map_id = region),
    color = "black", fill = "lightgray", size = 0.2
  ) +
  scale_color_viridis(discrete = FALSE, name = "FRic") +
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

## Climate change + extinctions
Map_SR_RCP45ext <- ggplot() +
  geom_point(
    data = SR_RCP45ext_map,
    aes(x = CenterLong, y = CenterLat, color = species_richness),
    size = 0.01,
    alpha = 0.3,
    inherit.aes = FALSE
  ) +
  geom_map(
    data = world, map = world,
    aes(x = long, y = lat, map_id = region),
    color = "black", fill = "lightgray", size = 0.2
  ) +
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

Map_FRic_RCP45ext <- ggplot() +
  geom_point(
    data = FD_RCP45ext,
    aes(x = CenterLong, y = CenterLat, color = fric),
    size = 0.01,
    alpha = 0.3,
    inherit.aes = FALSE
  ) +
  geom_map(
    data = world, map = world,
    aes(x = long, y = lat, map_id = region),
    color = "black", fill = "lightgray", size = 0.2
  ) +
  scale_color_viridis(discrete = FALSE, name = "FRic") +
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

# Calculate shifts under RCP 4.5 (climate change)
FD_shift_RCP45 <- FD_AqMap %>%
  inner_join(FD_RCP45, by = "Grid", suffix = c("_present", "_rcp45")) %>%
  mutate(
    SR_shift = sp_richn_rcp45 - sp_richn_present,
    FRic_shift = fric_rcp45 - fric_present
  ) %>%
  select(Grid, CenterLat_present, CenterLong_present, SR_shift, FRic_shift) %>%
  rename(CenterLat = CenterLat_present, CenterLong = CenterLong_present) %>% 
  na.omit()
save(FD_shift_RCP45,file="C:/Users/2022207/Dropbox/Jack's PhD/Chapter 3. Future shark FD/Data/R data/Maps/RCP45_shifts.RData")

# Calculate shifts under RCP 4.5 (climate change) + extinctions
FD_shift_RCP45ext <- FD_AqMap %>%
  inner_join(FD_RCP45ext, by = "Grid", suffix = c("_present", "_rcp45ext")) %>%
  mutate(
    SR_shift = sp_richn_rcp45ext - sp_richn_present,
    FRic_shift = fric_rcp45ext - fric_present
  ) %>%
  select(Grid, CenterLat_present, CenterLong_present, SR_shift, FRic_shift) %>%
  rename(CenterLat = CenterLat_present, CenterLong = CenterLong_present) %>% 
  na.omit()
save(FD_shift_RCP45ext,file="C:/Users/2022207/Dropbox/Jack's PhD/Chapter 3. Future shark FD/Data/R data/Maps/RCP45ext_shifts.RData")

# Calculate grid changes under both scenarios
sr_shift_RCP45 <- FD_shift_RCP45 %>%
  mutate(SR_change = case_when(
    SR_shift > 0 ~ "Increase",
    SR_shift < 0 ~ "Decrease",
    TRUE ~ "No change"
  )) %>%
  count(SR_change) %>%
  mutate(Proportion = n / sum(n))

fric_shift_RCP45 <- FD_shift_RCP45 %>%
  mutate(FRic_change = case_when(
    FRic_shift > 0 ~ "Increase",
    FRic_shift < 0 ~ "Decrease",
    TRUE ~ "No change"
  )) %>%
  count(FRic_change) %>%
  mutate(Proportion = n / sum(n))

sr_shift_RCP45ext <- FD_shift_RCP45ext %>%
  mutate(SR_change = case_when(
    SR_shift > 0 ~ "Increase",
    SR_shift < 0 ~ "Decrease",
    TRUE ~ "No change"
  )) %>%
  count(SR_change) %>%
  mutate(Proportion = n / sum(n))

fric_shift_RCP45ext <- FD_shift_RCP45ext %>%
  mutate(FRic_change = case_when(
    FRic_shift > 0 ~ "Increase",
    FRic_shift < 0 ~ "Decrease",
    TRUE ~ "No change"
  )) %>%
  count(FRic_change) %>%
  mutate(Proportion = n / sum(n))

# Test to see if FRic shift is significantly lower following extinctions
FRic_compare <- FD_shift_RCP45 %>%
  inner_join(FD_shift_RCP45ext, by = "Grid", suffix = c("_rcp45", "_rcp45ext"))

wilcox.test(FRic_compare$FRic_shift_rcp45, FRic_compare$FRic_shift_rcp45ext, 
            paired = TRUE, alternative = "greater")

# Calculate latitudinal gradient of SR and FRic shifts (part of Figure 1)
## Calculate median, 95% and 50% confidence intervals - species richness & FD separately
Lat_RCP45_FD <- FD_shift_RCP45 %>%
  mutate(CenterLat = as.numeric(CenterLat)) %>% 
  group_by(CenterLat) %>%  
  summarise(
    median_value_SR = median(SR_shift, na.rm = TRUE),
    lower_50CI_SR = quantile(SR_shift, 0.25, na.rm = TRUE),
    upper_50CI_SR = quantile(SR_shift, 0.75, na.rm = TRUE),
    median_value_FRic = median(FRic_shift, na.rm = TRUE),
    lower_50CI_FRic = quantile(FRic_shift, 0.25, na.rm = TRUE),
    upper_50CI_FRic = quantile(FRic_shift, 0.75, na.rm = TRUE),
 ) %>%
  arrange(CenterLat)

Lat_RCP45ext_FD <- FD_shift_RCP45ext %>%
  mutate(CenterLat = as.numeric(CenterLat)) %>% 
  group_by(CenterLat) %>%  
  summarise(
    median_value_SR = median(SR_shift, na.rm = TRUE),
    lower_50CI_SR = quantile(SR_shift, 0.25, na.rm = TRUE),
    upper_50CI_SR = quantile(SR_shift, 0.75, na.rm = TRUE),
    median_value_FRic = median(FRic_shift, na.rm = TRUE),
    lower_50CI_FRic = quantile(FRic_shift, 0.25, na.rm = TRUE),
    upper_50CI_FRic = quantile(FRic_shift, 0.75, na.rm = TRUE),
  ) %>%
  arrange(CenterLat)

Lat_RCP45_FD$Scenario <- "RCP 4.5"
Lat_RCP45ext_FD$Scenario <- "RCP 4.5 + Extinctions"
combined_RCP45_FD <- rbind(Lat_RCP45_FD, Lat_RCP45ext_FD)
save(combined_RCP45_FD,file="C:/Users/2022207/Dropbox/Jack's PhD/Chapter 3. Future shark FD/Data/R data/Maps/RCP45_latitudinal_gradients.RData")

## Plot SR changes
SR_RCP45 <- ggplot(combined_RCP45_FD, aes(x = CenterLat, y = median_value_SR, color = Scenario, fill = Scenario)) +
  geom_line(size = 1) + 
  geom_ribbon(aes(ymin = lower_50CI_SR, ymax = upper_50CI_SR), alpha = 0.2) + 
  labs(title = "SR latitudinal gradient", x = "Latitude", y = "SR shift") + 
  coord_flip() +  # Rotate the plot
  scale_x_continuous(limits = c(-60, 80), breaks = seq(-60, 80, by = 20)) + 
  scale_y_continuous(limits = c(-15, 15), breaks = seq(-15, 15, by = 5)) + 
  geom_hline(yintercept = 0, linetype = "dashed", color = "grey", size = 0.5) + 
  theme_minimal() +
  theme_bw() +
  theme(legend.position = "bottom", 
        axis.text = element_text(size = 8, color = "black"),
        axis.title = element_text(size = 10),
        panel.background = element_rect(fill = "white"),
        panel.grid = element_blank(),
        panel.border = element_rect(fill = "transparent"),
        plot.title = element_text(hjust = 0.5, size = 10)) + 
  scale_color_manual(
    name = NULL,
    values = c("RCP 4.5" = "#0072B5FF", "RCP 4.5 + Extinctions" = "#BC3C29FF"),
    labels = c("RCP 4.5" = "Climate change", "RCP 4.5 + Extinctions" = "Climate change + extinctions")
  ) +
  scale_fill_manual(
    name = NULL,
    values = c("RCP 4.5" = "#0072B5FF", "RCP 4.5 + Extinctions" = "#BC3C29FF"),
    labels = c("RCP 4.5" = "Climate change", "RCP 4.5 + Extinctions" = "Climate change + extinctions")
  )+
  guides(color = guide_legend(nrow = 2), fill = guide_legend(nrow = 2))
save(SR_RCP45,file="C:/Users/2022207/Dropbox/Jack's PhD/Chapter 3. Future shark FD/Data/R data/Maps/RCP45_latitudinal_gradients_SR.RData")

## Plot FRic changes
FRic_RCP45 <- ggplot(combined_RCP45_FD, aes(x = CenterLat, y = median_value_FRic, color = Scenario, fill = Scenario)) +
  geom_line(size = 1) + 
  geom_ribbon(aes(ymin = lower_50CI_FRic, ymax = upper_50CI_FRic), alpha = 0.2) + 
  labs(title = "FRic latitudinal gradient", x = "Latitude", y = "FRic shift") + 
  coord_flip() +  # Rotate the plot
  scale_x_continuous(limits = c(-60, 80), breaks = seq(-60, 80, by = 20)) + 
  scale_y_continuous(limits = c(-0.1, 0.2), breaks = seq(-0.1, 0.2, by = 0.1), labels = label_number(accuracy = 0.01)) + 
  geom_hline(yintercept = 0, linetype = "dashed", color = "grey", size = 0.5) + 
  theme_minimal() +
  theme_bw() +
  theme(legend.position = "bottom", 
        axis.text = element_text(size = 8, color = "black"),
        axis.title = element_text(size = 10),
        panel.background = element_rect(fill = "white"),
        panel.grid = element_blank(),
        panel.border = element_rect(fill = "transparent"),
        plot.title = element_text(hjust = 0.5, size = 10)) + 
  scale_color_manual(
    name = NULL,
    values = c("RCP 4.5" = "#0072B5FF", "RCP 4.5 + Extinctions" = "#BC3C29FF"),
    labels = c("RCP 4.5" = "Climate change", "RCP 4.5 + Extinctions" = "Climate change + extinctions")
  ) +
  scale_fill_manual(
    name = NULL,
    values = c("RCP 4.5" = "#0072B5FF", "RCP 4.5 + Extinctions" = "#BC3C29FF"),
    labels = c("RCP 4.5" = "Climate change", "RCP 4.5 + Extinctions" = "Climate change + extinctions")
  )+
  guides(color = guide_legend(nrow = 2), fill = guide_legend(nrow = 2))
save(FRic_RCP45,file="C:/Users/2022207/Dropbox/Jack's PhD/Chapter 3. Future shark FD/Data/R data/Maps/RCP45_latitudinal_gradients_FRic.RData")

# Plot and save all maps so far
## Set breaks and colours for plots
RCP45_colors <- c("#e31a1c", "#ff7f00", "#ffffff", "#6a3d9a", "#1f78b4")
SR_breaks <- c(-86,-43,0,1,37,74)
FRic_breaks <- c(-0.5,-0.25,0,0.0000001,0.16,0.33)

## Shift maps - Figure 1
load(file = "C:/Users/2022207/Dropbox/Jack's PhD/Chapter 3. Future shark FD/Data/R data/Spp_2100.RData")
load(file = "C:/Users/2022207/Dropbox/Jack's PhD/Chapter 3. Future shark FD/Data/R data/FRic_2100.RData")

## Climate change shifts
RCP45_SR_shift <- ggplot() +
  geom_point(
    data = FD_shift_RCP45,
    aes(CenterLong, CenterLat, color = SR_shift),
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
    limits = c(-86, 74),
    breaks = c(-86,-43,0,37,74),
    name = expression(Delta ~ SR)
  ) +
  labs(title = "SR (climate change)", x = "Longitude", y = "Latitude") +
  scale_x_continuous(breaks = seq(-180, 180, by = 45)) +
  scale_y_continuous(breaks = seq(-90, 90, by = 45)) +
  theme_minimal() +
  theme_bw() +
  theme(
    panel.grid = element_blank(),
    panel.border = element_rect(fill = "transparent"),
    legend.position = "",
    plot.title = element_text(hjust = 0.5, size = 10)
  )
save(RCP45_SR_shift,file="C:/Users/2022207/Dropbox/Jack's PhD/Chapter 3. Future shark FD/Data/R data/Maps/RCP45_SR_shift_map.RData")

RCP45_FRic_shift <- ggplot() +
  geom_point(
    data = FD_shift_RCP45,
    aes(CenterLong, CenterLat, color = FRic_shift),
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
    limits = c(-0.5, 0.33),
    breaks = c(-0.5,-0.25,0,0.16,0.33),
    labels = scales::label_number(accuracy = 0.01),
    name = expression(Delta ~ FRic)
  ) +
  labs(title = "FRic (climate change)", x = "Longitude", y = "Latitude") +
  scale_x_continuous(breaks = seq(-180, 180, by = 45)) +
  scale_y_continuous(breaks = seq(-90, 90, by = 45)) +
  theme_minimal() +
  theme_bw() +
  theme(
    panel.grid = element_blank(),
    panel.border = element_rect(fill = "transparent"),
    legend.position = "",
    plot.title = element_text(hjust = 0.5, size = 10)
  )
save(RCP45_FRic_shift,file="C:/Users/2022207/Dropbox/Jack's PhD/Chapter 3. Future shark FD/Data/R data/Maps/RCP45_FRic_shift_map.RData")

## Climate change + extinction shifts
RCP45ext_SR_shift <- ggplot() +
  geom_point(
    data = FD_shift_RCP45ext,
    aes(CenterLong, CenterLat, color = SR_shift),
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
    limits = c(-86, 74),
    breaks = c(-86,-43,0,37,74),
    name = expression(Delta ~ SR)
  ) +
  labs(title = "SR (climate change + extinctions)", x = "Longitude", y = "Latitude") +
  scale_x_continuous(breaks = seq(-180, 180, by = 45)) +
  scale_y_continuous(breaks = seq(-90, 90, by = 45)) +
  theme_minimal() +
  theme_bw() +
  theme(
    panel.grid = element_blank(),
    panel.border = element_rect(fill = "transparent"),
    legend.position = "bottom",
    plot.title = element_text(hjust = 0.5, size = 10)
  ) +
  theme(legend.key.width = unit(1, "cm"))
save(RCP45ext_SR_shift,file="C:/Users/2022207/Dropbox/Jack's PhD/Chapter 3. Future shark FD/Data/R data/Maps/RCP45ext_SR_shift_map.RData")

RCP45ext_FRic_shift <- ggplot() +
  geom_point(
    data = FD_shift_RCP45ext,
    aes(CenterLong, CenterLat, color = FRic_shift),
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
    limits = c(-0.5, 0.33),
    breaks = c(-0.5,-0.25,0,0.16,0.33),
    labels = scales::label_number(accuracy = 0.01),
    name = expression(Delta ~ FRic)
  ) +
  labs(title = "FRic (climate change + extinctions)", x = "Longitude", y = "Latitude") +
  scale_x_continuous(breaks = seq(-180, 180, by = 45)) +
  scale_y_continuous(breaks = seq(-90, 90, by = 45)) +
  theme_minimal() +
  theme_bw() +
  theme(
    panel.grid = element_blank(),
    panel.border = element_rect(fill = "transparent"),
    legend.position = "bottom",
    plot.title = element_text(hjust = 0.5, size = 10)
  ) +
  theme(legend.key.width = unit(1, "cm"))
save(RCP45ext_FRic_shift,file="C:/Users/2022207/Dropbox/Jack's PhD/Chapter 3. Future shark FD/Data/R data/Maps/RCP45ext_FRic_shift_map.RData")

# Plot Figure 1
row1 <- plot_grid(
  Spp_CritE_2100 + theme(plot.margin = margin(5, 5, 5, 5)),       
  RCP45_SR_shift + theme(plot.margin = margin(5, 5, 5, 5)),   
  RCP45_FRic_shift + theme(plot.margin = margin(5, 5, 5, 5)),  
  FRic_CritE_2100 + theme(plot.margin = margin(5, 5, 5, 5)),    
  labels = c("A", "B", "C", "D"),
  label_size = 12,
  nrow = 1,
  rel_widths = c(4, 8, 8, 4)
)

row2 <- plot_grid(
  SR_RCP45 + theme(plot.margin = margin(5, 5, 5, 5)),          
  RCP45ext_SR_shift + theme(plot.margin = margin(5, 5, 5, 5)),  
  RCP45ext_FRic_shift + theme(plot.margin = margin(5, 5, 5, 5)),        
  FRic_RCP45 + theme(plot.margin = margin(5, 5, 5, 5)),                   
  labels = c("E", "F", "G", "H"),
  label_size = 12,
  nrow = 1,
  rel_widths = c(4, 8, 8, 4)
)

# Combine both rows into full figure
full_figure <- plot_grid(
  row1,
  row2,
  ncol = 1,
  rel_heights = c(5, 6.2)
)

# Calculate FRic hotspots (climate change)
## Create top 2.5% loss and gain
lower_quantile_FRic <- quantile(FD_shift_RCP45$FRic_shift, probs = 0.025, na.rm = TRUE)
upper_quantile_FRic <- quantile(FD_shift_RCP45$FRic_shift, probs = 0.975, na.rm = TRUE)

lower_quantile_FRic_ext <- quantile(FD_shift_RCP45ext$FRic_shift, probs = 0.025, na.rm = TRUE)
upper_quantile_FRic_ext <- quantile(FD_shift_RCP45ext$FRic_shift, probs = 0.975, na.rm = TRUE)

## Classify hotspots
FD_hotspots_RCP45 <- FD_shift_RCP45 %>%
  mutate(FRic_Hotspot = case_when(
    FRic_shift <= lower_quantile_FRic ~ "Loss",
    FRic_shift >= upper_quantile_FRic ~ "Gain",
    TRUE ~ NA_character_
  ))

FD_hotspots_RCP45ext <- FD_shift_RCP45ext %>%
  mutate(FRic_Hotspot = case_when(
    FRic_shift <= lower_quantile_FRic ~ "Loss",
    FRic_shift >= upper_quantile_FRic ~ "Gain",
    TRUE ~ NA_character_
  ))

## Identify gain hotspots retained following extinctions
fric_gain_RCP45 <- FD_hotspots_RCP45 %>%
  filter(FRic_Hotspot == "Gain") %>%
  select(Grid)

fric_gain_RCP45ext <- FD_hotspots_RCP45ext %>%
  filter(FRic_Hotspot == "Gain") %>%
  select(Grid)

persistent_fric_gain <- intersect(fric_gain_RCP45$Grid, fric_gain_RCP45ext$Grid)

FD_hotspots_RCP45ext <- FD_hotspots_RCP45ext %>%
  mutate(FRic_Hotspot_mod = case_when(
    Grid %in% persistent_fric_gain ~ "Gain",
    FRic_Hotspot == "Loss" ~ "Loss",
    TRUE ~ NA_character_
  ))

# Plot hotspots (climate change - Figure 2a)
Map_FRic_Hotspots_RCP45 <- ggplot() +
  geom_point(data = FD_hotspots_RCP45 %>% filter(!is.na(FRic_Hotspot)),
             aes(CenterLong, CenterLat, color = FRic_Hotspot),
             size = 0.01, alpha = 0.8, pch = 16) +
  scale_color_manual(values = c("Gain" = "#1f78b4", "Loss" = "#e31a1c"),
                     name = "FRic Hotspot") +
  geom_map(data = world, map = world,
           aes(long, lat, map_id = region),
           color = "black", fill = "lightgray", size = 0.2) +
  labs(title = "FRic zones of change (climate change)", x = "Longitude", y = "Latitude") +
  scale_x_continuous(breaks = seq(-180, 180, by = 45)) +
  scale_y_continuous(breaks = seq(-90, 90, by = 45)) +
  theme_minimal() +
  theme_bw() +
  theme(
    panel.grid = element_blank(),
    panel.border = element_rect(fill = "transparent"),
    legend.position = "",
    plot.title = element_text(hjust = 0.5, size = 10)
  )
save(Map_FRic_Hotspots_RCP45,file="C:/Users/2022207/Dropbox/Jack's PhD/Chapter 3. Future shark FD/Data/R data/Maps/RCP45_hotspots_map.RData")

# Plot results (climate change + extinctions - Figure 2b)
Map_FRic_Hotspots_RCP45ext <- ggplot() +
  geom_point(data = FD_hotspots_RCP45ext %>% filter(!is.na(FRic_Hotspot_mod)),
             aes(CenterLong, CenterLat, color = FRic_Hotspot_mod),
             size = 0.01, alpha = 0.8, pch = 16) +
  scale_color_manual(values = c("Gain" = "#1f78b4", "Loss" = "#e31a1c"),
                     name = "FRic Hotspot") +
  geom_map(data = world, map = world,
           aes(long, lat, map_id = region),
           color = "black", fill = "lightgray", size = 0.2) +
  labs(title = "FRic zones of change (climate change + extinctions)", x = "Longitude", y = "Latitude") +
  scale_x_continuous(breaks = seq(-180, 180, by = 45)) +
  scale_y_continuous(breaks = seq(-90, 90, by = 45)) +
  theme_minimal() +
  theme_bw() +
  theme(
    panel.grid = element_blank(),
    panel.border = element_rect(fill = "transparent"),
    legend.position = "",
    plot.title = element_text(hjust = 0.5, size = 10)
  )
save(Map_FRic_Hotspots_RCP45ext,file="C:/Users/2022207/Dropbox/Jack's PhD/Chapter 3. Future shark FD/Data/R data/Maps/RCP45ext_hotspots_map.RData")

# Form Figure 2
## Make shared legend
legend_data <- data.frame(
  x = c(1, 2),
  y = c(1, 1),
  label = c("Gain (top 2.5%)", "Loss (bottom 2.5%)"),
  color = c("Gain", "Loss")
)

legend_plot <- ggplot(legend_data, aes(x, y, color = color)) +
  geom_point(size = 3.5) +
  scale_color_manual(values = c("Gain" = "#1f78b4", "Loss" = "#e31a1c"),
                     labels = legend_data$label,
                     name = NULL) +
  theme_void() +
  theme(legend.position = "bottom",
        legend.text = element_text(size = 11),
        legend.key = element_blank(),
        legend.spacing.x = unit(0.3, "cm"))

shared_legend <- get_legend(legend_plot)

## Combine maps with aligned legend at bottom
final_plot <- plot_grid(
  plot_grid(
    Map_FRic_Hotspots_RCP45,
    Map_FRic_Hotspots_RCP45ext,
    labels = c("A", "B"),
    label_size = 12,
    ncol = 2
  ),
  shared_legend,
  ncol = 1,
  rel_heights = c(1, 0.08)  # Adjust if more space is needed
)
