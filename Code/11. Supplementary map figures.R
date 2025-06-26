###############################################################################################################
# 11. Supplementary map figures
## This R code plots all maps of species and functional richness
## Also includes latitudinal gradients 
## it produces Figures S5-S7
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
load(file="C:/Users/Jack Cooper/Dropbox/Jack's PhD/Chapter 3. Future shark FD/Data/R data/Maps/FRic_map_present.RData")
load(file="C:/Users/Jack Cooper/Dropbox/Jack's PhD/Chapter 3. Future shark FD/Data/R data/Maps/FRic_map_RCP26.RData")
load(file="C:/Users/Jack Cooper/Dropbox/Jack's PhD/Chapter 3. Future shark FD/Data/R data/Maps/FRic_map_RCP26ext.RData")
load(file="C:/Users/Jack Cooper/Dropbox/Jack's PhD/Chapter 3. Future shark FD/Data/R data/Maps/FRic_map_RCP85.RData")
load(file="C:/Users/Jack Cooper/Dropbox/Jack's PhD/Chapter 3. Future shark FD/Data/R data/Maps/FRic_map_RCP85ext.RData")
load(file="C:/Users/Jack Cooper/Dropbox/Jack's PhD/Chapter 3. Future shark FD/Data/R data/Maps/SR_map_present.RData")
load(file="C:/Users/Jack Cooper/Dropbox/Jack's PhD/Chapter 3. Future shark FD/Data/R data/Maps/SR_map_RCP26.RData")
load(file="C:/Users/Jack Cooper/Dropbox/Jack's PhD/Chapter 3. Future shark FD/Data/R data/Maps/SR_map_RCP26ext.RData")
load(file="C:/Users/Jack Cooper/Dropbox/Jack's PhD/Chapter 3. Future shark FD/Data/R data/Maps/SR_map_RCP85.RData")
load(file="C:/Users/Jack Cooper/Dropbox/Jack's PhD/Chapter 3. Future shark FD/Data/R data/Maps/SR_map_RCP85ext.RData")

# Calculate shifts under RCP 2.6 (climate change)
FD_shift_RCP26 <- FD_AqMap %>%
  inner_join(FD_RCP26, by = "Grid", suffix = c("_present", "_RCP26")) %>%
  mutate(
    SR_shift = sp_richn_RCP26 - sp_richn_present,
    FRic_shift = fric_RCP26 - fric_present
  ) %>%
  select(Grid, CenterLat_present, CenterLong_present, SR_shift, FRic_shift) %>%
  rename(CenterLat = CenterLat_present, CenterLong = CenterLong_present) %>% 
  na.omit()
save(FD_shift_RCP26,file="C:/Users/Jack Cooper/Dropbox/Jack's PhD/Chapter 3. Future shark FD/Data/R data/Maps/RCP26_shifts.RData")

# Calculate shifts under RCP 2.6 (climate change) + extinctions
FD_shift_RCP26ext <- FD_AqMap %>%
  inner_join(FD_RCP26ext, by = "Grid", suffix = c("_present", "_RCP26ext")) %>%
  mutate(
    SR_shift = sp_richn_RCP26ext - sp_richn_present,
    FRic_shift = fric_RCP26ext - fric_present
  ) %>%
  select(Grid, CenterLat_present, CenterLong_present, SR_shift, FRic_shift) %>%
  rename(CenterLat = CenterLat_present, CenterLong = CenterLong_present) %>% 
  na.omit()
save(FD_shift_RCP26ext,file="C:/Users/Jack Cooper/Dropbox/Jack's PhD/Chapter 3. Future shark FD/Data/R data/Maps/RCP26ext_shifts.RData")

# Calculate shifts under RCP 8.5 (climate change)
FD_shift_RCP85 <- FD_AqMap %>%
  inner_join(FD_RCP85, by = "Grid", suffix = c("_present", "_RCP85")) %>%
  mutate(
    SR_shift = sp_richn_RCP85 - sp_richn_present,
    FRic_shift = fric_RCP85 - fric_present
  ) %>%
  select(Grid, CenterLat_present, CenterLong_present, SR_shift, FRic_shift) %>%
  rename(CenterLat = CenterLat_present, CenterLong = CenterLong_present) %>% 
  na.omit()
save(FD_shift_RCP85,file="C:/Users/Jack Cooper/Dropbox/Jack's PhD/Chapter 3. Future shark FD/Data/R data/Maps/RCP85_shifts.RData")

# Calculate shifts under RCP 8.5 (climate change) + extinctions
FD_shift_RCP85ext <- FD_AqMap %>%
  inner_join(FD_RCP85ext, by = "Grid", suffix = c("_present", "_RCP85ext")) %>%
  mutate(
    SR_shift = sp_richn_RCP85ext - sp_richn_present,
    FRic_shift = fric_RCP85ext - fric_present
  ) %>%
  select(Grid, CenterLat_present, CenterLong_present, SR_shift, FRic_shift) %>%
  rename(CenterLat = CenterLat_present, CenterLong = CenterLong_present) %>% 
  na.omit()
save(FD_shift_RCP85ext,file="C:/Users/Jack Cooper/Dropbox/Jack's PhD/Chapter 3. Future shark FD/Data/R data/Maps/RCP85ext_shifts.RData")

# Calculate grid changes under all scenarios
sr_shift_RCP26 <- FD_shift_RCP26 %>%
  mutate(SR_change = case_when(
    SR_shift > 0 ~ "Increase",
    SR_shift < 0 ~ "Decrease",
    TRUE ~ "No change"
  )) %>%
  count(SR_change) %>%
  mutate(Proportion = n / sum(n))

fric_shift_RCP26 <- FD_shift_RCP26 %>%
  mutate(FRic_change = case_when(
    FRic_shift > 0 ~ "Increase",
    FRic_shift < 0 ~ "Decrease",
    TRUE ~ "No change"
  )) %>%
  count(FRic_change) %>%
  mutate(Proportion = n / sum(n))

sr_shift_RCP26ext <- FD_shift_RCP26ext %>%
  mutate(SR_change = case_when(
    SR_shift > 0 ~ "Increase",
    SR_shift < 0 ~ "Decrease",
    TRUE ~ "No change"
  )) %>%
  count(SR_change) %>%
  mutate(Proportion = n / sum(n))

fric_shift_RCP26ext <- FD_shift_RCP26ext %>%
  mutate(FRic_change = case_when(
    FRic_shift > 0 ~ "Increase",
    FRic_shift < 0 ~ "Decrease",
    TRUE ~ "No change"
  )) %>%
  count(FRic_change) %>%
  mutate(Proportion = n / sum(n))

sr_shift_RCP85 <- FD_shift_RCP85 %>%
  mutate(SR_change = case_when(
    SR_shift > 0 ~ "Increase",
    SR_shift < 0 ~ "Decrease",
    TRUE ~ "No change"
  )) %>%
  count(SR_change) %>%
  mutate(Proportion = n / sum(n))

fric_shift_RCP85 <- FD_shift_RCP85 %>%
  mutate(FRic_change = case_when(
    FRic_shift > 0 ~ "Increase",
    FRic_shift < 0 ~ "Decrease",
    TRUE ~ "No change"
  )) %>%
  count(FRic_change) %>%
  mutate(Proportion = n / sum(n))

sr_shift_RCP85ext <- FD_shift_RCP85ext %>%
  mutate(SR_change = case_when(
    SR_shift > 0 ~ "Increase",
    SR_shift < 0 ~ "Decrease",
    TRUE ~ "No change"
  )) %>%
  count(SR_change) %>%
  mutate(Proportion = n / sum(n))

fric_shift_RCP85ext <- FD_shift_RCP85ext %>%
  mutate(FRic_change = case_when(
    FRic_shift > 0 ~ "Increase",
    FRic_shift < 0 ~ "Decrease",
    TRUE ~ "No change"
  )) %>%
  count(FRic_change) %>%
  mutate(Proportion = n / sum(n))

# Tests to see if FRic shift is significantly lower following extinctions
FRic_compare_RCP26 <- FD_shift_RCP26 %>%
  inner_join(FD_shift_RCP26ext, by = "Grid", suffix = c("_rcp26", "_rcp26ext"))

wilcox.test(FRic_compare_RCP26$FRic_shift_rcp26, FRic_compare_RCP26$FRic_shift_rcp26ext, 
            paired = TRUE, alternative = "greater")

FRic_compare_RCP85 <- FD_shift_RCP85 %>%
  inner_join(FD_shift_RCP85ext, by = "Grid", suffix = c("_rcp85", "_rcp85ext"))

wilcox.test(FRic_compare_RCP85$FRic_shift_rcp85, FRic_compare_RCP85$FRic_shift_rcp85ext, 
            paired = TRUE, alternative = "greater")

# Calculate latitudinal gradient of SR and FRic shifts
## Calculate median, 95% and 50% confidence intervals - RCP 2.6
Lat_RCP26_FD <- FD_shift_RCP26 %>%
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

Lat_RCP26ext_FD <- FD_shift_RCP26ext %>%
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

Lat_RCP26_FD$Scenario <- "RCP 2.6"
Lat_RCP26ext_FD$Scenario <- "RCP 2.6 + Extinctions"
combined_RCP26_FD <- rbind(Lat_RCP26_FD, Lat_RCP26ext_FD)
save(combined_RCP26_FD,file="C:/Users/Jack Cooper/Dropbox/Jack's PhD/Chapter 3. Future shark FD/Data/R data/Maps/RCP26_latitudinal_gradients.RData")

## Plot SR changes
SR_RCP26 <- ggplot(combined_RCP26_FD, aes(x = CenterLat, y = median_value_SR, color = Scenario, fill = Scenario)) +
  geom_line(size = 1) + 
  geom_ribbon(aes(ymin = lower_50CI_SR, ymax = upper_50CI_SR), alpha = 0.2) + 
  labs(x = "Latitude", y = "SR shift") + 
  coord_flip() +  # Rotate the plot
  scale_x_continuous(limits = c(-60, 80), breaks = seq(-60, 80, by = 20)) + 
  scale_y_continuous(limits = c(-15, 15), breaks = seq(-15, 15, by = 5)) + 
  geom_hline(yintercept = 0, linetype = "dashed", color = "grey", size = 0.5) + 
  theme_minimal() +
  theme_bw() +
  theme(legend.position = "top", 
        axis.text = element_text(size = 8, color = "black"),
        axis.title = element_text(size = 10),
        panel.background = element_rect(fill = "white"),
        panel.grid = element_blank(),
        panel.border = element_rect(fill = "transparent"),
        plot.title = element_text(hjust = 0.5)) + 
  scale_color_manual(
    name = NULL,
    values = c("RCP 2.6" = "#0072B5FF", "RCP 2.6 + Extinctions" = "#BC3C29FF"),
    labels = c("RCP 2.6" = "Climate change", "RCP 2.6 + Extinctions" = "Climate change + extinctions")
  ) +
  scale_fill_manual(
    name = NULL,
    values = c("RCP 2.6" = "#0072B5FF", "RCP 2.6 + Extinctions" = "#BC3C29FF"),
    labels = c("RCP 2.6" = "Climate change", "RCP 2.6 + Extinctions" = "Climate change + extinctions")
  )+
  guides(color = guide_legend(nrow = 2), fill = guide_legend(nrow = 2))
save(SR_RCP26,file="C:/Users/Jack Cooper/Dropbox/Jack's PhD/Chapter 3. Future shark FD/Data/R data/Maps/RCP26_latitudinal_gradients_SR.RData")

## Plot FRic changes
FRic_RCP26 <- ggplot(combined_RCP26_FD, aes(x = CenterLat, y = median_value_FRic, color = Scenario, fill = Scenario)) +
  geom_line(size = 1) + 
  geom_ribbon(aes(ymin = lower_50CI_FRic, ymax = upper_50CI_FRic), alpha = 0.2) + 
  labs(x = "Latitude", y = "FRic shift") + 
  coord_flip() +  # Rotate the plot
  scale_x_continuous(limits = c(-60, 80), breaks = seq(-60, 80, by = 20)) + 
  scale_y_continuous(limits = c(-0.1, 0.2), breaks = seq(-0.1, 0.2, by = 0.1), labels = label_number(accuracy = 0.01)) + 
  geom_hline(yintercept = 0, linetype = "dashed", color = "grey", size = 0.5) + 
  theme_minimal() +
  theme_bw() +
  theme(legend.position = "top", 
        axis.text = element_text(size = 8, color = "black"),
        axis.title = element_text(size = 10),
        panel.background = element_rect(fill = "white"),
        panel.grid = element_blank(),
        panel.border = element_rect(fill = "transparent"),
        plot.title = element_text(hjust = 0.5)) + 
  scale_color_manual(
    name = NULL,
    values = c("RCP 2.6" = "#0072B5FF", "RCP 2.6 + Extinctions" = "#BC3C29FF"),
    labels = c("RCP 2.6" = "Climate change", "RCP 2.6 + Extinctions" = "Climate change + extinctions")
  ) +
  scale_fill_manual(
    name = NULL,
    values = c("RCP 2.6" = "#0072B5FF", "RCP 2.6 + Extinctions" = "#BC3C29FF"),
    labels = c("RCP 2.6" = "Climate change", "RCP 2.6 + Extinctions" = "Climate change + extinctions")
  )+
  guides(color = guide_legend(nrow = 2), fill = guide_legend(nrow = 2))
save(FRic_RCP26,file="C:/Users/Jack Cooper/Dropbox/Jack's PhD/Chapter 3. Future shark FD/Data/R data/Maps/RCP26_latitudinal_gradients_FRic.RData")

## Calculate median, 95% and 50% confidence intervals - RCP 8.5
Lat_RCP85_FD <- FD_shift_RCP85 %>%
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

Lat_RCP85ext_FD <- FD_shift_RCP85ext %>%
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

Lat_RCP85_FD$Scenario <- "RCP 8.5"
Lat_RCP85ext_FD$Scenario <- "RCP 8.5 + Extinctions"
combined_RCP85_FD <- rbind(Lat_RCP85_FD, Lat_RCP85ext_FD)
save(combined_RCP85_FD,file="C:/Users/Jack Cooper/Dropbox/Jack's PhD/Chapter 3. Future shark FD/Data/R data/Maps/RCP85_latitudinal_gradients.RData")

## Plot SR changes
SR_RCP85 <- ggplot(combined_RCP85_FD, aes(x = CenterLat, y = median_value_SR, color = Scenario, fill = Scenario)) +
  geom_line(size = 1) + 
  geom_ribbon(aes(ymin = lower_50CI_SR, ymax = upper_50CI_SR), alpha = 0.2) + 
  labs(x = "Latitude", y = "SR shift") + 
  coord_flip() +  # Rotate the plot
  scale_x_continuous(limits = c(-60, 80), breaks = seq(-60, 80, by = 20)) + 
  #scale_y_continuous(limits = c(-15, 15), breaks = seq(-15, 15, by = 5)) + 
  geom_hline(yintercept = 0, linetype = "dashed", color = "grey", size = 0.5) + 
  theme_minimal() +
  theme_bw() +
  theme(legend.position = "top", 
        axis.text = element_text(size = 8, color = "black"),
        axis.title = element_text(size = 10),
        panel.background = element_rect(fill = "white"),
        panel.grid = element_blank(),
        panel.border = element_rect(fill = "transparent"),
        plot.title = element_text(hjust = 0.5)) + 
  scale_color_manual(
    name = NULL,
    values = c("RCP 8.5" = "#0072B5FF", "RCP 8.5 + Extinctions" = "#BC3C29FF"),
    labels = c("RCP 8.5" = "Climate change", "RCP 8.5 + Extinctions" = "Climate change + extinctions")
  ) +
  scale_fill_manual(
    name = NULL,
    values = c("RCP 8.5" = "#0072B5FF", "RCP 8.5 + Extinctions" = "#BC3C29FF"),
    labels = c("RCP 8.5" = "Climate change", "RCP 8.5 + Extinctions" = "Climate change + extinctions")
  )+
  guides(color = guide_legend(nrow = 2), fill = guide_legend(nrow = 2))
save(SR_RCP85,file="C:/Users/Jack Cooper/Dropbox/Jack's PhD/Chapter 3. Future shark FD/Data/R data/Maps/RCP85_latitudinal_gradients_SR.RData")

## Plot FRic changes
FRic_RCP85 <- ggplot(combined_RCP85_FD, aes(x = CenterLat, y = median_value_FRic, color = Scenario, fill = Scenario)) +
  geom_line(size = 1) + 
  geom_ribbon(aes(ymin = lower_50CI_FRic, ymax = upper_50CI_FRic), alpha = 0.2) + 
  labs(x = "Latitude", y = "FRic shift") + 
  coord_flip() +  # Rotate the plot
  scale_x_continuous(limits = c(-60, 80), breaks = seq(-60, 80, by = 20)) + 
  scale_y_continuous(limits = c(-0.1, 0.2), breaks = seq(-0.1, 0.2, by = 0.1), labels = label_number(accuracy = 0.01)) + 
  geom_hline(yintercept = 0, linetype = "dashed", color = "grey", size = 0.5) + 
  theme_minimal() +
  theme_bw() +
  theme(legend.position = "top", 
        axis.text = element_text(size = 8, color = "black"),
        axis.title = element_text(size = 10),
        panel.background = element_rect(fill = "white"),
        panel.grid = element_blank(),
        panel.border = element_rect(fill = "transparent"),
        plot.title = element_text(hjust = 0.5)) + 
  scale_color_manual(
    name = NULL,
    values = c("RCP 8.5" = "#0072B5FF", "RCP 8.5 + Extinctions" = "#BC3C29FF"),
    labels = c("RCP 8.5" = "Climate change", "RCP 8.5 + Extinctions" = "Climate change + extinctions")
  ) +
  scale_fill_manual(
    name = NULL,
    values = c("RCP 8.5" = "#0072B5FF", "RCP 8.5 + Extinctions" = "#BC3C29FF"),
    labels = c("RCP 8.5" = "Climate change", "RCP 8.5 + Extinctions" = "Climate change + extinctions")
  )+
  guides(color = guide_legend(nrow = 2), fill = guide_legend(nrow = 2))
save(FRic_RCP85,file="C:/Users/Jack Cooper/Dropbox/Jack's PhD/Chapter 3. Future shark FD/Data/R data/Maps/RCP85_latitudinal_gradients_FRic.RData")

# Plot and save all maps so far
## Set breaks and colours for plots
RCP_colors <- c("#e31a1c", "#ff7f00", "#ffffff", "#6a3d9a", "#1f78b4")
SR_breaks_RCP26 <- c(-77,-38,0,1,37,73)
FRic_breaks_RCP26 <- c(-0.45,-0.23,0,0.0000001,0.16,0.33)
SR_breaks_RCP85 <- c(-179,-80,0,1,54,107)
FRic_breaks_RCP85 <- c(-0.8,-0.4,0,0.0000001,0.25,0.5)

## Climate change shifts
RCP26_SR_shift <- ggplot() +
  geom_point(
    data = FD_shift_RCP26,
    aes(CenterLong, CenterLat, color = SR_shift),
    size = 0.01,
    alpha = 0.5,
    inherit.aes = FALSE
  ) +
  geom_map(data = world, map = world,
           aes(long, lat, map_id = region),
           color = "black", fill = "lightgray", size = 0.2) +
  scale_color_gradientn(
    colors = RCP_colors,
    values = scales::rescale(SR_breaks_RCP26, to = c(0, 1)),
    limits = c(-77, 73),
    breaks = c(-77,-38,0,37,73),
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
    legend.position = "right",
    plot.title = element_text(hjust = 0.5, size = 10)
  )
save(RCP26_SR_shift,file="C:/Users/Jack Cooper/Dropbox/Jack's PhD/Chapter 3. Future shark FD/Data/R data/Maps/RCP26_SR_shift_map.RData")

RCP26_FRic_shift <- ggplot() +
  geom_point(
    data = FD_shift_RCP26,
    aes(CenterLong, CenterLat, color = FRic_shift),
    size = 0.01,
    alpha = 0.5,
    inherit.aes = FALSE
  ) +
  geom_map(data = world, map = world,
           aes(long, lat, map_id = region),
           color = "black", fill = "lightgray", size = 0.2) +
  scale_color_gradientn(
    colors = RCP_colors,
    values = scales::rescale(FRic_breaks_RCP26, to = c(0, 1)),
    limits = c(-0.45, 0.33),
    breaks = c(-0.45,-0.23,0,0.16,0.33),
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
    legend.position = "right",
    plot.title = element_text(hjust = 0.5, size = 10)
  )
save(RCP26_FRic_shift,file="C:/Users/Jack Cooper/Dropbox/Jack's PhD/Chapter 3. Future shark FD/Data/R data/Maps/RCP26_FRic_shift_map.RData")

RCP85_SR_shift <- ggplot() +
  geom_point(
    data = FD_shift_RCP85,
    aes(CenterLong, CenterLat, color = SR_shift),
    size = 0.01,
    alpha = 0.5,
    inherit.aes = FALSE
  ) +
  geom_map(data = world, map = world,
           aes(long, lat, map_id = region),
           color = "black", fill = "lightgray", size = 0.2) +
  scale_color_gradientn(
    colors = RCP_colors,
    values = scales::rescale(SR_breaks_RCP85, to = c(0, 1)),
    limits = c(-179, 107),
    breaks = c(-179,-80,0,54,107),
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
    legend.position = "right",
    plot.title = element_text(hjust = 0.5, size = 10)
  )
save(RCP85_SR_shift,file="C:/Users/Jack Cooper/Dropbox/Jack's PhD/Chapter 3. Future shark FD/Data/R data/Maps/RCP85_SR_shift_map.RData")

RCP85_FRic_shift <- ggplot() +
  geom_point(
    data = FD_shift_RCP85,
    aes(CenterLong, CenterLat, color = FRic_shift),
    size = 0.01,
    alpha = 0.5,
    inherit.aes = FALSE
  ) +
  geom_map(data = world, map = world,
           aes(long, lat, map_id = region),
           color = "black", fill = "lightgray", size = 0.2) +
  scale_color_gradientn(
    colors = RCP_colors,
    values = scales::rescale(FRic_breaks_RCP85, to = c(0, 1)),
    limits = c(-0.8, 0.5),
    breaks = c(-0.8,-0.4,0,0.25,0.5),
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
    legend.position = "right",
    plot.title = element_text(hjust = 0.5, size = 10)
  )
save(RCP85_FRic_shift,file="C:/Users/Jack Cooper/Dropbox/Jack's PhD/Chapter 3. Future shark FD/Data/R data/Maps/RCP85_FRic_shift_map.RData")

## Climate change + extinction shifts
RCP26ext_SR_shift <- ggplot() +
  geom_point(
    data = FD_shift_RCP26ext,
    aes(CenterLong, CenterLat, color = SR_shift),
    size = 0.01,
    alpha = 0.5,
    inherit.aes = FALSE
  ) +
  geom_map(data = world, map = world,
           aes(long, lat, map_id = region),
           color = "black", fill = "lightgray", size = 0.2) +
  scale_color_gradientn(
    colors = RCP_colors,
    values = scales::rescale(SR_breaks_RCP26, to = c(0, 1)),
    limits = c(-77, 73),
    breaks = c(-77,-38,0,37,73),
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
    legend.position = "right",
    plot.title = element_text(hjust = 0.5, size = 10)
  ) 
save(RCP26ext_SR_shift,file="C:/Users/Jack Cooper/Dropbox/Jack's PhD/Chapter 3. Future shark FD/Data/R data/Maps/RCP26ext_SR_shift_map.RData")

RCP26ext_FRic_shift <- ggplot() +
  geom_point(
    data = FD_shift_RCP26ext,
    aes(CenterLong, CenterLat, color = FRic_shift),
    size = 0.01,
    alpha = 0.5,
    inherit.aes = FALSE
  ) +
  geom_map(data = world, map = world,
           aes(long, lat, map_id = region),
           color = "black", fill = "lightgray", size = 0.2) +
  scale_color_gradientn(
    colors = RCP_colors,
    values = scales::rescale(FRic_breaks_RCP26, to = c(0, 1)),
    limits = c(-0.45, 0.33),
    breaks = c(-0.45,-0.23,0,0.16,0.33),
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
    legend.position = "right",
    plot.title = element_text(hjust = 0.5, size = 10)
  ) 
save(RCP26ext_FRic_shift,file="C:/Users/Jack Cooper/Dropbox/Jack's PhD/Chapter 3. Future shark FD/Data/R data/Maps/RCP26ext_FRic_shift_map.RData")

RCP85ext_SR_shift <- ggplot() +
  geom_point(
    data = FD_shift_RCP85ext,
    aes(CenterLong, CenterLat, color = SR_shift),
    size = 0.01,
    alpha = 0.5,
    inherit.aes = FALSE
  ) +
  geom_map(data = world, map = world,
           aes(long, lat, map_id = region),
           color = "black", fill = "lightgray", size = 0.2) +
  scale_color_gradientn(
    colors = RCP_colors,
    values = scales::rescale(SR_breaks_RCP85, to = c(0, 1)),
    limits = c(-179, 107),
    breaks = c(-179,-80,0,54,107),
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
    legend.position = "right",
    plot.title = element_text(hjust = 0.5, size = 10)
  ) 
save(RCP85ext_SR_shift,file="C:/Users/Jack Cooper/Dropbox/Jack's PhD/Chapter 3. Future shark FD/Data/R data/Maps/RCP85ext_SR_shift_map.RData")

RCP85ext_FRic_shift <- ggplot() +
  geom_point(
    data = FD_shift_RCP85ext,
    aes(CenterLong, CenterLat, color = FRic_shift),
    size = 0.01,
    alpha = 0.5,
    inherit.aes = FALSE
  ) +
  geom_map(data = world, map = world,
           aes(long, lat, map_id = region),
           color = "black", fill = "lightgray", size = 0.2) +
  scale_color_gradientn(
    colors = RCP_colors,
    values = scales::rescale(FRic_breaks_RCP85, to = c(0, 1)),
    limits = c(-0.8, 0.5),
    breaks = c(-0.8,-0.4,0,0.25,0.5),
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
    legend.position = "right",
    plot.title = element_text(hjust = 0.5, size = 10)
  ) 
save(RCP85ext_FRic_shift,file="C:/Users/Jack Cooper/Dropbox/Jack's PhD/Chapter 3. Future shark FD/Data/R data/Maps/RCP85ext_FRic_shift_map.RData")

## Form Figures S5 and S6
Fig_S5 <- plot_grid(RCP26_SR_shift, RCP26ext_SR_shift, SR_RCP26,
                    RCP85_SR_shift, RCP85ext_SR_shift, SR_RCP85,
                    labels = c("A", "B", "C", "D", "E", "F"), 
                    label_size = 12, align = "hv", label_fontface = "bold",
                    nrow = 2,
                    rel_widths = c(1, 1, 0.5)) 

Fig_S6 <- plot_grid(RCP26_FRic_shift,RCP26ext_FRic_shift,FRic_RCP26,
                    RCP85_FRic_shift,RCP85ext_FRic_shift,FRic_RCP85,
                    labels = c("A", "B", "C", "D", "E", "F"), 
                    label_size = 12, align = "hv", label_fontface = "bold",
                    nrow = 2,
                    rel_widths = c(1, 1, 0.5))

# Hotspot analyses
# Calculate FRic hotspots (climate change)
## Create top 2.5% loss and gain
lower_quantile_FRic_RCP26 <- quantile(FD_shift_RCP26$FRic_shift, probs = 0.025, na.rm = TRUE)
upper_quantile_FRic_RCP26 <- quantile(FD_shift_RCP26$FRic_shift, probs = 0.975, na.rm = TRUE)

lower_quantile_FRic_RCP26ext <- quantile(FD_shift_RCP26ext$FRic_shift, probs = 0.025, na.rm = TRUE)
upper_quantile_FRic_RCP26ext <- quantile(FD_shift_RCP26ext$FRic_shift, probs = 0.975, na.rm = TRUE)

lower_quantile_FRic_RCP85 <- quantile(FD_shift_RCP85$FRic_shift, probs = 0.025, na.rm = TRUE)
upper_quantile_FRic_RCP85 <- quantile(FD_shift_RCP85$FRic_shift, probs = 0.975, na.rm = TRUE)

lower_quantile_FRic_RCP85ext <- quantile(FD_shift_RCP85ext$FRic_shift, probs = 0.025, na.rm = TRUE)
upper_quantile_FRic_RCP85ext <- quantile(FD_shift_RCP85ext$FRic_shift, probs = 0.975, na.rm = TRUE)

## Classify hotspots
FD_hotspots_RCP26 <- FD_shift_RCP26 %>%
  mutate(FRic_Hotspot = case_when(
    FRic_shift <= lower_quantile_FRic_RCP26 ~ "Loss",
    FRic_shift >= upper_quantile_FRic_RCP26 ~ "Gain",
    TRUE ~ NA_character_
  ))

FD_hotspots_RCP26ext <- FD_shift_RCP26ext %>%
  mutate(FRic_Hotspot = case_when(
    FRic_shift <= lower_quantile_FRic_RCP26ext ~ "Loss",
    FRic_shift >= upper_quantile_FRic_RCP26ext ~ "Gain",
    TRUE ~ NA_character_
  ))

FD_hotspots_RCP85 <- FD_shift_RCP85 %>%
  mutate(FRic_Hotspot = case_when(
    FRic_shift <= lower_quantile_FRic_RCP85 ~ "Loss",
    FRic_shift >= upper_quantile_FRic_RCP85 ~ "Gain",
    TRUE ~ NA_character_
  ))

FD_hotspots_RCP85ext <- FD_shift_RCP85ext %>%
  mutate(FRic_Hotspot = case_when(
    FRic_shift <= lower_quantile_FRic_RCP85ext ~ "Loss",
    FRic_shift >= upper_quantile_FRic_RCP85ext ~ "Gain",
    TRUE ~ NA_character_
  ))

## Identify gain hotspots retained following extinctions
fric_gain_RCP26 <- FD_hotspots_RCP26 %>%
  filter(FRic_Hotspot == "Gain") %>%
  select(Grid)

fric_gain_RCP26ext <- FD_hotspots_RCP26ext %>%
  filter(FRic_Hotspot == "Gain") %>%
  select(Grid)

persistent_fric_gain <- intersect(fric_gain_RCP26$Grid, fric_gain_RCP26ext$Grid)

FD_hotspots_RCP26ext <- FD_hotspots_RCP26ext %>%
  mutate(FRic_Hotspot_mod = case_when(
    Grid %in% persistent_fric_gain ~ "Gain",
    FRic_Hotspot == "Loss" ~ "Loss",
    TRUE ~ NA_character_
  ))

fric_gain_RCP85 <- FD_hotspots_RCP85 %>%
  filter(FRic_Hotspot == "Gain") %>%
  select(Grid)

fric_gain_RCP85ext <- FD_hotspots_RCP85ext %>%
  filter(FRic_Hotspot == "Gain") %>%
  select(Grid)

persistent_fric_gain <- intersect(fric_gain_RCP85$Grid, fric_gain_RCP85ext$Grid)

FD_hotspots_RCP85ext <- FD_hotspots_RCP85ext %>%
  mutate(FRic_Hotspot_mod = case_when(
    Grid %in% persistent_fric_gain ~ "Gain",
    FRic_Hotspot == "Loss" ~ "Loss",
    TRUE ~ NA_character_
  ))

# Plot hotspots (Figure S7)
Map_FRic_Hotspots_RCP26 <- ggplot() +
  geom_point(data = FD_hotspots_RCP26 %>% filter(!is.na(FRic_Hotspot)),
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
save(Map_FRic_Hotspots_RCP26, file = "C:/Users/Jack Cooper/Dropbox/Jack's PhD/Chapter 3. Future shark FD/Data/R data/Maps/RCP26_hotspots_map.RData")

Map_FRic_Hotspots_RCP26ext <- ggplot() +
  geom_point(data = FD_hotspots_RCP26ext %>% filter(!is.na(FRic_Hotspot_mod)),
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
save(Map_FRic_Hotspots_RCP26ext, file = "C:/Users/Jack Cooper/Dropbox/Jack's PhD/Chapter 3. Future shark FD/Data/R data/Maps/RCP26ext_hotspots_map.RData")

Map_FRic_Hotspots_RCP85 <- ggplot() +
  geom_point(data = FD_hotspots_RCP85 %>% filter(!is.na(FRic_Hotspot)),
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
save(Map_FRic_Hotspots_RCP85,file="C:/Users/Jack Cooper/Dropbox/Jack's PhD/Chapter 3. Future shark FD/Data/R data/Maps/RCP85_hotspots_map.RData")

Map_FRic_Hotspots_RCP85ext <- ggplot() +
  geom_point(data = FD_hotspots_RCP85ext %>% filter(!is.na(FRic_Hotspot_mod)),
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
save(Map_FRic_Hotspots_RCP85ext,file="C:/Users/Jack Cooper/Dropbox/Jack's PhD/Chapter 3. Future shark FD/Data/R data/Maps/RCP85ext_hotspots_map.RData")

# Form Figure
# Create legend
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

# Extract shared legend
shared_legend <- get_legend(legend_plot)

# Combine maps
selected_maps <- plot_grid(
  Map_FRic_Hotspots_RCP26,
  Map_FRic_Hotspots_RCP26ext,
  Map_FRic_Hotspots_RCP85,
  Map_FRic_Hotspots_RCP85ext,
  labels = c("A", "B", "C", "D"),
  label_size = 12,
  ncol = 2
)

# Final combined plot with legend
final_plot <- plot_grid(
  selected_maps,
  shared_legend,
  ncol = 1,
  rel_heights = c(1, 0.08)
)
