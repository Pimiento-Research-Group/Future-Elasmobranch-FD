###############################################################################################################
# 05. Future FUn-FSp scores
## This R code forms future functional spaces for FUn and FSp scores
## it produces produces Figure S6
###############################################################################################################

# Import packages
library(tidyverse)
library(tibble)
library(scales)
library(cowplot)
library(doBy)
library(reshape2)
library(RColorBrewer)
library(mFD)
library(ggsci)
library(data.table)
library(funrar)
library(viridis)

# Load data
load(file="C:/Users/2022207/Dropbox/Jack's PhD/Chapter 3. Future shark FD/Data/R data/Data.RData")

# Load mode IUCN_sim results
load(file = "C:/Users/2022207/Dropbox/Jack's PhD/Chapter 3. Future shark FD/Data/R data/2100_mode.RData")

# Make future functional space
# Occurrence data
Occ_Future <- data %>% 
  select(Species)

Occ_Future$Present <- 1

# Set up future scenarios
## 2100
Occ_Future <- Occ_Future %>%
  left_join(Spp_mode, by = "Species") %>%
  mutate(Future = ifelse(is.na(mode), 1, ifelse(mode <= 77, 0, 1))) %>%
  select(Species,Present,Future)

## 100 years
Occ_Future <- Occ_Future %>%
  left_join(Spp_mode, by = "Species") %>%
  mutate(Y100 = ifelse(is.na(mode), 1, ifelse(mode <= 100, 0, 1))) %>%
  select(Species,Present,Future,Y100)

## 200 years
Occ_Future <- Occ_Future %>%
  left_join(Spp_mode, by = "Species") %>%
  mutate(Y200 = ifelse(is.na(mode), 1, ifelse(mode <= 200, 0, 1))) %>%
  select(Species,Present,Future,Y100,Y200)

## 300 years
Occ_Future <- Occ_Future %>%
  left_join(Spp_mode, by = "Species") %>%
  mutate(Y300 = ifelse(is.na(mode), 1, ifelse(mode <= 300, 0, 1))) %>%
  select(Species,Present,Future,Y100,Y200,Y300)

## 400 years
Occ_Future <- Occ_Future %>%
  left_join(Spp_mode, by = "Species") %>%
  mutate(Y400 = ifelse(is.na(mode), 1, ifelse(mode <= 400, 0, 1))) %>%
  select(Species,Present,Future,Y100,Y200,Y300,Y400)

## 500 years
Occ_Future <- Occ_Future %>%
  left_join(Spp_mode, by = "Species") %>%
  mutate(Y500 = ifelse(is.na(mode), 1, ifelse(mode <= 500, 0, 1))) %>%
  select(Species,Present,Future,Y100,Y200,Y300,Y400,Y500)

baskets <- Occ_Future %>% 
  remove_rownames %>% 
  column_to_rownames(var="Species")

# Transpose to form final occurrence matrix
baskets_sharks_weights <- t(baskets)

# Species-trait matrix
traits <- data %>% 
  select(Species,habitat,vertical,terrestriality,thermo,feeding,diet,max.length)

sharks_traits <- traits %>% 
  remove_rownames %>% 
  column_to_rownames(var="Species")

# Form trait category matrix
sharks_traits_cat <- tibble(trait_name = c("habitat","vertical","terrestriality","thermo","feeding","diet","max.length"),
                            trait_type = c("N", "N", "O", "O", "O", "N", "Q"))

# Summarise dataset
sharks_traits_summ <- mFD::sp.tr.summary(
  tr_cat     = sharks_traits_cat,   
  sp_tr      = sharks_traits, 
  stop_if_NA = TRUE)

# Reveal number of states for each trait
sharks_traits_summ$tr_summary_list

# make assemblage data into a matrix & ensure all matrix values are numeric
baskets_sharks_weights <- data.matrix(baskets_sharks_weights, rownames.force = NA)
class(baskets_sharks_weights) <- "numeric"

#summarise dataset
asb_sp_sharks_summ <- mFD::asb.sp.summary(asb_sp_w = baskets_sharks_weights)

# retrieve species occurrences for all assemblages:
asb_sp_sharks_occ <- asb_sp_sharks_summ$asb_sp_occ

# Make sharks_traits matrices for each scenario (for FUn_FSp score plots)
## 2100
sharks_traits_2100.iucn <- sharks_traits %>% 
  rownames_to_column(var="Species")

sharks_traits_2100.iucn$Future <- Occ_Future$Future
sharks_traits_2100.iucn$iucn <- data$iucn

sharks_traits_2100.iucn <- sharks_traits_2100.iucn %>% 
  filter(Future == "1") %>% 
  select(-Future)

sharks_traits_2100 <- sharks_traits_2100.iucn %>% 
  remove_rownames %>% 
  column_to_rownames(var="Species") %>% 
  select(-iucn)

## 100 years
sharks_traits_Y100.iucn <- sharks_traits %>% 
  rownames_to_column(var="Species")

sharks_traits_Y100.iucn$Y100 <- Occ_Future$Y100
sharks_traits_Y100.iucn$iucn <- data$iucn

sharks_traits_Y100.iucn <- sharks_traits_Y100.iucn %>% 
  filter(Y100 == "1") %>% 
  select(-Y100)

sharks_traits_Y100 <- sharks_traits_Y100.iucn %>% 
  remove_rownames %>% 
  column_to_rownames(var="Species") %>% 
  select(-iucn)

## 200 years
sharks_traits_Y200.iucn <- sharks_traits %>% 
  rownames_to_column(var="Species")

sharks_traits_Y200.iucn$Y200 <- Occ_Future$Y200
sharks_traits_Y200.iucn$iucn <- data$iucn

sharks_traits_Y200.iucn <- sharks_traits_Y200.iucn %>% 
  filter(Y200 == "1") %>% 
  select(-Y200)

sharks_traits_Y200 <- sharks_traits_Y200.iucn %>% 
  remove_rownames %>% 
  column_to_rownames(var="Species") %>% 
  select(-iucn)

## 300 years
sharks_traits_Y300.iucn <- sharks_traits %>% 
  rownames_to_column(var="Species")

sharks_traits_Y300.iucn$Y300 <- Occ_Future$Y300
sharks_traits_Y300.iucn$iucn <- data$iucn

sharks_traits_Y300.iucn <- sharks_traits_Y300.iucn %>% 
  filter(Y300 == "1") %>% 
  select(-Y300)

sharks_traits_Y300 <- sharks_traits_Y300.iucn %>% 
  remove_rownames %>% 
  column_to_rownames(var="Species") %>% 
  select(-iucn)

## 400 years
sharks_traits_Y400.iucn <- sharks_traits %>% 
  rownames_to_column(var="Species")

sharks_traits_Y400.iucn$Y400 <- Occ_Future$Y400
sharks_traits_Y400.iucn$iucn <- data$iucn

sharks_traits_Y400.iucn <- sharks_traits_Y400.iucn %>% 
  filter(Y400 == "1") %>% 
  select(-Y400)

sharks_traits_Y400 <- sharks_traits_Y400.iucn %>% 
  remove_rownames %>% 
  column_to_rownames(var="Species") %>% 
  select(-iucn)

## 500 years
sharks_traits_Y500.iucn <- sharks_traits %>% 
  rownames_to_column(var="Species")

sharks_traits_Y500.iucn$Y500 <- Occ_Future$Y500
sharks_traits_Y500.iucn$iucn <- data$iucn

sharks_traits_Y500.iucn <- sharks_traits_Y500.iucn %>% 
  filter(Y500 == "1") %>% 
  select(-Y500)

sharks_traits_Y500 <- sharks_traits_Y500.iucn %>% 
  remove_rownames %>% 
  column_to_rownames(var="Species") %>% 
  select(-iucn)

# Future FUn and FSp score plots
## 2100
sp_dist_sharks_2100 <- mFD::funct.dist(
  sp_tr         = sharks_traits_2100,
  tr_cat        = sharks_traits_cat,
  metric        = "gower",
  scale_euclid  = "scale_center",
  ordinal_var   = "classic",
  weight_type   = "equal",
  stop_if_NA    = TRUE)

# Assess quality
fspaces_quality_sharks_2100 <- mFD::quality.fspaces(
  sp_dist             = sp_dist_sharks_2100,
  maxdim_pcoa         = 10,
  deviation_weighting = c("absolute","squared"),
  fdist_scaling       = FALSE,
  fdendro             = "average")

# Return coordinates of each axis
sp_faxes_coord_sharks_2100 <- fspaces_quality_sharks_2100$"details_fspaces"$"sp_pc_coord"

# Calculate per-species FUn and FSp
FUSE.future <- fuse(sp_dist        = sp_dist_sharks_2100, 
                    sp_faxes_coord = as.matrix(sp_faxes_coord_sharks_2100), 
                    nb_NN          = 5,  
                    GE             = sharks_traits_2100.iucn$iucn,
                    standGE        = TRUE)

FUSE.future$PC1 <- sp_faxes_coord_sharks_2100[, "PC1"]
FUSE.future$PC2 <- sp_faxes_coord_sharks_2100[, "PC2"]
FUSE.future$PC3 <- sp_faxes_coord_sharks_2100[, "PC3"]
FUSE.future$PC4 <- sp_faxes_coord_sharks_2100[, "PC4"]

# Form legends
color_scale_FUn <- scale_color_viridis(name = "FUn_std", guide = "none")
color_scale_FSp <- scale_color_viridis(name = "FSp_std", guide = "none")

#FUn
FUn_PC12_2100 <- ggplot(FUSE.future, aes(x = PC1, y = PC2, color = FUn_std)) +
  geom_point(shape = 16, size = 4) +
  labs(x = "PCoA1", y = "PCoA2") +
  color_scale_FUn +
  theme_minimal() +
  theme_bw() +
  theme(legend.position = "none") +
  theme(axis.text.y = element_text(size = 8, color = "black"),  
        axis.title = element_text(size = 12), panel.background = element_rect(fill = "white")) +
  theme(panel.grid = element_blank(), 
        panel.border = element_rect(fill = "transparent"))

# FSp
FSp_PC12_2100 <- ggplot(FUSE.future, aes(x = PC1, y = PC2, color = FSp_std)) +
  geom_point(shape = 16, size = 4) +
  labs(x = "PCoA1", y = "PCoA2") +
  color_scale_FSp +
  theme_minimal() +
  theme_bw() +
  theme(legend.position = "none") +
  theme(axis.text.y = element_text(size = 8, color = "black"),  
        axis.title = element_text(size = 12), panel.background = element_rect(fill = "white")) +
  theme(panel.grid = element_blank(), 
        panel.border = element_rect(fill = "transparent"))

## 100 years
sp_dist_sharks_Y100 <- mFD::funct.dist(
  sp_tr         = sharks_traits_Y100,
  tr_cat        = sharks_traits_cat,
  metric        = "gower",
  scale_euclid  = "scale_center",
  ordinal_var   = "classic",
  weight_type   = "equal",
  stop_if_NA    = TRUE)

# Assess quality
fspaces_quality_sharks_Y100 <- mFD::quality.fspaces(
  sp_dist             = sp_dist_sharks_Y100,
  maxdim_pcoa         = 10,
  deviation_weighting = c("absolute","squared"),
  fdist_scaling       = FALSE,
  fdendro             = "average")

# Return coordinates of each axis
sp_faxes_coord_sharks_Y100 <- fspaces_quality_sharks_Y100$"details_fspaces"$"sp_pc_coord"

# Calculate per-species FUn and FSp
FUSE.Y100 <- fuse(sp_dist        = sp_dist_sharks_Y100, 
                  sp_faxes_coord = as.matrix(sp_faxes_coord_sharks_Y100), 
                  nb_NN          = 5,  
                  GE             = sharks_traits_Y100.iucn$iucn,
                  standGE        = TRUE)

FUSE.Y100$PC1 <- sp_faxes_coord_sharks_Y100[, "PC1"]
FUSE.Y100$PC2 <- sp_faxes_coord_sharks_Y100[, "PC2"]
FUSE.Y100$PC3 <- sp_faxes_coord_sharks_Y100[, "PC3"]
FUSE.Y100$PC4 <- sp_faxes_coord_sharks_Y100[, "PC4"]

# Form legends
color_scale_FUn <- scale_color_viridis(name = "FUn_std", guide = "none")
color_scale_FSp <- scale_color_viridis(name = "FSp_std", guide = "none")

#FUn
FUn_PC12_100 <- ggplot(FUSE.Y100, aes(x = PC1, y = PC2, color = FUn_std)) +
  geom_point(shape = 16, size = 4) +
  labs(x = "PCoA1", y = "PCoA2") +
  color_scale_FUn +
  theme_minimal() +
  theme_bw() +
  theme(legend.position = "none") +
  theme(axis.text.y = element_text(size = 8, color = "black"),  
        axis.title = element_text(size = 12), panel.background = element_rect(fill = "white")) +
  theme(panel.grid = element_blank(), 
        panel.border = element_rect(fill = "transparent"))

# FSp
FSp_PC12_100 <- ggplot(FUSE.Y100, aes(x = PC1, y = PC2, color = FSp_std)) +
  geom_point(shape = 16, size = 4) +
  labs(x = "PCoA1", y = "PCoA2") +
  color_scale_FSp +
  theme_minimal() +
  theme_bw() +
  theme(legend.position = "none") +
  theme(axis.text.y = element_text(size = 8, color = "black"),  
        axis.title = element_text(size = 12), panel.background = element_rect(fill = "white")) +
  theme(panel.grid = element_blank(), 
        panel.border = element_rect(fill = "transparent"))

## 200 years
sp_dist_sharks_Y200 <- mFD::funct.dist(
  sp_tr         = sharks_traits_Y200,
  tr_cat        = sharks_traits_cat,
  metric        = "gower",
  scale_euclid  = "scale_center",
  ordinal_var   = "classic",
  weight_type   = "equal",
  stop_if_NA    = TRUE)

# Assess quality
fspaces_quality_sharks_Y200 <- mFD::quality.fspaces(
  sp_dist             = sp_dist_sharks_Y200,
  maxdim_pcoa         = 10,
  deviation_weighting = c("absolute","squared"),
  fdist_scaling       = FALSE,
  fdendro             = "average")

# Return coordinates of each axis
sp_faxes_coord_sharks_Y200 <- fspaces_quality_sharks_Y200$"details_fspaces"$"sp_pc_coord"

# Calculate per-species FUn and FSp
FUSE.Y200 <- fuse(sp_dist        = sp_dist_sharks_Y200, 
                  sp_faxes_coord = as.matrix(sp_faxes_coord_sharks_Y200), 
                  nb_NN          = 5,  
                  GE             = sharks_traits_Y200.iucn$iucn,
                  standGE        = TRUE)

FUSE.Y200$PC1 <- sp_faxes_coord_sharks_Y200[, "PC1"]
FUSE.Y200$PC2 <- sp_faxes_coord_sharks_Y200[, "PC2"]
FUSE.Y200$PC3 <- sp_faxes_coord_sharks_Y200[, "PC3"]
FUSE.Y200$PC4 <- sp_faxes_coord_sharks_Y200[, "PC4"]

# Form legends
color_scale_FUn <- scale_color_viridis(name = "FUn_std", guide = "none")
color_scale_FSp <- scale_color_viridis(name = "FSp_std", guide = "none")

#FUn
FUn_PC12_200 <- ggplot(FUSE.Y200, aes(x = PC1, y = PC2, color = FUn_std)) +
  geom_point(shape = 16, size = 4) +
  labs(x = "PCoA1", y = "PCoA2") +
  color_scale_FUn +
  theme_minimal() +
  theme_bw() +
  theme(legend.position = "none") +
  theme(axis.text.y = element_text(size = 8, color = "black"),  
        axis.title = element_text(size = 12), panel.background = element_rect(fill = "white")) +
  theme(panel.grid = element_blank(), 
        panel.border = element_rect(fill = "transparent"))

# FSp
FSp_PC12_200 <- ggplot(FUSE.Y200, aes(x = PC1, y = PC2, color = FSp_std)) +
  geom_point(shape = 16, size = 4) +
  labs(x = "PCoA1", y = "PCoA2") +
  color_scale_FSp +
  theme_minimal() +
  theme_bw() +
  theme(legend.position = "none") +
  theme(axis.text.y = element_text(size = 8, color = "black"),  
        axis.title = element_text(size = 12), panel.background = element_rect(fill = "white")) +
  theme(panel.grid = element_blank(), 
        panel.border = element_rect(fill = "transparent"))

## 300 years
sp_dist_sharks_Y300 <- mFD::funct.dist(
  sp_tr         = sharks_traits_Y300,
  tr_cat        = sharks_traits_cat,
  metric        = "gower",
  scale_euclid  = "scale_center",
  ordinal_var   = "classic",
  weight_type   = "equal",
  stop_if_NA    = TRUE)

# Assess quality
fspaces_quality_sharks_Y300 <- mFD::quality.fspaces(
  sp_dist             = sp_dist_sharks_Y300,
  maxdim_pcoa         = 10,
  deviation_weighting = c("absolute","squared"),
  fdist_scaling       = FALSE,
  fdendro             = "average")

# Return coordinates of each axis
sp_faxes_coord_sharks_Y300 <- fspaces_quality_sharks_Y300$"details_fspaces"$"sp_pc_coord"

# Calculate per-species FUn and FSp
FUSE.Y300 <- fuse(sp_dist        = sp_dist_sharks_Y300, 
                  sp_faxes_coord = as.matrix(sp_faxes_coord_sharks_Y300), 
                  nb_NN          = 5,  
                  GE             = sharks_traits_Y300.iucn$iucn,
                  standGE        = TRUE)

FUSE.Y300$PC1 <- sp_faxes_coord_sharks_Y300[, "PC1"]
FUSE.Y300$PC2 <- sp_faxes_coord_sharks_Y300[, "PC2"]
FUSE.Y300$PC3 <- sp_faxes_coord_sharks_Y300[, "PC3"]
FUSE.Y300$PC4 <- sp_faxes_coord_sharks_Y300[, "PC4"]

# Form legends
color_scale_FUn <- scale_color_viridis(name = "FUn_std", guide = "none")
color_scale_FSp <- scale_color_viridis(name = "FSp_std", guide = "none")

#FUn
FUn_PC12_300 <- ggplot(FUSE.Y300, aes(x = PC1, y = PC2, color = FUn_std)) +
  geom_point(shape = 16, size = 4) +
  labs(x = "PCoA1", y = "PCoA2") +
  color_scale_FUn +
  theme_minimal() +
  theme_bw() +
  theme(legend.position = "none") +
  theme(axis.text.y = element_text(size = 8, color = "black"),  
        axis.title = element_text(size = 12), panel.background = element_rect(fill = "white")) +
  theme(panel.grid = element_blank(), 
        panel.border = element_rect(fill = "transparent"))

# FSp
FSp_PC12_300 <- ggplot(FUSE.Y300, aes(x = PC1, y = PC2, color = FSp_std)) +
  geom_point(shape = 16, size = 4) +
  labs(x = "PCoA1", y = "PCoA2") +
  color_scale_FSp +
  theme_minimal() +
  theme_bw() +
  theme(legend.position = "none") +
  theme(axis.text.y = element_text(size = 8, color = "black"),  
        axis.title = element_text(size = 12), panel.background = element_rect(fill = "white")) +
  theme(panel.grid = element_blank(), 
        panel.border = element_rect(fill = "transparent"))

## 400 years
sp_dist_sharks_Y400 <- mFD::funct.dist(
  sp_tr         = sharks_traits_Y400,
  tr_cat        = sharks_traits_cat,
  metric        = "gower",
  scale_euclid  = "scale_center",
  ordinal_var   = "classic",
  weight_type   = "equal",
  stop_if_NA    = TRUE)

# Assess quality
fspaces_quality_sharks_Y400 <- mFD::quality.fspaces(
  sp_dist             = sp_dist_sharks_Y400,
  maxdim_pcoa         = 10,
  deviation_weighting = c("absolute","squared"),
  fdist_scaling       = FALSE,
  fdendro             = "average")

# Return coordinates of each axis
sp_faxes_coord_sharks_Y400 <- fspaces_quality_sharks_Y400$"details_fspaces"$"sp_pc_coord"

# Calculate per-species FUn and FSp
FUSE.Y400 <- fuse(sp_dist        = sp_dist_sharks_Y400, 
                  sp_faxes_coord = as.matrix(sp_faxes_coord_sharks_Y400), 
                  nb_NN          = 5,  
                  GE             = sharks_traits_Y400.iucn$iucn,
                  standGE        = TRUE)

FUSE.Y400$PC1 <- sp_faxes_coord_sharks_Y400[, "PC1"]
FUSE.Y400$PC2 <- sp_faxes_coord_sharks_Y400[, "PC2"]
FUSE.Y400$PC3 <- sp_faxes_coord_sharks_Y400[, "PC3"]
FUSE.Y400$PC4 <- sp_faxes_coord_sharks_Y400[, "PC4"]

# Form legends
color_scale_FUn <- scale_color_viridis(name = "FUn_std", guide = "none")
color_scale_FSp <- scale_color_viridis(name = "FSp_std", guide = "none")

#FUn
FUn_PC12_400 <- ggplot(FUSE.Y400, aes(x = PC1, y = PC2, color = FUn_std)) +
  geom_point(shape = 16, size = 4) +
  labs(x = "PCoA1", y = "PCoA2") +
  color_scale_FUn +
  theme_minimal() +
  theme_bw() +
  theme(legend.position = "none") +
  theme(axis.text.y = element_text(size = 8, color = "black"),  
        axis.title = element_text(size = 12), panel.background = element_rect(fill = "white")) +
  theme(panel.grid = element_blank(), 
        panel.border = element_rect(fill = "transparent"))

# FSp
FSp_PC12_400 <- ggplot(FUSE.Y400, aes(x = PC1, y = PC2, color = FSp_std)) +
  geom_point(shape = 16, size = 4) +
  labs(x = "PCoA1", y = "PCoA2") +
  color_scale_FSp +
  theme_minimal() +
  theme_bw() +
  theme(legend.position = "none") +
  theme(axis.text.y = element_text(size = 8, color = "black"),  
        axis.title = element_text(size = 12), panel.background = element_rect(fill = "white")) +
  theme(panel.grid = element_blank(), 
        panel.border = element_rect(fill = "transparent"))

## 500 years
sp_dist_sharks_Y500 <- mFD::funct.dist(
  sp_tr         = sharks_traits_Y500,
  tr_cat        = sharks_traits_cat,
  metric        = "gower",
  scale_euclid  = "scale_center",
  ordinal_var   = "classic",
  weight_type   = "equal",
  stop_if_NA    = TRUE)

# Assess quality
fspaces_quality_sharks_Y500 <- mFD::quality.fspaces(
  sp_dist             = sp_dist_sharks_Y500,
  maxdim_pcoa         = 10,
  deviation_weighting = c("absolute","squared"),
  fdist_scaling       = FALSE,
  fdendro             = "average")

# Return coordinates of each axis
sp_faxes_coord_sharks_Y500 <- fspaces_quality_sharks_Y500$"details_fspaces"$"sp_pc_coord"

# Calculate per-species FUn and FSp
FUSE.Y500 <- fuse(sp_dist        = sp_dist_sharks_Y500, 
                  sp_faxes_coord = as.matrix(sp_faxes_coord_sharks_Y500), 
                  nb_NN          = 5,  
                  GE             = sharks_traits_Y500.iucn$iucn,
                  standGE        = TRUE)

FUSE.Y500$PC1 <- sp_faxes_coord_sharks_Y500[, "PC1"]
FUSE.Y500$PC2 <- sp_faxes_coord_sharks_Y500[, "PC2"]
FUSE.Y500$PC3 <- sp_faxes_coord_sharks_Y500[, "PC3"]
FUSE.Y500$PC4 <- sp_faxes_coord_sharks_Y500[, "PC4"]

# Form legends
color_scale_FUn <- scale_color_viridis(name = "FUn_std", guide = "none")
color_scale_FSp <- scale_color_viridis(name = "FSp_std", guide = "none")

#FUn
FUn_PC12_500 <- ggplot(FUSE.Y500, aes(x = PC1, y = PC2, color = FUn_std)) +
  geom_point(shape = 16, size = 4) +
  labs(x = "PCoA1", y = "PCoA2") +
  color_scale_FUn +
  theme_minimal() +
  theme_bw() +
  theme(legend.position = "none") +
  theme(axis.text.y = element_text(size = 8, color = "black"),  
        axis.title = element_text(size = 12), panel.background = element_rect(fill = "white")) +
  theme(panel.grid = element_blank(), 
        panel.border = element_rect(fill = "transparent"))

# FSp
FSp_PC12_500 <- ggplot(FUSE.Y500, aes(x = PC1, y = PC2, color = FSp_std)) +
  geom_point(shape = 16, size = 4) +
  labs(x = "PCoA1", y = "PCoA2") +
  color_scale_FSp +
  theme_minimal() +
  theme_bw() +
  theme(legend.position = "none") +
  theme(axis.text.y = element_text(size = 8, color = "black"),  
        axis.title = element_text(size = 12), panel.background = element_rect(fill = "white")) +
  theme(panel.grid = element_blank(), 
        panel.border = element_rect(fill = "transparent"))

# Plot future scores - Produces Figure S6
FUn_Sp_scores <- plot_grid(FUn_PC12_2100,FSp_PC12_2100,
                           FUn_PC12_100,FSp_PC12_100,
                           FUn_PC12_200,FSp_PC12_200,
                           FUn_PC12_300,FSp_PC12_300,
                           FUn_PC12_400,FSp_PC12_400,
                           FUn_PC12_500,FSp_PC12_500,
                           labels = c("(a)","(b)","(c)","(d)","(e)","(f)","(g)","(h)",
                                      "(i)","(j)","(k)","(l)"),
                           label_size = 12,align = "hv", label_fontface = "bold", hjust = -0.15,  nrow=7)
