###############################################################################################################
# 03. Functional space
## This R code creates the functional space & assesses functional extent of AquaMaps
## it produces Figures S1-S4 and Tables S1-S3
###############################################################################################################

# Import packages
library(readxl)
library(tidyverse)
library(tibble)
library(scales)
library(cowplot)
library(doBy)
library(reshape2)
library(RColorBrewer)
library(mFD)
library(ggsci)
library(viridis)
library(data.table)
library(foreach)
library(parallel)
library(doParallel)
library(funrar)
library(ggplotify)

## Load Aquamaps and isolate species list
AqMap <- read_csv("C:/Users/2022207/Dropbox/Jack's PhD/Chapter 3. Future shark FD/Data/Aquamaps data/Elasmo_current.csv")

AqMap <- AqMap %>% 
  mutate(sci_name = paste(Genus, Species, sep = " "))

AqMapSpp <- AqMap %>% 
  select(sci_name) %>% 
  unique()

# Load data
data <- read_xlsx("C:/Users/2022207/Dropbox/Jack's PhD/Chapter 3. Future shark FD/Data/Trait data.xlsx", 4)

# Make all relevant traits (terrestriality, thermo, feeding), max.length & iucn numeric & the rest nominal
data$max.length <- as.numeric(data$max.length)
data$habitat <- as.factor(data$habitat)
data$vertical <- as.factor(data$vertical)
data$terrestriality <- ordered(data$terrestriality, 
                               levels=c("no","brackish","freshwater"))
data$thermo <- ordered(data$thermo, levels=c("meso","ecto"))
data$feeding <- ordered(data$feeding, levels=c("macropredator","filter feeder"))
data$diet <- as.factor(data$diet)
data$iucn <- as.numeric(data$iucn)

# Save data to load into future codes
save(data,file="C:/Users/2022207/Dropbox/Jack's PhD/Chapter 3. Future shark FD/Data/R data/Data.RData")

# Occurrence data
Occ <- data %>% 
  select(Species)

Occ$Present <- 1

# Make 2 new columns based on AquaMaps
Occ <- Occ %>%
  mutate(AquaMaps = as.integer(Species %in% AqMapSpp$sci_name)) %>% 
  mutate(AquaMapsAb = ifelse(Species %in% AqMapSpp$sci_name, 0, 1))

baskets <- Occ %>% 
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

## Functional diversity analyses
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

# summarise dataset
asb_sp_sharks_summ <- mFD::asb.sp.summary(asb_sp_w = baskets_sharks_weights)

# retrieve species occurrences for all assemblages:
asb_sp_sharks_occ <- asb_sp_sharks_summ$asb_sp_occ

# Construct trait distance matrix using species
sp_dist_sharks <- mFD::funct.dist(
  sp_tr         = sharks_traits,
  tr_cat        = sharks_traits_cat,
  metric        = "gower",
  scale_euclid  = "scale_center",
  ordinal_var   = "classic",
  weight_type   = "equal",
  stop_if_NA    = TRUE)

# Assess quality
fspaces_quality_sharks <- mFD::quality.fspaces(
  sp_dist             = sp_dist_sharks,
  maxdim_pcoa         = 10,
  deviation_weighting = c("absolute","squared"),
  fdist_scaling       = FALSE,
  fdendro             = "average")

# Check pocket
round(fspaces_quality_sharks$"quality_fspaces", 5)

# Plot quality - produces Figure S1
mFD::quality.fspaces.plot(
  fspaces_quality            = fspaces_quality_sharks,
  quality_metric             = "mad",
  fspaces_plot               = c("tree_average", "pcoa_2d", "pcoa_3d", 
                                 "pcoa_4d", "pcoa_5d"),
  name_file                  = NULL,
  range_dist                 = NULL,
  range_dev                  = NULL,
  range_qdev                 = NULL,
  gradient_deviation         = c(neg = "darkblue", nul = "grey80", pos = "darkred"),
  gradient_deviation_quality = c(low = "yellow", high = "red"),
  x_lab                      = "Trait-based distance")

# Check % of inertia accounted for by axes
eig <- fspaces_quality_sharks$details_fspaces$pc_eigenvalues %>% 
  as.data.frame()

eig$variance <- (eig$Eigenvalues/70.332073)*100

# Return coordinates of each axis
sp_faxes_coord_sharks <- fspaces_quality_sharks$"details_fspaces"$"sp_pc_coord"

# Correlate each axis against traits
sharks_tr_faxes <- mFD::traits.faxes.cor(
  sp_tr          = sharks_traits, 
  sp_faxes_coord = sp_faxes_coord_sharks[ , c("PC1", "PC2", "PC3", "PC4")], 
  plot           = TRUE)

# Print traits with significant effect - produces Table S2
sharks_tr_faxes$"tr_faxes_stat"[which(sharks_tr_faxes$"tr_faxes_stat"$"p.value" < 0.05), ]

# Return plots - produces Figure S2
sharks_tr_faxes$"tr_faxes_plot"

# Calculate FD metrics
alpha_fd_indices_sharks <- mFD::alpha.fd.multidim(
  sp_faxes_coord   = sp_faxes_coord_sharks[ , c("PC1", "PC2", "PC3", "PC4")],
  asb_sp_w         = baskets_sharks_weights,
  ind_vect         = c("fric", "fori", "fspe"),
  scaling          = TRUE,
  check_input      = TRUE,
  details_returned = TRUE)

fd_ind_values_sharks <- alpha_fd_indices_sharks$"functional_diversity_indices"

# Form dataframe
FDindices_taxon <- as.data.frame(fd_ind_values_sharks) %>% 
  tibble::rownames_to_column("Scenario")

# Plot functional space & identify highest & lowest contributors per axis
# Compute the range of functional axes:
range_sp_coord  <- range(sp_faxes_coord_sharks)

# Based on the range of species coordinates values, compute a nice range for functional axes:
range_faxes <- range_sp_coord +
  c(-1, 1) * (range_sp_coord[2] - range_sp_coord[1]) * 0.05
range_faxes

limits_PC1 <- c(-0.4, 0.4)
limits_PC2 <- c(-0.3, 0.5)
limits_PC3 <- c(-0.3, 0.4)
limits_PC4 <- c(-0.6, 0.4)

# get species coordinates along axes:
sp_faxes_coord_12 <- sp_faxes_coord_sharks[, c("PC1", "PC2")]
sp_faxes_coord_13 <- sp_faxes_coord_sharks[, c("PC1", "PC3")]
sp_faxes_coord_14 <- sp_faxes_coord_sharks[, c("PC1", "PC4")]
sp_faxes_coord_23 <- sp_faxes_coord_sharks[, c("PC2", "PC3")]
sp_faxes_coord_24 <- sp_faxes_coord_sharks[, c("PC2", "PC4")]
sp_faxes_coord_34 <- sp_faxes_coord_sharks[, c("PC3", "PC4")]

# Plot backgrounds:
## Done for PC1-PC2 & PC3-PC4 - can be replicated across any axis pairings
plot_12 <- mFD::background.plot(range_faxes = range_faxes, 
                                faxes_nm = c("PC1", "PC2"),
                                color_bg = "white")
plot_34 <- mFD::background.plot(range_faxes = range_faxes, 
                                faxes_nm = c("PC3", "PC4"),
                                color_bg = "white")

# Filter species to baskets based on Aquamaps
## PC1-PC2 Aquamaps baskets
sp_filter_AquaMaps12 <- mFD::sp.filter(asb_nm = c("AquaMaps"),
                                    sp_faxes_coord = sp_faxes_coord_12,
                                    asb_sp_w = baskets_sharks_weights)

sp_faxes_coord_AquaMaps12 <- sp_filter_AquaMaps12$`species coordinates`

sp_filter_AquaMapsAb12 <- mFD::sp.filter(asb_nm = c("AquaMapsAb"),
                                       sp_faxes_coord = sp_faxes_coord_12,
                                       asb_sp_w = baskets_sharks_weights)

sp_faxes_coord_AquaMapsAb12 <- sp_filter_AquaMapsAb12$`species coordinates`

## PC3-PC4 Aquamaps baskets
sp_filter_AquaMaps34 <- mFD::sp.filter(asb_nm = c("AquaMaps"),
                                       sp_faxes_coord = sp_faxes_coord_34,
                                       asb_sp_w = baskets_sharks_weights)

sp_faxes_coord_AquaMaps34 <- sp_filter_AquaMaps34$`species coordinates`

sp_filter_AquaMapsAb34 <- mFD::sp.filter(asb_nm = c("AquaMapsAb"),
                                         sp_faxes_coord = sp_faxes_coord_34,
                                         asb_sp_w = baskets_sharks_weights)

sp_faxes_coord_AquaMapsAb34 <- sp_filter_AquaMapsAb34$`species coordinates`

# Retrieve vertices coordinates along the functional axes:
vert_Aquamaps_12 <- mFD::vertices(sp_faxes_coord = sp_faxes_coord_AquaMaps12,  
                        order_2D = FALSE, 
                        check_input = TRUE)
vert_AquamapsAb_12 <- mFD::vertices(sp_faxes_coord = sp_faxes_coord_AquaMapsAb12,  
                                  order_2D = FALSE, 
                                  check_input = TRUE)
vert_Aquamaps_34 <- mFD::vertices(sp_faxes_coord = sp_faxes_coord_AquaMaps34,  
                                  order_2D = FALSE, 
                                  check_input = TRUE)
vert_AquamapsAb_34 <- mFD::vertices(sp_faxes_coord = sp_faxes_coord_AquaMapsAb34,  
                                    order_2D = FALSE, 
                                    check_input = TRUE)

# Plot spaces
plot_sp_vert12 <- mFD::fric.plot(ggplot_bg = plot_12,
                                 asb_sp_coord2D = list("AquaMaps" = sp_faxes_coord_AquaMaps12,
                                                   "AquaMapsAb" = sp_faxes_coord_AquaMapsAb12),
                                 asb_vertices_nD = list("AquaMaps" = vert_Aquamaps_12,
                                                    "AquaMapsAb" = vert_AquamapsAb_12),
                                 plot_sp = TRUE,
                                 color_ch = NA,
                                 fill_ch = c("AquaMaps" = "#91bfdb",
                                             "AquaMapsAb" = "#fc8d59"),
                                 alpha_ch =  c("AquaMaps" = 0.4,
                                               "AquaMapsAb" = 0.4),
                                 shape_sp = c("AquaMaps" = 21,
                                              "AquaMapsAb" = 24),
                                 size_sp = c("AquaMaps" = 1.5,
                                             "AquaMapsAb" = 1.5),
                                 color_sp = c("AquaMaps" = "#91bfdb",
                                              "AquaMapsAb" = "#fc8d59"),
                                 fill_sp = c("AquaMaps" = "#91bfdb",
                                             "AquaMapsAb" = "#fc8d59"),
                                 shape_vert = c("AquaMaps" = 21,
                                                "AquaMapsAb" = 24),
                                 size_vert = c("AquaMaps" = 1.5,
                                               "AquaMapsAb" = 1.5),
                                 color_vert = c("AquaMaps" = "#91bfdb",
                                                "AquaMapsAb" = "#fc8d59"),
                                 fill_vert = c("AquaMaps" = "#91bfdb",
                                               "AquaMapsAb" = "#fc8d59"))+
  labs(x = "PCoA1", y = "PCoA2")+
  coord_fixed(ratio = 1, xlim = limits_PC1, ylim = limits_PC2)+
  annotate("point", x = -0.3346027, y = 0.13891530, color = "black", size = 3) +
  annotate("point", x = 0.3753394, y = 0.089805027, color = "black", size = 3) +
  annotate("point", x = -0.06238338, y = -0.1725660, color = "black", size = 3) +
  annotate("point", x = -0.023538397, y = 0.4653613, color = "black", size = 3)

plot_sp_vert34 <- mFD::fric.plot(ggplot_bg = plot_34,
                                 asb_sp_coord2D = list("AquaMaps" = sp_faxes_coord_AquaMaps34,
                                                       "AquaMapsAb" = sp_faxes_coord_AquaMapsAb34),
                                 asb_vertices_nD = list("AquaMaps" = vert_Aquamaps_34,
                                                        "AquaMapsAb" = vert_AquamapsAb_34),
                                 plot_sp = TRUE,
                                 color_ch = NA,
                                 fill_ch = c("AquaMaps" = "#91bfdb",
                                             "AquaMapsAb" = "#fc8d59"),
                                 alpha_ch =  c("AquaMaps" = 0.4,
                                               "AquaMapsAb" = 0.4),
                                 shape_sp = c("AquaMaps" = 21,
                                              "AquaMapsAb" = 24),
                                 size_sp = c("AquaMaps" = 1.5,
                                             "AquaMapsAb" = 1.5),
                                 color_sp = c("AquaMaps" = "#91bfdb",
                                              "AquaMapsAb" = "#fc8d59"),
                                 fill_sp = c("AquaMaps" = "#91bfdb",
                                             "AquaMapsAb" = "#fc8d59"),
                                 shape_vert = c("AquaMaps" = 21,
                                                "AquaMapsAb" = 24),
                                 size_vert = c("AquaMaps" = 1.5,
                                               "AquaMapsAb" = 1.5),
                                 color_vert = c("AquaMaps" = "#91bfdb",
                                                "AquaMapsAb" = "#fc8d59"),
                                 fill_vert = c("AquaMaps" = "#91bfdb",
                                               "AquaMapsAb" = "#fc8d59"))+
  labs(x = "PCoA3", y = "PCoA4")+
  coord_fixed(ratio = 1, xlim = limits_PC3, ylim = limits_PC4)+
  annotate("point", x = -0.23857812, y = -0.044301832, color = "black", size = 3) +
  annotate("point", x = 0.3845090, y = -0.57501448, color = "black", size = 3) +
  annotate("point", x = 0.3845090, y = -0.57501448, color = "black", size = 3) +
  annotate("point", x = 0.26992461, y = 0.28880660, color = "black", size = 3)

# Plot space 
Fig_S3a <- plot_grid(plot_sp_vert12,plot_sp_vert13,plot_sp_vert14,
                    plot_sp_vert23,plot_sp_vert24,plot_sp_vert34,
          labels = c("(a)","(b)","(c)","(d)","(e)","(f)"),
          label_size = 12,align = "hv", label_fontface = "bold", hjust = -0.15,  nrow=2)

# FUSE calculation
FUSE.present <- fuse(sp_dist        = sp_dist_sharks, 
                    sp_faxes_coord = as.matrix(sp_faxes_coord_sharks), 
                    nb_NN          = 5,  
                    GE             = data$iucn,
                    standGE        = TRUE)

# Status colours
status_colors <- c("LC" = "#008000", "NT" = "#ADFF2F", "VU" = "#FFFF00", "EN" = "#FFA500", "CR" = "#FF0000")

# Set up df
FUSE.present$IUCN <- data$iucn
FUSE.present$IUCN <- as.factor(FUSE.present$IUCN)
FUSE_sum <- FUSE.present %>% 
  mutate(IUCN = recode(IUCN, 
                       "0" = "LC", 
                       "1" = "NT",
                       "2" = "VU", 
                       "3" = "EN", 
                       "4" = "CR"))
FUSE_sum$IUCN <- ordered(FUSE_sum$IUCN, levels=c("LC", "NT", "VU", "EN", "CR"))
FUSE_sum$AquaMaps <- Occ$AquaMaps
FUSE_sum <- FUSE_sum %>% 
  rownames_to_column(var = "Species")

# Top 25 contributors
FUSE_rank <- FUSE_sum %>%
  select(Species, FUSE, IUCN, AquaMaps) %>%
  arrange(desc(FUSE)) %>%
  head(25)
FUSE_rank$Species <- factor(FUSE_rank$Species, levels = FUSE_rank$Species[order(FUSE_rank$FUSE, decreasing = TRUE)])

# What % of top 20, 50 & 100 FUSE species are present in AquaMaps? - Produces Table S3
FUSE_rank_20 <- FUSE_sum %>%
  select(Species, FUSE, IUCN, AquaMaps) %>%
  arrange(desc(FUSE)) %>%
  head(20)

percentage_FUSE_20 <- FUSE_rank_20 %>%
  filter(AquaMaps == 1) %>%
  nrow() / nrow(FUSE_rank_20) * 100     # 18/20 = 90%

FUSE_rank_50 <- FUSE_sum %>%
  select(Species, FUSE, IUCN, AquaMaps) %>%
  arrange(desc(FUSE)) %>%
  head(50)

percentage_FUSE_50 <- FUSE_rank_50 %>%
  filter(AquaMaps == 1) %>%
  nrow() / nrow(FUSE_rank_50) * 100     # 41/50 = 82%

FUSE_rank_100 <- FUSE_sum %>%
  select(Species, FUSE, IUCN, AquaMaps) %>%
  arrange(desc(FUSE)) %>%
  head(100)

percentage_FUSE_100 <- FUSE_rank_100 %>%
  filter(AquaMaps == 1) %>%
  nrow() / nrow(FUSE_rank_100) * 100     # 72/100 = 72%

# Create lollipop plot for FUSE
plot_FUSE_rank_lollipop <- ggplot(FUSE_rank, aes(x = Species, y = FUSE, fill = IUCN, color = "darkgrey", shape = factor(AquaMaps))) +
  geom_segment(aes(xend = Species, yend = 0), size = 1.5) +
  geom_point(size = 3) +
  coord_flip() +
  labs(x = "", y = "FUSE") +
  scale_fill_manual(values = status_colors, guide = guide_legend(title = "Status")) +
  scale_color_manual(values = status_colors, guide = guide_legend(title = "Status")) +
  scale_shape_manual(values = c("1" = 21, "0" = 24), guide = guide_legend(title = "AquaMaps")) + 
  theme_minimal() +
  theme_bw() +
  theme(legend.position = "right") +
  theme(axis.text.y = element_text(size = 8, color = "black", face = "italic"),  # Set y-axis labels to italic
        axis.title = element_text(size = 12), panel.background = element_rect(fill = "white")) +
  theme(panel.grid = element_blank(), 
        panel.border = element_rect(fill = "transparent"))


# Plot Figure S3
Fig_S3 <- plot_grid(Fig_S3a,plot_FUSE_rank_lollipop,
                    labels = c("","(g)"),
                    label_size = 12,align = "hv", label_fontface = "bold", hjust = -0.15,  nrow=1)

# Plot FUn and FSp scores
## Present - produces Figure S6a
FUSE.present$PC1 <- sp_faxes_coord_sharks[, "PC1"]
FUSE.present$PC2 <- sp_faxes_coord_sharks[, "PC2"]
FUSE.present$PC3 <- sp_faxes_coord_sharks[, "PC3"]
FUSE.present$PC4 <- sp_faxes_coord_sharks[, "PC4"]

color_scale_FUn <- scale_color_viridis(name = "FUn_std", guide = "none")
color_scale_FSp <- scale_color_viridis(name = "FSp_std", guide = "none")

#FUn - done for PC1-PC2 but can be replicated across axis pairings
FUn_legend <- ggplot(FUSE.present, aes(x = 0, y = FUn_std, color = FUn_std)) +
  geom_point() +
  scale_color_viridis(name = "FUn_std") +
  theme_void() +
  theme(legend.position = "right")

FUNleg <- get_legend(FUn_legend)

legend_FUn <- as.ggplot(FUNleg)

FUn_PC12 <- ggplot(FUSE.present, aes(x = PC1, y = PC2, color = FUn_std)) +
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

# FSp  - done for PC1-PC2 but can be replicated across axis pairings
FSp_legend <- ggplot(FUSE.present, aes(x = 0, y = FSp_std, color = FSp_std)) +
  geom_point() +
  scale_color_viridis(name = "FSp_std") +
  theme_void() +
  theme(legend.position = "right")

FSPleg <- get_legend(FSp_legend)

legend_FSp <- as.ggplot(FSPleg)

FSp_PC12 <- ggplot(FUSE.present, aes(x = PC1, y = PC2, color = FSp_std)) +
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

# Plot overall functional space 
## Done for PC1-PC2 & PC3-PC4 - can be replicated across any axis pairings
# Plot background:
plot_12 <- mFD::background.plot(range_faxes = range_faxes, 
                                faxes_nm = c("PC1", "PC2"),
                                color_bg = "white")
plot_34 <- mFD::background.plot(range_faxes = range_faxes, 
                                faxes_nm = c("PC3", "PC4"),
                                color_bg = "white")

# Retrieve vertices coordinates along the functional axes:
vert12 <- mFD::vertices(sp_faxes_coord = sp_faxes_coord_12,  
                        order_2D = FALSE, 
                        check_input = TRUE)
vert34 <- mFD::vertices(sp_faxes_coord = sp_faxes_coord_34,  
                        order_2D = FALSE, 
                        check_input = TRUE)

# Plot vertices species
plot_sp_vert12 <- mFD::pool.plot(ggplot_bg = plot_12,
                                 sp_coord2D = sp_faxes_coord_12,
                                 vertices_nD = vert12,
                                 plot_pool = TRUE,
                                 color_pool = "darkgrey",
                                 fill_pool = "darkgrey",
                                 alpha_ch =  0.8,
                                 color_ch = "black",
                                 fill_ch = "white",
                                 shape_pool = 21,
                                 size_pool = 1.5,
                                 shape_vert = 16,
                                 size_vert = 1,
                                 color_vert = "darkgrey",
                                 fill_vert = "darkgrey")+
  coord_fixed(ratio = 1, xlim = limits_PC1, ylim = limits_PC2)+
  annotate("point", x = -0.3346027, y = 0.13891530, color = "black", size = 3) +
  annotate("point", x = 0.3753394, y = 0.089805027, color = "black", size = 3) +
  annotate("point", x = -0.06238338, y = -0.1725660, color = "black", size = 3) +
  annotate("point", x = -0.023538397, y = 0.4653613, color = "black", size = 3)+
  labs(x = "PCoA1", y = "PCoA2")

plot_sp_vert34 <- mFD::pool.plot(ggplot_bg = plot_34,
                                 sp_coord2D = sp_faxes_coord_34,
                                 vertices_nD = vert34,
                                 plot_pool = TRUE,
                                 color_pool = "darkgrey",
                                 fill_pool = "darkgrey",
                                 alpha_ch =  1.5,
                                 color_ch = "black",
                                 fill_ch = "white",
                                 shape_pool = 21,
                                 size_pool = 1.5,
                                 shape_vert = 16,
                                 size_vert = 1,
                                 color_vert = "darkgrey",
                                 fill_vert = "darkgrey")+
  labs(x = "PCoA3", y = "PCoA4")+
  coord_fixed(ratio = 1, xlim = limits_PC3, ylim = limits_PC4)+
  annotate("point", x = -0.23857812, y = -0.044301832, color = "black", size = 3) +
  annotate("point", x = 0.3845090, y = -0.57501448, color = "black", size = 3) +
  annotate("point", x = 0.3845090, y = -0.57501448, color = "black", size = 3) +
  annotate("point", x = 0.26992461, y = 0.28880660, color = "black", size = 3)

# Produce full space - produces Figure S4
Space <- plot_grid(plot_sp_vert12,plot_sp_vert34,
                    labels = c("(c)","(d)"),
                    label_size = 12,align = "hv", label_fontface = "bold", hjust = -0.15,  nrow=1)