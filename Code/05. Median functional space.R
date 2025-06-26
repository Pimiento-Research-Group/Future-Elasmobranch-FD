###############################################################################################################
# 05. Median functional space
## This R code forms future functional spaces for our extinction scenarios
## it produces Figure S2
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
library(data.table)
library(geometry)

# Load data
load(file="C:/Users/2022207/Dropbox/Jack's PhD/Chapter 3. Future shark FD/Data/R data/Data.RData")
# Load median results for spatial results
load(file = "C:/Users/2022207/Dropbox/Jack's PhD/Chapter 3. Future shark FD/Analyses/Final paper R codes/data/Crit_E/CritE_median.RData")
# Load IUCN status for AT scenario
load(file = "C:/Users/2022207/Dropbox/Jack's PhD/Chapter 3. Future shark FD/Data/R data/Species_list_IUCN.RData")

# Load and correct synonyms
Synonyms <- read_xlsx("C:/Users/2022207/Dropbox/Jack's PhD/Chapter 3. Future shark FD/Analyses/DatasetSynonyms.xlsx")

Spp_median.CritE <- Spp_median.CritE %>%
  left_join(Synonyms, by = c("Species" = "iucnsim_name")) %>%
  mutate(Species = ifelse(!is.na(Accepted_name), Accepted_name, Species)) %>%
  select(-c(Accepted_name,Marked_iucnsim,Aquamaps_name,Marked_aquamaps)) %>% 
  distinct(Species, .keep_all = TRUE)

# Source functions
source("C:/Users/2022207/Dropbox/Jack's PhD/Chapter 3. Future shark FD/Analyses/R code/Final pipeline R code/Functions/get_indicator_function 2.R")
source("C:/Users/2022207/Dropbox/Jack's PhD/Chapter 3. Future shark FD/Analyses/R code/Final pipeline R code/Functions/fonction_FRIC_Global_full.R")

# Make future functional space
# Occurrence data
Occ_Future <- data %>% 
  select(Species,iucn)

Occ_Future$Present <- 1

# Set up future scenarios
## 2100
Occ_Future <- Occ_Future %>%
  left_join(Spp_median.CritE, by = "Species") %>%
  mutate(Future = ifelse(is.na(median), 1, ifelse(median <= 77, 0, 1))) %>%
  select(Species,Present,Future,iucn)

baskets <- Occ_Future %>%
  select(Species,Present,Future) %>% 
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

# make assemblage data into a matrix & ensure all matrix values are numeric
baskets_sharks_weights <- data.matrix(baskets_sharks_weights, rownames.force = NA)
class(baskets_sharks_weights) <- "numeric"

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

# Return coordinates of each axis
sp_faxes_coord_sharks <- fspaces_quality_sharks$"details_fspaces"$"sp_pc_coord"

# Transform dataframe from presence-absence matrix in a list
test <- apply(baskets_sharks_weights,1,function(x){colnames(baskets_sharks_weights)[which(x==1)]})

# Calculate FUn
lala <-lapply(test, function(y){get_FU_region (Mat_dist=as.matrix(sp_dist_sharks), nb_NN=5, Selected_sp=y)})

# Calculate FD metrics from function
pcoa <- list(li = sp_faxes_coord_sharks)
ax <- c("PC1", "PC2", "PC3", "PC4")
coord_d <- pcoa$li[, ax]

lala1<-lapply(test, function(y){ get_FV_Sp (ax=c(1:4), pcoa=pcoa, Selected_sp=y)})

# Calculate all FD metrics per scenario
FD.metrics <- data.frame(
  Present = rbind(lala1$Present$data$RS,
                  lala1$Present$data$Fric_r,
                  lala1$Present$specialization[1],
                  lala$Present$FU_Region_Mean),
  Future = rbind(lala1$Future$data$RS,
                 lala1$Future$data$Fric_r,
                 lala1$Future$specialization[1],
                 lala$Future$FU_Region_Mean)
)

rownames(FD.metrics) <- c("sp_richn","fric","fsp","fun")
colnames(FD.metrics) <- c("Present","Future")

FDindices_taxon.CritE <- t(FD.metrics) %>% as.data.frame()

FDindices_taxon.CritE <- FDindices_taxon.CritE %>% 
  rownames_to_column(var = "Scenario")

FDindices_taxon.CritE$Scenario <- ordered(FDindices_taxon.CritE$Scenario, 
                                    levels=c("Present","Future"))

# Save future FD results for plotting
save(FDindices_taxon.CritE, file = "C:/Users/2022207/Dropbox/Jack's PhD/Chapter 3. Future shark FD/Analyses/Final paper R codes/data/Crit_E/CritE_FD.RData")

# Form future functional spaces
limits_PC1 <- c(-0.4, 0.4)
limits_PC2 <- c(-0.3, 0.5)
limits_PC3 <- c(-0.3, 0.4)
limits_PC4 <- c(-0.6, 0.3)

# mFD FRic for plotting (identical to results above)
# Calculate FD metrics
alpha_fd_indices_sharks <- mFD::alpha.fd.multidim(
  sp_faxes_coord   = sp_faxes_coord_sharks[ , c("PC1", "PC2", "PC3", "PC4")],
  asb_sp_w         = baskets_sharks_weights,
  ind_vect         = c("fric"),
  scaling          = TRUE,
  check_input      = TRUE,
  details_returned = TRUE)

fd_ind_values_sharks <- alpha_fd_indices_sharks$"functional_diversity_indices" %>% as.data.frame()

## Plot
Present_space <- mFD::alpha.multidim.plot(
  output_alpha_fd_multidim = alpha_fd_indices_sharks,
  plot_asb_nm              = "Present",
  ind_nm                   = c("fric"),
  faxes                    = NULL,
  faxes_nm                 = NULL,
  range_faxes              = c(NA, NA),
  color_bg                 = NA,
  shape_sp                 = c(pool = 21, asb1 = 21),
  size_sp                  = c(pool = 1, asb1 = 1),
  color_sp                 = c(pool = "#F8766D", asb1 = "#F8766D"),
  color_vert               = c(pool = "#F8766D", asb1 = "#F8766D"),
  fill_sp                  = c(pool = "#F8766D", asb1 = "#F8766D"), 
  fill_vert                = c(pool = "#F8766D", asb1 = "#F8766D"),
  color_ch                 = c(pool = 'grey50', asb1 = "#F8766D"),  # Set dark grey outline for convex hull
  fill_ch                  = c(pool = "white", asb1 = "#F8766D"),
  alpha_ch                 = c(pool = 0.1, asb1 = 0.1),
  shape_centroid_fdis      = c(asb1 = 22),
  shape_centroid_fdiv      = c(asb1 = 22),
  shape_centroid_fspe      = 23,
  color_centroid_fspe      = "black",
  size_sp_nm               = 3, 
  color_sp_nm              = "black",
  plot_sp_nm               = NULL,
  fontface_sp_nm           = "plain",
  save_file                = FALSE,
  check_input              = TRUE)

Present_a <- Present_space$fric$PC1_PC2 +
  labs(x = "PCoA1", y = "PCoA2")+
  coord_fixed(ratio = 1, xlim = limits_PC1, ylim = limits_PC2)

Present_b <- Present_space$fric$PC3_PC4 +
  labs(x = "PCoA3", y = "PCoA4")+
  coord_fixed(ratio = 1, xlim = limits_PC3, ylim = limits_PC4)

Future_space <- mFD::alpha.multidim.plot(
  output_alpha_fd_multidim = alpha_fd_indices_sharks,
  plot_asb_nm              = "Future",
  ind_nm                   = c("fric"),
  faxes                    = NULL,
  faxes_nm                 = NULL,
  range_faxes              = c(NA, NA),
  color_bg                 = NA,
  shape_sp                 = c(pool = 24, asb1 = 21),
  size_sp                  = c(pool = 1, asb1 = 1),
  color_sp                 = c(pool = "black", asb1 = "#F8766D"),
  color_vert               = c(pool = "black", asb1 = "#F8766D"),
  fill_sp                  = c(pool = "black", asb1 = "#F8766D"),  
  fill_vert                = c(pool = "black", asb1 = "#F8766D"),
  color_ch                 = c(pool = 'grey50', asb1 = "#F8766D"),  # Set dark grey outline for convex hull
  fill_ch                  = c(pool = "white", asb1 = "#F8766D"),
  alpha_ch                 = c(pool = 1, asb1 = 0.1),
  shape_centroid_fdis      = c(asb1 = 22),
  shape_centroid_fdiv      = c(asb1 = 22),
  shape_centroid_fspe      = 23,
  color_centroid_fspe      = "black",
  size_sp_nm               = 3, 
  color_sp_nm              = "black",
  plot_sp_nm               = NULL,
  fontface_sp_nm           = "plain",
  save_file                = FALSE,
  check_input              = TRUE)

Future_2100a <- Future_space$fric$PC1_PC2 +
  labs(x = "PCoA1", y = "PCoA2")+
  coord_fixed(ratio = 1, xlim = limits_PC1, ylim = limits_PC2)

Future_2100b <- Future_space$fric$PC3_PC4 +
  labs(x = "PCoA3", y = "PCoA4")+
  coord_fixed(ratio = 1, xlim = limits_PC3, ylim = limits_PC4)

# Plot all future spaces
Sup_Fig <- plot_grid(Present_a,Present_b,
                     Future_2100a,Future_2100b,
                     labels = c("(a)","(b)","(c)","(d)"),
                     label_size = 12,align = "hv", label_fontface = "bold", hjust = -0.15,  nrow=2)

# Create a plot to extract the legend
legend_plot <- ggplot() +
  geom_point(aes(x = 1, y = 1, shape = "extant"), color = "#F8766D", size = 3) +
  geom_point(aes(x = 1, y = 1, shape = "extinct"), color = "black", size = 3) +
  scale_shape_manual(name = NULL, values = c(extant = 16, extinct = 17), labels = c("extant", "extinct")) +
  guides(shape = guide_legend(override.aes = list(color = c("#F8766D", "black")))) +
  theme_void() +
  theme(legend.position = "top", legend.direction = "horizontal", legend.title = element_blank())

# Extract the legend
legend <- get_legend(legend_plot)

# Produce final figure - produces supplement figure
final_sup_fig <- plot_grid(legend, Sup_Fig, ncol = 1, rel_heights = c(0.1, 1))
