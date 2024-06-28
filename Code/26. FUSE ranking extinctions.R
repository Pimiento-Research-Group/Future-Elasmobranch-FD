###############################################################################################################
# 26. Core list
## This R code identifies the extinct species most important to avoid forecast changes in 2100
## it produces Table S10
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
library(plotly)
library(gridExtra)
library(ggpubr)

# Load data
load(file="C:/Users/2022207/Dropbox/Jack's PhD/Chapter 3. Future shark FD/Data/R data/Data.RData")

# Load mode IUCN_sim results
load(file = "C:/Users/2022207/Dropbox/Jack's PhD/Chapter 3. Future shark FD/Data/R data/2100_mode.RData")

# Occurrence data
Occ <- data %>% 
  select(Species) 

# Join with extinction times for the current iteration
Occ <- Occ %>%
  left_join(Spp_mode, by = "Species")

Occ <- Occ %>% 
  mutate(Present = 1,
         Future_2100 = ifelse(is.na(mode), 1, ifelse(mode <= 77, 0, 1))) %>% 
  select(Species,Present,Future_2100)

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

# Form trait category matrix - weighting all equal
sharks_traits_cat <- tibble(trait_name = c("habitat","vertical","terrestriality","thermo","feeding","diet","max.length"),
                            trait_type = c("N", "N", "O", "O", "O", "N", "Q"))

## Functional diversity analyses
# Summarise dataset
sharks_traits_summ <- mFD::sp.tr.summary(
  tr_cat     = sharks_traits_cat,   
  sp_tr      = sharks_traits, 
  stop_if_NA = TRUE)

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

# Return coordinates of each axis
sp_faxes_coord_sharks <- fspaces_quality_sharks$"details_fspaces"$"sp_pc_coord"

# Calculate FD metrics
alpha_fd_indices_sharks <- mFD::alpha.fd.multidim(
  sp_faxes_coord   = sp_faxes_coord_sharks[ , c("PC1", "PC2", "PC3", "PC4")],
  asb_sp_w         = baskets_sharks_weights,
  ind_vect         = c("fric", "fori", "fspe"),
  scaling          = TRUE,
  check_input      = TRUE,
  details_returned = TRUE)

fd_ind_values_sharks <- alpha_fd_indices_sharks$"functional_diversity_indices"

## Form dataframe
FDindices_taxon <- as.data.frame(fd_ind_values_sharks) %>% 
  tibble::rownames_to_column("Scenario")

# Calculate extinction probability in 75 years (based on Davis et al. 2018)
data <- data %>%
  mutate(ext_prob_75 = recode(iucn, "0"=0.0007497188203, "1"=0.007471945181, "2"=0.075979, "3"=0.566901,"4"=0.994476))

# FUSE
FUSE.present <- fuse(sp_dist        = sp_dist_sharks, 
                     sp_faxes_coord = as.matrix(sp_faxes_coord_sharks), 
                     nb_NN          = 5,  
                     GE             = data$ext_prob_75,
                     standGE        = TRUE)

# Status colours & extinction colours (blue for extant, red for extinct)
status_colors <- c("LC" = "#008000", "NT" = "#ADFF2F", "VU" = "#FFFF00", "EN" = "#FFA500", "CR" = "#FF0000")
extinct_colours <- c("Extant" = "#4292c6", "Extinct" = "darkred")

# Set up df
FUSE.present$Ext_prob <- data$ext_prob_75
FUSE.present$Ext_prob <- as.factor(FUSE.present$Ext_prob)
FUSE_sum <- FUSE.present %>% 
  mutate(Ext_prob = recode(Ext_prob, 
                           "0.0007497188203" = "LC", 
                           "0.007471945181" = "NT",
                           "0.075979" = "VU", 
                           "0.566901" = "EN", 
                           "0.994476" = "CR"))
FUSE_sum$Ext_prob <- ordered(FUSE_sum$Ext_prob, levels=c("LC", "NT", "VU", "EN", "CR"))
FUSE_sum$stat_2100 <- Occ$Future_2100
FUSE_sum <- FUSE_sum %>% 
  mutate(stat_2100 = recode(stat_2100, 
                            "1" = "Extant", 
                            "0" = "Extinct"))
FUSE_sum <- FUSE_sum %>% 
  rownames_to_column(var = "Species")

# Top FUSE contributors
# Plot
FUSE_rank <- FUSE_sum %>%
  select(Species, FUSE, Ext_prob, stat_2100) %>%
  arrange(desc(FUSE)) %>% 
  head(50)
FUSE_rank$Species <- factor(FUSE_rank$Species, levels = FUSE_rank$Species[order(FUSE_rank$FUSE, decreasing = TRUE)])

# Create lollipop plot for FUSE
plot_FUSE_rank_lollipop <- ggplot(FUSE_rank, aes(x = Species, y = FUSE, fill = stat_2100, color = stat_2100)) +
  geom_segment(aes(xend = Species, yend = 0), size = 1.5) +
  geom_point(size = 3) +
  coord_flip() +
  labs(x = "", y = "FUSE") +
  scale_fill_manual(values = extinct_colours, guide = guide_legend(title = "Status")) +
  scale_color_manual(values = extinct_colours, guide = guide_legend(title = "Status")) +
  theme_minimal() +
  theme_bw() +
  theme(legend.position = "right") +
  theme(axis.text.y = element_text(size = 8, color = "black", face = "italic"),  
        axis.title = element_text(size = 12), panel.background = element_rect(fill = "white")) +
  theme(panel.grid = element_blank(), 
        panel.border = element_rect(fill = "transparent"))

# 90th quantile test for FUSE species extinct vs extant in 2100
## Test that extinct FUSE outliers significantly differ from the same distribution of extant sharks
# Extract FUSE values for extinct and extant species
FUSE_extinct <- FUSE_rank$FUSE[FUSE_rank$stat_2100 == "Extinct"]
FUSE_extant <- FUSE_rank$FUSE[FUSE_rank$stat_2100 == "Extant"]

# Compute 90th quantiles for each group
quantile_extinct <- quantile(FUSE_extinct, probs = 0.9)
quantile_extant <- quantile(FUSE_extant, probs = 0.9)

# Define a function to compute the difference in 90th quantiles
compute_diff <- function(x, y) {
  quantile(x, probs = 0.9) - quantile(y, probs = 0.9)
}

# Compute observed difference
observed_diff <- compute_diff(FUSE_extinct, FUSE_extant)

# Generate permutations of Status variable
n_permutations <- 5000
permuted_diffs <- replicate(n_permutations, {
  status_permuted <- sample(FUSE_rank$stat_2100)
  FUSE_extinct_permuted <- FUSE_rank$FUSE[status_permuted == "Extinct"]
  FUSE_extant_permuted <- FUSE_rank$FUSE[status_permuted == "Extant"]
  compute_diff(FUSE_extinct_permuted, FUSE_extant_permuted)
})

# Compute p-value
p_value <- mean(abs(permuted_diffs) >= abs(observed_diff))

# Output results
cat("Observed Difference:", observed_diff, "\n")
cat("P-value:", p_value, "\n")