###############################################################################################################
# 25. Extinction GLMs
## This R code assesses ranked FUSE species against simulated extinction
## it produces Figures 4, S10 and S11 and Table S11
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
library(visreg)
library(lme4)

# Load data
load(file="C:/Users/2022207/Dropbox/Jack's PhD/Chapter 3. Future shark FD/Data/R data/Data.RData")

# Load mode IUCN_sim results
load(file = "C:/Users/2022207/Dropbox/Jack's PhD/Chapter 3. Future shark FD/Data/R data/2100_mode.RData")

# Make future functional space
# Occurrence data
Occ <- data %>% 
  select(Species)

Occ$Present <- 1

Occ_Future <- Occ %>%
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

# Form occurrences through time
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

# Reveal number of states for each trait
sharks_traits_summ$tr_summary_list

# make assemblage data into a matrix & ensure all matrix values are numeric
baskets_sharks_weights <- data.matrix(baskets_sharks_weights, rownames.force = NA)
class(baskets_sharks_weights) <- "numeric"

#summarise dataset
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

# Calculate functional uniqueness per species
Funct_Uni <- uniqueness(baskets_sharks_weights, as.matrix(sp_dist_sharks)) %>% 
  as.data.frame()
Funct_Uni$Present <- baskets$Present
Funct_Uni$Future <- baskets$Future

# Calculate functional uniqueness per scenario
Present_FUn <- mean(Funct_Uni$Ui)
Future_FUn_df <- Funct_Uni %>% 
  filter(Future == "1") %>% 
  select(Ui)
Future_FUn <- mean(Future_FUn_df$Ui)

# Add uniqueness to final FD results data frame
Fun <- data.frame(
  Fun = c(Present_FUn,Future_FUn),
  row.names = c("Present","Future")
)

FDindices_taxon$fun <- Fun$Fun
FDindices_taxon$Scenario <- ordered(FDindices_taxon$Scenario, levels=c("Present","Future"))

# FUSE calculation
FUSE.present <- fuse(sp_dist        = sp_dist_sharks, 
                     sp_faxes_coord = as.matrix(sp_faxes_coord_sharks), 
                     nb_NN          = 5,  
                     GE             = data$iucn,
                     standGE        = TRUE)

# Status colours
status_colors <- c("LC" = "#008000", "NT" = "#ADFF2F", "VU" = "#FFFF00", "EN" = "#FFA500", "CR" = "#FF0000")

# Set up df
FUSE.present$Order <- data$Order
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
FUSE_sum <- FUSE_sum %>% 
  rownames_to_column(var = "Species")

# Add scenarios
FUSE_sum$Present_ext <- as.numeric(!baskets$Present)
FUSE_sum$Future_ext <- as.numeric(!baskets$Future)

# Generealised mixed effect model including Order 
FUSE_sum$Rank_FUSE <- rank(-FUSE_sum$FUSE)
FUSE_sum$Order <- as.factor(FUSE_sum$Order)
FUSE_Future_GLMM <- glmer(Future_ext ~ Rank_FUSE + (Rank_FUSE | Order), data = FUSE_sum, family = binomial)
summary(FUSE_Future_GLMM)
visreg(FUSE_Future_GLMM, scale = "response", main = "Extinction~FUSE in 2100", 
       xlab = "FUSE (ranked)", ylab = "Extinction probability (2100)", 
       line.par= list(col = "#AA1638"))

# Plot partial pooling/random trend by Order - produces Figure 4b
v <- visreg(FUSE_Future_GLMM, "Rank_FUSE",by="Order", plot = FALSE, scale = "response")
subOrder <- sample(levels(FUSE_sum$Order),12)
vv <- subset(v, Order %in% subOrder)
plot(vv,xlab="FUSE (ranked)",ylab="Extinction probability (2100)", layout=c(4,3),
     line.par= list(col = "#AA1638"))

# Predict the extinction probabilities using the GLMM
predictions <- predict(FUSE_Future_GLMM, type = "response")

# Combine predictions with the original data to get the Order information
data_with_predictions <- data.frame(FUSE_sum, Predicted_Prob = predictions)

# Filter the data to include only the selected subOrder
filtered_data <- subset(data_with_predictions, Order %in% subOrder)

# Calculate the average extinction probability per order
avg_extinction <- aggregate(Predicted_Prob ~ Order, data = filtered_data, FUN = mean)

# plots of 12 orders - raw data
plot_visreg_for_order <- function(order_name) {
  # Filter data for the specified order
  order_data <- FUSE_sum %>% 
    filter(Order == order_name)
  
  # Check if there is enough variation in the predictor variable
  if(length(unique(order_data$Rank_FUSE)) > 1) {
    # Fit the GLM for the order
    order_glm <- glm(Future_ext ~ Rank_FUSE, data = order_data, family = binomial)
    
    # Plot using visreg
    visreg(order_glm, scale = "response", main = order_name, xlab = "FUSE (ranked)", ylab = "Extinction probability")
  } else {
    # Calculate average extinction probability for the order
    avg_ext_prob <- mean(order_data$Future_ext)
    
    # Create a plot with a horizontal line representing the average extinction probability
    plot(NULL, xlim = c(0, 1), ylim = c(0, 1), xlab = "FUSE (ranked)", ylab = "Extinction probability", main = order_name)
    abline(h = avg_ext_prob, col = "#AA1638")
  }
}

# List of orders to plot
orders <- c("Carcharhiniformes", "Lamniformes", "Orectolobiformes", "Heterodontiformes", 
            "Squaliformes", "Pristiophoriformes", "Squatiniformes", "Hexanchiformes", 
            "Myliobatiformes", "Rhinopristiformes", "Torpediniformes", "Rajiformes")

# Set up plotting area
par(mfrow = c(3, 4))

# Plot for each order with raw data - produces Figure S10
for (order in orders) {
  plot_visreg_for_order(order)
}

# Reset plotting area
par(mfrow = c(1, 1))

### 100 years in future
baskets_100 <- Occ_Future %>% 
  filter(Future == "1") %>% 
  select(Species,Future,Y100) %>% 
  remove_rownames %>% 
  column_to_rownames(var="Species")

# Transpose to form final occurrence matrix
baskets_sharks_weights_100 <- t(baskets_100)

# Species-trait matrix
traits <- data %>% 
  select(Species,Order,habitat,vertical,terrestriality,thermo,feeding,diet,max.length,iucn)

sharks_traits <- traits %>% 
  remove_rownames %>% 
  column_to_rownames(var="Species")

sharks_traits_100filter <- sharks_traits %>%
  filter(row.names(.) %in% row.names(baskets_100))

sharks_traits_100 <- sharks_traits_100filter %>% 
  select(-c(iucn,Order))

# Form trait category matrix - weighting all equal
sharks_traits_cat <- tibble(trait_name = c("habitat","vertical","terrestriality","thermo","feeding","diet","max.length"),
                            trait_type = c("N", "N", "O", "O", "O", "N", "Q"))

# Summarise dataset
sharks_traits_summ <- mFD::sp.tr.summary(
  tr_cat     = sharks_traits_cat,   
  sp_tr      = sharks_traits_100, 
  stop_if_NA = TRUE)

# Reveal number of states for each trait
sharks_traits_summ$tr_summary_list

# make assemblage data into a matrix & ensure all matrix values are numeric
baskets_sharks_weights_100 <- data.matrix(baskets_sharks_weights_100, rownames.force = NA)
class(baskets_sharks_weights_100) <- "numeric"

#summarise dataset
asb_sp_sharks_summ <- mFD::asb.sp.summary(asb_sp_w = baskets_sharks_weights_100)

# retrieve species occurrences for all assemblages:
asb_sp_sharks_occ <- asb_sp_sharks_summ$asb_sp_occ

# Construct trait distance matrix using species
sp_dist_sharks_100 <- mFD::funct.dist(
  sp_tr         = sharks_traits_100,
  tr_cat        = sharks_traits_cat,
  metric        = "gower",
  scale_euclid  = "scale_center",
  ordinal_var   = "classic",
  weight_type   = "equal",
  stop_if_NA    = TRUE)

# Assess quality
fspaces_quality_sharks_100 <- mFD::quality.fspaces(
  sp_dist             = sp_dist_sharks_100,
  maxdim_pcoa         = 10,
  deviation_weighting = c("absolute","squared"),
  fdist_scaling       = FALSE,
  fdendro             = "average")

# Return coordinates of each axis
sp_faxes_coord_sharks_100 <- fspaces_quality_sharks_100$"details_fspaces"$"sp_pc_coord"

# Calculate FD metrics
alpha_fd_indices_sharks_100 <- mFD::alpha.fd.multidim(
  sp_faxes_coord   = sp_faxes_coord_sharks_100[ , c("PC1", "PC2", "PC3", "PC4")],
  asb_sp_w         = baskets_sharks_weights_100,
  ind_vect         = c("fric", "fori", "fspe"),
  scaling          = TRUE,
  check_input      = TRUE,
  details_returned = TRUE)

fd_ind_values_sharks_100 <- alpha_fd_indices_sharks_100$"functional_diversity_indices"

## Form dataframe
FDindices_taxon_100 <- as.data.frame(fd_ind_values_sharks_100) %>% 
  tibble::rownames_to_column("Scenario")

# Calculate functional uniqueness per species
Funct_Uni_100 <- uniqueness(baskets_sharks_weights_100, as.matrix(sp_dist_sharks_100)) %>% 
  as.data.frame()
Funct_Uni_100$Future <- baskets_100$Future
Funct_Uni_100$Y100 <- baskets_100$Y100

# Calculate functional uniqueness per scenario
Future_FUn_100df <- Funct_Uni_100 %>% 
  filter(Future == "1") %>% 
  select(Ui)
Future_FUn <- mean(Future_FUn_100df$Ui)

Y100_FUn_100df <- Funct_Uni_100 %>% 
  filter(Y100 == "1") %>% 
  select(Ui)
Y100_FUn <- mean(Y100_FUn_100df$Ui)

# Add uniqueness to final FD results data frame
Fun_100 <- data.frame(
  Fun = c(Future_FUn,Y100_FUn),
  row.names = c("Future","Y100")
)

FDindices_taxon_100$fun <- Fun_100$Fun
FDindices_taxon_100$Scenario <- ordered(FDindices_taxon_100$Scenario, levels=c("Future","Y100"))

# FUSE calculation
FUSE.present_100 <- fuse(sp_dist        = sp_dist_sharks_100, 
                         sp_faxes_coord = as.matrix(sp_faxes_coord_sharks_100), 
                         nb_NN          = 5,  
                         GE             = sharks_traits_100filter$iucn,
                         standGE        = TRUE)

# Status colours
status_colors <- c("LC" = "#008000", "NT" = "#ADFF2F", "VU" = "#FFFF00", "EN" = "#FFA500", "CR" = "#FF0000")

# Set up df
FUSE.present_100$Order <- sharks_traits_100filter$Order
FUSE.present_100$IUCN <- sharks_traits_100filter$iucn
FUSE.present_100$IUCN <- as.factor(FUSE.present_100$IUCN)
FUSE_sum_100 <- FUSE.present_100 %>% 
  mutate(IUCN = recode(IUCN, 
                       "0" = "LC", 
                       "1" = "NT",
                       "2" = "VU", 
                       "3" = "EN", 
                       "4" = "CR"))
FUSE_sum_100$IUCN <- ordered(FUSE_sum_100$IUCN, levels=c("LC", "NT", "VU", "EN", "CR"))
FUSE_sum_100 <- FUSE_sum_100 %>% 
  rownames_to_column(var = "Species")

# Add scenarios
FUSE_sum_100$Future_ext <- as.numeric(!baskets_100$Future)
FUSE_sum_100$Y100_ext <- as.numeric(!baskets_100$Y100)

## GLM
FUSE_sum_100$Rank_FUSE <- rank(-FUSE_sum_100$FUSE)
FUSE_sum_100$Order <- as.factor(FUSE_sum_100$Order)
FUSE_Y100_GLMM <- glmer(Y100_ext ~ Rank_FUSE + (Rank_FUSE | Order), data = FUSE_sum_100, family = binomial)
summary(FUSE_Y100_GLMM)

### 200 years in future
baskets_200 <- Occ_Future %>% 
  filter(Y100 == "1") %>% 
  select(Species,Y100,Y200) %>% 
  remove_rownames %>% 
  column_to_rownames(var="Species")

# Transpose to form final occurrence matrix
baskets_sharks_weights_200 <- t(baskets_200)

# Species-trait matrix
traits <- data %>% 
  select(Species,Order,habitat,vertical,terrestriality,thermo,feeding,diet,max.length,iucn)

sharks_traits <- traits %>% 
  remove_rownames %>% 
  column_to_rownames(var="Species")

sharks_traits_200filter <- sharks_traits %>%
  filter(row.names(.) %in% row.names(baskets_200))

sharks_traits_200 <- sharks_traits_200filter %>% 
  select(-c(iucn,Order))

# Form trait category matrix - weighting all equal
sharks_traits_cat <- tibble(trait_name = c("habitat","vertical","terrestriality","thermo","feeding","diet","max.length"),
                            trait_type = c("N", "N", "O", "O", "O", "N", "Q"))

# Summarise dataset
sharks_traits_summ <- mFD::sp.tr.summary(
  tr_cat     = sharks_traits_cat,   
  sp_tr      = sharks_traits_200, 
  stop_if_NA = TRUE)

# Reveal number of states for each trait
sharks_traits_summ$tr_summary_list

# make assemblage data into a matrix & ensure all matrix values are numeric
baskets_sharks_weights_200 <- data.matrix(baskets_sharks_weights_200, rownames.force = NA)
class(baskets_sharks_weights_200) <- "numeric"

#summarise dataset
asb_sp_sharks_summ <- mFD::asb.sp.summary(asb_sp_w = baskets_sharks_weights_200)

# retrieve species occurrences for all assemblages:
asb_sp_sharks_occ <- asb_sp_sharks_summ$asb_sp_occ

# Construct trait distance matrix using species
sp_dist_sharks_200 <- mFD::funct.dist(
  sp_tr         = sharks_traits_200,
  tr_cat        = sharks_traits_cat,
  metric        = "gower",
  scale_euclid  = "scale_center",
  ordinal_var   = "classic",
  weight_type   = "equal",
  stop_if_NA    = TRUE)

# Assess quality
fspaces_quality_sharks_200 <- mFD::quality.fspaces(
  sp_dist             = sp_dist_sharks_200,
  maxdim_pcoa         = 10,
  deviation_weighting = c("absolute","squared"),
  fdist_scaling       = FALSE,
  fdendro             = "average")

# Return coordinates of each axis
sp_faxes_coord_sharks_200 <- fspaces_quality_sharks_200$"details_fspaces"$"sp_pc_coord"

# Calculate FD metrics
alpha_fd_indices_sharks_200 <- mFD::alpha.fd.multidim(
  sp_faxes_coord   = sp_faxes_coord_sharks_200[ , c("PC1", "PC2", "PC3", "PC4")],
  asb_sp_w         = baskets_sharks_weights_200,
  ind_vect         = c("fric", "fori", "fspe"),
  scaling          = TRUE,
  check_input      = TRUE,
  details_returned = TRUE)

fd_ind_values_sharks_200 <- alpha_fd_indices_sharks_200$"functional_diversity_indices"

## Form dataframe
FDindices_taxon_200 <- as.data.frame(fd_ind_values_sharks_200) %>% 
  tibble::rownames_to_column("Scenario")

# Calculate functional uniqueness per species
Funct_Uni_200 <- uniqueness(baskets_sharks_weights_200, as.matrix(sp_dist_sharks_200)) %>% 
  as.data.frame()
Funct_Uni_200$Y100 <- baskets_200$Y100
Funct_Uni_200$Y200 <- baskets_200$Y200

# Calculate functional uniqueness per scenario
Y100_FUn_200df <- Funct_Uni_200 %>% 
  filter(Y100 == "1") %>% 
  select(Ui)
Y100_FUn <- mean(Y100_FUn_200df$Ui)

Y200_FUn_200df <- Funct_Uni_200 %>% 
  filter(Y200 == "1") %>% 
  select(Ui)
Y200_FUn <- mean(Y200_FUn_200df$Ui)

# Add uniqueness to final FD results data frame
Fun_200 <- data.frame(
  Fun = c(Y100_FUn,Y200_FUn),
  row.names = c("Y100","Y200")
)

FDindices_taxon_200$fun <- Fun_200$Fun
FDindices_taxon_200$Scenario <- ordered(FDindices_taxon_200$Scenario, levels=c("Y100","Y200"))

# FUSE calculation
FUSE.present_200 <- fuse(sp_dist        = sp_dist_sharks_200, 
                         sp_faxes_coord = as.matrix(sp_faxes_coord_sharks_200), 
                         nb_NN          = 5,  
                         GE             = sharks_traits_200filter$iucn,
                         standGE        = TRUE)

# Status colours
status_colors <- c("LC" = "#008000", "NT" = "#ADFF2F", "VU" = "#FFFF00", "EN" = "#FFA500", "CR" = "#FF0000")

# Set up df
FUSE.present_200$Order <- sharks_traits_200filter$Order
FUSE.present_200$IUCN <- sharks_traits_200filter$iucn
FUSE.present_200$IUCN <- as.factor(FUSE.present_200$IUCN)
FUSE_sum_200 <- FUSE.present_200 %>% 
  mutate(IUCN = recode(IUCN, 
                       "0" = "LC", 
                       "1" = "NT",
                       "2" = "VU", 
                       "3" = "EN", 
                       "4" = "CR"))
FUSE_sum_200$IUCN <- ordered(FUSE_sum_200$IUCN, levels=c("LC", "NT", "VU", "EN", "CR"))
FUSE_sum_200 <- FUSE_sum_200 %>% 
  rownames_to_column(var = "Species")

# Add scenarios
FUSE_sum_200$Y100_ext <- as.numeric(!baskets_200$Y100)
FUSE_sum_200$Y200_ext <- as.numeric(!baskets_200$Y200)

## GLMs
FUSE_sum_200$Rank_FUSE <- rank(-FUSE_sum_200$FUSE)
FUSE_sum_200$Order <- as.factor(FUSE_sum_200$Order)
FUSE_Y200_GLMM <- glmer(Y200_ext ~ Rank_FUSE + (Rank_FUSE | Order), data = FUSE_sum_200, family = binomial)
summary(FUSE_Y200_GLMM)

### 300 years in future
baskets_300 <- Occ_Future %>% 
  filter(Y200 == "1") %>% 
  select(Species,Y200,Y300) %>% 
  remove_rownames %>% 
  column_to_rownames(var="Species")

# Transpose to form final occurrence matrix
baskets_sharks_weights_300 <- t(baskets_300)

# Species-trait matrix
traits <- data %>% 
  select(Species,Order,habitat,vertical,terrestriality,thermo,feeding,diet,max.length,iucn)

sharks_traits <- traits %>% 
  remove_rownames %>% 
  column_to_rownames(var="Species")

sharks_traits_300filter <- sharks_traits %>%
  filter(row.names(.) %in% row.names(baskets_300))

sharks_traits_300 <- sharks_traits_300filter %>% 
  select(-c(iucn,Order))

# Form trait category matrix - weighting all equal
sharks_traits_cat <- tibble(trait_name = c("habitat","vertical","terrestriality","thermo","feeding","diet","max.length"),
                            trait_type = c("N", "N", "O", "O", "O", "N", "Q"))

# Summarise dataset
sharks_traits_summ <- mFD::sp.tr.summary(
  tr_cat     = sharks_traits_cat,   
  sp_tr      = sharks_traits_300, 
  stop_if_NA = TRUE)

# Reveal number of states for each trait
sharks_traits_summ$tr_summary_list

# make assemblage data into a matrix & ensure all matrix values are numeric
baskets_sharks_weights_300 <- data.matrix(baskets_sharks_weights_300, rownames.force = NA)
class(baskets_sharks_weights_300) <- "numeric"

#summarise dataset
asb_sp_sharks_summ <- mFD::asb.sp.summary(asb_sp_w = baskets_sharks_weights_300)

# retrieve species occurrences for all assemblages:
asb_sp_sharks_occ <- asb_sp_sharks_summ$asb_sp_occ

# Construct trait distance matrix using species
sp_dist_sharks_300 <- mFD::funct.dist(
  sp_tr         = sharks_traits_300,
  tr_cat        = sharks_traits_cat,
  metric        = "gower",
  scale_euclid  = "scale_center",
  ordinal_var   = "classic",
  weight_type   = "equal",
  stop_if_NA    = TRUE)

# Assess quality
fspaces_quality_sharks_300 <- mFD::quality.fspaces(
  sp_dist             = sp_dist_sharks_300,
  maxdim_pcoa         = 10,
  deviation_weighting = c("absolute","squared"),
  fdist_scaling       = FALSE,
  fdendro             = "average")

# Return coordinates of each axis
sp_faxes_coord_sharks_300 <- fspaces_quality_sharks_300$"details_fspaces"$"sp_pc_coord"

# Calculate FD metrics
alpha_fd_indices_sharks_300 <- mFD::alpha.fd.multidim(
  sp_faxes_coord   = sp_faxes_coord_sharks_300[ , c("PC1", "PC2", "PC3", "PC4")],
  asb_sp_w         = baskets_sharks_weights_300,
  ind_vect         = c("fric", "fori", "fspe"),
  scaling          = TRUE,
  check_input      = TRUE,
  details_returned = TRUE)

fd_ind_values_sharks_300 <- alpha_fd_indices_sharks_300$"functional_diversity_indices"

## Form dataframe
FDindices_taxon_300 <- as.data.frame(fd_ind_values_sharks_300) %>% 
  tibble::rownames_to_column("Scenario")

# Calculate functional uniqueness per species
Funct_Uni_300 <- uniqueness(baskets_sharks_weights_300, as.matrix(sp_dist_sharks_300)) %>% 
  as.data.frame()
Funct_Uni_300$Y200 <- baskets_300$Y200
Funct_Uni_300$Y300 <- baskets_300$Y300

# Calculate functional uniqueness per scenario
Y200_FUn_300df <- Funct_Uni_300 %>% 
  filter(Y200 == "1") %>% 
  select(Ui)
Y200_FUn <- mean(Y200_FUn_300df$Ui)

Y300_FUn_300df <- Funct_Uni_300 %>% 
  filter(Y300 == "1") %>% 
  select(Ui)
Y300_FUn <- mean(Y300_FUn_300df$Ui)

# Add uniqueness to final FD results data frame
Fun_300 <- data.frame(
  Fun = c(Y200_FUn,Y300_FUn),
  row.names = c("Y200","Y300")
)

FDindices_taxon_300$fun <- Fun_300$Fun
FDindices_taxon_300$Scenario <- ordered(FDindices_taxon_300$Scenario, levels=c("Y200","Y300"))

# FUSE calculation
FUSE.present_300 <- fuse(sp_dist        = sp_dist_sharks_300, 
                         sp_faxes_coord = as.matrix(sp_faxes_coord_sharks_300), 
                         nb_NN          = 5,  
                         GE             = sharks_traits_300filter$iucn,
                         standGE        = TRUE)

# Status colours
status_colors <- c("LC" = "#008000", "NT" = "#ADFF2F", "VU" = "#FFFF00", "EN" = "#FFA500", "CR" = "#FF0000")

# Set up df
FUSE.present_300$Order <- sharks_traits_300filter$Order
FUSE.present_300$IUCN <- sharks_traits_300filter$iucn
FUSE.present_300$IUCN <- as.factor(FUSE.present_300$IUCN)
FUSE_sum_300 <- FUSE.present_300 %>% 
  mutate(IUCN = recode(IUCN, 
                       "0" = "LC", 
                       "1" = "NT",
                       "2" = "VU", 
                       "3" = "EN", 
                       "4" = "CR"))
FUSE_sum_300$IUCN <- ordered(FUSE_sum_300$IUCN, levels=c("LC", "NT", "VU", "EN", "CR"))
FUSE_sum_300 <- FUSE_sum_300 %>% 
  rownames_to_column(var = "Species")

# Add scenarios
FUSE_sum_300$Y200_ext <- as.numeric(!baskets_300$Y200)
FUSE_sum_300$Y300_ext <- as.numeric(!baskets_300$Y300)

## GLMs
FUSE_sum_300$Rank_FUSE <- rank(-FUSE_sum_300$FUSE)
FUSE_sum_300$Order <- as.factor(FUSE_sum_300$Order)
FUSE_Y300_GLMM <- glmer(Y300_ext ~ Rank_FUSE + (Rank_FUSE | Order), data = FUSE_sum_300, family = binomial)
summary(FUSE_Y300_GLMM)

### 400 years in future
baskets_400 <- Occ_Future %>% 
  filter(Y300 == "1") %>% 
  select(Species,Y300,Y400) %>% 
  remove_rownames %>% 
  column_to_rownames(var="Species")

# Transpose to form final occurrence matrix
baskets_sharks_weights_400 <- t(baskets_400)

# Species-trait matrix
traits <- data %>% 
  select(Species,Order,habitat,vertical,terrestriality,thermo,feeding,diet,max.length,iucn)

sharks_traits <- traits %>% 
  remove_rownames %>% 
  column_to_rownames(var="Species")

sharks_traits_400filter <- sharks_traits %>%
  filter(row.names(.) %in% row.names(baskets_400))

sharks_traits_400 <- sharks_traits_400filter %>% 
  select(-c(iucn,Order))

# Form trait category matrix - weighting all equal
sharks_traits_cat <- tibble(trait_name = c("habitat","vertical","terrestriality","thermo","feeding","diet","max.length"),
                            trait_type = c("N", "N", "O", "O", "O", "N", "Q"))

# Summarise dataset
sharks_traits_summ <- mFD::sp.tr.summary(
  tr_cat     = sharks_traits_cat,   
  sp_tr      = sharks_traits_400, 
  stop_if_NA = TRUE)

# Reveal number of states for each trait
sharks_traits_summ$tr_summary_list

# make assemblage data into a matrix & ensure all matrix values are numeric
baskets_sharks_weights_400 <- data.matrix(baskets_sharks_weights_400, rownames.force = NA)
class(baskets_sharks_weights_400) <- "numeric"

#summarise dataset
asb_sp_sharks_summ <- mFD::asb.sp.summary(asb_sp_w = baskets_sharks_weights_400)

# retrieve species occurrences for all assemblages:
asb_sp_sharks_occ <- asb_sp_sharks_summ$asb_sp_occ

# Construct trait distance matrix using species
sp_dist_sharks_400 <- mFD::funct.dist(
  sp_tr         = sharks_traits_400,
  tr_cat        = sharks_traits_cat,
  metric        = "gower",
  scale_euclid  = "scale_center",
  ordinal_var   = "classic",
  weight_type   = "equal",
  stop_if_NA    = TRUE)

# Assess quality
fspaces_quality_sharks_400 <- mFD::quality.fspaces(
  sp_dist             = sp_dist_sharks_400,
  maxdim_pcoa         = 10,
  deviation_weighting = c("absolute","squared"),
  fdist_scaling       = FALSE,
  fdendro             = "average")

# Return coordinates of each axis
sp_faxes_coord_sharks_400 <- fspaces_quality_sharks_400$"details_fspaces"$"sp_pc_coord"

# Calculate FD metrics
alpha_fd_indices_sharks_400 <- mFD::alpha.fd.multidim(
  sp_faxes_coord   = sp_faxes_coord_sharks_400[ , c("PC1", "PC2", "PC3", "PC4")],
  asb_sp_w         = baskets_sharks_weights_400,
  ind_vect         = c("fric", "fori", "fspe"),
  scaling          = TRUE,
  check_input      = TRUE,
  details_returned = TRUE)

fd_ind_values_sharks_400 <- alpha_fd_indices_sharks_400$"functional_diversity_indices"

## Form dataframe
FDindices_taxon_400 <- as.data.frame(fd_ind_values_sharks_400) %>% 
  tibble::rownames_to_column("Scenario")

# Calculate functional uniqueness per species
Funct_Uni_400 <- uniqueness(baskets_sharks_weights_400, as.matrix(sp_dist_sharks_400)) %>% 
  as.data.frame()
Funct_Uni_400$Y300 <- baskets_400$Y300
Funct_Uni_400$Y400 <- baskets_400$Y400

# Calculate functional uniqueness per scenario
Y300_FUn_400df <- Funct_Uni_400 %>% 
  filter(Y300 == "1") %>% 
  select(Ui)
Y300_FUn <- mean(Y300_FUn_400df$Ui)

Y400_FUn_400df <- Funct_Uni_400 %>% 
  filter(Y400 == "1") %>% 
  select(Ui)
Y400_FUn <- mean(Y400_FUn_400df$Ui)

# Add uniqueness to final FD results data frame
Fun_400 <- data.frame(
  Fun = c(Y300_FUn,Y400_FUn),
  row.names = c("Y300","Y400")
)

FDindices_taxon_400$fun <- Fun_400$Fun
FDindices_taxon_400$Scenario <- ordered(FDindices_taxon_400$Scenario, levels=c("Y300","Y400"))

# FUSE calculation
FUSE.present_400 <- fuse(sp_dist        = sp_dist_sharks_400, 
                         sp_faxes_coord = as.matrix(sp_faxes_coord_sharks_400), 
                         nb_NN          = 5,  
                         GE             = sharks_traits_400filter$iucn,
                         standGE        = TRUE)

# Status colours
status_colors <- c("LC" = "#008000", "NT" = "#ADFF2F", "VU" = "#FFFF00", "EN" = "#FFA500", "CR" = "#FF0000")

# Set up df
FUSE.present_400$Order <- sharks_traits_400filter$Order
FUSE.present_400$IUCN <- sharks_traits_400filter$iucn
FUSE.present_400$IUCN <- as.factor(FUSE.present_400$IUCN)
FUSE_sum_400 <- FUSE.present_400 %>% 
  mutate(IUCN = recode(IUCN, 
                       "0" = "LC", 
                       "1" = "NT",
                       "2" = "VU", 
                       "3" = "EN", 
                       "4" = "CR"))
FUSE_sum_400$IUCN <- ordered(FUSE_sum_400$IUCN, levels=c("LC", "NT", "VU", "EN", "CR"))
FUSE_sum_400 <- FUSE_sum_400 %>% 
  rownames_to_column(var = "Species")

# Add scenarios
FUSE_sum_400$Y300_ext <- as.numeric(!baskets_400$Y300)
FUSE_sum_400$Y400_ext <- as.numeric(!baskets_400$Y400)

## GLMs
FUSE_sum_400$Rank_FUSE <- rank(-FUSE_sum_400$FUSE)
FUSE_sum_400$Order <- as.factor(FUSE_sum_400$Order)
FUSE_Y400_GLMM <- glmer(Y400_ext ~ Rank_FUSE + (Rank_FUSE | Order), data = FUSE_sum_400, family = binomial)
summary(FUSE_Y400_GLMM)

### 500 years in future
baskets_500 <- Occ_Future %>% 
  filter(Y400 == "1") %>% 
  select(Species,Y400,Y500) %>% 
  remove_rownames %>% 
  column_to_rownames(var="Species")

# Transpose to form final occurrence matrix
baskets_sharks_weights_500 <- t(baskets_500)

# Species-trait matrix
traits <- data %>% 
  select(Species,Order,habitat,vertical,terrestriality,thermo,feeding,diet,max.length,iucn)

sharks_traits <- traits %>% 
  remove_rownames %>% 
  column_to_rownames(var="Species")

sharks_traits_500filter <- sharks_traits %>%
  filter(row.names(.) %in% row.names(baskets_500))

sharks_traits_500 <- sharks_traits_500filter %>% 
  select(-c(iucn,Order))

# Form trait category matrix - weighting all equal
sharks_traits_cat <- tibble(trait_name = c("habitat","vertical","terrestriality","thermo","feeding","diet","max.length"),
                            trait_type = c("N", "N", "O", "O", "O", "N", "Q"))

# Summarise dataset
sharks_traits_summ <- mFD::sp.tr.summary(
  tr_cat     = sharks_traits_cat,   
  sp_tr      = sharks_traits_500, 
  stop_if_NA = TRUE)

# Reveal number of states for each trait
sharks_traits_summ$tr_summary_list

# make assemblage data into a matrix & ensure all matrix values are numeric
baskets_sharks_weights_500 <- data.matrix(baskets_sharks_weights_500, rownames.force = NA)
class(baskets_sharks_weights_500) <- "numeric"

#summarise dataset
asb_sp_sharks_summ <- mFD::asb.sp.summary(asb_sp_w = baskets_sharks_weights_500)

# retrieve species occurrences for all assemblages:
asb_sp_sharks_occ <- asb_sp_sharks_summ$asb_sp_occ

# Construct trait distance matrix using species
sp_dist_sharks_500 <- mFD::funct.dist(
  sp_tr         = sharks_traits_500,
  tr_cat        = sharks_traits_cat,
  metric        = "gower",
  scale_euclid  = "scale_center",
  ordinal_var   = "classic",
  weight_type   = "equal",
  stop_if_NA    = TRUE)

# Assess quality
fspaces_quality_sharks_500 <- mFD::quality.fspaces(
  sp_dist             = sp_dist_sharks_500,
  maxdim_pcoa         = 10,
  deviation_weighting = c("absolute","squared"),
  fdist_scaling       = FALSE,
  fdendro             = "average")

# Return coordinates of each axis
sp_faxes_coord_sharks_500 <- fspaces_quality_sharks_500$"details_fspaces"$"sp_pc_coord"

# Calculate FD metrics
alpha_fd_indices_sharks_500 <- mFD::alpha.fd.multidim(
  sp_faxes_coord   = sp_faxes_coord_sharks_500[ , c("PC1", "PC2", "PC3", "PC4")],
  asb_sp_w         = baskets_sharks_weights_500,
  ind_vect         = c("fric", "fori", "fspe"),
  scaling          = TRUE,
  check_input      = TRUE,
  details_returned = TRUE)

fd_ind_values_sharks_500 <- alpha_fd_indices_sharks_500$"functional_diversity_indices"

## Form dataframe
FDindices_taxon_500 <- as.data.frame(fd_ind_values_sharks_500) %>% 
  tibble::rownames_to_column("Scenario")

# Calculate functional uniqueness per species
Funct_Uni_500 <- uniqueness(baskets_sharks_weights_500, as.matrix(sp_dist_sharks_500)) %>% 
  as.data.frame()
Funct_Uni_500$Y400 <- baskets_500$Y400
Funct_Uni_500$Y500 <- baskets_500$Y500

# Calculate functional uniqueness per scenario
Y400_FUn_500df <- Funct_Uni_500 %>% 
  filter(Y400 == "1") %>% 
  select(Ui)
Y400_FUn <- mean(Y400_FUn_500df$Ui)

Y500_FUn_500df <- Funct_Uni_500 %>% 
  filter(Y500 == "1") %>% 
  select(Ui)
Y500_FUn <- mean(Y500_FUn_500df$Ui)

# Add uniqueness to final FD results data frame
Fun_500 <- data.frame(
  Fun = c(Y400_FUn,Y500_FUn),
  row.names = c("Y400","Y500")
)

FDindices_taxon_500$fun <- Fun_500$Fun
FDindices_taxon_500$Scenario <- ordered(FDindices_taxon_500$Scenario, levels=c("Y400","Y500"))

# FUSE calculation
FUSE.present_500 <- fuse(sp_dist        = sp_dist_sharks_500, 
                         sp_faxes_coord = as.matrix(sp_faxes_coord_sharks_500), 
                         nb_NN          = 5,  
                         GE             = sharks_traits_500filter$iucn,
                         standGE        = TRUE)

# Status colours
status_colors <- c("LC" = "#008000", "NT" = "#ADFF2F", "VU" = "#FFFF00", "EN" = "#FFA500", "CR" = "#FF0000")

# Set up df
FUSE.present_500$Order <- sharks_traits_500filter$Order
FUSE.present_500$IUCN <- sharks_traits_500filter$iucn
FUSE.present_500$IUCN <- as.factor(FUSE.present_500$IUCN)
FUSE_sum_500 <- FUSE.present_500 %>% 
  mutate(IUCN = recode(IUCN, 
                       "0" = "LC", 
                       "1" = "NT",
                       "2" = "VU", 
                       "3" = "EN", 
                       "4" = "CR"))
FUSE_sum_500$IUCN <- ordered(FUSE_sum_500$IUCN, levels=c("LC", "NT", "VU", "EN", "CR"))
FUSE_sum_500 <- FUSE_sum_500 %>% 
  rownames_to_column(var = "Species")

# Add scenarios
FUSE_sum_500$Y400_ext <- as.numeric(!baskets_500$Y400)
FUSE_sum_500$Y500_ext <- as.numeric(!baskets_500$Y500)

## GLMs
FUSE_sum_500$Rank_FUSE <- rank(-FUSE_sum_500$FUSE)
FUSE_sum_500$Order <- as.factor(FUSE_sum_500$Order)
FUSE_Y500_GLMM <- glmer(Y500_ext ~ Rank_FUSE + (Rank_FUSE | Order), data = FUSE_sum_500, family = binomial)
summary(FUSE_Y500_GLMM)

# AICs
AIC(FUSE_Future_GLMM,FUSE_Y100_GLMM,FUSE_Y200_GLMM,FUSE_Y300_GLMM,FUSE_Y400_GLMM,FUSE_Y500_GLMM)

# Plot FUSE figure - Figure 4
Fig_4a <- visreg(FUSE_Future_GLMM, scale = "response", main = "Extinction~FUSE in 2100", 
                 xlab = "FUSE (ranked)", ylab = "Extinction probability (2100)", 
                 line.par= list(col = "#AA1638"))
Fig_4b <- plot(vv,xlab="FUSE (ranked)",ylab="Extinction probability (2100)", layout=c(4,3),
               line.par= list(col = "#AA1638"))

# Plot FUSE all other futures - produces Figure S11
par(mfrow=c(5,1))
visreg(FUSE_Y100_GLMM, scale = "response", main = "Extinction~Ranked FUSE in 100 years", 
       xlab = "FUSE (ranking)", ylab = "Extinction probability", ylim = c(0,0.5), line.par= list(col = "#AA1638"))
visreg(FUSE_Y200_GLMM, scale = "response", main = "Extinction~Ranked FUSE in 200 years", 
       xlab = "FUSE (ranking)", ylab = "Extinction probability", ylim = c(0,0.5),line.par= list(col = "#AA1638"))
visreg(FUSE_Y300_GLMM, scale = "response", main = "Extinction~Ranked FUSE in 300 years", 
       xlab = "FUSE (ranking)", ylab = "Extinction probability", ylim = c(0,0.5), line.par= list(col = "#AA1638"))
visreg(FUSE_Y400_GLMM, scale = "response", main = "Extinction~Ranked FUSE in 400 years", 
       xlab = "FUSE (ranking)", ylab = "Extinction probability", ylim = c(0,0.5), line.par= list(col = "#AA1638"))
visreg(FUSE_Y500_GLMM, scale = "response", main = "Extinction~Ranked FUSE in 500 years", 
       xlab = "FUSE (ranking)", ylab = "Extinction probability", ylim = c(0,0.5), line.par= list(col = "#AA1638"))
par(mfrow=c(1,1))
