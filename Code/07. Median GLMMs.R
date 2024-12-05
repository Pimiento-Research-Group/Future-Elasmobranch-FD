###############################################################################################################
# 07. Median GLMMs
## This R code assesses the relationship between extinction probability and ranked FUS scores
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
library(visreg)
library(lme4)

# Load data
load(file="C:/Users/2022207/Dropbox/Jack's PhD/Chapter 3. Future shark FD/Data/R data/Data.RData")

# Load median IUCN_sim results
### Crit E EX
load(file = "C:/Users/2022207/Dropbox/Jack's PhD/Chapter 3. Future shark FD/Analyses/Final paper R codes/data/Crit_E/CritE_median.RData")

# Load synonyms and correct in Crit E data
Synonyms <- read_xlsx("C:/Users/2022207/Dropbox/Jack's PhD/Chapter 3. Future shark FD/Analyses/DatasetSynonyms.xlsx")

Spp_median.CritE <- Spp_median.CritE %>%
  left_join(Synonyms, by = c("Species" = "iucnsim_name")) %>%
  mutate(Species = ifelse(!is.na(Accepted_name), Accepted_name, Species)) %>%
  select(-c(Accepted_name,Marked_iucnsim,Aquamaps_name,Marked_aquamaps)) %>% 
  distinct(Species, .keep_all = TRUE)

# Make future functional space
# Occurrence data - under Crit E EX mode here but can be replicated for all 3
Occ_Future <- data %>% 
  select(Species,iucn)

Occ_Future$Present <- 1

# Set up future scenarios
## 2100
Occ_Future <- Occ_Future %>%
  left_join(Spp_median.CritE, by = "Species") %>%
  mutate(Future = ifelse(is.na(median), 1, ifelse(median <= 77, 0, 1))) %>%
  select(Species,Present,Future,median,iucn)

## IUCN AT
Occ_Future <- Occ_Future %>%
  mutate(AT = ifelse(iucn %in% c(2, 3, 4), 0, ifelse(iucn %in% c(0, 1), 1, NA))) %>%
  select(Species,Present,Future,AT)

baskets <- Occ_Future %>%
  select(Species,Present,Future,AT) %>% 
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
  mutate(FUS = FUn_std + FSp_std/2)
FUSE_sum <- FUSE_sum %>% 
  rownames_to_column(var = "Species")

# Add scenarios
FUSE_sum$Present_ext <- as.numeric(!baskets$Present)
FUSE_sum$Future_ext <- as.numeric(!baskets$Future)
FUSE_sum$Rank_FUS <- rank(-FUSE_sum$FUS)
FUSE_sum$Rank_FUn <- rank(-FUSE_sum$FUn_std)
FUSE_sum$Rank_FSp <- rank(-FUSE_sum$FSp_std)

# Generealised mixed effect model including Order
FUSE_sum$Order <- as.factor(FUSE_sum$Order)
FUS_Future_GLMM <- glmer(Future_ext ~ Rank_FUS + (Rank_FUS | Order), data = FUSE_sum, family = binomial)
summary(FUS_Future_GLMM)
visreg(FUS_Future_GLMM, scale = "response", main = "Extinction~FUS in 2100", 
       xlab = "FUS (ranked)", ylab = "Extinction probability (2100)", 
       line.par= list(col = "#AA1638"))

# Plot partial pooling/random trend by Order
v <- visreg(FUS_Future_GLMM, "Rank_FUS",by="Order", plot = FALSE, scale = "response")
subOrder <- sample(levels(FUSE_sum$Order),12)
vv <- subset(v, Order %in% subOrder)
plot(vv,xlab="FUS (ranked)",ylab="Extinction probability (2100)", layout=c(4,3),
     line.par= list(col = "#AA1638"))

# Predict the extinction probabilities using the GLMM
predictions <- predict(FUS_Future_GLMM, type = "response")

# Combine predictions with the original data to get the Order information
data_with_predictions <- data.frame(FUSE_sum, Predicted_Prob = predictions)

# Filter the data to include only the selected subOrder
filtered_data <- subset(data_with_predictions, Order %in% subOrder)

# Calculate the average extinction probability overall and per order
overall_avg_extinction <- mean(data_with_predictions$Predicted_Prob)
avg_extinction <- aggregate(Predicted_Prob ~ Order, data = filtered_data, FUN = mean)

# plots of 12 orders - raw data
glm_results <- data.frame(
  Order = character(),
  Estimate = numeric(),
  Std_Error = numeric(),
  Z_value = numeric(),
  Pr_Z = numeric(),
  stringsAsFactors = FALSE
)

# List of orders to plot in the specified order
shark_orders <- c("Carcharhiniformes", "Lamniformes", "Orectolobiformes", "Heterodontiformes", 
                  "Squaliformes", "Pristiophoriformes", "Squatiniformes", "Hexanchiformes", 
                  "Myliobatiformes", "Rhinopristiformes", "Torpediniformes", "Rajiformes")

# Set up the plotting layout to 3 rows and 4 columns
par(mfrow = c(3, 4), mar = c(4, 4, 2, 1))  # Adjust margins as needed

# Function to fit glm, plot using visreg, and store results for each order
plot_visreg_and_save_results <- function(order_name) {
  # Filter data for the specified order
  order_data <- FUSE_sum %>% 
    filter(Order == order_name)
  
  # Check if there is enough variation in the predictor variable
  if(length(unique(order_data$Rank_FUS)) > 1) {
    # Fit the GLM for the order
    order_glm <- glm(Future_ext ~ Rank_FUS, data = order_data, family = binomial)
    
    # Extract the coefficients summary
    glm_summary <- summary(order_glm)$coefficients
    
    # Extract relevant values from the summary (for the Rank_FUS predictor)
    rank_fus_results <- glm_summary["Rank_FUS", ]
    
    # Store the results in the dataframe
    glm_results <<- rbind(
      glm_results,
      data.frame(
        Order = order_name,
        Estimate = rank_fus_results["Estimate"],
        Std_Error = rank_fus_results["Std. Error"],
        Z_value = rank_fus_results["z value"],
        Pr_Z = rank_fus_results["Pr(>|z|)"]
      )
    )
    
    # Plot using visreg
    visreg(order_glm, scale = "response", main = order_name, xlab = "FUS (ranked)", ylab = "Extinction probability (2100)")
  } else {
    # Plot a flat line showing the average extinction probability if not enough variation
    avg_ext_prob <- mean(order_data$Future_ext)
    plot(NULL, xlim = c(0, 1), ylim = c(0, 1), xlab = "FUS (ranked)", ylab = "Extinction probability", main = order_name)
    abline(h = avg_ext_prob, col = "#AA1638")
  }
}

# Apply the function to the orders in the specified order
for(order_name in shark_orders) {
  if(order_name %in% unique(FUSE_sum$Order)) {
    plot_visreg_and_save_results(order_name)
  }
}

# Print the GLM results for all orders
print(glm_results)

# Reset graphics parameters to default
par(mfrow = c(1, 1))

# Identify top 10 FUS and presence/absence in each scenario
FUSE_sum$Present <- baskets$Present
FUSE_sum$Future <- baskets$Future
FUSE_sum$AT <- baskets$AT

Top25 <- FUSE_sum %>%
  left_join(Spp_median.CritE %>% select(Species, median), by = "Species") %>%
  mutate(Year = round(median + 2023)) %>%
  select(Species, FUS, Rank_FUS, Future, AT, median, Year)

### IUCN AT
baskets_AT <- Occ_Future %>% 
  filter(Future == "1") %>% 
  select(Species,Future,AT) %>% 
  remove_rownames %>% 
  column_to_rownames(var="Species")

# Transpose to form final occurrence matrix
baskets_sharks_weights_AT <- t(baskets_AT)

# Species-trait matrix
traits <- data %>% 
  select(Species,Order,habitat,vertical,terrestriality,thermo,feeding,diet,max.length,iucn)

sharks_traits <- traits %>% 
  remove_rownames %>% 
  column_to_rownames(var="Species")

sharks_traits_ATfilter <- sharks_traits %>%
  filter(row.names(.) %in% row.names(baskets_AT))

sharks_traits_AT <- sharks_traits_ATfilter %>% 
  select(-c(iucn,Order))

# Form trait category matrix - weighting all equal
sharks_traits_cat <- tibble(trait_name = c("habitat","vertical","terrestriality","thermo","feeding","diet","max.length"),
                            trait_type = c("N", "N", "O", "O", "O", "N", "Q"))

# Summarise dataset
sharks_traits_summ <- mFD::sp.tr.summary(
  tr_cat     = sharks_traits_cat,   
  sp_tr      = sharks_traits_AT, 
  stop_if_NA = TRUE)

# Reveal number of states for each trait
sharks_traits_summ$tr_summary_list

# make assemblage data into a matrix & ensure all matrix values are numeric
baskets_sharks_weights_AT <- data.matrix(baskets_sharks_weights_AT, rownames.force = NA)
class(baskets_sharks_weights_AT) <- "numeric"

#summarise dataset
asb_sp_sharks_summ <- mFD::asb.sp.summary(asb_sp_w = baskets_sharks_weights_AT)

# retrieve species occurrences for all assemblages:
asb_sp_sharks_occ <- asb_sp_sharks_summ$asb_sp_occ

# Construct trait distance matrix using species
sp_dist_sharks_AT <- mFD::funct.dist(
  sp_tr         = sharks_traits_AT,
  tr_cat        = sharks_traits_cat,
  metric        = "gower",
  scale_euclid  = "scale_center",
  ordinal_var   = "classic",
  weight_type   = "equal",
  stop_if_NA    = TRUE)

# Assess quality
fspaces_quality_sharks_AT <- mFD::quality.fspaces(
  sp_dist             = sp_dist_sharks_AT,
  maxdim_pcoa         = 10,
  deviation_weighting = c("absolute","squared"),
  fdist_scaling       = FALSE,
  fdendro             = "average")

# Return coordinates of each axis
sp_faxes_coord_sharks_AT <- fspaces_quality_sharks_AT$"details_fspaces"$"sp_pc_coord"

# FUSE calculation
FUSE.present_AT <- fuse(sp_dist        = sp_dist_sharks_AT, 
                         sp_faxes_coord = as.matrix(sp_faxes_coord_sharks_AT), 
                         nb_NN          = 5,  
                         GE             = sharks_traits_ATfilter$iucn,
                         standGE        = TRUE)

# Status colours
status_colors <- c("LC" = "#008000", "NT" = "#ADFF2F", "VU" = "#FFFF00", "EN" = "#FFA500", "CR" = "#FF0000")

# Set up df
FUSE.present_AT$Order <- sharks_traits_ATfilter$Order
FUSE.present_AT$IUCN <- sharks_traits_ATfilter$iucn
FUSE.present_AT$IUCN <- as.factor(FUSE.present_AT$IUCN)
FUSE_sum_AT <- FUSE.present_AT %>% 
  mutate(IUCN = recode(IUCN, 
                       "0" = "LC", 
                       "1" = "NT",
                       "2" = "VU", 
                       "3" = "EN", 
                       "4" = "CR"))
FUSE_sum_AT$IUCN <- ordered(FUSE_sum_AT$IUCN, levels=c("LC", "NT", "VU", "EN", "CR"))

FUSE_sum_AT <- FUSE_sum_AT %>% 
  mutate(FUS = FUn_std + FSp_std/2)
FUSE_sum_AT <- FUSE_sum_AT %>% 
  rownames_to_column(var = "Species")

# Add scenarios
FUSE_sum_AT$Future_ext <- as.numeric(!baskets_AT$Future)
FUSE_sum_AT$AT_ext <- as.numeric(!baskets_AT$AT)

## GLM
FUSE_sum_AT$Rank_FUS <- rank(-FUSE_sum_AT$FUS)
FUSE_sum$Rank_FUn <- rank(-FUSE_sum$FUn_std)
FUSE_sum$Rank_FSp <- rank(-FUSE_sum$FSp_std)
FUSE_sum_AT$Order <- as.factor(FUSE_sum_AT$Order)
FUSE_AT_GLMM <- glmer(AT_ext ~ Rank_FUS + (Rank_FUS | Order), data = FUSE_sum_AT, family = binomial)
summary(FUSE_AT_GLMM)

visreg(FUSE_AT_GLMM, scale = "response", main = "Extinction~FUS in IUCN AT", 
       xlab = "FUS (ranking)", ylab = "Extinction probability",  line.par= list(col = "#AA1638"))

# Create an empty dataframe to store the results
glm_results_AT <- data.frame(
  Order = character(),
  Estimate = numeric(),
  Std_Error = numeric(),
  Z_value = numeric(),
  Pr_Z = numeric(),
  stringsAsFactors = FALSE
)

# Function to fit glm, plot using visreg, and store results for each order in the AT scenario
plot_visreg_and_save_results_AT <- function(order_name) {
  # Filter data for the specified order
  order_data <- FUSE_sum_AT %>% 
    filter(Order == order_name)
  
  # Check if there is enough variation in the predictor variable
  if(length(unique(order_data$Rank_FUS)) > 1) {
    # Fit the GLM for the order
    order_glm <- glm(AT_ext ~ Rank_FUS, data = order_data, family = binomial)
    
    # Extract the coefficients summary
    glm_summary <- summary(order_glm)$coefficients
    
    # Extract relevant values from the summary (for the Rank_FUS predictor)
    rank_fus_results <- glm_summary["Rank_FUS", ]
    
    # Store the results in the dataframe
    glm_results_AT <<- rbind(
      glm_results_AT,
      data.frame(
        Order = order_name,
        Estimate = rank_fus_results["Estimate"],
        Std_Error = rank_fus_results["Std. Error"],
        Z_value = rank_fus_results["z value"],
        Pr_Z = rank_fus_results["Pr(>|z|)"]
      )
    )
    
    # Plot using visreg
    visreg(order_glm, scale = "response", main = paste(order_name), 
           xlab = "FUS (ranked)", ylab = "Extinction probability (IUCN AT)")
  } else {
    # Plot a flat line showing the average extinction probability if not enough variation
    avg_ext_prob <- mean(order_data$AT_ext)
    plot(NULL, xlim = c(0, 1), ylim = c(0, 1), xlab = "FUS (ranked)", ylab = "Extinction probability (IUCN AT)", main = order_name)
    abline(h = avg_ext_prob, col = "#AA1638")
  }
}

# List of orders to plot
shark_orders <- c("Carcharhiniformes", "Lamniformes", "Orectolobiformes", "Heterodontiformes", 
                  "Squaliformes", "Pristiophoriformes", "Squatiniformes", "Hexanchiformes", 
                  "Myliobatiformes", "Rhinopristiformes", "Torpediniformes", "Rajiformes")

# Filter FUSE_sum_AT for the 12 shark orders
FUSE_sum_AT <- FUSE_sum_AT %>% filter(Order %in% shark_orders)

# Set up the plotting area for a 3x4 grid (3 rows, 4 columns)
par(mfrow = c(3, 4))

# Loop through the list of orders and plot for each
for (order in shark_orders) {
  plot_visreg_and_save_results_AT(order)
}

# Reset the plotting area back to a single plot (1x1 grid)
par(mfrow = c(1, 1))

# Print the GLM results for all 12 shark orders
print(glm_results_AT)

# Plot FUS figure - Figure 2
par(mfrow=c(1,2))
Fig_2a <- visreg(FUS_Future_GLMM, scale = "response", main = "", 
                 xlab = "FUS (ranked)", ylab = "Extinction probability (2100)", 
                 line.par= list(col = "#AA1638"))

Fig_2c <- visreg(FUSE_AT_GLMM, scale = "response", main = "", 
       xlab = "FUS (ranking)", ylab = "Extinction probability (IUCN AT)", line.par= list(col = "#AA1638"))
par(mfrow=c(1,1))
