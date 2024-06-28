###############################################################################################################
# 21. Trait and environmental GLMs
## This R code uses GLMs to analyse relationships between traits & extinction
## It also analyses relationship between occurrence probability and environmental data in AquaMaps
## it produces Figure S7 and Table S5
###############################################################################################################

# Import packages
library(tidyverse)
library(tibble)
library(scales)
library(cowplot)
library(doBy)
library(reshape2)
library(data.table)
library(visreg)

# Load data 
load(file="C:/Users/2022207/Dropbox/Jack's PhD/Chapter 3. Future shark FD/Data/R data/Data.RData")

# Load mode IUCN_sim results
load(file = "C:/Users/2022207/Dropbox/Jack's PhD/Chapter 3. Future shark FD/Data/R data/2100_mode.RData")

# Form species trait matrix & add IUCN status + extinction status
Occ <- data %>% 
  select(Species)

Occ$Present <- 1

Occ_Future <- Occ %>%
  left_join(Spp_mode, by = "Species") %>%
  mutate(Future = ifelse(is.na(mode), 1, ifelse(mode <= 77, 0, 1))) %>%
  select(Species,Present,Future)

## Swap values in Occ_Future dataframe to set 1 as "extinct"
Occ_Future_swapped <- Occ_Future %>%
  mutate(
    Future_ext = ifelse(Future == 0, 1, 0))

baskets <- Occ_Future_swapped %>% 
  select(Species,Future_ext) %>% 
  remove_rownames()

traits <- data %>% 
  select(Species,habitat,vertical,terrestriality,thermo,feeding,diet,max.length,iucn)

# Join traits and baskets by the Species column
merged_data <- left_join(traits, baskets, by = "Species") %>% 
  remove_rownames() %>% 
  column_to_rownames("Species")

# Run GLMS of traits against extinction in 2100 - produces Figure S7 and Table S5
par(mfrow=c(4,2))
## Habitat
Hab <- glm(Future_ext ~ habitat, data = merged_data, family = "binomial")
Hab_drop<- drop1(Hab, test = "Chi")
Hab_vis <- visreg(Hab, scale = "response", main = "Extinction~Habitat", ylab = "Extinction probability")

## Vertical position
VP <- glm(Future_ext ~ vertical, data = merged_data, family = "binomial")
VP_drop<- drop1(VP, test = "Chi")
VP_vis <- visreg(VP, scale = "response", main = "Extinction~Vertical position", ylab = "Extinction probability")

## Terrestriality
Terrest <- glm(Future_ext ~ terrestriality, data = merged_data, family = "binomial")
T_drop<- drop1(Terrest, test = "Chi")
T_vis <- visreg(Terrest, scale = "response", main = "Extinction~Terrestriality", ylab = "Extinction probability")

## Thermoregulation
Therm <- glm(Future_ext ~ thermo, data = merged_data, family = "binomial")
Therm_drop<- drop1(Therm, test = "Chi")
Thermo_vis <- visreg(Therm, scale = "response", main = "Extinction~Thermoregulation", ylab = "Extinction probability")

## Feeding mechanism
FM <- glm(Future_ext ~ feeding, data = merged_data, family = "binomial")
FM_drop<- drop1(FM, test = "Chi")
FM_vis <- visreg(FM, scale = "response", main = "Exinction~Feeding mechanism", ylab = "Extinction probability")

## Diet
Di <- glm(Future_ext ~ diet, data = merged_data, family = "binomial")
D_drop<- drop1(Di, test = "Chi")
D_vis <- visreg(Di, scale = "response", main = "Extinction~Diet", ylab = "Extinction probability")

## Maximum body size
BS <- glm(Future_ext ~ max.length, data = merged_data, family = "binomial")
BS_drop<- drop1(BS, test = "Chi")
BS_vis <- visreg(BS, scale = "response", main = "Extinction~Max body size", ylab = "Extinction probability")
par(mfrow=c(1,1))

#################################################################################
# Relationships between environmental data and occurrence probabilities - Aquamaps data
## RCP26
RCP26 <- read_csv("C:/Users/2022207/Dropbox/Jack's PhD/Chapter 3. Future shark FD/Data/Aquamaps data/Elasmo_rcp26_2100.csv")
RCP26_environ <- read_csv("C:/Users/2022207/Dropbox/Jack's PhD/Chapter 3. Future shark FD/Data/Aquamaps data/hcaf_with_rcp26/hcaf_rcp26.csv")

merged_RCP26 <- merge(RCP26, RCP26_environ, by.x = "csquarecode", by.y = "CsquareCode", all.x = TRUE)

### Mean depth
MeanDepth_model <- glm(probability ~ DepthMean, data = merged_RCP26, family = binomial)
summary(MeanDepth_model)

### Min depth
MinDepth_model <- glm(probability ~ DepthMin, data = merged_RCP26, family = binomial)
summary(MinDepth_model)

### Max depth
MaxDepth_model <- glm(probability ~ DepthMax, data = merged_RCP26, family = binomial)
summary(MaxDepth_model)

### Mean surface temp
MeanSST_model <- glm(probability ~ SST2100c, data = merged_RCP26, family = binomial)
summary(MeanSST_model)

### Mean bottom surface temp
MeanSBT_model <- glm(probability ~ SBT2100c, data = merged_RCP26, family = binomial)
summary(MeanSBT_model)

### Mean surface salinity
MeanSalinity_model <- glm(probability ~ Salinity2100c, data = merged_RCP26, family = binomial)
summary(MeanSalinity_model)

### Mean bottom salinity
MeanBSalinity_model <- glm(probability ~ SalinityB2100c, data = merged_RCP26, family = binomial)
summary(MeanBSalinity_model)

### Mean primary productivity
PrimProd_model <- glm(probability ~ PrimProd2100c, data = merged_RCP26, family = binomial)
summary(PrimProd_model)

### Mean sea ice concentration
Ice_model <- glm(probability ~ IceCon2100c, data = merged_RCP26, family = binomial)
summary(Ice_model)

### Mean dissolved oxygen
Oxy_model <- glm(probability ~ OxyMean, data = merged_RCP26, family = binomial)
summary(Oxy_model)

### Mean bottom dissolved oxygen
OxyB_model <- glm(probability ~ OxyB2100c, data = merged_RCP26, family = binomial)
summary(OxyB_model)

AIC(MeanDepth_model,MinDepth_model,MaxDepth_model,MeanSST_model,MeanSBT_model,MeanSalinity_model,MeanBSalinity_model,
    PrimProd_model,Ice_model,Oxy_model,OxyB_model)

# Plot - produces supplementary figure x
par(mfrow=c(4,3))
visreg(MeanDepth_model, scale = "response", xlab = "Mean Depth", ylab = "P(occurrence)")
visreg(MinDepth_model, scale = "response", xlab = "Min Depth", ylab = "P(occurrence)")
visreg(MaxDepth_model, scale = "response", xlab = "Max Depth", ylab = "P(occurrence)")
visreg(MeanSST_model, scale = "response", xlab = "SST (2100)", ylab = "P(occurrence)")
visreg(MeanSBT_model, scale = "response", xlab = "SBT (2100)", ylab = "P(occurrence)")
visreg(MeanSalinity_model, scale = "response", xlab = "Salinity (2100)", ylab = "P(occurrence)")
visreg(MeanBSalinity_model, scale = "response", xlab = "Bottom Salinity (2100)", ylab = "P(occurrence)")
visreg(PrimProd_model, scale = "response", xlab = "Primary productivity (2100)", ylab = "P(occurrence)")
visreg(Ice_model, scale = "response", xlab = "ice Concentration (2100)", ylab = "P(occurrence)")
visreg(Oxy_model, scale = "response", xlab = "Mean Dissolved Oxygen", ylab = "P(occurrence)")
visreg(OxyB_model, scale = "response", xlab = "Bottom Oxygen (2100)", ylab = "P(occurrence)")

par(mfrow=c(1,1))

##RCP45
RCP45 <- read_csv("C:/Users/2022207/Dropbox/Jack's PhD/Chapter 3. Future shark FD/Data/Aquamaps data/Elasmo_rcp45_2100.csv")
RCP45_environ <- read_csv("C:/Users/2022207/Dropbox/Jack's PhD/Chapter 3. Future shark FD/Data/Aquamaps data/hcaf_with_rcp45/hcaf_rcp45.csv")

merged_RCP45 <- merge(RCP45, RCP45_environ, by = "CsquareCode", all.x = TRUE)

# Correlations between environmental factors (prior to 2100)
cor_df_RCP45 <- merged_RCP45 %>% 
  select(DepthMin,DepthMax,DepthMean,SSTAnMean,SBTAnMean,SalinityMean,SalinityBMean,PrimProdMean,IceConAnn,OxyMean,OxyBMean)
cor_matrix <- cor(cor_df_RCP45)
corrplot(cor_matrix)

### Mean depth
MeanDepth_model <- glm(Probability ~ DepthMean, data = merged_RCP45, family = binomial)
summary(MeanDepth_model)

### Min depth
MinDepth_model <- glm(Probability ~ DepthMin, data = merged_RCP45, family = binomial)
summary(MinDepth_model)

### Max depth
MaxDepth_model <- glm(Probability ~ DepthMax, data = merged_RCP45, family = binomial)
summary(MaxDepth_model)

### Mean surface temp
MeanSST_model <- glm(Probability ~ SST2100c, data = merged_RCP45, family = binomial)
summary(MeanSST_model)

### Mean bottom surface temp
MeanSBT_model <- glm(Probability ~ SBT2100c, data = merged_RCP45, family = binomial)
summary(MeanSBT_model)

### Mean surface salinity
MeanSalinity_model <- glm(Probability ~ Salinity2100c, data = merged_RCP45, family = binomial)
summary(MeanSalinity_model)

### Mean bottom salinity
MeanBSalinity_model <- glm(Probability ~ SalinityB2100c, data = merged_RCP45, family = binomial)
summary(MeanBSalinity_model)

### Mean primary productivity
PrimProd_model <- glm(Probability ~ PrimProd2100c, data = merged_RCP45, family = binomial)
summary(PrimProd_model)

### Mean sea ice concentration
Ice_model <- glm(Probability ~ IceCon2100c, data = merged_RCP45, family = binomial)
summary(Ice_model)

### Mean dissolved oxygen
Oxy_model <- glm(Probability ~ OxyMean, data = merged_RCP45, family = binomial)
summary(Oxy_model)

### Mean bottom dissolved oxygen
OxyB_model <- glm(Probability ~ OxyB2100c, data = merged_RCP45, family = binomial)
summary(OxyB_model)

AIC(MeanDepth_model,MinDepth_model,MaxDepth_model,MeanSST_model,MeanSBT_model,MeanSalinity_model,MeanBSalinity_model,
    PrimProd_model,Ice_model,Oxy_model,OxyB_model)

# Plot - produces supplementary figure x
par(mfrow=c(4,3))
visreg(MeanDepth_model, scale = "response", xlab = "Mean Depth", ylab = "P(occurrence)")
visreg(MinDepth_model, scale = "response", xlab = "Min Depth", ylab = "P(occurrence)")
visreg(MaxDepth_model, scale = "response", xlab = "Max Depth", ylab = "P(occurrence)")
visreg(MeanSST_model, scale = "response", xlab = "SST (2100)", ylab = "P(occurrence)")
visreg(MeanSBT_model, scale = "response", xlab = "SBT (2100)", ylab = "P(occurrence)")
visreg(MeanSalinity_model, scale = "response", xlab = "Salinity (2100)", ylab = "P(occurrence)")
visreg(MeanBSalinity_model, scale = "response", xlab = "Bottom Salinity (2100)", ylab = "P(occurrence)")
visreg(PrimProd_model, scale = "response", xlab = "Primary productivity (2100)", ylab = "P(occurrence)")
visreg(Ice_model, scale = "response", xlab = "ice Concentration (2100)", ylab = "P(occurrence)")
visreg(Oxy_model, scale = "response", xlab = "Mean Dissolved Oxygen", ylab = "P(occurrence)")
visreg(OxyB_model, scale = "response", xlab = "Bottom Oxygen (2100)", ylab = "P(occurrence)")

par(mfrow=c(1,1))

##RCP85
RCP85 <- read_csv("C:/Users/2022207/Dropbox/Jack's PhD/Chapter 3. Future shark FD/Data/Aquamaps data/Elasmo_rcp85_2100.csv") 
RCP85_environ <- read_csv("C:/Users/2022207/Dropbox/Jack's PhD/Chapter 3. Future shark FD/Data/Aquamaps data/hcaf_with_rcp85/hcaf_with_rcp85.csv")

merged_RCP85 <- merge(RCP85, RCP85_environ, by = "CsquareCode", all.x = TRUE)

# Correlations between environmental factors (prior to 2100)
cor_df_RCP85 <- merged_RCP85 %>% 
  select(DepthMin,DepthMax,DepthMean,SSTAnMean,SBTAnMean,SalinityMean,SalinityBMean,PrimProdMean,IceConAnn,OxyMean,OxyBMean)
cor_matrix <- cor(cor_df_RCP85)
corrplot(cor_matrix)

### Mean depth
MeanDepth_model <- glm(Probability ~ DepthMean, data = merged_RCP85, family = binomial)
summary(MeanDepth_model)

### Min depth
MinDepth_model <- glm(Probability ~ DepthMin, data = merged_RCP85, family = binomial)
summary(MinDepth_model)

### Max depth
MaxDepth_model <- glm(Probability ~ DepthMax, data = merged_RCP85, family = binomial)
summary(MaxDepth_model)

### Mean surface temp
MeanSST_model <- glm(Probability ~ SST2100c, data = merged_RCP85, family = binomial)
summary(MeanSST_model)

### Mean bottom surface temp
MeanSBT_model <- glm(Probability ~ SBT2100c, data = merged_RCP85, family = binomial)
summary(MeanSBT_model)

### Mean surface salinity
MeanSalinity_model <- glm(Probability ~ Salinity2100c, data = merged_RCP85, family = binomial)
summary(MeanSalinity_model)

### Mean bottom salinity
MeanBSalinity_model <- glm(Probability ~ SalinityB2100c, data = merged_RCP85, family = binomial)
summary(MeanBSalinity_model)

### Mean primary productivity
PrimProd_model <- glm(Probability ~ PrimProd2100c, data = merged_RCP85, family = binomial)
summary(PrimProd_model)

### Mean sea ice concentration
Ice_model <- glm(Probability ~ IceCon2100c, data = merged_RCP85, family = binomial)
summary(Ice_model)

### Mean dissolved oxygen
Oxy_model <- glm(Probability ~ OxyMean, data = merged_RCP85, family = binomial)
summary(Oxy_model)

### Mean bottom dissolved oxygen
OxyB_model <- glm(Probability ~ OxyB2100c, data = merged_RCP85, family = binomial)
summary(OxyB_model)

AIC(MeanDepth_model,MinDepth_model,MaxDepth_model,MeanSST_model,MeanSBT_model,MeanSalinity_model,MeanBSalinity_model,
    PrimProd_model,Ice_model,Oxy_model,OxyB_model)

# Plot 
par(mfrow=c(4,3))
visreg(MeanDepth_model, scale = "response", xlab = "Mean Depth", ylab = "P(occurrence)")
visreg(MinDepth_model, scale = "response", xlab = "Min Depth", ylab = "P(occurrence)")
visreg(MaxDepth_model, scale = "response", xlab = "Max Depth", ylab = "P(occurrence)")
visreg(MeanSST_model, scale = "response", xlab = "SST (2100)", ylab = "P(occurrence)")
visreg(MeanSBT_model, scale = "response", xlab = "SBT (2100)", ylab = "P(occurrence)")
visreg(MeanSalinity_model, scale = "response", xlab = "Salinity (2100)", ylab = "P(occurrence)")
visreg(MeanBSalinity_model, scale = "response", xlab = "Bottom Salinity (2100)", ylab = "P(occurrence)")
visreg(PrimProd_model, scale = "response", xlab = "Primary productivity (2100)", ylab = "P(occurrence)")
visreg(Ice_model, scale = "response", xlab = "ice Concentration (2100)", ylab = "P(occurrence)")
visreg(Oxy_model, scale = "response", xlab = "Mean Dissolved Oxygen", ylab = "P(occurrence)")
visreg(OxyB_model, scale = "response", xlab = "Bottom Oxygen (2100)", ylab = "P(occurrence)")

par(mfrow=c(1,1))