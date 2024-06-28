###############################################################################################################
# 23. Environmental drivers of RCP 4.5 FD
## This R code analyses the relationship of environmental data in 2100 from RCP 4.5 AquaMaps 
## to projected functional diversity
## it produces part of Tables S6-S9
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

## Environmental drivers of FD after extinctions
# Load RCP 4.5 with extinctions
load(file = "C:/Users/2022207/Dropbox/Jack's PhD/Chapter 3. Future shark FD/Data/R data/FD_RCP45_extinctions_map_df.RData")

# Load RCP 4.5 environmental data
RCP45_environ <- read_csv("C:/Users/2022207/Dropbox/Jack's PhD/Chapter 3. Future shark FD/Data/Aquamaps data/hcaf_with_rcp45/hcaf_rcp45.csv")

# Add environmental data to FD dataframe by grid
## Merge dataframes by grid & isolate environmental data by 2100
FD_RCP45_full <- merge(final_map_df_RCP45_extinctions, RCP45_environ, by.x = "Grid", by.y = "CsquareCode", all.x = TRUE)

FD_RCP45 <- FD_RCP45_full %>% 
  select(Grid,CenterLat.x,CenterLong.x,sp_richn,fric,fori,fspe,FUn,
         DepthMean,DepthMin,DepthMax,SST2100c,SBT2100c,Salinity2100c,SalinityB2100c,
         PrimProd2100c,IceCon2100c,OxyMean,OxyB2100c)

# Relationship between each environmental factor and functional diversity metrics across grids
## Species richness
par(mfrow=c(4,3))
MeanDepth_Spp <- glm(sp_richn ~ DepthMean, data = FD_RCP45, family = poisson)
visreg(MeanDepth_Spp, scale = "response", main = "", xlab = "Mean Depth", ylab = "P(Sp_Richn)")
MeanDepth_Spp_drop<- drop1(MeanDepth_Spp, test = "Chi")

MinDepth_Spp <- glm(sp_richn ~ DepthMin, data = FD_RCP45, family = poisson)
visreg(MinDepth_Spp, scale = "response", main = "", xlab = "Min Depth", ylab = "P(Sp_Richn)")
MinDepth_Spp_drop<- drop1(MinDepth_Spp, test = "Chi")

MaxDepth_Spp <- glm(sp_richn ~ DepthMax, data = FD_RCP45, family = poisson)
visreg(MaxDepth_Spp, scale = "response", main = "", xlab = "Max Depth", ylab = "P(Sp_Richn)")
MaxDepth_Spp_drop<- drop1(MaxDepth_Spp, test = "Chi")

SST2100c_Spp <- glm(sp_richn ~ SST2100c, data = FD_RCP45, family = poisson)
visreg(SST2100c_Spp, scale = "response", main = "", xlab = "SST", ylab = "P(Sp_Richn)")
SST2100c_Spp_drop<- drop1(SST2100c_Spp, test = "Chi")

SBT2100c_Spp <- glm(sp_richn ~ SBT2100c, data = FD_RCP45, family = poisson)
visreg(SBT2100c_Spp, scale = "response", main = "", xlab = "SBT", ylab = "P(Sp_Richn)")
SBT2100c_Spp_drop<- drop1(SBT2100c_Spp, test = "Chi")

Salinity2100c_Spp <- glm(sp_richn ~ Salinity2100c, data = FD_RCP45, family = poisson)
visreg(Salinity2100c_Spp, scale = "response", main = "", xlab = "Salinity", ylab = "P(Sp_Richn)")
Salinity2100c_Spp_drop<- drop1(Salinity2100c_Spp, test = "Chi")

SalinityB2100c_Spp <- glm(sp_richn ~ SalinityB2100c, data = FD_RCP45, family = poisson)
visreg(SalinityB2100c_Spp, scale = "response", main = "", xlab = "Bottom Salinity", ylab = "P(Sp_Richn)")
SalinityB2100c_Spp_drop<- drop1(SalinityB2100c_Spp, test = "Chi")

PrimProd2100c_Spp <- glm(sp_richn ~ PrimProd2100c, data = FD_RCP45, family = poisson)
visreg(PrimProd2100c_Spp, scale = "response", main = "", xlab = "Primary productivity", ylab = "P(Sp_Richn)")
PrimProd2100c_Spp_drop<- drop1(PrimProd2100c_Spp, test = "Chi")

IceCon2100c_Spp <- glm(sp_richn ~ IceCon2100c, data = FD_RCP45, family = poisson)
visreg(IceCon2100c_Spp, scale = "response", main = "", xlab = "Ice concentration", ylab = "P(Sp_Richn)")
IceCon2100c_Spp_drop<- drop1(IceCon2100c_Spp, test = "Chi")

OxyMean_Spp <- glm(sp_richn ~ OxyMean, data = FD_RCP45, family = poisson)
visreg(OxyMean_Spp, scale = "response", main = "", xlab = "Mean dissolved oxygen", ylab = "P(Sp_Richn)")
OxyMean_Spp_drop<- drop1(OxyMean_Spp, test = "Chi")

OxyB2100c_Spp <- glm(sp_richn ~ OxyB2100c, data = FD_RCP45, family = poisson)
visreg(OxyB2100c_Spp, scale = "response", main = "", xlab = "Dissolved bottom oxygen", ylab = "P(Sp_Richn)")
OxyB2100c_Spp_drop<- drop1(OxyB2100c_Spp, test = "Chi")
par(mfrow=c(1,1))

AIC(MeanDepth_Spp,MinDepth_Spp,MaxDepth_Spp,SST2100c_Spp,SBT2100c_Spp,Salinity2100c_Spp,SalinityB2100c_Spp,
    PrimProd2100c_Spp,IceCon2100c_Spp,OxyMean_Spp,OxyB2100c_Spp)

## Functional richness
par(mfrow=c(4,3))
MeanDepth_FRic <- glm(fric ~ DepthMean, data = FD_RCP45, family = binomial)
visreg(MeanDepth_FRic, scale = "response", main = "", xlab = "Mean Depth", ylab = "P(FRic)")
MeanDepth_FRic_drop<- drop1(MeanDepth_FRic, test = "Chi")

MinDepth_FRic <- glm(fric ~ DepthMin, data = FD_RCP45, family = binomial)
visreg(MinDepth_FRic, scale = "response", main = "", xlab = "Min Depth", ylab = "P(FRic)")
MinDepth_FRic_drop<- drop1(MinDepth_FRic, test = "Chi")

MaxDepth_FRic <- glm(fric ~ DepthMax, data = FD_RCP45, family = binomial)
visreg(MaxDepth_FRic, scale = "response", main = "", xlab = "Max Depth", ylab = "P(FRic)")
MaxDepth_FRic_drop<- drop1(MaxDepth_FRic, test = "Chi")

SST2100c_FRic <- glm(fric ~ SST2100c, data = FD_RCP45, family = binomial)
visreg(SST2100c_FRic, scale = "response", main = "", xlab = "SST", ylab = "P(FRic)")
SST2100c_FRic_drop<- drop1(SST2100c_FRic, test = "Chi")

SBT2100c_FRic <- glm(fric ~ SBT2100c, data = FD_RCP45, family = binomial)
visreg(SBT2100c_FRic, scale = "response", main = "", xlab = "SBT", ylab = "P(FRic)")
SBT2100c_FRic_drop<- drop1(SBT2100c_FRic, test = "Chi")

Salinity2100c_FRic <- glm(fric ~ Salinity2100c, data = FD_RCP45, family = binomial)
visreg(Salinity2100c_FRic, scale = "response", main = "", xlab = "Salinity", ylab = "P(FRic)")
Salinity2100c_FRic_drop<- drop1(Salinity2100c_FRic, test = "Chi")

SalinityB2100c_FRic <- glm(fric ~ SalinityB2100c, data = FD_RCP45, family = binomial)
visreg(SalinityB2100c_FRic, scale = "response", main = "", xlab = "Bottom Salinity", ylab = "P(FRic)")
SalinityB2100c_FRic_drop<- drop1(SalinityB2100c_FRic, test = "Chi")

PrimProd2100c_FRic <- glm(fric ~ PrimProd2100c, data = FD_RCP45, family = binomial)
visreg(PrimProd2100c_FRic, scale = "response", main = "", xlab = "Primary productivity", ylab = "P(FRic)")
PrimProd2100c_FRic_drop<- drop1(PrimProd2100c_FRic, test = "Chi")

IceCon2100c_FRic <- glm(fric ~ IceCon2100c, data = FD_RCP45, family = binomial)
visreg(IceCon2100c_FRic, scale = "response", main = "", xlab = "Ice concentration", ylab = "P(FRic)")
IceCon2100c_FRic_drop<- drop1(IceCon2100c_FRic, test = "Chi")

OxyMean_FRic <- glm(fric ~ OxyMean, data = FD_RCP45, family = binomial)
visreg(OxyMean_FRic, scale = "response", main = "", xlab = "Mean dissolved oxygen", ylab = "P(FRic)")
OxyMean_FRic_drop<- drop1(OxyMean_FRic, test = "Chi")

OxyB2100c_FRic <- glm(fric ~ OxyB2100c, data = FD_RCP45, family = binomial)
visreg(OxyB2100c_FRic, scale = "response", main = "", xlab = "Dissolved bottom oxygen", ylab = "P(FRic)")
OxyB2100c_FRic_drop<- drop1(OxyB2100c_FRic, test = "Chi")
par(mfrow=c(1,1))

AIC(MeanDepth_FRic,MinDepth_FRic,MaxDepth_FRic,SST2100c_FRic,SBT2100c_FRic,Salinity2100c_FRic,SalinityB2100c_FRic,
    PrimProd2100c_FRic,IceCon2100c_FRic,OxyMean_FRic,OxyB2100c_FRic)

## Functional uniqueness
par(mfrow=c(4,3))
MeanDepth_FUn <- glm(FUn ~ DepthMean, data = FD_RCP45, family = binomial)
visreg(MeanDepth_FUn, scale = "response", main = "", xlab = "Mean Depth", ylab = "P(FUn)")
MeanDepth_FUn_drop<- drop1(MeanDepth_FUn, test = "Chi")

MinDepth_FUn <- glm(FUn ~ DepthMin, data = FD_RCP45, family = binomial)
visreg(MinDepth_FUn, scale = "response", main = "", xlab = "Min Depth", ylab = "P(FUn)")
MinDepth_FUn_drop<- drop1(MinDepth_FUn, test = "Chi")

MaxDepth_FUn <- glm(FUn ~ DepthMax, data = FD_RCP45, family = binomial)
visreg(MaxDepth_FUn, scale = "response", main = "", xlab = "Max Depth", ylab = "P(FUn)")
MaxDepth_FUn_drop<- drop1(MaxDepth_FUn, test = "Chi")

SST2100c_FUn <- glm(FUn ~ SST2100c, data = FD_RCP45, family = binomial)
visreg(SST2100c_FUn, scale = "response", main = "", xlab = "SST", ylab = "P(FUn)")
SST2100c_FUn_drop<- drop1(SST2100c_FUn, test = "Chi")

SBT2100c_FUn <- glm(FUn ~ SBT2100c, data = FD_RCP45, family = binomial)
visreg(SBT2100c_FUn, scale = "response", main = "", xlab = "SBT", ylab = "P(FUn)")
SBT2100c_FUn_drop<- drop1(SBT2100c_FUn, test = "Chi")

Salinity2100c_FUn <- glm(FUn ~ Salinity2100c, data = FD_RCP45, family = binomial)
visreg(Salinity2100c_FUn, scale = "response", main = "", xlab = "Salinity", ylab = "P(FUn)")
Salinity2100c_FUn_drop<- drop1(Salinity2100c_FUn, test = "Chi")

SalinityB2100c_FUn <- glm(FUn ~ SalinityB2100c, data = FD_RCP45, family = binomial)
visreg(SalinityB2100c_FUn, scale = "response", main = "", xlab = "Bottom Salinity", ylab = "P(FUn)")
SalinityB2100c_FUn_drop<- drop1(SalinityB2100c_FUn, test = "Chi")

PrimProd2100c_FUn <- glm(FUn ~ PrimProd2100c, data = FD_RCP45, family = binomial)
visreg(PrimProd2100c_FUn, scale = "response", main = "", xlab = "Primary productivity", ylab = "P(FUn)")
PrimProd2100c_FUn_drop<- drop1(PrimProd2100c_FUn, test = "Chi")

IceCon2100c_FUn <- glm(FUn ~ IceCon2100c, data = FD_RCP45, family = binomial)
visreg(IceCon2100c_FUn, scale = "response", main = "", xlab = "Ice concentration", ylab = "P(FUn)")
IceCon2100c_FUn_drop<- drop1(IceCon2100c_FUn, test = "Chi")

OxyMean_FUn <- glm(FUn ~ OxyMean, data = FD_RCP45, family = binomial)
visreg(OxyMean_FUn, scale = "response", main = "", xlab = "Mean dissolved oxygen", ylab = "P(FUn)")
OxyMean_FUn_drop<- drop1(OxyMean_FUn, test = "Chi")

OxyB2100c_FUn <- glm(FUn ~ OxyB2100c, data = FD_RCP45, family = binomial)
visreg(OxyB2100c_FUn, scale = "response", main = "", xlab = "Dissolved bottom oxygen", ylab = "P(FUn)")
OxyB2100c_FUn_drop<- drop1(OxyB2100c_FUn, test = "Chi")
par(mfrow=c(1,1))

AIC(MeanDepth_FUn,MinDepth_FUn,MaxDepth_FUn,SST2100c_FUn,SBT2100c_FUn,Salinity2100c_FUn,SalinityB2100c_FUn,
    PrimProd2100c_FUn,IceCon2100c_FUn,OxyMean_FUn,OxyB2100c_FUn)

## Functional specialisation
par(mfrow=c(4,3))
MeanDepth_FSp <- glm(fspe ~ DepthMean, data = FD_RCP45, family = binomial)
visreg(MeanDepth_FSp, scale = "response", main = "", xlab = "Mean Depth", ylab = "P(FSp)")
MeanDepth_FSp_drop<- drop1(MeanDepth_FSp, test = "Chi")

MinDepth_FSp <- glm(fspe ~ DepthMin, data = FD_RCP45, family = binomial)
visreg(MinDepth_FSp, scale = "response", main = "", xlab = "Min Depth", ylab = "P(FSp)")
MinDepth_FSp_drop<- drop1(MinDepth_FSp, test = "Chi")

MaxDepth_FSp <- glm(fspe ~ DepthMax, data = FD_RCP45, family = binomial)
visreg(MaxDepth_FSp, scale = "response", main = "", xlab = "Max Depth", ylab = "P(FSp)")
MaxDepth_FSp_drop<- drop1(MaxDepth_FSp, test = "Chi")

SST2100c_FSp <- glm(fspe ~ SST2100c, data = FD_RCP45, family = binomial)
visreg(SST2100c_FSp, scale = "response", main = "", xlab = "SST", ylab = "P(FSp)")
SST2100c_FSp_drop<- drop1(SST2100c_FSp, test = "Chi")

SBT2100c_FSp <- glm(fspe ~ SBT2100c, data = FD_RCP45, family = binomial)
visreg(SBT2100c_FSp, scale = "response", main = "", xlab = "SBT", ylab = "P(FSp)")
SBT2100c_FSp_drop<- drop1(SBT2100c_FSp, test = "Chi")

Salinity2100c_FSp <- glm(fspe ~ Salinity2100c, data = FD_RCP45, family = binomial)
visreg(Salinity2100c_FSp, scale = "response", main = "", xlab = "Salinity", ylab = "P(FSp)")
Salinity2100c_FSp_drop<- drop1(Salinity2100c_FSp, test = "Chi")

SalinityB2100c_FSp <- glm(fspe ~ SalinityB2100c, data = FD_RCP45, family = binomial)
visreg(SalinityB2100c_FSp, scale = "response", main = "", xlab = "Bottom Salinity", ylab = "P(FSp)")
SalinityB2100c_FSp_drop<- drop1(SalinityB2100c_FSp, test = "Chi")

PrimProd2100c_FSp <- glm(fspe ~ PrimProd2100c, data = FD_RCP45, family = binomial)
visreg(PrimProd2100c_FSp, scale = "response", main = "", xlab = "Primary productivity", ylab = "P(FSp)")
PrimProd2100c_FSp_drop<- drop1(PrimProd2100c_FSp, test = "Chi")

IceCon2100c_FSp <- glm(fspe ~ IceCon2100c, data = FD_RCP45, family = binomial)
visreg(IceCon2100c_FSp, scale = "response", main = "", xlab = "Ice concentration", ylab = "P(FSp)")
IceCon2100c_FSp_drop<- drop1(IceCon2100c_FSp, test = "Chi")

OxyMean_FSp <- glm(fspe ~ OxyMean, data = FD_RCP45, family = binomial)
visreg(OxyMean_FSp, scale = "response", main = "", xlab = "Mean dissolved oxygen", ylab = "P(FSp)")
OxyMean_FSp_drop<- drop1(OxyMean_FSp, test = "Chi")

OxyB2100c_FSp <- glm(fspe ~ OxyB2100c, data = FD_RCP45, family = binomial)
visreg(OxyB2100c_FSp, scale = "response", main = "", xlab = "Dissolved bottom oxygen", ylab = "P(FSp)")
OxyB2100c_FSp_drop<- drop1(OxyB2100c_FSp, test = "Chi")
par(mfrow=c(1,1))

AIC(MeanDepth_FSp,MinDepth_FSp,MaxDepth_FSp,SST2100c_FSp,SBT2100c_FSp,Salinity2100c_FSp,SalinityB2100c_FSp,
    PrimProd2100c_FSp,IceCon2100c_FSp,OxyMean_FSp,OxyB2100c_FSp)