###############################################################################################################
# 06. Null model
## This R code calculates null distribution of functional diversity analyses
## This code is designed to be run in the supercomputer; but is run in 2 halves (i.e., 1:500)
## it produces Rdata files for loading null models into Figure 1
###############################################################################################################

# Import packages
library(tidyverse)
library(ade4)
library(tibble)
library(scales)
library(cowplot)
library(doBy)
library(reshape2)
library(RColorBrewer)
library(mFD)
library(ggsci)
library(data.table)
library(foreach)
library(parallel)
library(doParallel)
library(funrar)

# Load data
load(file="C:/Users/2022207/Dropbox/Jack's PhD/Chapter 3. Future shark FD/Data/R data/Data.RData")

# Null model
registerDoParallel(cores = 10)
getDoParWorkers()
res.null<-NULL
res.null<- lapply(1:1000, function(x){
  # Species-trait matrix
  sharks_traits_null <- data %>% 
    select(Species,habitat,vertical,terrestriality,thermo,feeding,diet,max.length) %>% 
    column_to_rownames(var = "Species")
  
  # Sample species
  Resamp <- data %>% 
    select(Species)
  
  Resamp_2100 <- Resamp %>% 
    sample_n(1034, replace = FALSE)
  
  Resamp_Y100 <- Resamp %>% 
    sample_n(1004, replace = FALSE)
  
  Resamp_Y200 <- Resamp %>% 
    sample_n(887, replace = FALSE)
  
  Resamp_Y300 <- Resamp %>% 
    sample_n(786, replace = FALSE)
  
  Resamp_Y400 <- Resamp %>% 
    sample_n(654, replace = FALSE)
  
  Resamp_Y500 <- Resamp %>% 
    sample_n(554, replace = FALSE)
  
  Resamp$Present <- 1
  
  Resamp <- Resamp %>% 
    mutate(Future = Species) %>% 
    mutate(Y100 = Species) %>% 
    mutate(Y200 = Species) %>% 
    mutate(Y300 = Species) %>% 
    mutate(Y400 = Species) %>% 
    mutate(Y500 = Species) %>% 
    rowwise() %>% 
    mutate(Future = replace(Future, !Future %in% Resamp_2100$Species, 0)) %>% 
    mutate(Future = replace(Future,Future!=0,1)) %>% 
    mutate(Y100 = replace(Y100, !Y100 %in% Resamp_Y100$Species, 0)) %>% 
    mutate(Y100 = replace(Y100,Y100!=0,1)) %>% 
    mutate(Y200 = replace(Y200, !Y200 %in% Resamp_Y200$Species, 0)) %>% 
    mutate(Y200 = replace(Y200,Y200!=0,1)) %>% 
    mutate(Y300 = replace(Y300, !Y300 %in% Resamp_Y300$Species, 0)) %>% 
    mutate(Y300 = replace(Y300,Y300!=0,1)) %>% 
    mutate(Y400 = replace(Y400, !Y400 %in% Resamp_Y400$Species, 0)) %>% 
    mutate(Y400 = replace(Y400,Y400!=0,1)) %>% 
    mutate(Y500 = replace(Y500, !Y500 %in% Resamp_Y500$Species, 0)) %>% 
    mutate(Y500 = replace(Y500,Y500!=0,1))
  
  # Form final occurrence matrix
  baskets_null <- Resamp %>%
    column_to_rownames(var = "Species")
  
  # Transpose to form final occurrence matrix
  baskets_sharks_weights_null <- t(baskets_null)

 # Form trait category matrix
  sharks_traits_cat_null <- tibble(trait_name = c("habitat","vertical","terrestriality","thermo","feeding","diet","max.length"),
                                   trait_type = c("N", "N", "O", "O", "O", "N", "Q"))
  
  # Summarise dataset
  sharks_traits_summ_null <- mFD::sp.tr.summary(
    tr_cat     = sharks_traits_cat_null,   
    sp_tr      = sharks_traits_null, 
    stop_if_NA = TRUE)
  
  # Reveal number of states for each trait
  sharks_traits_summ_null$tr_summary_list
  
  # make assemblage data into a matrix & ensure all matrix values are numeric
  baskets_sharks_weights_null <- data.matrix(baskets_sharks_weights_null, rownames.force = NA)
  class(baskets_sharks_weights_null) <- "numeric"
  
  #summarise dataset
  asb_sp_sharks_summ_null <- mFD::asb.sp.summary(asb_sp_w = baskets_sharks_weights_null)
  
  # retrieve species occurrences for all assemblages:
  asb_sp_sharks_occ_null <- asb_sp_sharks_summ_null$asb_sp_occ
  
  # Construct trait distance matrix using species
  sp_dist_sharks_null <- mFD::funct.dist(
    sp_tr         = sharks_traits_null,
    tr_cat        = sharks_traits_cat_null,
    metric        = "gower",
    scale_euclid  = "scale_center",
    ordinal_var   = "classic",
    weight_type   = "equal",
    stop_if_NA    = TRUE)
  
  # Assess quality
  fspaces_quality_sharks_null <- mFD::quality.fspaces(
    sp_dist             = sp_dist_sharks_null,
    maxdim_pcoa         = 10,
    deviation_weighting = c("absolute","squared"),
    fdist_scaling       = FALSE,
    fdendro             = "average")
  
  # Return coordinates of each axis
  sp_faxes_coord_sharks_null <- fspaces_quality_sharks_null$"details_fspaces"$"sp_pc_coord"
  
  # Calculate FD metrics
  alpha_fd_indices_sharks_null <- mFD::alpha.fd.multidim(
    sp_faxes_coord   = sp_faxes_coord_sharks_null[ , c("PC1", "PC2", "PC3", "PC4")],
    asb_sp_w         = baskets_sharks_weights_null,
    ind_vect         = c("fric", "fori", "fspe"),
    scaling          = TRUE,
    check_input      = TRUE,
    details_returned = TRUE)
  
  fd_ind_values_sharks_null <- alpha_fd_indices_sharks_null$"functional_diversity_indices"
  
  ## Form dataframe
  FDindices_taxon_null <- as.data.frame(fd_ind_values_sharks_null) %>% 
    tibble::rownames_to_column("Scenario")
  
  # Calculate functional uniqueness per species
  Funct_Uni_null <- uniqueness(baskets_sharks_weights_null, as.matrix(sp_dist_sharks_null)) %>% 
    as.data.frame()
  Funct_Uni_null$Present <- baskets_null$Present
  Funct_Uni_null$Future <- baskets_null$Future
  Funct_Uni_null$Y100 <- baskets_null$Y100
  Funct_Uni_null$Y200 <- baskets_null$Y200
  Funct_Uni_null$Y300 <- baskets_null$Y300
  Funct_Uni_null$Y400 <- baskets_null$Y400
  Funct_Uni_null$Y500 <- baskets_null$Y500
  
  # Calculate functional uniqueness per scenario
  Present_FUn_null <- mean(Funct_Uni_null$Ui)
  Future_FUn_null_df <- Funct_Uni_null %>% 
    filter(Future == "1") %>% 
    select(Ui)
  Future_FUn_null <- mean(Future_FUn_null_df$Ui)
  
  Y100_FUn_null_df <- Funct_Uni_null%>% 
    filter(Y100 == "1") %>% 
    select(Ui)
  Y100_FUn_null <- mean(Y100_FUn_null_df$Ui)
  
  Y200_FUn_null_df <- Funct_Uni_null %>% 
    filter(Y200 == "1") %>% 
    select(Ui)
  Y200_FUn_null <- mean(Y200_FUn_null_df$Ui)
  
  Y300_FUn_null_df <- Funct_Uni_null %>% 
    filter(Y300 == "1") %>% 
    select(Ui)
  Y300_FUn_null <- mean(Y300_FUn_null_df$Ui)
  
  Y400_FUn_null_df <- Funct_Uni_null %>% 
    filter(Y400 == "1") %>% 
    select(Ui)
  Y400_FUn_null <- mean(Y400_FUn_null_df$Ui)
  
  Y500_FUn_null_df <- Funct_Uni_null %>% 
    filter(Y500 == "1") %>% 
    select(Ui)
  Y500_FUn_null <- mean(Y500_FUn_null_df$Ui)
  
  # Add uniqueness to final FD results data frame
  Fun_null <- data.frame(
    Fun = c(Present_FUn_null,Future_FUn_null,Y100_FUn_null,
            Y200_FUn_null,Y300_FUn_null,Y400_FUn_null,Y500_FUn_null),
    row.names = c("Present","Future","Y100","Y200","Y300","Y400","Y500")
  )
  
  FDindices_taxon_null$fun <- Fun_null$Fun
  FDindices_taxon_null$Scenario <- ordered(FDindices_taxon_null$Scenario, 
                                           levels=c("Present","Future","Y100","Y200","Y300","Y400","Y500"))
  
  # Form list to merge datasets
  FDind.var_null = list(FDindices_taxon_null)
  
  # Output
  FDind.var_null %>% 
    reduce(inner_join, by = "Scenario")
  
  #close loop
})

stopImplicitCluster()

# Merge lists into 1 dataframe
res_df_null <- res.null %>% 
  bind_rows()

# Format data
Null_FDmetrics_taxonvar<- res_df_null %>% 
  select(Scenario:fun)
save(Null_FDmetrics_taxonvar, file = "C:/Users/2022207/Dropbox/Jack's PhD/Chapter 3. Future shark FD/Data/R data/Future_null_metrics.RData")

# Melt data
FDmetrics_null_long_taxonvar<- melt(Null_FDmetrics_taxonvar, id.vars= "Scenario")
save(FDmetrics_null_long_taxonvar, file = "C:/Users/2022207/Dropbox/Jack's PhD/Chapter 3. Future shark FD/Data/R data/Future_null_long_metrics.RData")

# Form dataframe of mean, median and standard deviation of all FD metrics
## Spp_richness + Median and SD values
Future_Null_var <- Null_FDmetrics_taxonvar %>% 
  group_by(Scenario) %>%
  summarise(Sp_mean = mean(sp_richn),
            Sp_med = median(sp_richn),
            Sp_sd = sd(sp_richn),
            FRic_mean = mean(fric),
            FRic_med = median(fric),
            FRic_sd = sd(fric),
            Fori_mean = mean(fori),
            Fori_med = median(fori),
            Fori_sd = sd(fori),
            Fspe_mean = mean(fspe),
            Fspe_med = median(fspe),
            Fspe_sd = sd(fspe),
            Fun_mean = mean(fun),
            Fun_med = median(fun),
            Fun_sd = sd(fun))

Future_Null_var$Scenario <- ordered(Future_Null_var$Scenario, levels=c("Present","Future","Y100","Y200","Y300","Y400","Y500"))

# Save iteration data
save(Future_Null_var, file = "C:/Users/2022207/Dropbox/Jack's PhD/Chapter 3. Future shark FD/Data/R data/Future_Mean_Null_metrics.RData")
