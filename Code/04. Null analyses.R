###############################################################################################################
# 04. Null analyses
## This R code forms the null model based on simulated extinction times
## based on Crit E EX mode
###############################################################################################################

# Import packages
library(tidyverse)
library(tibble)
library(reshape2)
library(mFD)
library(data.table)
library(geometry)
library(purrr)
library(doParallel)

# Source functions
source("C:/Users/2022207/Dropbox/Jack's PhD/Chapter 3. Future shark FD/Analyses/R code/Final pipeline R code/Functions/get_indicator_function 2.R")
source("C:/Users/2022207/Dropbox/Jack's PhD/Chapter 3. Future shark FD/Analyses/R code/Final pipeline R code/Functions/fonction_FRIC_Global_full.R")

# Load data
load(file="C:/Users/2022207/Dropbox/Jack's PhD/Chapter 3. Future shark FD/Data/R data/Data.RData")
# Load median results for spatial results
load(file = "C:/Users/2022207/Dropbox/Jack's PhD/Chapter 3. Future shark FD/Analyses/Final paper R codes/data/Crit_E/CritE_median.RData")
# Load IUCN status for AT scenario
load(file = "C:/Users/2022207/Dropbox/Jack's PhD/Chapter 3. Future shark FD/Data/R data/Species_list_IUCN.RData")
# Load extinction times corrected for synonyms
load(file = "C:/Users/2022207/Dropbox/Jack's PhD/Chapter 3. Future shark FD/Analyses/Final paper R codes/data/Crit_E/Extinction_times.RData")

# Parallelised loop
registerDoParallel(cores = 5)
getDoParWorkers()
res.null<-NULL
res.null<- lapply(1:10, function(x){
  sharks_traits_null <- data %>% 
    select(Species,habitat,vertical,terrestriality,thermo,feeding,diet,max.length) %>% 
    column_to_rownames(var = "Species")
  
  # Calculate the actual number of species lost for `Future` and `AT` scenarios in the simulation
  sim_col <- paste0("sim_", x)
  
  # Identify extinct species (extinction times <= 77) for the current simulation
  ext_species <- dat_ext %>%
    select(Species, !!sym(sim_col)) %>%
    filter(!!sym(sim_col) <= 77) %>%
    pull(Species)  # Get list of species marked as extinct
  
  # Set up present and 2100 scenarios
  Resamp <- data %>% 
    select(Species) %>% 
    mutate(Present = 1)
  
  Resamp_2100 <- data %>%
    select(Species) %>%
    mutate(
      Future = ifelse(Species %in% ext_species, 0, 1)  # Mark extinct species as 0, others as 1
    )
  
  # Randomise extinctions
  num_extinct_future <- sum(Resamp_2100$Future == 0)
  Resamp_2100 <- Resamp_2100 %>%
    mutate(Future = sample(c(rep(0, num_extinct_future), rep(1, nrow(Resamp) - num_extinct_future)), size = nrow(Resamp)))
  
  Resamp_AT <- Resamp %>% 
    sample_n(682, replace = FALSE) %>%
    mutate(AT = 1)
  
  null_data <- Resamp %>%
    left_join(Resamp_2100, by = "Species")
  
  null_data <- null_data %>%
    left_join(Resamp_AT %>% select(Species,AT), by = "Species") %>%
    mutate(AT = ifelse(is.na(AT), 0, 1)) %>% 
    select(Species,Present,Future,AT)
  
  # Form final occurrence matrix
  baskets_null <- null_data %>%
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
  
  # make assemblage data into a matrix & ensure all matrix values are numeric
  baskets_sharks_weights_null <- data.matrix(baskets_sharks_weights_null, rownames.force = NA)
  class(baskets_sharks_weights_null) <- "numeric"
  
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
  
  # Transform dataframe from presence-absence matrix in a list
  test_null <- apply(baskets_sharks_weights_null,1,function(x){colnames(baskets_sharks_weights_null)[which(x==1)]})
  
  # Calculate FUn
  lala <-lapply(test_null, function(y){get_FU_region (Mat_dist=as.matrix(sp_dist_sharks_null), nb_NN=5, Selected_sp=y)})
  
  # Calculate FD metrics from function
  pcoa <- list(li = sp_faxes_coord_sharks_null)
  ax <- c("PC1", "PC2", "PC3", "PC4")
  coord_d <- pcoa$li[, ax]
  
  lala1<-lapply(test_null, function(y){ get_FV_Sp (ax=c(1:4), pcoa=pcoa, Selected_sp=y)})
  
  # Calculate all FD metrics per scenario
  FD.metrics_null <- data.frame(
    Present = rbind(lala1$Present$data$RS,
                    lala1$Present$data$Fric_r,
                    lala1$Present$specialization[1],
                    lala$Present$FU_Region_Mean),
    Future = rbind(lala1$Future$data$RS,
                   lala1$Future$data$Fric_r,
                   lala1$Future$specialization[1],
                   lala$Future$FU_Region_Mean),
    AT = rbind(lala1$AT$data$RS,
               lala1$AT$data$Fric_r,
               lala1$AT$specialization[1],
               lala$AT$FU_Region_Mean)
  )
  
  rownames(FD.metrics_null) <- c("sp_richn","fric","fsp","fun")
  colnames(FD.metrics_null) <- c("Present","Future","AT")
  
  FDindices_taxon_null <- t(FD.metrics_null) %>% as.data.frame()
  
  FDindices_taxon_null <- FDindices_taxon_null %>% 
    rownames_to_column(var = "Scenario")
  
  FDindices_taxon_null$Scenario <- ordered(FDindices_taxon_null$Scenario, 
                                           levels=c("Present","Future","AT"))
  
  # Form list to merge datasets
  FDind.var_null = list(FDindices_taxon_null)
  
  # Output
  FDind.var_null %>% 
    reduce(inner_join, by = "Scenario")
  
  
  # Close loop
})

stopImplicitCluster()

# Merge lists into 1 dataframe
res_df_null <- res.null %>% 
  bind_rows()

# Format and save results
Null_FDmetrics_taxonvar.CritE<- res_df_null %>% 
  select(Scenario:fun)
save(Null_FDmetrics_taxonvar.CritE, file = "C:/Users/2022207/Dropbox/Jack's PhD/Chapter 3. Future shark FD/Analyses/Final paper R codes/data/Crit_E/Null_wide.RData")

# Melt data
FDmetrics_null_long_taxonvar.CritE<- melt(Null_FDmetrics_taxonvar.CritE, id.vars= "Scenario")
save(FDmetrics_null_long_taxonvar.CritE, file = "C:/Users/2022207/Dropbox/Jack's PhD/Chapter 3. Future shark FD/Analyses/Final paper R codes/data/Crit_E/Null_long.RData")
