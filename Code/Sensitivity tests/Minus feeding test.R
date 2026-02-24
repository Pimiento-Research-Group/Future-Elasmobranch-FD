###############################################################################################################
# Sensitivity tests 5: Minus feeding
## This R code conducts senstivity tests based on individual trait removal
## This code produces all sensitivity tests for supplementary material
###############################################################################################################

# Import packages
library(tidyverse)
library(devtools)
library(rredlist)
library(reticulate)
library(iucnsim)
library(tibble)
library(here)
library(readxl)
library(scales)
library(cowplot)
library(doBy)
library(reshape2)
library(RColorBrewer)
library(mFD)
library(ggsci)
library(data.table)
library(geometry)
library(viridis)
library(purrr)
library(foreach)
library(progress)
library(parallel)
library(doParallel)

# Source functions
source("C:/Users/Jack Cooper/Dropbox/Jack's PhD/Chapter 3. Future shark FD/Analyses/R code/Final pipeline R code/Functions/get_indicator_function 2.R")
source("C:/Users/Jack Cooper/Dropbox/Jack's PhD/Chapter 3. Future shark FD/Analyses/R code/Final pipeline R code/Functions/fonction_FRIC_Global_full.R")

#Import data
load(file="C:/Users/Jack Cooper/Documents/iucn/Data.RData")

# Load median results for spatial results
load(file = "C:/Users/Jack Cooper/Documents/iucn/iucn_sim/iucn_data/CritE_median_D24.RData")

# Load synonyms for matching up names in loop
Synonyms <- read_xlsx("C:/Users/Jack Cooper/Dropbox/Jack's PhD/Chapter 3. Future shark FD/Analyses/DatasetSynonyms.xlsx")

# Load extinction times
load(file = "C:/Users/Jack Cooper/Documents/iucn/iucn_sim/iucn_data/Extinction_times.RData")

## Aquamaps coverage
# Load Aquamaps and isolate species list
AqMap <- read_csv("C:/Users/Jack Cooper/Dropbox/Jack's PhD/Chapter 3. Future shark FD/Data/Aquamaps data/Elasmo_current.csv")

AqMap <- AqMap %>% 
  mutate(sci_name = paste(Genus, Species, sep = " "))

AqMapSpp <- AqMap %>% 
  select(sci_name) %>% 
  unique()

# Correct synonyms
AqMap <- AqMap %>%
  left_join(Synonyms, by = c("sci_name" = "Aquamaps_name")) %>%
  mutate(sci_name = ifelse(!is.na(Accepted_name), Accepted_name, sci_name)) %>%
  select(-c(Accepted_name,Marked_iucnsim,iucnsim_name,Marked_aquamaps)) %>% 
  distinct(sci_name, .keep_all = TRUE)

AqMapSpp <- AqMapSpp %>%
  left_join(Synonyms, by = c("sci_name" = "Aquamaps_name")) %>%
  mutate(sci_name = ifelse(!is.na(Accepted_name), Accepted_name, sci_name)) %>%
  select(-c(Accepted_name,Marked_iucnsim,iucnsim_name,Marked_aquamaps)) %>% 
  distinct(sci_name, .keep_all = TRUE)

# Conduct coverage analyses
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
  select(Species,habitat,vertical,terrestriality,thermo,
         diet,max.length,deepwater,body.plan,reproduction)

sharks_traits <- traits %>% 
  remove_rownames %>% 
  column_to_rownames(var="Species")

# Form trait category matrix
sharks_traits_cat <- tibble(trait_name = c("habitat","vertical","terrestriality","thermo",
                                           "diet","max.length","deepwater","body.plan","reproduction"),
                            trait_type = c("N", "N", "O", "O", "N", "Q", "N", "N", "N"))

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

# Return coordinates of each axis
sp_faxes_coord_sharks <- fspaces_quality_sharks$"details_fspaces"$"sp_pc_coord"

# Calculate FD metrics
alpha_fd_indices_sharks <- mFD::alpha.fd.multidim(
  sp_faxes_coord   = sp_faxes_coord_sharks[ , c("PC1", "PC2", "PC3", "PC4")],
  asb_sp_w         = baskets_sharks_weights,
  ind_vect         = c("fric"),
  scaling          = TRUE,
  check_input      = TRUE,
  details_returned = TRUE)

fd_ind_values_sharks <- alpha_fd_indices_sharks$"functional_diversity_indices"

# Form dataframe
FDindices_taxon <- as.data.frame(fd_ind_values_sharks) %>% 
  tibble::rownames_to_column("Scenario")

## Simulation sensitivity tests
# Parallelised loop for simulations boxplot
registerDoParallel(cores = 8)
getDoParWorkers()
res.TaxonVar<-NULL
res.TaxonVar<- lapply(1:10000, function(x){
  # Identify the current simulation column dynamically
  sim_col <- paste0("sim_", x)
  
  # Identify extinct species (extinction times <= 77) for the current simulation
  ext_species <- dat_ext %>%
    select(Species, !!sym(sim_col)) %>%
    filter(!!sym(sim_col) <= 77) %>%
    pull(Species)  # Get list of species marked as extinct
  
  # Set up present and 2100 scenarios
  Occ_Future <- data %>%
    select(Species,Order) %>%
    mutate(
      Present = 1,  # Set all to present
      Future = ifelse(Species %in% ext_species, 0, 1)  # Mark extinct species as 0, others as 1
    )
  
  baskets <- Occ_Future %>%
    select(Species,Present,Future) %>% 
    remove_rownames %>% 
    column_to_rownames(var="Species")
  
  # Transpose to form final occurrence matrix
  baskets_sharks_weights <- t(baskets)
  
  # Species-trait matrix
  traits <- data %>% 
    select(Species,habitat,vertical,terrestriality,thermo,
           diet,max.length,deepwater,body.plan,reproduction)
  
  sharks_traits <- traits %>% 
    remove_rownames %>% 
    column_to_rownames(var="Species")
  
  # Form trait category matrix
  sharks_traits_cat <- tibble(trait_name = c("habitat","vertical","terrestriality","thermo",
                                             "diet","max.length","deepwater","body.plan","reproduction"),
                              trait_type = c("N", "N", "O", "O", "N", "Q", "N", "N", "N"))
  
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
  # Form list to merge datasets
  FDind.var = list(FDindices_taxon.CritE)
  
  # Output
  FDind.var %>% 
    reduce(inner_join, by = "Scenario")
  
  # Close loop
})

stopImplicitCluster()

# Merge lists into 1 dataframe
res_df <- res.TaxonVar %>% 
  bind_rows()

# Format and save results
FDmetrics_taxonvar.feeding<- res_df %>% 
  select(Scenario:fun)
save(FDmetrics_taxonvar.feeding,file = "C:/Users/Jack Cooper/Documents/iucn/iucn_sim/Robustness tests/Minus_FM/CritE_full_buffer.RData")

FDmetrics_long_TaxonVar.feeding<- melt(FDmetrics_taxonvar.feeding, id.vars= "Scenario")
save(FDmetrics_long_TaxonVar.feeding,file = "C:/Users/Jack Cooper/Documents/iucn/iucn_sim/Robustness tests/Minus_FM/CritE_full_long_buffer.RData")
