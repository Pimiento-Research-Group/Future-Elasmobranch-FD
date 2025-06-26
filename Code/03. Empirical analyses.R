###############################################################################################################
# 03. Empirical analyses
## This R code quantifies functional diversity for each extinction scenario
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
source("C:/Users/2022207/Dropbox/Jack's PhD/Chapter 3. Future shark FD/Analyses/R code/Final pipeline R code/Functions/get_indicator_function 2.R")
source("C:/Users/2022207/Dropbox/Jack's PhD/Chapter 3. Future shark FD/Analyses/R code/Final pipeline R code/Functions/fonction_FRIC_Global_full.R")

# Load data
load(file="C:/Users/2022207/Dropbox/Jack's PhD/Chapter 3. Future shark FD/Data/R data/Data.RData")
# Load median results for spatial results
load(file = "C:/Users/2022207/Dropbox/Jack's PhD/Chapter 3. Future shark FD/Analyses/Final paper R codes/data/Crit_E/CritE_median.RData")
# Load IUCN status for AT scenario
load(file = "C:/Users/2022207/Dropbox/Jack's PhD/Chapter 3. Future shark FD/Data/R data/Species_list_IUCN.RData")
# Load synonyms for matching up names in loop
Synonyms <- read_xlsx("C:/Users/2022207/Dropbox/Jack's PhD/Chapter 3. Future shark FD/Analyses/DatasetSynonyms.xlsx")

# Produce a dataframe to determine which species are extinct in the future
dat.ext <- read_delim("C:/Users/2022207/Dropbox/Jack's PhD/Chapter 3. Future shark FD/Analyses/Final paper R codes/data/Crit_E/iucn_sim/future_simulations_future_status/te_all_species.txt")

# Make dat.ext a dataframe & tidy
ext.df <- dat.ext %>% as.data.frame()

ext.df1 <- rbind(colnames(ext.df), ext.df)

ext.df1[1, ] <- lapply(ext.df1[1, ], function(x) ifelse(is.na(as.numeric(x)), NA, as.numeric(x)))

ext.df1[1, 1] <- "Acroteriobatus annulatus"
names(ext.df1)[1] <- "Species"
simulation_names <- paste0("sim_", seq(1, ncol(ext.df1) - 1))
names(ext.df1)[2:ncol(ext.df1)] <- simulation_names
ext.df1 <- as.data.frame(ext.df1)

# Correct for synonyms
ext.df1 <- ext.df1 %>% 
  left_join(Synonyms, by = c("Species" = "iucnsim_name")) %>%
  mutate(Species = ifelse(!is.na(Accepted_name), Accepted_name, Species)) %>%
  select(-c(Accepted_name,Marked_iucnsim,Aquamaps_name,Marked_aquamaps)) %>% 
  distinct(Species, .keep_all = TRUE)

# Show that you can calculate no of extinctions per simulation (sim 1 as example)
ext.df1 %>% 
  as_tibble() %>% 
  select(1,2) %>% 
  mutate(sim_1 = as.double(sim_1)) %>% 
  filter(sim_1 <= 77)

# Calculate no. of extinctions per simulation (per column; to be checked in loop) and calculate mean/median
ext.df1[ , -1] <- lapply(ext.df1[ , -1], as.numeric)
extinction_counts <- sapply(ext.df1[ , -1], function(col) sum(col <+ 77, na.rm = TRUE))
ext.df1 <- rbind(ext.df1, c("extinctions", extinction_counts))
extinction_mean <- mean(extinction_counts)
extinction_median <- median(extinction_counts)

# Set up data for looping
dat_ext <- ext.df1 %>% 
  as_tibble() %>% 
  mutate(across(-Species, as.double)) 

# Save dat_ext for supercomputer
save(dat_ext,file = "C:/Users/2022207/Dropbox/Jack's PhD/Chapter 3. Future shark FD/Analyses/Final paper R codes/data/Crit_E/Extinction_times.RData")

# Parallelised loop
registerDoParallel(cores = 5)
getDoParWorkers()
res.TaxonVar<-NULL
res.TaxonVar<- lapply(1:10, function(x){
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
  
  # Set up IUCN AT scenario in Occ_Future (including IUCN buffer)
  iucn <- as_tibble(species_iucn_updated_fixed) %>%
    rename(iucn = category) %>%
    dplyr::select(Species, iucn) %>%
    mutate(iucn=replace(iucn, is.na(iucn), "NE")) %>%
    mutate(iucn=replace(iucn, iucn=="CR",4))%>%
    mutate(iucn=replace(iucn, iucn=="EN",3))%>%  
    mutate(iucn=replace(iucn, iucn=="VU",2))%>%
    mutate(iucn=replace(iucn, iucn=="NT",1))%>%
    mutate(iucn=replace(iucn, iucn=="LC",0))%>%
    mutate(iucn=replace(iucn, iucn=="DD",NA))%>%
    mutate(iucn=replace(iucn, iucn=="NE",NA))
  iucn$iucn <- as.numeric(iucn$iucn)
  
  Occ_Future <- Occ_Future %>%
    left_join(iucn, by = "Species")
  
  # Infer IUCN - based on Pimiento et al. 2020
  iucn_distribution <- Occ_Future %>%
    filter(!is.na(iucn)) %>%
    group_by(Order, iucn) %>%
    tally() %>%
    group_by(Order) %>%
    mutate(prop = n / sum(n))
  
  sample_iucn <- function(order) {
    dist <- iucn_distribution %>%
      filter(Order == order) %>%
      select(iucn, prop)
    
    sample(dist$iucn, 1, prob = dist$prop)
  }
  
  Occ_Future <- Occ_Future %>%
    rowwise() %>%
    mutate(iucn = ifelse(is.na(iucn), sample_iucn(Order), iucn)) %>%
    ungroup()
  
  Occ_Future <- Occ_Future %>%
    mutate(AT = ifelse(iucn %in% c(2, 3, 4), 0, ifelse(iucn %in% c(0, 1), 1, NA))) %>%
    select(Species, Present, Future, AT)
  
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
                   lala$Future$FU_Region_Mean),
    AT = rbind(lala1$AT$data$RS,
               lala1$AT$data$Fric_r,
               lala1$AT$specialization[1],
               lala$AT$FU_Region_Mean)
  )
  
  rownames(FD.metrics) <- c("sp_richn","fric","fsp","fun")
  colnames(FD.metrics) <- c("Present","Future","AT")
  
  FDindices_taxon.CritE <- t(FD.metrics) %>% as.data.frame()
  
  FDindices_taxon.CritE <- FDindices_taxon.CritE %>% 
    rownames_to_column(var = "Scenario")
  
  FDindices_taxon.CritE$Scenario <- ordered(FDindices_taxon.CritE$Scenario, 
                                            levels=c("Present","Future","AT"))
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
FDmetrics_taxonvar.CritE<- res_df %>% 
  select(Scenario:fun)
save(FDmetrics_taxonvar.CritE,file = "C:/Users/2022207/Dropbox/Jack's PhD/Chapter 3. Future shark FD/Analyses/Final paper R codes/data/Crit_E/CritE_full_buffer.RData")

FDmetrics_long_TaxonVar.CritE<- melt(FDmetrics_taxonvar.CritE, id.vars= "Scenario")
save(FDmetrics_long_TaxonVar.CritE,file = "C:/Users/2022207/Dropbox/Jack's PhD/Chapter 3. Future shark FD/Analyses/Final paper R codes/data/Crit_E/CritE_full_long_buffer.RData")
