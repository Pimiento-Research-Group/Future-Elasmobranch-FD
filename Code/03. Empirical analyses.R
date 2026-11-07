###############################################################################################################
# 03. Empirical analyses
## This R code quantifies functional diversity for each extinction scenario
## This code is based on Crit E EX mode from iucnsim
## it produces Table S1
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
library(DescTools)
library(parallel)
library(doParallel)

# Source functions
source("C:/Users/2022207/Dropbox/Jack's PhD/Chapter 3. Future shark FD/Analyses/R code/Final pipeline R code/Functions/get_indicator_function 2.R")
source("C:/Users/2022207/Dropbox/Jack's PhD/Chapter 3. Future shark FD/Analyses/R code/Final pipeline R code/Functions/fonction_FRIC_Global_full.R")

# Load & save trait data
data <- read_xlsx(
  path = "C:/Users/Jack Cooper/Dropbox/Jack's PhD (1)/Chapter 3. Future shark FD/Data/Trait data.xlsx",
  sheet = 4)

data <- data %>%
  mutate(
    # Nominal factors
    habitat      = as.factor(habitat),
    vertical     = as.factor(vertical),
    diet         = as.factor(diet),
    deepwater    = as.factor(deepwater),
    body.plan    = as.factor(body.plan),
    reproduction = as.factor(reproduction),
    # Quantitative
    max.length   = as.numeric(max.length),
    iucn   = as.numeric(iucn),
    # Ordinal factors (with defined order)
    terrestriality = factor(terrestriality,
                            levels = c("no", "brackish", "freshwater"),
                            ordered = TRUE),
    thermo = factor(thermo,
                    levels = c("meso", "ecto"),
                    ordered = TRUE),
    feeding = factor(feeding,
                     levels = c("macropredator", "filter feeder"),
                     ordered = TRUE)
  )

save(data, file="C:/Users/Jack Cooper/Documents/iucn/Data.RData")

# Check for trait associations - produces Table S1
cor_traits <- data %>% 
  select(habitat,vertical,terrestriality,thermo,feeding,diet,max.length,
         deepwater,body.plan,reproduction)

# Function for correlation ratio (eta)
eta <- function(num, cat) {
  
  keep <- complete.cases(num, cat)
  num <- num[keep]
  cat <- as.factor(cat[keep])
  
  grand_mean <- mean(num)
  
  ss_between <- sum(
    tapply(num, cat, length) *
      (tapply(num, cat, mean) - grand_mean)^2
  )
  
  ss_total <- sum((num - grand_mean)^2)
  
  sqrt(ss_between / ss_total)
}

# Function choosing association measure according to variable type
assoc <- function(x, y) {
  
  # Numeric x numeric
  if (is.numeric(x) && is.numeric(y)) {
    return(abs(cor(x, y, method = "spearman", use = "complete.obs")))
  }
  
  # Numeric x categorical
  if (is.numeric(x) && !is.numeric(y)) {
    return(eta(x, y))
  }
  
  # Categorical x numeric
  if (!is.numeric(x) && is.numeric(y)) {
    return(eta(y, x))
  }
  
  # Categorical x categorical
  if (!is.numeric(x) && !is.numeric(y)) {
    return(
      CramerV(
        table(x, y),
        conf.level = NA
      )
    )
  }
}

# Create empty matrix
trait_assoc <- matrix(
  NA,
  nrow = ncol(cor_traits),
  ncol = ncol(cor_traits),
  dimnames = list(names(cor_traits), names(cor_traits))
)

# Calculate all pairwise associations
for (i in seq_along(cor_traits)) {
  for (j in seq_along(cor_traits)) {
    
    if (i == j) {
      trait_assoc[i, j] <- 1
    } else {
      trait_assoc[i, j] <- assoc(
        cor_traits[[i]],
        cor_traits[[j]]
      )
    }
  }
}

# View matrix
round(trait_assoc, 2)

# Load median results for spatial results
load(file = "C:/Users/Jack Cooper/Documents/iucn/iucn_sim/iucn_data/CritE_median_D24.RData")
# Load synonyms for matching up names in loop
Synonyms <- read_xlsx("C:/Users/Jack Cooper/Dropbox/Jack's PhD/Chapter 3. Future shark FD/Analyses/DatasetSynonyms.xlsx")

# Produce a dataframe to determine which species are extinct in the future
dat.ext <- read_delim("C:/Users/Jack Cooper/Documents/iucn/iucn_sim/iucn_data/future_simulations_future_status_D24_CritE/te_all_species.txt")

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

# Save corrected synonym data to be used in null analyses
save(dat_ext,file = "C:/Users/Jack Cooper/Documents/iucn/iucn_sim/iucn_data/Extinction_times.RData")

# Set up functional matrices for calculations
# Species-trait matrix
traits <- data %>% 
  select(Species,habitat,vertical,terrestriality,thermo,feeding,
         diet,max.length,deepwater,body.plan,reproduction)

sharks_traits <- traits %>% 
  remove_rownames %>% 
  column_to_rownames(var="Species")

# Form trait category matrix
sharks_traits_cat <- tibble(trait_name = c("habitat","vertical","terrestriality","thermo","feeding",
                                           "diet","max.length","deepwater","body.plan","reproduction"),
                            trait_type = c("N", "N", "O", "O", "O", "N", "Q", "N", "N", "N"),
                            trait_weight = c(0.5,1,1,1,0.5,0.5,1,0.5,1,1))

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

pcoa <- list(li = sp_faxes_coord_sharks)
ax <- c("PC1", "PC2", "PC3", "PC4")
coord_d <- pcoa$li[, ax]

# Parallelised loop
registerDoParallel(cores = 8)
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
  
  baskets <- Occ_Future %>%
    select(Species,Present,Future) %>% 
    remove_rownames %>% 
    column_to_rownames(var="Species")
  
  # Transpose to form final occurrence matrix
  baskets_sharks_weights <- t(baskets)
  
  # make assemblage data into a matrix & ensure all matrix values are numeric
  baskets_sharks_weights <- data.matrix(baskets_sharks_weights, rownames.force = NA)
  class(baskets_sharks_weights) <- "numeric"
  
  
  # Transform dataframe from presence-absence matrix in a list
  test <- apply(baskets_sharks_weights,1,function(x){colnames(baskets_sharks_weights)[which(x==1)]})

  lala1<-lapply(test, function(y){ get_FV_Sp (ax=c(1:4), pcoa=pcoa, Selected_sp=y)})
  
  # Calculate all FD metrics per scenario
  FD.metrics <- data.frame(
    Present = rbind(lala1$Present$data$RS,
                    lala1$Present$data$Fric_r),
    Future = rbind(lala1$Future$data$RS,
                   lala1$Future$data$Fric_r)
  )
  
  rownames(FD.metrics) <- c("sp_richn","fric")
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
FDmetrics_taxonvar.CritE<- res_df %>% 
  select(Scenario:fric)
save(FDmetrics_taxonvar.CritE,file = "C:/Users/Jack Cooper/Documents/iucn/iucn_sim/Crit E analyses/CritE_full_buffer.RData")

FDmetrics_long_TaxonVar.CritE<- melt(FDmetrics_taxonvar.CritE, id.vars= "Scenario")
save(FDmetrics_long_TaxonVar.CritE,file = "C:/Users/Jack Cooper/Documents/iucn/iucn_sim/Crit E analyses/CritE_full_long_buffer.RData")
