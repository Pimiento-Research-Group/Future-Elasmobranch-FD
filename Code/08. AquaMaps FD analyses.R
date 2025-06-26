###############################################################################################################
# 08. AquaMaps FD analyses
## This R code analyses global functional diversity in the present
## it produces data necessary to form all map figures
###############################################################################################################

# Import packages
library(readxl)
library(tidyverse)
library(tibble)
library(scales)
library(cowplot)
library(doBy)
library(reshape2)
library(mFD)
library(data.table)
library(progress)
library(geometry)

# Source functions
source("C:/Users/2022207/Dropbox/Jack's PhD/Chapter 3. Future shark FD/Analyses/R code/Final pipeline R code/Functions/get_indicator_function 2.R")
source("C:/Users/2022207/Dropbox/Jack's PhD/Chapter 3. Future shark FD/Analyses/R code/Final pipeline R code/Functions/fonction_FRIC_Global_full.R")

## 1. Load all data
# Load data
load(file="C:/Users/2022207/Dropbox/Jack's PhD/Chapter 3. Future shark FD/Data/R data/Data.RData")

# Trait type matrix
sharks_traits_cat <- tibble(trait_name = c("habitat","vertical","terrestriality","thermo","feeding","diet","max.length"),
                            trait_type = c("N", "N", "O", "O", "O", "N", "Q"))

# Load median IUCN_sim results - Crit E
load(file = "C:/Users/2022207/Dropbox/Jack's PhD/Chapter 3. Future shark FD/Analyses/Final paper R codes/data/Crit_E/CritE_median.RData")

# Load synonyms and correct
Synonyms <- read_xlsx("C:/Users/2022207/Dropbox/Jack's PhD/Chapter 3. Future shark FD/Analyses/DatasetSynonyms.xlsx")

Spp_median.CritE <- Spp_median.CritE %>%
  left_join(Synonyms, by = c("Species" = "iucnsim_name")) %>%
  mutate(Species = ifelse(!is.na(Accepted_name), Accepted_name, Species)) %>%
  select(-c(Accepted_name,Marked_iucnsim,Aquamaps_name,Marked_aquamaps)) %>% 
  distinct(Species, .keep_all = TRUE)

# Load Present-day AquaMaps data & correct for synonymns
AqMap <- read_csv("C:/Users/2022207/Dropbox/Jack's PhD/Chapter 3. Future shark FD/Data/Aquamaps data/Elasmo_current.csv")

AqMap <- AqMap %>% 
  mutate(sci_name = paste(Genus, Species, sep = " "))

AqMap <- AqMap %>%
  left_join(Synonyms, by = c("sci_name" = "Aquamaps_name")) %>%
  mutate(Species = ifelse(!is.na(Accepted_name), Accepted_name, sci_name)) %>%  
  select(-c(Accepted_name, Marked_iucnsim, iucnsim_name, Marked_aquamaps, sci_name)) 

AqMap_filtered <- AqMap %>% 
  filter(Probability >= 0.3)

AqMap_filtered <- AqMap_filtered %>%
  filter(Species %in% data$Species)

AqMap_Occ <- AqMap_filtered %>% 
  select(Species,SpeciesID,SpecCode,CsquareCode,CenterLat,CenterLong,Probability)

save(AqMap_Occ,file="C:/Users/2022207/Dropbox/Jack's PhD/Chapter 3. Future shark FD/Data/R data/Maps/AqMap_Occ.RData")

# Load RCP 4.5 AquaMaps data & correct for synonymns
RCP45 <- read_csv("C:/Users/2022207/Dropbox/Jack's PhD/Chapter 3. Future shark FD/Data/Aquamaps data/Elasmo_rcp45_2100.csv") 

RCP45 <- RCP45 %>% 
  mutate(sci_name = paste(Genus, Species, sep = " "))

RCP45 <- RCP45 %>%
  left_join(Synonyms, by = c("sci_name" = "Aquamaps_name")) %>%
  mutate(Species = ifelse(!is.na(Accepted_name), Accepted_name, sci_name)) %>%
  select(-c(Accepted_name, Marked_iucnsim, iucnsim_name, Marked_aquamaps, sci_name))
 
RCP45_filtered <- RCP45 %>% 
  filter (Probability >= 0.3)

RCP45_filtered <- RCP45_filtered %>%
  filter(Species %in% data$Species)

RCP45_Occ <- RCP45_filtered %>% 
  select(Species,SpeciesID,SpecCode,CsquareCode,CenterLat,CenterLong,Probability)

## 2. Form species-trait matrix
# Form species-trait matrix using species present in AquaMaps
traits <- data %>% 
  select(Species,habitat,vertical,terrestriality,thermo,feeding,diet,max.length)

traits_filtered <- traits %>%
  filter(Species %in% AqMap_Occ$Species)

sharks_traits <- traits_filtered %>% 
  remove_rownames %>% 
  column_to_rownames(var="Species")
save(sharks_traits,file="C:/Users/2022207/Dropbox/Jack's PhD/Chapter 3. Future shark FD/Data/R data/Maps/AquaMaps_trait matrix.RData")

## 3. Run PCoA analyses
# Summarise dataset
sharks_traits_summ <- mFD::sp.tr.summary(
  tr_cat     = sharks_traits_cat,   
  sp_tr      = sharks_traits, 
  stop_if_NA = TRUE)

# Construct trait distance matrix using species
sp_dist_sharks <- mFD::funct.dist(
  sp_tr         = sharks_traits,
  tr_cat        = sharks_traits_cat,
  metric        = "gower",
  scale_euclid  = "scale_center",
  ordinal_var   = "classic",
  weight_type   = "equal",
  stop_if_NA    = TRUE)
save(sp_dist_sharks,file="C:/Users/2022207/Dropbox/Jack's PhD/Chapter 3. Future shark FD/Data/R data/Maps/AquaMaps_sp_dist_matrix.RData")

# Assess quality
fspaces_quality_sharks <- mFD::quality.fspaces(
  sp_dist             = sp_dist_sharks,
  maxdim_pcoa         = 10,
  deviation_weighting = c("absolute","squared"),
  fdist_scaling       = FALSE,
  fdendro             = "average")

# Return coordinates of each axis
sp_faxes_coord_sharksAqMap <- fspaces_quality_sharks$"details_fspaces"$"sp_pc_coord"
save(sp_faxes_coord_sharksAqMap,file="C:/Users/2022207/Dropbox/Jack's PhD/Chapter 3. Future shark FD/Data/R data/Maps/AquaMaps_pcoa_coordinates.RData")

# Form final inputs for FD analyses
pcoa <- list(li = sp_faxes_coord_sharksAqMap)
ax <- c("PC1", "PC2", "PC3", "PC4")
coord_d <- pcoa$li[, ax]

## 4. Form species-grid matrices (present, RCP 4.5 & climate + extinctions)
# Form grids to be saved and added to matrices for maps
grids <- AqMap_Occ %>%
  distinct(CsquareCode, CenterLat, CenterLong) %>% 
  rename(Grid = CsquareCode)
save(grids,file="C:/Users/2022207/Dropbox/Jack's PhD/Chapter 3. Future shark FD/Data/R data/Maps/Grids.RData")

# Species-grid matrix: Aquamaps present day
baskets_AqMap <- AqMap_Occ %>%
  distinct(CsquareCode, Species) %>%
  mutate(presence = 1) %>%
  pivot_wider(names_from = Species, values_from = presence, values_fill = 0) %>%
  rename(Grid = CsquareCode)

# Filter for grids with >5 species for FD analyses
baskets_AqMaps_weights <- baskets_AqMap %>%
  filter(rowSums(select(., -Grid)) >= 5) %>% 
  column_to_rownames("Grid")

baskets_sharks_weightsAqMap <- data.matrix(baskets_AqMaps_weights, rownames.force = NA)
class(baskets_sharks_weightsAqMap) <- "numeric"
save(baskets_sharks_weightsAqMap,file="C:/Users/2022207/Dropbox/Jack's PhD/Chapter 3. Future shark FD/Data/R data/Maps/AquaMaps_species_assemblage.RData")

# Species-grid matrix: Aquamaps climate change (RCP 4.5)
baskets_RCP45 <- RCP45_Occ %>%
  distinct(CsquareCode, Species) %>%
  mutate(presence = 1) %>%
  pivot_wider(names_from = Species, values_from = presence, values_fill = 0) %>%
  rename(Grid = CsquareCode)

# Filter for grids with >5 species for FD analyses
baskets_RCP45_weights <- baskets_RCP45 %>%
  filter(rowSums(select(., -Grid)) >= 5) %>% 
  column_to_rownames("Grid")

baskets_sharks_weightsRCP45 <- data.matrix(baskets_RCP45_weights, rownames.force = NA)
class(baskets_sharks_weightsRCP45) <- "numeric"

# Species-grid matrix: Aquamaps climate change + extinctions (RCP 4.5)
# Remove species & reform matrix
RCP45ext_Occ <- RCP45_Occ %>%
  semi_join(Spp_median.CritE %>% filter(median >= 77), by = "Species")

baskets_RCP45ext <- RCP45ext_Occ %>%
  distinct(CsquareCode, Species) %>%
  mutate(presence = 1) %>%
  pivot_wider(names_from = Species, values_from = presence, values_fill = 0) %>%
  rename(Grid = CsquareCode)

# Filter for grids with >5 species for FD analyses
baskets_RCP45ext_weights <- baskets_RCP45ext %>%
  filter(rowSums(select(., -Grid)) >= 5) %>%
  column_to_rownames("Grid")

baskets_sharks_weightsRCP45ext <- data.matrix(baskets_RCP45ext_weights, rownames.force = NA)
class(baskets_sharks_weightsRCP45ext) <- "numeric"

## 5. Run FD analyses (x3; 1 = Aquamaps; 2 = RCP 4.5; 3 = RCP 4.5 + extinctions)
# Step 1: Create the list of species present in each grid per scenario with progress bar
pb_test1 <- progress_bar$new(
  format = "  Processing test [:bar] :percent in :elapsed",
  total = nrow(baskets_sharks_weightsAqMap), clear = FALSE, width = 60
)

pb_test2 <- progress_bar$new(
  format = "  Processing test [:bar] :percent in :elapsed",
  total = nrow(baskets_sharks_weightsRCP45), clear = FALSE, width = 60
)

pb_test3 <- progress_bar$new(
  format = "  Processing test [:bar] :percent in :elapsed",
  total = nrow(baskets_sharks_weightsRCP45ext), clear = FALSE, width = 60
)

test1 <- apply(baskets_sharks_weightsAqMap, 1, function(x) {
  pb_test1$tick()  # Update progress bar for test
  colnames(baskets_sharks_weightsAqMap)[which(x == 1)]
})

test2 <- apply(baskets_sharks_weightsRCP45, 1, function(x) {
  pb_test2$tick()  # Update progress bar for test
  colnames(baskets_sharks_weightsRCP45)[which(x == 1)]
})

test3 <- apply(baskets_sharks_weightsRCP45ext, 1, function(x) {
  pb_test3$tick()  # Update progress bar for test
  colnames(baskets_sharks_weightsRCP45ext)[which(x == 1)]
})

grids.summary1 <- as.data.frame(lengths(test1))
min(grids.summary1$`lengths(test1)`)
grids.summary2 <- as.data.frame(lengths(test2))
min(grids.summary2$`lengths(test2)`)
grids.summary3 <- as.data.frame(lengths(test3))
min(grids.summary3$`lengths(test3)`)

# Step 2: Calculate FUn with progress tracking
pb_fun1 <- progress_bar$new(
  format = "  Calculating FU_region [:bar] :percent in :elapsed",
  total = length(test1), clear = FALSE, width = 60
)

pb_fun2 <- progress_bar$new(
  format = "  Calculating FU_region [:bar] :percent in :elapsed",
  total = length(test2), clear = FALSE, width = 60
)

pb_fun3 <- progress_bar$new(
  format = "  Calculating FU_region [:bar] :percent in :elapsed",
  total = length(test3), clear = FALSE, width = 60
)

sp_dist_sharks_mat <- as.matrix(sp_dist_sharks)

lala1 <- lapply(test1, function(y) {
  if (!pb_fun1$finished) pb_fun1$tick()
  if (length(y) >= 2) {
    get_FU_region(Mat_dist = sp_dist_sharks_mat, nb_NN = 5, Selected_sp = y)
  } else {
    NULL
  }
})

lala2 <- lapply(test2, function(y) {
  if (!pb_fun2$finished) pb_fun2$tick()
  selected <- intersect(y, rownames(sp_dist_sharks_mat))
  if (length(selected) >= 2) {
    get_FU_region(Mat_dist = sp_dist_sharks_mat, nb_NN = 5, Selected_sp = selected)
  } else {
    NULL
  }
})

lala3 <- lapply(test3, function(y) {
  if (!pb_fun3$finished) pb_fun3$tick()
  selected <- intersect(y, rownames(sp_dist_sharks_mat))
  if (length(selected) >= 2) {
    get_FU_region(Mat_dist = sp_dist_sharks_mat, nb_NN = 5, Selected_sp = selected)
  } else {
    NULL
  }
})

# Step 3: Calculate FV Sp with progress tracking
pb_fvsp1 <- progress_bar$new(
  format = "  Calculating FV_Sp [:bar] :percent in :elapsed",
  total = length(test1), clear = FALSE, width = 60
)

pb_fvsp2 <- progress_bar$new(
  format = "  Calculating FV_Sp [:bar] :percent in :elapsed",
  total = length(test2), clear = FALSE, width = 60
)

pb_fvsp3 <- progress_bar$new(
  format = "  Calculating FV_Sp [:bar] :percent in :elapsed",
  total = length(test3), clear = FALSE, width = 60
)

# Use a for-loop to have more control over the progress bar and tracking
lala4 <- vector("list", length(test1))
lala5 <- vector("list", length(test2))
lala6 <- vector("list", length(test3))

for (i in seq_along(test1)) {
  result <- tryCatch({
    get_FV_Sp(ax = c(1:4), pcoa = pcoa, Selected_sp = test1[[i]])
  })
  
  lala4[[i]] <- result
  
  pb_fvsp1$tick()  # Update progress bar after processing each element
}

for (i in seq_along(test2)) {
  result2 <- tryCatch({
    get_FV_Sp(ax = c(1:4), pcoa = pcoa, Selected_sp = test2[[i]])
  }, error = function(e) {
    warning(paste("Error in get_FV_Sp for species:", paste(test2[[i]], collapse = ",")))
    return(NULL)
  })
  
  lala5[[i]] <- result2
  
  pb_fvsp2$tick()  # Update progress bar after processing each element
}

for (i in seq_along(test3)) {
  result3 <- tryCatch({
    get_FV_Sp(ax = c(1:4), pcoa = pcoa, Selected_sp = test3[[i]])
  }, error = function(e) {
    warning(paste("Error in get_FV_Sp for species:", paste(test3[[i]], collapse = ",")))
    return(NULL)
  })
  
  lala6[[i]] <- result3
  
  pb_fvsp3$tick()  # Update progress bar after processing each element
}

# Step 4: Extract FD metrics per cell with progress tracking
pb_metrics1 <- progress_bar$new(
  format = "  Extracting FD metrics [:bar] :percent in :elapsed",
  total = length(test1) * 4, clear = FALSE, width = 60
)

pb_metrics2 <- progress_bar$new(
  format = "  Extracting FD metrics [:bar] :percent in :elapsed",
  total = length(test2) * 4, clear = FALSE, width = 60
)

pb_metrics3 <- progress_bar$new(
  format = "  Extracting FD metrics [:bar] :percent in :elapsed",
  total = length(test3) * 4, clear = FALSE, width = 60
)

## Present-day
sp_richn_AqMap <- unlist(lapply(lala4, function(k) {
  pb_metrics1$tick()  # Update progress bar for sp_richn extraction
  return(k$data$RS)
}))

fric_AqMap <- unlist(lapply(lala4, function(k) {
  pb_metrics1$tick()  # Update progress bar for fric extraction
  return(k$data$Fric_r)
}))

fsp_AqMap <- unlist(lapply(lala4, function(k) {
  pb_metrics1$tick()  # Update progress bar for fsp extraction
  return(k$specialization[1])
}))

fun_AqMap <- unlist(lapply(lala1, function(k) {
  pb_metrics1$tick()  # Update progress bar for fun extraction
  return(k$FU_Region_Mean)
}))

## Climate change (RCP 4.5)
sp_richn_RCP45 <- unlist(lapply(lala5, function(k) {
  pb_metrics2$tick()  # Update progress bar for sp_richn extraction
  return(k$data$RS)
}))

fric_RCP45 <- unlist(lapply(lala5, function(k) {
  pb_metrics2$tick()  # Update progress bar for fric extraction
  return(k$data$Fric_r)
}))

fsp_RCP45 <- unlist(lapply(lala5, function(k) {
  pb_metrics2$tick()  # Update progress bar for fsp extraction
  return(k$specialization[1])
}))

fun_RCP45 <- unlist(lapply(lala2, function(k) {
  pb_metrics2$tick()  # Update progress bar for fun extraction
  return(k$FU_Region_Mean)
}))

## Climate change (RCP 4.5) + extinctions
sp_richn_RCP45ext <- unlist(lapply(lala6, function(k) {
  pb_metrics3$tick()  # Update progress bar for sp_richn extraction
  return(k$data$RS)
}))

fric_RCP45ext <- unlist(lapply(lala6, function(k) {
  pb_metrics3$tick()  # Update progress bar for fric extraction
  return(k$data$Fric_r)
}))

fsp_RCP45ext <- unlist(lapply(lala6, function(k) {
  pb_metrics3$tick()  # Update progress bar for fsp extraction
  return(k$specialization[1])
}))

fun_RCP45ext <- unlist(lapply(lala3, function(k) {
  pb_metrics3$tick()  # Update progress bar for fun extraction
  return(k$FU_Region_Mean)
}))

# Step 5: Precaution; removing any missing data so all 5 metrics have the same number of entries
valid_indices1 <- !is.na(sp_richn_AqMap) & !is.na(fric_AqMap) & !is.na(fsp_AqMap) & !is.na(fun_AqMap)
valid_indices2 <- !is.na(sp_richn_RCP45) & !is.na(fric_RCP45) & !is.na(fsp_RCP45) & !is.na(fun_RCP45)
valid_indices3 <- !is.na(sp_richn_RCP45ext) & !is.na(fric_RCP45ext) & !is.na(fsp_RCP45ext) & !is.na(fun_RCP45ext)

sp_richn_AqMap <- sp_richn_AqMap[valid_indices1]
fric_AqMap <- fric_AqMap[valid_indices1]
fsp_AqMap <- fsp_AqMap[valid_indices1]
fun_AqMap <- fun_AqMap[valid_indices1]
grid_names1 <- rownames(grids.summary1)[valid_indices1]

sp_richn_RCP45 <- sp_richn_RCP45[valid_indices2]
fric_RCP45 <- fric_RCP45[valid_indices2]
fsp_RCP45 <- fsp_RCP45[valid_indices2]
fun_RCP45 <- fun_RCP45[valid_indices2]
grid_names2 <- rownames(grids.summary2)[valid_indices2]

sp_richn_RCP45ext <- sp_richn_RCP45ext[valid_indices3]
fric_RCP45ext <- fric_RCP45ext[valid_indices3]
fsp_RCP45ext <- fsp_RCP45ext[valid_indices3]
fun_RCP45ext <- fun_RCP45ext[valid_indices3]
grid_names3 <- rownames(grids.summary3)[valid_indices3]

# Step 6: Combine the list of FD metrics into a single data frame for each scenario
FD_AqMap <- data.frame(
  Grid = grid_names1, 
  sp_richn = sp_richn_AqMap,
  fric = fric_AqMap,
  fsp = fsp_AqMap,
  fun = fun_AqMap,
  row.names = NULL
)

FD_RCP45 <- data.frame(
  Grid = grid_names2, 
  sp_richn = sp_richn_RCP45,
  fric = fric_RCP45,
  fsp = fsp_RCP45,
  fun = fun_RCP45,
  row.names = NULL
)

FD_RCP45ext <- data.frame(
  Grid = grid_names3, 
  sp_richn = sp_richn_RCP45ext,
  fric = fric_RCP45ext,
  fsp = fsp_RCP45ext,
  fun = fun_RCP45ext,
  row.names = NULL
)

# Join spatial coordinates to each FD dataframe
FD_AqMap <- FD_AqMap %>%
  left_join(grids, by = "Grid")

FD_RCP45 <- FD_RCP45 %>%
  left_join(grids, by = "Grid")

FD_RCP45ext <- FD_RCP45ext %>%
  left_join(grids, by = "Grid")

# Save all 3 dataframes for shift calculations and mapping
save(FD_AqMap,file="C:/Users/2022207/Dropbox/Jack's PhD/Chapter 3. Future shark FD/Data/R data/Maps/FRic_map_present.RData")
save(FD_RCP45,file="C:/Users/2022207/Dropbox/Jack's PhD/Chapter 3. Future shark FD/Data/R data/Maps/FRic_map_RCP45.RData")
save(FD_RCP45ext,file="C:/Users/2022207/Dropbox/Jack's PhD/Chapter 3. Future shark FD/Data/R data/Maps/FRic_map_RCP45ext.RData")

# Create individual species richness maps and save
SR_AqMap <- AqMap_Occ %>%
  select(CsquareCode, CenterLat, CenterLong, Species) %>%
  rename(Grid = CsquareCode)

SR_AqMap_map <- SR_AqMap %>%
  distinct(Grid, Species, .keep_all = TRUE) %>% 
  group_by(Grid, CenterLat, CenterLong) %>%
  summarise(species_richness = n(), .groups = "drop")

save(SR_AqMap_map,file="C:/Users/2022207/Dropbox/Jack's PhD/Chapter 3. Future shark FD/Data/R data/Maps/SR_map_present.RData")

SR_RCP45 <- RCP45_Occ %>%
  select(CsquareCode, CenterLat, CenterLong, Species) %>%
  rename(Grid = CsquareCode)

SR_RCP45_map <- SR_RCP45 %>%
  distinct(Grid, Species, .keep_all = TRUE) %>% 
  group_by(Grid, CenterLat, CenterLong) %>%
  summarise(species_richness = n(), .groups = "drop")

save(SR_RCP45_map,file="C:/Users/2022207/Dropbox/Jack's PhD/Chapter 3. Future shark FD/Data/R data/Maps/SR_map_RCP45.RData")

SR_RCP45ext <- RCP45ext_Occ %>%
  select(CsquareCode, CenterLat, CenterLong, Species) %>%
  rename(Grid = CsquareCode)

SR_RCP45ext_map <- SR_RCP45ext %>%
  distinct(Grid, Species, .keep_all = TRUE) %>% 
  group_by(Grid, CenterLat, CenterLong) %>%
  summarise(species_richness = n(), .groups = "drop")

save(SR_RCP45ext_map,file="C:/Users/2022207/Dropbox/Jack's PhD/Chapter 3. Future shark FD/Data/R data/Maps/SR_map_RCP45ext.RData")