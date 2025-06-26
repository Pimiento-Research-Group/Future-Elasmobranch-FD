###############################################################################################################
# 11. AquaMaps supplementary FD analyses
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

# Load AquaMaps present
load(file="C:/Users/2022207/Dropbox/Jack's PhD/Chapter 3. Future shark FD/Data/R data/Maps/AqMap_Occ.RData")

# Load RCP 2.6 AquaMaps data & correct for synonyms
RCP26 <- read_csv("C:/Users/2022207/Dropbox/Jack's PhD/Chapter 3. Future shark FD/Data/Aquamaps data/Elasmo_rcp26_2100.csv") 

RCP26 <- RCP26 %>% 
  mutate(sci_name = paste(Genus, Species, sep = " "))

RCP26 <- RCP26 %>%
  left_join(Synonyms, by = c("sci_name" = "Aquamaps_name")) %>%
  mutate(Species = ifelse(!is.na(Accepted_name), Accepted_name, sci_name)) %>%
  select(-c(Accepted_name, Marked_iucnsim, iucnsim_name, Marked_aquamaps, sci_name))

RCP26_filtered <- RCP26 %>% 
  filter (probability >= 0.3)

RCP26_filtered <- RCP26_filtered %>%
  filter(Species %in% data$Species)

RCP26_Occ <- RCP26_filtered %>% 
  select(Species,speciesid,SpecCode,csquarecode,centerlat,centerlong,probability) %>% 
  rename(
    SpeciesID = speciesid,
    CsquareCode = csquarecode,
    CenterLat = centerlat,
    CenterLong = centerlong,
    Probability = probability
  )

# Load RCP 8.5 AquaMaps data & correct for synonymns
RCP85 <- read_csv("C:/Users/2022207/Dropbox/Jack's PhD/Chapter 3. Future shark FD/Data/Aquamaps data/Elasmo_rcp85_2100.csv") 

RCP85 <- RCP85 %>% 
  mutate(sci_name = paste(Genus, Species, sep = " "))

RCP85 <- RCP85 %>%
  left_join(Synonyms, by = c("sci_name" = "Aquamaps_name")) %>%
  mutate(Species = ifelse(!is.na(Accepted_name), Accepted_name, sci_name)) %>%
  select(-c(Accepted_name, Marked_iucnsim, iucnsim_name, Marked_aquamaps, sci_name))

RCP85_filtered <- RCP85 %>% 
  filter (Probability >= 0.3)

RCP85_filtered <- RCP85_filtered %>%
  filter(Species %in% data$Species)

RCP85_Occ <- RCP85_filtered %>% 
  select(Species,SpeciesID,SpecCode,CsquareCode,CenterLat,CenterLong,Probability)

## 2. Load AquaMaps grids, species-trait matrix, species-assemblage matrix & PCoA coorindates
load(file="C:/Users/2022207/Dropbox/Jack's PhD/Chapter 3. Future shark FD/Data/R data/Maps/AquaMaps_trait matrix.RData")
load(file="C:/Users/2022207/Dropbox/Jack's PhD/Chapter 3. Future shark FD/Data/R data/Maps/AquaMaps_pcoa_coordinates.RData")
load(file="C:/Users/2022207/Dropbox/Jack's PhD/Chapter 3. Future shark FD/Data/R data/Maps/AquaMaps_species_assemblage.RData")
load(file="C:/Users/2022207/Dropbox/Jack's PhD/Chapter 3. Future shark FD/Data/R data/Maps/Grids.RData")
load(file="C:/Users/2022207/Dropbox/Jack's PhD/Chapter 3. Future shark FD/Data/R data/Maps/AquaMaps_sp_dist_matrix.RData")

# Form final inputs for FD analyses
pcoa <- list(li = sp_faxes_coord_sharksAqMap)
ax <- c("PC1", "PC2", "PC3", "PC4")
coord_d <- pcoa$li[, ax]

## 3. Form RCP 2.6 and RCP 8.5 species assemblage matrices (with & without extinctions)
# Species-grid matrix: Aquamaps climate change (RCP 2.6)
baskets_RCP26 <- RCP26_Occ %>%
  distinct(CsquareCode, Species) %>%
  mutate(presence = 1) %>%
  pivot_wider(names_from = Species, values_from = presence, values_fill = 0) %>%
  rename(Grid = CsquareCode)

# Filter for grids with >5 species for FD analyses
baskets_RCP26_weights <- baskets_RCP26 %>%
  filter(rowSums(select(., -Grid)) >= 5) %>% 
  column_to_rownames("Grid")

baskets_sharks_weightsRCP26 <- data.matrix(baskets_RCP26_weights, rownames.force = NA)
class(baskets_sharks_weightsRCP26) <- "numeric"

# Species-grid matrix: Aquamaps climate change + extinctions (RCP 2.6)
# Remove species & reform matrix
RCP26ext_Occ <- RCP26_Occ %>%
  semi_join(Spp_median.CritE %>% filter(median >= 77), by = "Species")

baskets_RCP26ext <- RCP26ext_Occ %>%
  distinct(CsquareCode, Species) %>%
  mutate(presence = 1) %>%
  pivot_wider(names_from = Species, values_from = presence, values_fill = 0) %>%
  rename(Grid = CsquareCode)

# Filter for grids with >5 species for FD analyses
baskets_RCP26ext_weights <- baskets_RCP26ext %>%
  filter(rowSums(select(., -Grid)) >= 5) %>%
  column_to_rownames("Grid")

baskets_sharks_weightsRCP26ext <- data.matrix(baskets_RCP26ext_weights, rownames.force = NA)
class(baskets_sharks_weightsRCP26ext) <- "numeric"

# Species-grid matrix: Aquamaps climate change (RCP 8.5)
baskets_RCP85 <- RCP85_Occ %>%
  distinct(CsquareCode, Species) %>%
  mutate(presence = 1) %>%
  pivot_wider(names_from = Species, values_from = presence, values_fill = 0) %>%
  rename(Grid = CsquareCode)

# Filter for grids with >5 species for FD analyses
baskets_RCP85_weights <- baskets_RCP85 %>%
  filter(rowSums(select(., -Grid)) >= 5) %>% 
  column_to_rownames("Grid")

baskets_sharks_weightsRCP85 <- data.matrix(baskets_RCP85_weights, rownames.force = NA)
class(baskets_sharks_weightsRCP85) <- "numeric"

# Species-grid matrix: Aquamaps climate change + extinctions (RCP 8.5)
# Remove species & reform matrix
RCP85ext_Occ <- RCP85_Occ %>%
  semi_join(Spp_median.CritE %>% filter(median >= 77), by = "Species")

baskets_RCP85ext <- RCP85ext_Occ %>%
  distinct(CsquareCode, Species) %>%
  mutate(presence = 1) %>%
  pivot_wider(names_from = Species, values_from = presence, values_fill = 0) %>%
  rename(Grid = CsquareCode)

# Filter for grids with >5 species for FD analyses
baskets_RCP85ext_weights <- baskets_RCP85ext %>%
  filter(rowSums(select(., -Grid)) >= 5) %>%
  column_to_rownames("Grid")

baskets_sharks_weightsRCP85ext <- data.matrix(baskets_RCP85ext_weights, rownames.force = NA)
class(baskets_sharks_weightsRCP85ext) <- "numeric"

## 4. Run FD analyses (n = 4)
# Load AquaMaps present-day results - no need to redo these
load(file="C:/Users/2022207/Dropbox/Jack's PhD/Chapter 3. Future shark FD/Data/R data/Maps/FRic_map_present.RData")
load(file="C:/Users/2022207/Dropbox/Jack's PhD/Chapter 3. Future shark FD/Data/R data/Maps/SR_map_present.RData")

# Step 1: Create the list of species present in each grid per scenario with progress bar
pb_test1 <- progress_bar$new(
  format = "  Processing test [:bar] :percent in :elapsed",
  total = nrow(baskets_sharks_weightsRCP26), clear = FALSE, width = 60
)

pb_test2 <- progress_bar$new(
  format = "  Processing test [:bar] :percent in :elapsed",
  total = nrow(baskets_sharks_weightsRCP26ext), clear = FALSE, width = 60
)

pb_test3 <- progress_bar$new(
  format = "  Processing test [:bar] :percent in :elapsed",
  total = nrow(baskets_sharks_weightsRCP85), clear = FALSE, width = 60
)

pb_test4 <- progress_bar$new(
  format = "  Processing test [:bar] :percent in :elapsed",
  total = nrow(baskets_sharks_weightsRCP85ext), clear = FALSE, width = 60
)

test1 <- apply(baskets_sharks_weightsRCP26, 1, function(x) {
  pb_test1$tick()  # Update progress bar for test
  colnames(baskets_sharks_weightsRCP26)[which(x == 1)]
})

test2 <- apply(baskets_sharks_weightsRCP26ext, 1, function(x) {
  pb_test2$tick()  # Update progress bar for test
  colnames(baskets_sharks_weightsRCP26ext)[which(x == 1)]
})

test3 <- apply(baskets_sharks_weightsRCP85, 1, function(x) {
  pb_test3$tick()  # Update progress bar for test
  colnames(baskets_sharks_weightsRCP85ext)[which(x == 1)]
})

test4 <- apply(baskets_sharks_weightsRCP85ext, 1, function(x) {
  pb_test4$tick()  # Update progress bar for test
  colnames(baskets_sharks_weightsRCP85ext)[which(x == 1)]
})

grids.summary1 <- as.data.frame(lengths(test1))
min(grids.summary1$`lengths(test1)`)
grids.summary2 <- as.data.frame(lengths(test2))
min(grids.summary2$`lengths(test2)`)
grids.summary3 <- as.data.frame(lengths(test3))
min(grids.summary3$`lengths(test3)`)
grids.summary4 <- as.data.frame(lengths(test4))
min(grids.summary4$`lengths(test4)`)

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

pb_fun4 <- progress_bar$new(
  format = "  Calculating FU_region [:bar] :percent in :elapsed",
  total = length(test4), clear = FALSE, width = 60
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

lala4 <- lapply(test4, function(y) {
  if (!pb_fun4$finished) pb_fun4$tick()
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

pb_fvsp4 <- progress_bar$new(
  format = "  Calculating FV_Sp [:bar] :percent in :elapsed",
  total = length(test4), clear = FALSE, width = 60
)

# Use a for-loop to have more control over the progress bar and tracking
lala5 <- vector("list", length(test1))
lala6 <- vector("list", length(test2))
lala7 <- vector("list", length(test3))
lala8 <- vector("list", length(test4))

for (i in seq_along(test1)) {
  result <- tryCatch({
    get_FV_Sp(ax = c(1:4), pcoa = pcoa, Selected_sp = test1[[i]])
  })
  
  lala5[[i]] <- result
  
  pb_fvsp1$tick()  # Update progress bar after processing each element
}

for (i in seq_along(test2)) {
  result2 <- tryCatch({
    get_FV_Sp(ax = c(1:4), pcoa = pcoa, Selected_sp = test2[[i]])
  }, error = function(e) {
    warning(paste("Error in get_FV_Sp for species:", paste(test2[[i]], collapse = ",")))
    return(NULL)
  })
  
  lala6[[i]] <- result2
  
  pb_fvsp2$tick()  # Update progress bar after processing each element
}

for (i in seq_along(test3)) {
  result3 <- tryCatch({
    get_FV_Sp(ax = c(1:4), pcoa = pcoa, Selected_sp = test3[[i]])
  }, error = function(e) {
    warning(paste("Error in get_FV_Sp for species:", paste(test3[[i]], collapse = ",")))
    return(NULL)
  })
  
  lala7[[i]] <- result3
  
  pb_fvsp3$tick()  # Update progress bar after processing each element
}

for (i in seq_along(test4)) {
  result4 <- tryCatch({
    get_FV_Sp(ax = c(1:4), pcoa = pcoa, Selected_sp = test4[[i]])
  }, error = function(e) {
    warning(paste("Error in get_FV_Sp for species:", paste(test4[[i]], collapse = ",")))
    return(NULL)
  })
  
  lala8[[i]] <- result4
  
  pb_fvsp4$tick()  # Update progress bar after processing each element
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

pb_metrics4 <- progress_bar$new(
  format = "  Extracting FD metrics [:bar] :percent in :elapsed",
  total = length(test4) * 4, clear = FALSE, width = 60
)

## Climate change (RCP 2.6)
sp_richn_RCP26 <- unlist(lapply(lala5, function(k) {
  pb_metrics1$tick()  # Update progress bar for sp_richn extraction
  return(k$data$RS)
}))

fric_RCP26 <- unlist(lapply(lala5, function(k) {
  pb_metrics1$tick()  # Update progress bar for fric extraction
  return(k$data$Fric_r)
}))

fsp_RCP26 <- unlist(lapply(lala5, function(k) {
  pb_metrics1$tick()  # Update progress bar for fsp extraction
  return(k$specialization[1])
}))

fun_RCP26 <- unlist(lapply(lala1, function(k) {
  pb_metrics1$tick()  # Update progress bar for fun extraction
  return(k$FU_Region_Mean)
}))

## Climate change (RCP 2.6) + extinctions
sp_richn_RCP26ext <- unlist(lapply(lala6, function(k) {
  pb_metrics2$tick()  # Update progress bar for sp_richn extraction
  return(k$data$RS)
}))

fric_RCP26ext <- unlist(lapply(lala6, function(k) {
  pb_metrics2$tick()  # Update progress bar for fric extraction
  return(k$data$Fric_r)
}))

fsp_RCP26ext <- unlist(lapply(lala6, function(k) {
  pb_metrics2$tick()  # Update progress bar for fsp extraction
  return(k$specialization[1])
}))

fun_RCP26ext <- unlist(lapply(lala2, function(k) {
  pb_metrics2$tick()  # Update progress bar for fun extraction
  return(k$FU_Region_Mean)
}))

## Climate change (RCP 8.5)
sp_richn_RCP85 <- unlist(lapply(lala7, function(k) {
  pb_metrics3$tick()  # Update progress bar for sp_richn extraction
  return(k$data$RS)
}))

fric_RCP85 <- unlist(lapply(lala7, function(k) {
  pb_metrics3$tick()  # Update progress bar for fric extraction
  return(k$data$Fric_r)
}))

fsp_RCP85 <- unlist(lapply(lala7, function(k) {
  pb_metrics3$tick()  # Update progress bar for fsp extraction
  return(k$specialization[1])
}))

fun_RCP85 <- unlist(lapply(lala3, function(k) {
  pb_metrics3$tick()  # Update progress bar for fun extraction
  return(k$FU_Region_Mean)
}))

## Climate change (RCP 8.5) + extinctions
sp_richn_RCP85ext <- unlist(lapply(lala8, function(k) {
  pb_metrics4$tick()  # Update progress bar for sp_richn extraction
  return(k$data$RS)
}))

fric_RCP85ext <- unlist(lapply(lala8, function(k) {
  pb_metrics4$tick()  # Update progress bar for fric extraction
  return(k$data$Fric_r)
}))

fsp_RCP85ext <- unlist(lapply(lala8, function(k) {
  pb_metrics4$tick()  # Update progress bar for fsp extraction
  return(k$specialization[1])
}))

fun_RCP85ext <- unlist(lapply(lala4, function(k) {
  pb_metrics4$tick()  # Update progress bar for fun extraction
  return(k$FU_Region_Mean)
}))

# Step 5: Precaution; removing any missing data so all 5 metrics have the same number of entries
valid_indices1 <- !is.na(sp_richn_RCP26) & !is.na(fric_RCP26) & !is.na(fsp_RCP26) & !is.na(fun_RCP26)
valid_indices2 <- !is.na(sp_richn_RCP26ext) & !is.na(fric_RCP26ext) & !is.na(fsp_RCP26ext) & !is.na(fun_RCP26ext)
valid_indices3 <- !is.na(sp_richn_RCP85) & !is.na(fric_RCP85) & !is.na(fsp_RCP85) & !is.na(fun_RCP85)
valid_indices4 <- !is.na(sp_richn_RCP85ext) & !is.na(fric_RCP85ext) & !is.na(fsp_RCP85ext) & !is.na(fun_RCP85ext)

sp_richn_RCP26 <- sp_richn_RCP26[valid_indices1]
fric_RCP26 <- fric_RCP26[valid_indices1]
fsp_RCP26 <- fsp_RCP26[valid_indices1]
fun_RCP26 <- fun_RCP26[valid_indices1]
grid_names1 <- rownames(grids.summary1)[valid_indices1]

sp_richn_RCP26ext <- sp_richn_RCP26ext[valid_indices2]
fric_RCP26ext <- fric_RCP26ext[valid_indices2]
fsp_RCP26ext <- fsp_RCP26ext[valid_indices2]
fun_RCP26ext <- fun_RCP26ext[valid_indices2]
grid_names2 <- rownames(grids.summary2)[valid_indices2]

sp_richn_RCP85 <- sp_richn_RCP85[valid_indices3]
fric_RCP85 <- fric_RCP85[valid_indices3]
fsp_RCP85 <- fsp_RCP85[valid_indices3]
fun_RCP85 <- fun_RCP85[valid_indices3]
grid_names3 <- rownames(grids.summary3)[valid_indices3]

sp_richn_RCP85ext <- sp_richn_RCP85ext[valid_indices4]
fric_RCP85ext <- fric_RCP85ext[valid_indices4]
fsp_RCP85ext <- fsp_RCP85ext[valid_indices4]
fun_RCP85ext <- fun_RCP85ext[valid_indices4]
grid_names4 <- rownames(grids.summary4)[valid_indices4]

# Step 6: Combine the list of FD metrics into a single data frame for each scenario
FD_RCP26 <- data.frame(
  Grid = grid_names1, 
  sp_richn = sp_richn_RCP26,
  fric = fric_RCP26,
  fsp = fsp_RCP26,
  fun = fun_RCP26,
  row.names = NULL
)

FD_RCP26ext <- data.frame(
  Grid = grid_names2, 
  sp_richn = sp_richn_RCP26ext,
  fric = fric_RCP26ext,
  fsp = fsp_RCP26ext,
  fun = fun_RCP26ext,
  row.names = NULL
)

FD_RCP85 <- data.frame(
  Grid = grid_names3, 
  sp_richn = sp_richn_RCP85,
  fric = fric_RCP85,
  fsp = fsp_RCP85,
  fun = fun_RCP85,
  row.names = NULL
)

FD_RCP85ext <- data.frame(
  Grid = grid_names4, 
  sp_richn = sp_richn_RCP85ext,
  fric = fric_RCP85ext,
  fsp = fsp_RCP85ext,
  fun = fun_RCP85ext,
  row.names = NULL
)

# Join spatial coordinates to each FD dataframe
FD_RCP26 <- FD_RCP26 %>%
  left_join(grids, by = "Grid")

FD_RCP26ext <- FD_RCP26ext %>%
  left_join(grids, by = "Grid")

FD_RCP85 <- FD_RCP85 %>%
  left_join(grids, by = "Grid")

FD_RCP85ext <- FD_RCP85ext %>%
  left_join(grids, by = "Grid")

# Save all dataframes for shift calculations and mapping
save(FD_RCP26,file="C:/Users/2022207/Dropbox/Jack's PhD/Chapter 3. Future shark FD/Data/R data/Maps/FRic_map_RCP26.RData")
save(FD_RCP26ext,file="C:/Users/2022207/Dropbox/Jack's PhD/Chapter 3. Future shark FD/Data/R data/Maps/FRic_map_RCP26ext.RData")
save(FD_RCP85,file="C:/Users/2022207/Dropbox/Jack's PhD/Chapter 3. Future shark FD/Data/R data/Maps/FRic_map_RCP85.RData")
save(FD_RCP85ext,file="C:/Users/2022207/Dropbox/Jack's PhD/Chapter 3. Future shark FD/Data/R data/Maps/FRic_map_RCP85ext.RData")

# Create individual species richness maps and save
SR_RCP26 <- RCP26_Occ %>%
  select(CsquareCode, CenterLat, CenterLong, Species) %>%
  rename(Grid = CsquareCode)

SR_RCP26_map <- SR_RCP26 %>%
  distinct(Grid, Species, .keep_all = TRUE) %>% 
  group_by(Grid, CenterLat, CenterLong) %>%
  summarise(species_richness = n(), .groups = "drop")

save(SR_RCP26_map,file="C:/Users/2022207/Dropbox/Jack's PhD/Chapter 3. Future shark FD/Data/R data/Maps/SR_map_RCP26.RData")

SR_RCP26ext <- RCP26ext_Occ %>%
  select(CsquareCode, CenterLat, CenterLong, Species) %>%
  rename(Grid = CsquareCode)

SR_RCP26ext_map <- SR_RCP26ext %>%
  distinct(Grid, Species, .keep_all = TRUE) %>% 
  group_by(Grid, CenterLat, CenterLong) %>%
  summarise(species_richness = n(), .groups = "drop")

save(SR_RCP26ext_map,file="C:/Users/2022207/Dropbox/Jack's PhD/Chapter 3. Future shark FD/Data/R data/Maps/SR_map_RCP26ext.RData")

SR_RCP85 <- RCP85_Occ %>%
  select(CsquareCode, CenterLat, CenterLong, Species) %>%
  rename(Grid = CsquareCode)

SR_RCP85_map <- SR_RCP85 %>%
  distinct(Grid, Species, .keep_all = TRUE) %>% 
  group_by(Grid, CenterLat, CenterLong) %>%
  summarise(species_richness = n(), .groups = "drop")

save(SR_RCP85_map,file="C:/Users/2022207/Dropbox/Jack's PhD/Chapter 3. Future shark FD/Data/R data/Maps/SR_map_RCP85.RData")

SR_RCP85ext <- RCP85ext_Occ %>%
  select(CsquareCode, CenterLat, CenterLong, Species) %>%
  rename(Grid = CsquareCode)

SR_RCP85ext_map <- SR_RCP85ext %>%
  distinct(Grid, Species, .keep_all = TRUE) %>% 
  group_by(Grid, CenterLat, CenterLong) %>%
  summarise(species_richness = n(), .groups = "drop")

save(SR_RCP85ext_map,file="C:/Users/2022207/Dropbox/Jack's PhD/Chapter 3. Future shark FD/Data/R data/Maps/SR_map_RCP85ext.RData")