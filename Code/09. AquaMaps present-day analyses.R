###############################################################################################################
# 09. AquaMaps present-day analyses
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

# Load data
load(file="C:/Users/2022207/Dropbox/Jack's PhD/Chapter 3. Future shark FD/Data/R data/Data.RData")

# Trait type matrix
sharks_traits_cat <- tibble(trait_name = c("habitat","vertical","terrestriality","thermo","feeding","diet","max.length"),
                            trait_type = c("N", "N", "O", "O", "O", "N", "Q"))

# Load Present-day AquaMaps data
AqMap <- read_csv("C:/Users/2022207/Dropbox/Jack's PhD/Chapter 3. Future shark FD/Data/Aquamaps data/Elasmo_current.csv")

AqMap <- AqMap %>% 
  mutate(sci_name = paste(Genus, Species, sep = " "))

# Load synonyms and correct
Synonyms <- read_xlsx("C:/Users/2022207/Dropbox/Jack's PhD/Chapter 3. Future shark FD/Analyses/DatasetSynonyms.xlsx")

AqMap <- AqMap %>%
  left_join(Synonyms, by = c("sci_name" = "Aquamaps_name")) %>%
  mutate(sci_name = ifelse(!is.na(Accepted_name), Accepted_name, sci_name)) %>%
  select(-c(Accepted_name,Marked_iucnsim,iucnsim_name,Marked_aquamaps)) 

AqMap_filtered <- AqMap %>% 
  filter(Probability >= 0.3)

AqMap_Occ <- AqMap_filtered %>% 
  select(sci_name,SpeciesID,SpecCode,CsquareCode,CenterLat,CenterLong,Probability) %>% 
  rename(Species = sci_name)

# Remove species in data not seen in AquaMaps
traits <- data %>% 
  select(Species,habitat,vertical,terrestriality,thermo,feeding,diet,max.length)

traits_filtered <- traits %>%
  filter(Species %in% AqMap_Occ$Species)

sharks_traits <- traits_filtered %>% 
  remove_rownames %>% 
  column_to_rownames(var="Species")

## Calculate abundance per grid cell
AqMap_Occ_abund <- AqMap_filtered %>%
  select(sci_name, SpeciesID, SpecCode, CsquareCode, CenterLat, CenterLong, Probability) %>%
  rename(Species = sci_name) %>% 
  group_by(CsquareCode) %>%
  summarize(total_species = n_distinct(Species))

# Put together
AqMap_Occ_filtered <- inner_join(AqMap_Occ, AqMap_Occ_abund, by = "CsquareCode")

# Form occurrence matrix with grid cells as assemblages
## Keep coordinates until right before transposition so map can be formed later based on them 
Occ_start <- AqMap_Occ_filtered %>%
  group_by(Species) %>% 
  distinct(SpeciesID, SpecCode, CsquareCode, CenterLat, CenterLong, Probability)
Occ_start <- Occ_start[order(Occ_start$Species),]
Occ_start$Taxa_dup<-paste(Occ_start$Species, Occ_start$CsquareCode, sep="+")

# Cast into wide-format data
wp.abun <- Occ_start %>% 
  dcast(Species ~ CsquareCode, value.var = "Taxa_dup")

# Ensure species match between trait matrix and occurrence matrix
wp.abun <- wp.abun %>%
  filter(Species %in% rownames(sharks_traits))

baskets <- wp.abun %>%
  mutate(across(-Species, ~ifelse(is.na(.) | . == "0", "0", "1")))

# Move species into row names
baskets_filtered <- baskets %>% 
  remove_rownames %>% 
  column_to_rownames(var="Species")

# Calculate species richness, add coordinates & save map
baskets_Spp <- t(baskets_filtered) %>% 
  as.data.frame() %>%
  mutate_all(as.numeric) %>% 
  mutate(Spp_Richn = rowSums(.))

baskets_Spp_richn <- baskets_Spp %>%
  rownames_to_column(var="Grid") %>% 
  select(Grid,Spp_Richn)

Occ_map <- Occ_start %>% 
  ungroup() %>%  
  select(CsquareCode, CenterLat, CenterLong) %>% 
  group_by(CsquareCode)

Occ_map <- Occ_map %>%
  rename(Grid = CsquareCode)

Occ_map <- Occ_map %>%
  distinct() %>%
  ungroup()

# Join occurrence and map; save map to plot species richness
AqMap_Occ_map <- inner_join(Occ_map, baskets_Spp_richn, by = "Grid")
save(AqMap_Occ_map, file = "C:/Users/2022207/Dropbox/Jack's PhD/Chapter 3. Future shark FD/Analyses/Final paper R codes/data/AquaMaps/FD_AqMap_Species_richness_map.RData")

# Remove grids with < 5 species for FD analyses 
baskets_sharks_weights <- baskets_Spp %>% 
  filter(Spp_Richn >= 5) %>% 
  select(-Spp_Richn)

## Functional diversity analyses
# Summarise dataset
sharks_traits_summ <- mFD::sp.tr.summary(
  tr_cat     = sharks_traits_cat,   
  sp_tr      = sharks_traits, 
  stop_if_NA = TRUE)

# make assemblage data into a matrix & ensure all matrix values are numeric
baskets_sharks_weightsAqMap <- data.matrix(baskets_sharks_weights, rownames.force = NA)
class(baskets_sharks_weightsAqMap) <- "numeric"

# Remove any species absent from all assemblages (if none, will have no effect)
absent_species <- colnames(baskets_sharks_weightsAqMap)[colSums(baskets_sharks_weightsAqMap) == 0]
baskets_sharks_weightsAqMap <- baskets_sharks_weightsAqMap[, colSums(baskets_sharks_weightsAqMap) > 0]
sharks_traits <- sharks_traits[row.names(sharks_traits) %in% colnames(baskets_sharks_weightsAqMap), ]

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
sp_faxes_coord_sharksAqMap <- fspaces_quality_sharks$"details_fspaces"$"sp_pc_coord"

# Calculate FD metrics
pcoa <- list(li = sp_faxes_coord_sharksAqMap)
ax <- c("PC1", "PC2", "PC3", "PC4")
coord_d <- pcoa$li[, ax]

# Step 1: Create the list of species present in each grid with progress bar
pb_test <- progress_bar$new(
  format = "  Processing test [:bar] :percent in :elapsed",
  total = nrow(baskets_sharks_weights), clear = FALSE, width = 60
)

test <- apply(baskets_sharks_weightsAqMap, 1, function(x) {
  pb_test$tick()  # Update progress bar for test
  colnames(baskets_sharks_weightsAqMap)[which(x == 1)]
})

grids.summary <- as.data.frame(lengths(test))
min(grids.summary$`lengths(test)`)

# Step 2: Calculate FUn with progress tracking
pb_fun <- progress_bar$new(
  format = "  Calculating FU_region [:bar] :percent in :elapsed",
  total = length(test), clear = FALSE, width = 60
)

lala <- lapply(test, function(y) {
  pb_fun$tick()  # Update progress bar for FU_region calculation
  get_FU_region(Mat_dist = as.matrix(sp_dist_sharks), nb_NN = 5, Selected_sp = y)
})

# Step 3: Calculate FV Sp with progress tracking
pb_fvsp <- progress_bar$new(
  format = "  Calculating FV_Sp [:bar] :percent in :elapsed",
  total = length(test), clear = FALSE, width = 60
)

# Use a for-loop to have more control over the progress bar and tracking
lala1 <- vector("list", length(test))

for (i in seq_along(test)) {
  result <- tryCatch({
    get_FV_Sp(ax = c(1:4), pcoa = pcoa, Selected_sp = test[[i]])
  }, error = function(e) {
    warning(paste("Error in get_FV_Sp for species:", paste(test[[i]], collapse = ",")))
    return(NULL)
  })
  
  lala1[[i]] <- result
  
  pb_fvsp$tick()  # Update progress bar after processing each element
}

# Step 4: Extract FD metrics per cell with progress tracking
pb_metrics <- progress_bar$new(
  format = "  Extracting FD metrics [:bar] :percent in :elapsed",
  total = length(test) * 4, clear = FALSE, width = 60
)

sp_richn <- unlist(lapply(lala1, function(k) {
  pb_metrics$tick()  # Update progress bar for sp_richn extraction
  return(k$data$RS)
}))

fric <- unlist(lapply(lala1, function(k) {
  pb_metrics$tick()  # Update progress bar for fric extraction
  return(k$data$Fric_r)
}))

fsp <- unlist(lapply(lala1, function(k) {
  pb_metrics$tick()  # Update progress bar for fsp extraction
  return(k$specialization[1])
}))

fun <- unlist(lapply(lala, function(k) {
  pb_metrics$tick()  # Update progress bar for fun extraction
  return(k$FU_Region_Mean)
}))

# Step 5: Precaution; removing any missing data so all 5 metrics have the same number of entries
valid_indices <- !is.na(sp_richn) & !is.na(fric) & !is.na(fsp) & !is.na(fun)

sp_richn <- sp_richn[valid_indices]
fric <- fric[valid_indices]
fsp <- fsp[valid_indices]
fun <- fun[valid_indices]
grid_names <- rownames(grids.summary)[valid_indices]

# Step 6: Combine the list of FD metrics into a single data frame
fd_ind_values_sharks_AqMap <- data.frame(
  Grid = grid_names, 
  sp_richn = sp_richn,
  fric = fric,
  fsp = fsp,
  fun = fun,
  row.names = NULL
)

# Add map coordinates
## Isolate cell and coordinates to be included in map results
Occ_map <- Occ_start %>% 
  ungroup() %>%  
  select(CsquareCode, CenterLat, CenterLong) %>% 
  group_by(CsquareCode)

## Remove grids not present in FD dataframe; then check order is the same & add coordinates
Occ_map <- Occ_map %>%
  rename(Grid = CsquareCode)

filtered_Occ_map <- inner_join(Occ_map, fd_ind_values_sharks_AqMap, by = "Grid")

final_map_df <- filtered_Occ_map %>%
  distinct(Grid, .keep_all = TRUE)

final_map_df_AqMap <- final_map_df %>% as.data.frame()
# Save final dataframe to be used in plotting
save(final_map_df_AqMap, file = "C:/Users/2022207/Dropbox/Jack's PhD/Chapter 3. Future shark FD/Analyses/Final paper R codes/data/AquaMaps/FD_AqMap_map_df.RData")
