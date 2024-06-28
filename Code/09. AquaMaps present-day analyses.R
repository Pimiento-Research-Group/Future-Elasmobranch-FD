###############################################################################################################
# 09. AquaMaps present-day analyses
## This R code analyses global functional diversity in the present
## it produces data necessary to form all map figures
###############################################################################################################

# Import packages
library(tidyverse)
library(tibble)
library(scales)
library(cowplot)
library(doBy)
library(reshape2)
library(mFD)
library(data.table)
library(funrar)

# Load data
load(file="C:/Users/2022207/Dropbox/Jack's PhD/Chapter 3. Future shark FD/Data/R data/Data.RData")

# Trait type matrix
sharks_traits_cat <- tibble(trait_name = c("habitat","vertical","terrestriality","thermo","feeding","diet","max.length"),
                            trait_type = c("N", "N", "O", "O", "O", "N", "Q"))

# Load Present-day AquaMaps data
AqMap <- read_csv("C:/Users/2022207/Dropbox/Jack's PhD/Chapter 3. Future shark FD/Data/Aquamaps data/Elasmo_current.csv")

AqMap_filtered <- AqMap %>% 
  filter(Probability >= 0.3)

# Isolate species  
AqMap_filtered <- AqMap_filtered %>% 
  mutate(sci_name = paste(Genus, Species, sep = " "))

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
save(AqMap_Occ_map, file = "C:/Users/2022207/Dropbox/Jack's PhD/Chapter 3. Future shark FD/Data/R data/FD_AqMap_Species_richness_map.RData")

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

# Reveal number of states for each trait
sharks_traits_summ$tr_summary_list

# make assemblage data into a matrix & ensure all matrix values are numeric
baskets_sharks_weightsAqMap <- data.matrix(baskets_sharks_weights, rownames.force = NA)
class(baskets_sharks_weightsAqMap) <- "numeric"

# Remove any species absent from all assemblages (if none, will have no effect)
absent_species <- colnames(baskets_sharks_weightsAqMap)[colSums(baskets_sharks_weightsAqMap) == 0]
baskets_sharks_weightsAqMap <- baskets_sharks_weightsAqMap[, colSums(baskets_sharks_weightsAqMap) > 0]
sharks_traits <- sharks_traits[row.names(sharks_traits) %in% colnames(baskets_sharks_weightsAqMap), ]

# summarise dataset
asb_sp_sharks_summ <- mFD::asb.sp.summary(asb_sp_w = baskets_sharks_weightsAqMap)

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
sp_faxes_coord_sharksAqMap <- fspaces_quality_sharks$"details_fspaces"$"sp_pc_coord"

# Save data for supercomputer
save(sp_faxes_coord_sharksAqMap, file = "C:/Users/2022207/Dropbox/Jack's PhD/Chapter 3. Future shark FD/Data/R data/AqMapext_coordinates.RData")
save(baskets_sharks_weightsAqMap, file = "C:/Users/2022207/Dropbox/Jack's PhD/Chapter 3. Future shark FD/Data/R data/AqMapext_occurrences.RData")

# FD metric calculations (run in supercomputer)
alpha_fd_indices_sharksAqMap <- mFD::alpha.fd.multidim(
  sp_faxes_coord   = sp_faxes_coord_sharksAqMap[ , c("PC1", "PC2", "PC3", "PC4")],
  asb_sp_w         = baskets_sharks_weightsAqMap,
  ind_vect         = c("fric", "fori", "fspe"),
  scaling          = TRUE,
  check_input      = TRUE,
  details_returned = TRUE)

fd_ind_values_sharks_AqMap <- alpha_fd_indices_sharksAqMap$"functional_diversity_indices"

# Load results of supercomputer analysis
load(file = "C:/Users/2022207/Dropbox/Jack's PhD/Chapter 3. Future shark FD/Data/R data/AqMapext_FDmetrics.RData")

# Form dataframe
FDindices_taxon_AqMap <- as.data.frame(fd_ind_values_sharks_AqMap) %>% 
  tibble::rownames_to_column("Grid")

# Calculate functional uniqueness per grid cell
## Per species
Funct_Uni_AqMap <- uniqueness(baskets_sharks_weightsAqMap, as.matrix(sp_dist_sharks)) %>% 
  as.data.frame()
## Per grid cell
baskets_updated <- baskets %>%
  mutate(across(-Species, ~ifelse(. == "1", Funct_Uni_AqMap[["Ui"]], .))) 

baskets_FUn <- baskets_updated %>% 
  remove_rownames %>% 
  column_to_rownames(var="Species")

FUn_baskets <- t(baskets_FUn) 

FUn_baskets<- data.matrix(FUn_baskets, rownames.force = NA)
class(FUn_baskets) <- "numeric"
FUn_baskets <- as.data.frame(FUn_baskets)
FUn_baskets$FUn <- rowMeans(FUn_baskets, na.rm = TRUE)

AqMap_FUn <- FUn_baskets %>% 
  select(FUn)

# Check FUn grids match main dataframe order & add to FD dataframe
AqMap_FUn <- AqMap_FUn %>% 
  rownames_to_column("Grid")

FDindices_taxon_AqMap_FD <- inner_join(FDindices_taxon_AqMap, AqMap_FUn, by = "Grid")

# Save dataframe
save(FDindices_taxon_AqMap_FD, file = "C:/Users/2022207/Dropbox/Jack's PhD/Chapter 3. Future shark FD/Data/R data/FD_AqMap.RData")

# Add map coordinates
## Isolate cell and coordinates to be included in map results
Occ_map <- Occ_start %>% 
  ungroup() %>%  
  select(CsquareCode, CenterLat, CenterLong) %>% 
  group_by(CsquareCode)

## Remove grids not present in FD dataframe; then check order is the same & add coordinates
Occ_map <- Occ_map %>%
  rename(Grid = CsquareCode)

filtered_Occ_map <- inner_join(Occ_map, FDindices_taxon_AqMap_FD, by = "Grid")

final_map_df <- filtered_Occ_map %>%
  distinct(Grid, .keep_all = TRUE)

final_map_df_AqMap <- final_map_df %>% as.data.frame()
# Save final dataframe to be used in plotting
save(final_map_df_AqMap, file = "C:/Users/2022207/Dropbox/Jack's PhD/Chapter 3. Future shark FD/Data/R data/FD_AqMap_map_df.RData")
