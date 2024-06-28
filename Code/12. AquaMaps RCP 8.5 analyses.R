###############################################################################################################
# 12. AquaMaps RCP 8.5 analyses
## This R code analyses global functional diversity under the RCP 8.5 scenario
## it produces data necessary to form map figures
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

# Load mode IUCN_sim results
load(file = "C:/Users/2022207/Dropbox/Jack's PhD/Chapter 3. Future shark FD/Data/R data/2100_mode.RData")

## RCP85 without extinctions
RCP85 <- read_csv("C:/Users/2022207/Dropbox/Jack's PhD/Chapter 3. Future shark FD/Data/Aquamaps data/Elasmo_rcp85_2100.csv") 

RCP85_filtered <- RCP85 %>% 
  filter (Probability >= 0.3)

# Isolate species 
RCP85_filtered <- RCP85_filtered %>% 
  mutate(sci_name = paste(Genus, Species, sep = " "))

RCP85_Occ <- RCP85_filtered %>% 
  select(sci_name,SpeciesID,SpecCode,CsquareCode,CenterLat,CenterLong,Probability) %>% 
  rename(Species = sci_name)

# Form species-trait matrix & remove species in data not seen in AquaMaps
traits <- data %>% 
  select(Species,habitat,vertical,terrestriality,thermo,feeding,diet,max.length)

traits_filtered <- traits %>%
  filter(Species %in% RCP85_Occ$Species)

sharks_traits_2100 <- traits_filtered %>% 
  remove_rownames %>% 
  column_to_rownames(var="Species")

## Calculate abundance per grid cell
RCP85_Occ_abund <- RCP85_filtered %>%
  select(sci_name, SpeciesID, SpecCode, CsquareCode, CenterLat, CenterLong, Probability) %>%
  rename(Species = sci_name) %>% 
  group_by(CsquareCode) %>%
  summarize(total_species = n_distinct(Species))

# Put together
RCP85_Occ_filtered <- inner_join(RCP85_Occ, RCP85_Occ_abund, by = "CsquareCode")

# Form occurrence matrix with grid cells as assemblages
## Keep coordinates until right before transposition so map can be formed later based on them 
Occ_start <- RCP85_Occ_filtered %>%
  group_by(Species) %>% 
  distinct(SpeciesID, SpecCode, CsquareCode, CenterLat, CenterLong, Probability)
Occ_start <- Occ_start[order(Occ_start$Species),]
Occ_start$Taxa_dup<-paste(Occ_start$Species, Occ_start$CsquareCode, sep="+")

# Cast into wide-format data
wp.abun <- Occ_start %>% 
  dcast(Species ~ CsquareCode, value.var = "Taxa_dup")

wp.abun <- wp.abun %>%
  filter(Species %in% rownames(sharks_traits_2100))

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

RCP85_Occ_map <- inner_join(Occ_map, baskets_Spp_richn, by = "Grid")
save(RCP85_Occ_map, file = "C:/Users/2022207/Dropbox/Jack's PhD/Chapter 3. Future shark FD/Data/R data/FD_RCP85_Species_richness_map.RData")

# Remove grids with < 5 species for FD analyses 
baskets_sharks_weights <- baskets_Spp %>% 
  filter(Spp_Richn >= 5) %>%
  select(-Spp_Richn)

# Summarise dataset
sharks_traits_summ <- mFD::sp.tr.summary(
  tr_cat     = sharks_traits_cat,   
  sp_tr      = sharks_traits_2100, 
  stop_if_NA = TRUE)

# Reveal number of states for each trait
sharks_traits_summ$tr_summary_list

# make assemblage data into a matrix & ensure all matrix values are numeric
baskets_sharks_weightsRCP85 <- data.matrix(baskets_sharks_weights, rownames.force = NA)
class(baskets_sharks_weightsRCP85) <- "numeric"

# Remove any species absent from all assemblages (if none, will have no effect)
absent_species <- colnames(baskets_sharks_weightsRCP85)[colSums(baskets_sharks_weightsRCP85) == 0]
baskets_sharks_weightsRCP85 <- baskets_sharks_weightsRCP85[, colSums(baskets_sharks_weightsRCP85) > 0]
sharks_traits_2100 <- sharks_traits_2100[row.names(sharks_traits_2100) %in% colnames(baskets_sharks_weightsRCP85), ]

# summarise dataset
asb_sp_sharks_summ <- mFD::asb.sp.summary(asb_sp_w = baskets_sharks_weightsRCP85)

# retrieve species occurrences for all assemblages:
asb_sp_sharks_occ <- asb_sp_sharks_summ$asb_sp_occ

# Construct trait distance matrix using species
sp_dist_sharks <- mFD::funct.dist(
  sp_tr         = sharks_traits_2100,
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
sp_faxes_coord_sharksRCP85 <- fspaces_quality_sharks$"details_fspaces"$"sp_pc_coord"

# Save for supercomputer
save(sp_faxes_coord_sharksRCP85, file = "C:/Users/2022207/Dropbox/Jack's PhD/Chapter 3. Future shark FD/Data/R data/RCP85_coordinates.RData")
save(baskets_sharks_weightsRCP85, file = "C:/Users/2022207/Dropbox/Jack's PhD/Chapter 3. Future shark FD/Data/R data/RCP85_occurrences.RData")

# Calculate FD metrics (used in supercomputer)
alpha_fd_indices_sharksRCP85 <- mFD::alpha.fd.multidim(
  sp_faxes_coord   = sp_faxes_coord_sharksRCP85[ , c("PC1", "PC2", "PC3", "PC4")],
  asb_sp_w         = baskets_sharks_weightsRCP85,
  ind_vect         = c("fric", "fori", "fspe"),
  scaling          = TRUE,
  check_input      = TRUE,
  details_returned = TRUE)

# Load supercomputer results
load(file = "C:/Users/2022207/Dropbox/Jack's PhD/Chapter 3. Future shark FD/Data/R data/RCP85_FDmetrics.RData")

## Form dataframe
FDindices_taxon_RCP85 <- as.data.frame(fd_ind_values_sharks_RCP85) %>% 
  tibble::rownames_to_column("Grid")

# Calculate functional uniqueness per grid cell
## Per species
Funct_Uni_RCP85 <- uniqueness(baskets_sharks_weightsRCP85, as.matrix(sp_dist_sharks)) %>% 
  as.data.frame()
## Per grid cell
baskets_updated <- baskets %>%
  mutate(across(-Species, ~ifelse(. == "1", Funct_Uni_RCP85[["Ui"]], .))) 

baskets_FUn <- baskets_updated %>% 
  remove_rownames %>% 
  column_to_rownames(var="Species")

FUn_baskets <- t(baskets_FUn) 

FUn_baskets<- data.matrix(FUn_baskets, rownames.force = NA)
class(FUn_baskets) <- "numeric"
FUn_baskets <- as.data.frame(FUn_baskets)
FUn_baskets$FUn <- rowMeans(FUn_baskets, na.rm = TRUE)

RCP85_FUn <- FUn_baskets %>% 
  select(FUn)

# Check FUn grids match main dataframe order & add to FD dataframe
RCP85_FUn <- RCP85_FUn %>% 
  rownames_to_column("Grid")

FDindices_taxon_RCP85_FD <- inner_join(FDindices_taxon_RCP85, RCP85_FUn, by = "Grid")

# Save dataframe
save(FDindices_taxon_RCP85_FD, file = "C:/Users/2022207/Dropbox/Jack's PhD/Chapter 3. Future shark FD/Data/R data/FD_RCP85.RData")

# Add map coordinates
## Isolate cell and coordinates to be included in map results
Occ_map <- Occ_start %>% 
  ungroup() %>%  
  select(CsquareCode, CenterLat, CenterLong) %>% 
  group_by(CsquareCode)

## Remove grids not present in FD dataframe; then check order is the same & add coordinates
Occ_map <- Occ_map %>%
  rename(Grid = CsquareCode)

filtered_Occ_map <- inner_join(Occ_map, FDindices_taxon_RCP85_FD, by = "Grid")

final_map_df <- filtered_Occ_map %>%
  distinct(Grid, .keep_all = TRUE)

final_map_df_RCP85 <- final_map_df %>% as.data.frame()
# Save final dataframe to be used in plotting
save(final_map_df_RCP85, file = "C:/Users/2022207/Dropbox/Jack's PhD/Chapter 3. Future shark FD/Data/R data/FD_RCP85_map_df.RData")


# Now re-do with sharks deemed to be extinct by 2100
sharks_traits_filtered <- traits_filtered %>%
  semi_join(Spp_mode %>% filter(mode >= 77), by = "Species") 

sharks_traits_2100 <- sharks_traits_filtered %>% 
  remove_rownames %>% 
  column_to_rownames(var="Species")

## Calculate abundance per grid cell
RCP85_Occ_abund <- RCP85_filtered %>%
  select(sci_name, SpeciesID, SpecCode, CsquareCode, CenterLat, CenterLong, Probability) %>%
  rename(Species = sci_name) %>% 
  group_by(CsquareCode) %>%
  summarize(total_species = n_distinct(Species))

# Put together
RCP85_Occ_filtered <- inner_join(RCP85_Occ, RCP85_Occ_abund, by = "CsquareCode")

# Form occurrence matrix with grid cells as assemblages
## Keep coordinates until right before transposition so map can be formed later based on them 
Occ_start <- RCP85_Occ_filtered %>%
  group_by(Species) %>% 
  distinct(SpeciesID, SpecCode, CsquareCode, CenterLat, CenterLong, Probability)
Occ_start <- Occ_start[order(Occ_start$Species),]
Occ_start$Taxa_dup<-paste(Occ_start$Species, Occ_start$CsquareCode, sep="+")

# Cast into wide-format data
wp.abun <- Occ_start %>% 
  dcast(Species ~ CsquareCode, value.var = "Taxa_dup")

# Remove species going extinct in 2100 so that it matches species-trait matrix
wp.abun <- wp.abun %>%
  filter(Species %in% rownames(sharks_traits_2100))

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

RCP85ext_Occ_map <- inner_join(Occ_map, baskets_Spp_richn, by = "Grid")
save(RCP85ext_Occ_map, file = "C:/Users/2022207/Dropbox/Jack's PhD/Chapter 3. Future shark FD/Data/R data/FD_RCP85_extinctions_Species_richness_map.RData")

# Remove grids with < 5 species for FD analyses 
baskets_sharks_weights <- baskets_Spp %>% 
  filter(Spp_Richn >= 5) %>%
  select(-Spp_Richn)

# Summarise dataset
sharks_traits_summ <- mFD::sp.tr.summary(
  tr_cat     = sharks_traits_cat,   
  sp_tr      = sharks_traits_2100, 
  stop_if_NA = TRUE)

# Reveal number of states for each trait
sharks_traits_summ$tr_summary_list

# make assemblage data into a matrix & ensure all matrix values are numeric
baskets_sharks_weightsRCP85_ext <- data.matrix(baskets_sharks_weights, rownames.force = NA)
class(baskets_sharks_weightsRCP85_ext) <- "numeric"

# Remove any species absent from all assemblages (if none, will have no effect)
absent_species <- colnames(baskets_sharks_weightsRCP85_ext)[colSums(baskets_sharks_weightsRCP85_ext) == 0]
baskets_sharks_weightsRCP85_ext <- baskets_sharks_weightsRCP85_ext[, colSums(baskets_sharks_weightsRCP85_ext) > 0]
sharks_traits_2100 <- sharks_traits_2100[row.names(sharks_traits_2100) %in% colnames(baskets_sharks_weightsRCP85_ext), ]

# summarise dataset
asb_sp_sharks_summ <- mFD::asb.sp.summary(asb_sp_w = baskets_sharks_weightsRCP85_ext)

# retrieve species occurrences for all assemblages:
asb_sp_sharks_occ <- asb_sp_sharks_summ$asb_sp_occ

# Construct trait distance matrix using species
sp_dist_sharks <- mFD::funct.dist(
  sp_tr         = sharks_traits_2100,
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
sp_faxes_coord_sharksRCP85_ext <- fspaces_quality_sharks$"details_fspaces"$"sp_pc_coord"

# Save for supercomputer
save(sp_faxes_coord_sharksRCP85_ext, file = "C:/Users/2022207/Dropbox/Jack's PhD/Chapter 3. Future shark FD/Data/R data/RCP85ext_coordinates.RData")
save(baskets_sharks_weightsRCP85_ext, file = "C:/Users/2022207/Dropbox/Jack's PhD/Chapter 3. Future shark FD/Data/R data/RCP85ext_occurrences.RData")

# Calculate FD metrics (used in supercomputer)
alpha_fd_indices_sharksRCP85_ext <- mFD::alpha.fd.multidim(
  sp_faxes_coord   = sp_faxes_coord_sharksRCP85_ext[ , c("PC1", "PC2", "PC3", "PC4")],
  asb_sp_w         = baskets_sharks_weightsRCP85_ext,
  ind_vect         = c("fric", "fori", "fspe"),
  scaling          = TRUE,
  check_input      = TRUE,
  details_returned = TRUE)

# Load supercomputer result
load(file = "C:/Users/2022207/Dropbox/Jack's PhD/Chapter 3. Future shark FD/Data/R data/RCP85ext_FDmetrics.RData")

## Form dataframe
FDindices_taxon_RCP85_extinctions <- as.data.frame(fd_ind_values_sharks_RCP85_extinctions) %>% 
  tibble::rownames_to_column("Grid")

# Calculate functional uniqueness per grid cell
## Per species
Funct_Uni_RCP85_extinctions <- uniqueness(baskets_sharks_weightsRCP85_ext, as.matrix(sp_dist_sharks)) %>% 
  as.data.frame()
## Per grid cell
baskets_updated <- baskets %>%
  mutate(across(-Species, ~ifelse(. == "1", Funct_Uni_RCP85_extinctions[["Ui"]], .))) 

baskets_FUn <- baskets_updated %>% 
  remove_rownames %>% 
  column_to_rownames(var="Species")

FUn_baskets <- t(baskets_FUn) 

FUn_baskets<- data.matrix(FUn_baskets, rownames.force = NA)
class(FUn_baskets) <- "numeric"
FUn_baskets <- as.data.frame(FUn_baskets)
FUn_baskets$FUn <- rowMeans(FUn_baskets, na.rm = TRUE)

RCP85_extinctions_FUn <- FUn_baskets %>% 
  select(FUn)

# Check FUn grids match main dataframe order & add to FD dataframe
RCP85_extinctions_FUn <- RCP85_extinctions_FUn %>% 
  rownames_to_column("Grid")

FDindices_taxon_RCP85_extinctions_FD <- inner_join(FDindices_taxon_RCP85_extinctions, RCP85_extinctions_FUn, by = "Grid")

# Save dataframe
save(FDindices_taxon_RCP85_extinctions_FD, file = "C:/Users/2022207/Dropbox/Jack's PhD/Chapter 3. Future shark FD/Data/R data/FD_RCP85_extinctions.RData")

# Add map coordinates
## Isolate cell and coordinates to be included in map results
Occ_map <- Occ_start %>% 
  ungroup() %>%  
  select(CsquareCode, CenterLat, CenterLong) %>% 
  group_by(CsquareCode)

## Remove grids not present in FD dataframe; then check order is the same & add coordinates
Occ_map <- Occ_map %>%
  rename(Grid = CsquareCode)

filtered_Occ_map <- inner_join(Occ_map, FDindices_taxon_RCP85_extinctions_FD, by = "Grid")

final_map_df <- filtered_Occ_map %>%
  distinct(Grid, .keep_all = TRUE)

final_map_df_RCP85_extinctions <- final_map_df %>% as.data.frame()
# Save final dataframe to be used in plotting
save(final_map_df_RCP85_extinctions, file = "C:/Users/2022207/Dropbox/Jack's PhD/Chapter 3. Future shark FD/Data/R data/FD_RCP85_extinctions_map_df.RData")
