###############################################################################################################
# 10. AquaMaps RCP 2.6 analyses
## This R code analyses global functional diversity under the RCP 2.6 scenario
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

## RCP 2.6 without extinctions
RCP26 <- read_csv("C:/Users/2022207/Dropbox/Jack's PhD/Chapter 3. Future shark FD/Data/Aquamaps data/Elasmo_rcp26_2100.csv") 

RCP26_filtered <- RCP26 %>% 
  filter (probability >= 0.3)

# Form species-trait matrix from AquaMaps data 
RCP26_filtered <- RCP26_filtered %>% 
  mutate(sci_name = paste(Genus, Species, sep = " "))

RCP26_Occ <- RCP26_filtered %>% 
  select(sci_name,speciesid,SpecCode,csquarecode,centerlat,centerlong,probability) %>% 
  rename(Species = sci_name)

# Form species-trait matrix & remove species in data not seen in AquaMaps
traits <- data %>% 
  select(Species,habitat,vertical,terrestriality,thermo,feeding,diet,max.length)

traits_filtered <- traits %>%
  filter(Species %in% RCP26_Occ$Species)

sharks_traits_2100 <- traits_filtered %>% 
  remove_rownames %>% 
  column_to_rownames(var="Species")

## Calculate abundance per grid cell
RCP26_Occ_abund <- RCP26_filtered %>%
  select(sci_name, speciesid, SpecCode, csquarecode, centerlat, centerlong, probability) %>%
  rename(Species = sci_name) %>% 
  group_by(csquarecode) %>%
  summarize(total_species = n_distinct(Species))

# Put together
RCP26_Occ_filtered <- inner_join(RCP26_Occ, RCP26_Occ_abund, by = "csquarecode")

# Form occurrence matrix with grid cells as assemblages
## Keep coordinates until right before transposition so map can be formed later based on them 
Occ_start <- RCP26_Occ_filtered %>%
  group_by(Species) %>% 
  distinct(speciesid, SpecCode, csquarecode, centerlat, centerlong, probability)
Occ_start <- Occ_start[order(Occ_start$Species),]
Occ_start$Taxa_dup<-paste(Occ_start$Species, Occ_start$csquarecode, sep="+")

# Cast into wide-format data
wp.abun <- Occ_start %>% 
  dcast(Species ~ csquarecode, value.var = "Taxa_dup")

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
  select(csquarecode, centerlat, centerlong) %>% 
  group_by(csquarecode)

Occ_map <- Occ_map %>%
  rename(Grid = csquarecode)

Occ_map <- Occ_map %>%
  distinct() %>%
  ungroup()

# Join occurrence and map; save map to plot species richness
RCP26_Occ_map <- inner_join(Occ_map, baskets_Spp_richn, by = "Grid")
save(RCP26_Occ_map, file = "C:/Users/2022207/Dropbox/Jack's PhD/Chapter 3. Future shark FD/Data/R data/FD_RCP26_Species_richness_map.RData")

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
baskets_sharks_weightsRCP26 <- data.matrix(baskets_sharks_weights, rownames.force = NA)
class(baskets_sharks_weightsRCP26) <- "numeric"

# Remove any species absent from all assemblages (if none, will have no effect)
absent_species <- colnames(baskets_sharks_weightsRCP26)[colSums(baskets_sharks_weightsRCP26) == 0]
baskets_sharks_weightsRCP26 <- baskets_sharks_weightsRCP26[, colSums(baskets_sharks_weightsRCP26) > 0]
sharks_traits_2100 <- sharks_traits_2100[row.names(sharks_traits_2100) %in% colnames(baskets_sharks_weightsRCP26), ]

# summarise dataset
asb_sp_sharks_summ <- mFD::asb.sp.summary(asb_sp_w = baskets_sharks_weightsRCP26)

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
sp_faxes_coord_sharksRCP26 <- fspaces_quality_sharks$"details_fspaces"$"sp_pc_coord"

# Save coordinates and occurrence matrix for supercomputer
save(sp_faxes_coord_sharksRCP26, file = "C:/Users/2022207/Dropbox/Jack's PhD/Chapter 3. Future shark FD/Data/R data/RCP26_coordinates.RData")
save(baskets_sharks_weightsRCP26, file = "C:/Users/2022207/Dropbox/Jack's PhD/Chapter 3. Future shark FD/Data/R data/RCP26_occurrences.RData")

# Calculate FD metrics (used in supercomputer)
alpha_fd_indices_sharksRCP26 <- mFD::alpha.fd.multidim(
  sp_faxes_coord   = sp_faxes_coord_sharksRCP26[ , c("PC1", "PC2", "PC3", "PC4")],
  asb_sp_w         = baskets_sharks_weightsRCP26,
  ind_vect         = c("fric", "fori", "fspe"),
  scaling          = TRUE,
  check_input      = TRUE,
  details_returned = TRUE)

# Load supercomputer results
load(file = "C:/Users/2022207/Dropbox/Jack's PhD/Chapter 3. Future shark FD/Data/R data/RCP26_FDmetrics.RData")

## Form dataframe
FDindices_taxon_RCP26 <- as.data.frame(fd_ind_values_sharks_RCP26) %>% 
  tibble::rownames_to_column("Grid")

# Calculate functional uniqueness per grid cell
## Per species
Funct_Uni_RCP26 <- uniqueness(baskets_sharks_weightsRCP26, as.matrix(sp_dist_sharks)) %>% 
  as.data.frame()
## Per grid cell
baskets_updated <- baskets %>%
  mutate(across(-Species, ~ifelse(. == "1", Funct_Uni_RCP26[["Ui"]], .))) 

baskets_FUn <- baskets_updated %>% 
  remove_rownames %>% 
  column_to_rownames(var="Species")

FUn_baskets <- t(baskets_FUn) 

FUn_baskets<- data.matrix(FUn_baskets, rownames.force = NA)
class(FUn_baskets) <- "numeric"
FUn_baskets <- as.data.frame(FUn_baskets)
FUn_baskets$FUn <- rowMeans(FUn_baskets, na.rm = TRUE)

RCP26_FUn <- FUn_baskets %>% 
  select(FUn)

# Check FUn grids match main dataframe order & add to FD dataframe
RCP26_FUn <- RCP26_FUn %>% 
  rownames_to_column("Grid")

FDindices_taxon_RCP26_FD <- inner_join(FDindices_taxon_RCP26, RCP26_FUn, by = "Grid")

# Save dataframe
save(FDindices_taxon_RCP26_FD, file = "C:/Users/2022207/Dropbox/Jack's PhD/Chapter 3. Future shark FD/Data/R data/FD_RCP26.RData")

# Add map coordinates
## Isolate cell and coordinates to be included in map results
Occ_map <- Occ_start %>% 
  ungroup() %>%  
  select(csquarecode, centerlat, centerlong) %>% 
  group_by(csquarecode)

## Remove grids not present in FD dataframe; then check order is the same & add coordinates
Occ_map <- Occ_map %>%
  rename(Grid = csquarecode)

filtered_Occ_map <- inner_join(Occ_map, FDindices_taxon_RCP26_FD, by = "Grid")

final_map_df <- filtered_Occ_map %>%
  distinct(Grid, .keep_all = TRUE)

final_map_df_RCP26 <- final_map_df %>% as.data.frame()
# Save final dataframe to be used in plotting
save(final_map_df_RCP26, file = "C:/Users/2022207/Dropbox/Jack's PhD/Chapter 3. Future shark FD/Data/R data/FD_RCP26_map_df.RData")


## RCP2.6 with extinctions
# Now re-do with sharks deemed to be extinct by 2100
sharks_traits_filtered <- traits_filtered %>%
  semi_join(Spp_mode %>% filter(mode >= 77), by = "Species") 

sharks_traits_2100 <- sharks_traits_filtered %>% 
  remove_rownames %>% 
  column_to_rownames(var="Species")

## Calculate abundance per grid cell
RCP26_Occ_abund <- RCP26_filtered %>%
  select(sci_name, speciesid, SpecCode, csquarecode, centerlat, centerlong, probability) %>%
  rename(Species = sci_name) %>% 
  group_by(csquarecode) %>%
  summarize(total_species = n_distinct(Species))

# Put together
RCP26_Occ_filtered <- inner_join(RCP26_Occ, RCP26_Occ_abund, by = "csquarecode")

# Form occurrence matrix with grid cells as assemblages
Occ_start <- RCP26_Occ_filtered %>%
  group_by(Species) %>% 
  distinct(speciesid, SpecCode, csquarecode, centerlat, centerlong, probability)
Occ_start <- Occ_start[order(Occ_start$Species),]
Occ_start$Taxa_dup<-paste(Occ_start$Species, Occ_start$csquarecode, sep="+")

# Cast into wide-format data
wp.abun <- Occ_start %>% 
  dcast(Species ~ csquarecode, value.var = "Taxa_dup")

# Ensure species match between trait matrix and occurrence matrix
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
  select(csquarecode, centerlat, centerlong) %>% 
  group_by(csquarecode)

Occ_map <- Occ_map %>%
  rename(Grid = csquarecode)

Occ_map <- Occ_map %>%
  distinct() %>%
  ungroup()

RCP26_ext_Occ_map <- inner_join(Occ_map, baskets_Spp_richn, by = "Grid")
save(RCP26_ext_Occ_map, file = "C:/Users/2022207/Dropbox/Jack's PhD/Chapter 3. Future shark FD/Data/R data/FD_RCP26_extinctions_Species_richness_map.RData")

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
baskets_sharks_weightsRCP26_ext <- data.matrix(baskets_sharks_weights, rownames.force = NA)
class(baskets_sharks_weightsRCP26_ext) <- "numeric"

# Remove any species absent from all assemblages (if none, will have no effect)
absent_species <- colnames(baskets_sharks_weightsRCP26_ext)[colSums(baskets_sharks_weightsRCP26_ext) == 0]
baskets_sharks_weightsRCP26_ext <- baskets_sharks_weightsRCP26_ext[, colSums(baskets_sharks_weightsRCP26_ext) > 0]
sharks_traits_2100 <- sharks_traits_2100[row.names(sharks_traits_2100) %in% colnames(baskets_sharks_weightsRCP26_ext), ]

# summarise dataset
asb_sp_sharks_summ <- mFD::asb.sp.summary(asb_sp_w = baskets_sharks_weightsRCP26_ext)

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
sp_faxes_coord_sharksRCP26_ext <- fspaces_quality_sharks$"details_fspaces"$"sp_pc_coord"

# Save coordinates and occurrence matrix for supercomputer
save(sp_faxes_coord_sharksRCP26_ext, file = "C:/Users/2022207/Dropbox/Jack's PhD/Chapter 3. Future shark FD/Data/R data/RCP26ext_coordinates.RData")
save(baskets_sharks_weightsRCP26_ext, file = "C:/Users/2022207/Dropbox/Jack's PhD/Chapter 3. Future shark FD/Data/R data/RCP26ext_occurrences.RData")

# Calculate FD metrics (used in supercomputer)
alpha_fd_indices_sharksRCP26_ext <- mFD::alpha.fd.multidim(
  sp_faxes_coord   = sp_faxes_coord_sharksRCP26_ext[ , c("PC1", "PC2", "PC3", "PC4")],
  asb_sp_w         = baskets_sharks_weightsRCP26_ext,
  ind_vect         = c("fric", "fori", "fspe"),
  scaling          = TRUE,
  check_input      = TRUE,
  details_returned = TRUE)

# Load supercomputer results
load(file = "C:/Users/2022207/Dropbox/Jack's PhD/Chapter 3. Future shark FD/Data/R data/RCP26ext_FDmetrics.RData")

## Form dataframe
FDindices_taxon_RCP26_extinctions <- as.data.frame(fd_ind_values_sharks_RCP26_extinctions) %>% 
  tibble::rownames_to_column("Grid")

# Calculate functional uniqueness per grid cell
## Per species
Funct_Uni_RCP26_extinctions <- uniqueness(baskets_sharks_weightsRCP26_ext, as.matrix(sp_dist_sharks)) %>% 
  as.data.frame()
## Per grid cell
baskets_updated <- baskets %>%
  mutate(across(-Species, ~ifelse(. == "1", Funct_Uni_RCP26_extinctions[["Ui"]], .))) 

baskets_FUn <- baskets_updated %>% 
  remove_rownames %>% 
  column_to_rownames(var="Species")

FUn_baskets <- t(baskets_FUn) 

FUn_baskets<- data.matrix(FUn_baskets, rownames.force = NA)
class(FUn_baskets) <- "numeric"
FUn_baskets <- as.data.frame(FUn_baskets)
FUn_baskets$FUn <- rowMeans(FUn_baskets, na.rm = TRUE)

RCP26_extinctions_FUn <- FUn_baskets %>% 
  select(FUn)

# Check FUn grids match main dataframe order & add to FD dataframe
RCP26_extinctions_FUn <- RCP26_extinctions_FUn %>% 
  rownames_to_column("Grid")

FDindices_taxon_RCP26_extinctions_FD <- inner_join(FDindices_taxon_RCP26_extinctions, RCP26_extinctions_FUn, by = "Grid")

# Save dataframe
save(FDindices_taxon_RCP26_extinctions_FD, file = "C:/Users/2022207/Dropbox/Jack's PhD/Chapter 3. Future shark FD/Data/R data/FD_RCP26_extinctions.RData")

# Add map coordinates
## Isolate cell and coordinates to be included in map results
Occ_map <- Occ_start %>% 
  ungroup() %>%  
  select(csquarecode, centerlat, centerlong) %>% 
  group_by(csquarecode)

## Remove grids not present in FD dataframe; then check order is the same & add coordinates
Occ_map <- Occ_map %>%
  rename(Grid = csquarecode)

filtered_Occ_map <- inner_join(Occ_map, FDindices_taxon_RCP26_extinctions_FD, by = "Grid")

final_map_df <- filtered_Occ_map %>%
  distinct(Grid, .keep_all = TRUE)

final_map_df_RCP26_extinctions <- final_map_df %>% as.data.frame()
# Save final dataframe to be used in plotting
save(final_map_df_RCP26_extinctions, file = "C:/Users/2022207/Dropbox/Jack's PhD/Chapter 3. Future shark FD/Data/R data/FD_RCP26_extinctions_map_df.RData")
