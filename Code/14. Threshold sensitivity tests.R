###############################################################################################################
# 14. Spatial sensitivity tests to threshold
## This R code tests the spatial results to sensitivity against thresholds
## it produces Tables S8 and S9
###############################################################################################################

# 1. Packages
library(readxl)
library(tidyverse)
library(tibble)
library(mFD)
library(progress)
library(geometry)

# Source functions
source("C:/Users/2022207/Dropbox/Jack's PhD/Chapter 3. Future shark FD/Analyses/R code/Final pipeline R code/Functions/get_indicator_function 2.R")
source("C:/Users/2022207/Dropbox/Jack's PhD/Chapter 3. Future shark FD/Analyses/R code/Final pipeline R code/Functions/fonction_FRIC_Global_full.R")

# 2. Load primary analysis (RCP 4.5; threshold 0.3)
load(file="C:/Users/Jack Cooper/Documents/iucn/iucn_sim/Maps/FRic_map_present.RData")
load(file="C:/Users/Jack Cooper/Documents/iucn/iucn_sim/Maps/FRic_map_RCP45.RData")
load(file="C:/Users/Jack Cooper/Documents/iucn/iucn_sim/Maps/FRic_map_RCP45ext.RData")
load(file="C:/Users/Jack Cooper/Documents/iucn/iucn_sim/Maps/SR_map_present.RData")
load(file="C:/Users/Jack Cooper/Documents/iucn/iucn_sim/Maps/SR_map_RCP45.RData")
load(file="C:/Users/Jack Cooper/Documents/iucn/iucn_sim/Maps/SR_map_RCP45ext.RData")

# 3. Load data required for threshold sensitivity analyses
## Load main species/trait dataset
load(file="C:/Users/Jack Cooper/Documents/iucn/Data.RData")

# Load median extinction times
load(file="C:/Users/Jack Cooper/Documents/iucn/iucn_sim/iucn_data/CritE_median_D24.RData")

# Load synonyms
Synonyms <- read_xlsx(
  "C:/Users/Jack Cooper/Dropbox/Jack's PhD/Chapter 3. Future shark FD/Analyses/DatasetSynonyms.xlsx"
)

# Trait type matrix
sharks_traits_cat <- tibble(trait_name = c("habitat","vertical","terrestriality","thermo","feeding",
                                           "diet","max.length","deepwater","body.plan","reproduction"),
                            trait_type = c("N", "N", "O", "O", "O", "N", "Q", "N", "N", "N"),
                            trait_weight = c(0.5,1,1,1,0.5,0.5,1,0.5,1,1))

# 4. Primary analysis — threshold 0.3
## Calculate shifts under RCP 4.5 (climate change)
FD_shift_RCP45 <- FD_AqMap %>%
  inner_join(
    FD_RCP45,
    by = "Grid",
    suffix = c("_present", "_rcp45")
  ) %>%
  mutate(
    SR_shift = sp_richn_rcp45 - sp_richn_present,
    FRic_shift = fric_rcp45 - fric_present
  ) %>%
  select(
    Grid,
    CenterLat_present,
    CenterLong_present,
    SR_shift,
    FRic_shift
  ) %>%
  rename(
    CenterLat = CenterLat_present,
    CenterLong = CenterLong_present
  ) %>% 
  na.omit()

# Calculate shifts under RCP 4.5 + extinctions
FD_shift_RCP45ext <- FD_AqMap %>%
  inner_join(
    FD_RCP45ext,
    by = "Grid",
    suffix = c("_present", "_rcp45ext")
  ) %>%
  mutate(
    SR_shift = sp_richn_rcp45ext - sp_richn_present,
    FRic_shift = fric_rcp45ext - fric_present
  ) %>%
  select(
    Grid,
    CenterLat_present,
    CenterLong_present,
    SR_shift,
    FRic_shift
  ) %>%
  rename(
    CenterLat = CenterLat_present,
    CenterLong = CenterLong_present
  ) %>% 
  na.omit()

# Calculate grid changes under both scenarios
sr_shift_RCP45 <- FD_shift_RCP45 %>%
  mutate(
    SR_change = case_when(
      SR_shift > 0 ~ "Increase",
      SR_shift < 0 ~ "Decrease",
      TRUE ~ "No change"
    )
  ) %>%
  count(SR_change) %>%
  mutate(Proportion = n / sum(n))

fric_shift_RCP45 <- FD_shift_RCP45 %>%
  mutate(
    FRic_change = case_when(
      FRic_shift > 0 ~ "Increase",
      FRic_shift < 0 ~ "Decrease",
      TRUE ~ "No change"
    )
  ) %>%
  count(FRic_change) %>%
  mutate(Proportion = n / sum(n))

sr_shift_RCP45ext <- FD_shift_RCP45ext %>%
  mutate(
    SR_change = case_when(
      SR_shift > 0 ~ "Increase",
      SR_shift < 0 ~ "Decrease",
      TRUE ~ "No change"
    )
  ) %>%
  count(SR_change) %>%
  mutate(Proportion = n / sum(n))

fric_shift_RCP45ext <- FD_shift_RCP45ext %>%
  mutate(
    FRic_change = case_when(
      FRic_shift > 0 ~ "Increase",
      FRic_shift < 0 ~ "Decrease",
      TRUE ~ "No change"
    )
  ) %>%
  count(FRic_change) %>%
  mutate(Proportion = n / sum(n))

# 5. THRESHOLD 0.2
## Load Present-day AquaMaps data & correct for synonyms
AqMap_02 <- read_csv(
  "C:/Users/Jack Cooper/Dropbox/Jack's PhD/Chapter 3. Future shark FD/Data/Aquamaps data/Elasmo_current.csv"
)

AqMap_02 <- AqMap_02 %>% 
  mutate(
    sci_name = paste(Genus, Species, sep = " ")
  )

AqMap_02 <- AqMap_02 %>%
  left_join(
    Synonyms,
    by = c("sci_name" = "Aquamaps_name")
  ) %>%
  mutate(
    Species = ifelse(!is.na(Accepted_name), Accepted_name, sci_name)
  ) %>%  
  select(
    -c(
      Accepted_name,
      Marked_iucnsim,
      iucnsim_name,
      Marked_aquamaps,
      sci_name
    )
  )

# Apply 0.2 occurrence threshold
AqMap_02_filtered <- AqMap_02 %>% 
  filter(Probability >= 0.2)

AqMap_02_filtered <- AqMap_02_filtered %>%
  filter(Species %in% data$Species)

AqMap_02_Occ <- AqMap_02_filtered %>% 
  select(
    Species,
    SpeciesID,
    SpecCode,
    CsquareCode,
    CenterLat,
    CenterLong,
    Probability
  )

# Load RCP 4.5 AquaMaps data & correct for synonyms
RCP45_02 <- read_csv(
  "C:/Users/Jack Cooper/Dropbox/Jack's PhD/Chapter 3. Future shark FD/Data/Aquamaps data/Elasmo_rcp45_2100.csv"
)

RCP45_02 <- RCP45_02 %>% 
  mutate(
    sci_name = paste(Genus, Species, sep = " ")
  )

RCP45_02 <- RCP45_02 %>%
  left_join(
    Synonyms,
    by = c("sci_name" = "Aquamaps_name")
  ) %>%
  mutate(
    Species = ifelse(!is.na(Accepted_name), Accepted_name, sci_name)
  ) %>%  
  select(
    -c(
      Accepted_name,
      Marked_iucnsim,
      iucnsim_name,
      Marked_aquamaps,
      sci_name
    )
  )

# Apply 0.2 occurrence threshold
RCP45_02_filtered <- RCP45_02 %>% 
  filter(Probability >= 0.2)

RCP45_02_filtered <- RCP45_02_filtered %>%
  filter(Species %in% data$Species)

RCP45_02_Occ <- RCP45_02_filtered %>% 
  select(
    Species,
    SpeciesID,
    SpecCode,
    CsquareCode,
    CenterLat,
    CenterLong,
    Probability
  )

# 6. Threshold-specific trait matrix — 0.2
traits_02 <- data %>% 
  select(
    Species,habitat,vertical,terrestriality,thermo,feeding,
    diet,max.length,deepwater,body.plan,reproduction
  )

traits_02_filtered <- traits_02 %>%
  filter(
    Species %in% unique(
      c(
        AqMap_02_Occ$Species,
        RCP45_02_Occ$Species
      )
    )
  )

sharks_traits_02 <- traits_02_filtered %>% 
  remove_rownames %>% 
  column_to_rownames(var = "Species")

# 7. Functional space — 0.2
sp_dist_sharks_02 <- mFD::funct.dist(
  sp_tr = sharks_traits_02,
  tr_cat = sharks_traits_cat,
  metric = "gower",
  scale_euclid = "scale_center",
  ordinal_var = "classic",
  weight_type = "equal",
  stop_if_NA = TRUE
)

fspaces_quality_sharks_02 <- mFD::quality.fspaces(
  sp_dist = sp_dist_sharks_02,
  maxdim_pcoa = 10,
  deviation_weighting = c("absolute", "squared"),
  fdist_scaling = FALSE,
  fdendro = "average"
)

sp_faxes_coord_sharksAqMap_02 <-
  fspaces_quality_sharks_02$details_fspaces$sp_pc_coord

pcoa_02 <- list(
  li = sp_faxes_coord_sharksAqMap_02
)

# 8. Species-grid matrices — 0.2
## Present day
baskets_AqMap_02 <- AqMap_02_Occ %>%
  distinct(CsquareCode, Species) %>%
  mutate(presence = 1) %>%
  pivot_wider(
    names_from = Species,
    values_from = presence,
    values_fill = 0
  ) %>%
  rename(Grid = CsquareCode)

baskets_AqMaps_weights_02 <- baskets_AqMap_02 %>%
  filter(
    rowSums(select(., -Grid)) >= 5
  ) %>% 
  column_to_rownames("Grid")

baskets_sharks_weightsAqMap_02 <- data.matrix(
  baskets_AqMaps_weights_02,
  rownames.force = NA
)

class(baskets_sharks_weightsAqMap_02) <- "numeric"

# Climate change
baskets_RCP45_02 <- RCP45_02_Occ %>%
  distinct(CsquareCode, Species) %>%
  mutate(presence = 1) %>%
  pivot_wider(
    names_from = Species,
    values_from = presence,
    values_fill = 0
  ) %>%
  rename(Grid = CsquareCode)

baskets_RCP45_weights_02 <- baskets_RCP45_02 %>%
  filter(
    rowSums(select(., -Grid)) >= 5
  ) %>% 
  column_to_rownames("Grid")

baskets_sharks_weightsRCP45_02 <- data.matrix(
  baskets_RCP45_weights_02,
  rownames.force = NA
)

class(baskets_sharks_weightsRCP45_02) <- "numeric"

# Climate change + extinctions
RCP45ext_02_Occ <- RCP45_02_Occ %>%
  semi_join(
    Spp_median.CritE %>%
      filter(median >= 77),
    by = "Species"
  )

baskets_RCP45ext_02 <- RCP45ext_02_Occ %>%
  distinct(CsquareCode, Species) %>%
  mutate(presence = 1) %>%
  pivot_wider(
    names_from = Species,
    values_from = presence,
    values_fill = 0
  ) %>%
  rename(Grid = CsquareCode)

baskets_RCP45ext_weights_02 <- baskets_RCP45ext_02 %>%
  filter(
    rowSums(select(., -Grid)) >= 5
  ) %>% 
  column_to_rownames("Grid")

baskets_sharks_weightsRCP45ext_02 <- data.matrix(
  baskets_RCP45ext_weights_02,
  rownames.force = NA
)

class(baskets_sharks_weightsRCP45ext_02) <- "numeric"

# 9. Species lists — 0.2
test1_02 <- apply(
  baskets_sharks_weightsAqMap_02,
  1,
  function(x)
    colnames(baskets_sharks_weightsAqMap_02)[which(x == 1)]
)

test2_02 <- apply(
  baskets_sharks_weightsRCP45_02,
  1,
  function(x)
    colnames(baskets_sharks_weightsRCP45_02)[which(x == 1)]
)

test3_02 <- apply(
  baskets_sharks_weightsRCP45ext_02,
  1,
  function(x)
    colnames(baskets_sharks_weightsRCP45ext_02)[which(x == 1)]
)

# 10. Calculate FV_Sp — 0.2
pb_fvsp1_02 <- progress_bar$new(
  format = "  Calculating FV_Sp present 0.2 [:bar] :percent in :elapsed",
  total = length(test1_02),
  clear = FALSE,
  width = 60
)

pb_fvsp2_02 <- progress_bar$new(
  format = "  Calculating FV_Sp climate 0.2 [:bar] :percent in :elapsed",
  total = length(test2_02),
  clear = FALSE,
  width = 60
)

pb_fvsp3_02 <- progress_bar$new(
  format = "  Calculating FV_Sp climate + extinction 0.2 [:bar] :percent in :elapsed",
  total = length(test3_02),
  clear = FALSE,
  width = 60
)


lala4_02 <- vector("list", length(test1_02))
lala5_02 <- vector("list", length(test2_02))
lala6_02 <- vector("list", length(test3_02))


for (i in seq_along(test1_02)) {
  
  result <- get_FV_Sp(
    ax = c(1:4),
    pcoa = pcoa_02,
    Selected_sp = test1_02[[i]]
  )
  
  lala4_02[[i]] <- result
  
  pb_fvsp1_02$tick()
}


for (i in seq_along(test2_02)) {
  
  result2 <- get_FV_Sp(
    ax = c(1:4),
    pcoa = pcoa_02,
    Selected_sp = test2_02[[i]]
  )
  
  lala5_02[[i]] <- result2
  
  pb_fvsp2_02$tick()
}


for (i in seq_along(test3_02)) {
  
  result3 <- get_FV_Sp(
    ax = c(1:4),
    pcoa = pcoa_02,
    Selected_sp = test3_02[[i]]
  )
  
  lala6_02[[i]] <- result3
  
  pb_fvsp3_02$tick()
}

# 11. Extract FRic — 0.2
fric_AqMap_02 <- sapply(
  seq_along(lala4_02),
  function(i) {
    x <- lala4_02[[i]]$data$Fric_r
    
    if (length(x) == 0) {
      NA_real_
    } else {
      as.numeric(x)
    }
  }
)

names(fric_AqMap_02) <- names(test1_02)

fric_RCP45_02 <- sapply(
  seq_along(lala5_02),
  function(i) {
    x <- lala5_02[[i]]$data$Fric_r
    
    if (length(x) == 0) {
      NA_real_
    } else {
      as.numeric(x)
    }
  }
)

names(fric_RCP45_02) <- names(test2_02)

fric_RCP45ext_02 <- sapply(
  seq_along(lala6_02),
  function(i) {
    x <- lala6_02[[i]]$data$Fric_r
    
    if (length(x) == 0) {
      NA_real_
    } else {
      as.numeric(x)
    }
  }
)

names(fric_RCP45ext_02) <- names(test3_02)

# 12. Create FRic + SR + coordinate dataframes — 0.2
## Present-day grid information
grid_info_AqMap_02 <- AqMap_02_Occ %>%
  distinct(
    CsquareCode,
    CenterLat,
    CenterLong
  )

# RCP45 grid information
grid_info_RCP45_02 <- RCP45_02_Occ %>%
  distinct(
    CsquareCode,
    CenterLat,
    CenterLong
  )

# RCP45 + extinction grid information
grid_info_RCP45ext_02 <- RCP45ext_02_Occ %>%
  distinct(
    CsquareCode,
    CenterLat,
    CenterLong
  )

# Present day
FD_AqMap_02 <- data.frame(
  Grid = names(fric_AqMap_02),
  fric = as.numeric(fric_AqMap_02),
  row.names = NULL
) %>%
  left_join(
    baskets_AqMaps_weights_02 %>%
      rownames_to_column("Grid") %>%
      mutate(
        sp_richn = rowSums(
          select(., -Grid)
        )
      ) %>%
      select(Grid, sp_richn),
    by = "Grid"
  ) %>%
  left_join(
    grid_info_AqMap_02 %>%
      rename(Grid = CsquareCode),
    by = "Grid"
  )

save(FD_AqMap_02,file="C:/Users/Jack Cooper/Documents/iucn/iucn_sim/Maps/AqMap_02.RData")

# RCP45
FD_RCP45_02 <- data.frame(
  Grid = names(fric_RCP45_02),
  fric = as.numeric(fric_RCP45_02),
  row.names = NULL
) %>%
  left_join(
    baskets_RCP45_weights_02 %>%
      rownames_to_column("Grid") %>%
      mutate(
        sp_richn = rowSums(
          select(., -Grid)
        )
      ) %>%
      select(Grid, sp_richn),
    by = "Grid"
  ) %>%
  left_join(
    grid_info_RCP45_02 %>%
      rename(Grid = CsquareCode),
    by = "Grid"
  )

save(FD_RCP45_02,file="C:/Users/Jack Cooper/Documents/iucn/iucn_sim/Maps/RCP45_02.RData")

# RCP45 + extinction
FD_RCP45ext_02 <- data.frame(
  Grid = names(fric_RCP45ext_02),
  fric = as.numeric(fric_RCP45ext_02),
  row.names = NULL
) %>%
  left_join(
    baskets_RCP45ext_weights_02 %>%
      rownames_to_column("Grid") %>%
      mutate(
        sp_richn = rowSums(
          select(., -Grid)
        )
      ) %>%
      select(Grid, sp_richn),
    by = "Grid"
  ) %>%
  left_join(
    grid_info_RCP45ext_02 %>%
      rename(Grid = CsquareCode),
    by = "Grid"
  )

save(FD_RCP45ext_02,file="C:/Users/Jack Cooper/Documents/iucn/iucn_sim/Maps/RCP45ext_02.RData")

# 13. Calculate shifts — 0.2
FD_shift_RCP45_02 <- FD_AqMap_02 %>%
  inner_join(
    FD_RCP45_02,
    by = "Grid",
    suffix = c("_present", "_rcp45")
  ) %>%
  mutate(
    SR_shift = sp_richn_rcp45 - sp_richn_present,
    FRic_shift = fric_rcp45 - fric_present
  ) %>%
  select(
    Grid,
    CenterLat_present,
    CenterLong_present,
    SR_shift,
    FRic_shift
  ) %>%
  rename(
    CenterLat = CenterLat_present,
    CenterLong = CenterLong_present
  ) %>%
  na.omit()

FD_shift_RCP45ext_02 <- FD_AqMap_02 %>%
  inner_join(
    FD_RCP45ext_02,
    by = "Grid",
    suffix = c("_present", "_rcp45ext")
  ) %>%
  mutate(
    SR_shift = sp_richn_rcp45ext - sp_richn_present,
    FRic_shift = fric_rcp45ext - fric_present
  ) %>%
  select(
    Grid,
    CenterLat_present,
    CenterLong_present,
    SR_shift,
    FRic_shift
  ) %>%
  rename(
    CenterLat = CenterLat_present,
    CenterLong = CenterLong_present
  ) %>%
  na.omit()

# 14. Calculate grid changes — 0.2
sr_shift_RCP45_02 <- FD_shift_RCP45_02 %>%
  mutate(
    SR_change = case_when(
      SR_shift > 0 ~ "Increase",
      SR_shift < 0 ~ "Decrease",
      TRUE ~ "No change"
    )
  ) %>%
  count(SR_change) %>%
  mutate(Proportion = n / sum(n))


fric_shift_RCP45_02 <- FD_shift_RCP45_02 %>%
  mutate(
    FRic_change = case_when(
      FRic_shift > 0 ~ "Increase",
      FRic_shift < 0 ~ "Decrease",
      TRUE ~ "No change"
    )
  ) %>%
  count(FRic_change) %>%
  mutate(Proportion = n / sum(n))

sr_shift_RCP45ext_02 <- FD_shift_RCP45ext_02 %>%
  mutate(
    SR_change = case_when(
      SR_shift > 0 ~ "Increase",
      SR_shift < 0 ~ "Decrease",
      TRUE ~ "No change"
    )
  ) %>%
  count(SR_change) %>%
  mutate(Proportion = n / sum(n))

fric_shift_RCP45ext_02 <- FD_shift_RCP45ext_02 %>%
  mutate(
    FRic_change = case_when(
      FRic_shift > 0 ~ "Increase",
      FRic_shift < 0 ~ "Decrease",
      TRUE ~ "No change"
    )
  ) %>%
  count(FRic_change) %>%
  mutate(Proportion = n / sum(n))

# 15. THRESHOLD 0.5
## Load Present-day AquaMaps data & correct for synonyms
AqMap_05 <- read_csv(
  "C:/Users/Jack Cooper/Dropbox/Jack's PhD/Chapter 3. Future shark FD/Data/Aquamaps data/Elasmo_current.csv"
)

AqMap_05 <- AqMap_05 %>% 
  mutate(
    sci_name = paste(Genus, Species, sep = " ")
  )

AqMap_05 <- AqMap_05 %>%
  left_join(
    Synonyms,
    by = c("sci_name" = "Aquamaps_name")
  ) %>%
  mutate(
    Species = ifelse(!is.na(Accepted_name), Accepted_name, sci_name)
  ) %>%  
  select(
    -c(
      Accepted_name,
      Marked_iucnsim,
      iucnsim_name,
      Marked_aquamaps,
      sci_name
    )
  )

# Apply 0.5 occurrence threshold
AqMap_05_filtered <- AqMap_05 %>% 
  filter(Probability >= 0.5)

AqMap_05_filtered <- AqMap_05_filtered %>%
  filter(Species %in% data$Species)

AqMap_05_Occ <- AqMap_05_filtered %>% 
  select(
    Species,
    SpeciesID,
    SpecCode,
    CsquareCode,
    CenterLat,
    CenterLong,
    Probability
  )

# Load RCP 4.5 AquaMaps data & correct for synonyms
RCP45_05 <- read_csv(
  "C:/Users/Jack Cooper/Dropbox/Jack's PhD/Chapter 3. Future shark FD/Data/Aquamaps data/Elasmo_rcp45_2100.csv"
)

RCP45_05 <- RCP45_05 %>% 
  mutate(
    sci_name = paste(Genus, Species, sep = " ")
  )

RCP45_05 <- RCP45_05 %>%
  left_join(
    Synonyms,
    by = c("sci_name" = "Aquamaps_name")
  ) %>%
  mutate(
    Species = ifelse(!is.na(Accepted_name), Accepted_name, sci_name)
  ) %>%  
  select(
    -c(
      Accepted_name,
      Marked_iucnsim,
      iucnsim_name,
      Marked_aquamaps,
      sci_name
    )
  )

# Apply 0.5 occurrence threshold
RCP45_05_filtered <- RCP45_05 %>% 
  filter(Probability >= 0.5)

RCP45_05_filtered <- RCP45_05_filtered %>%
  filter(Species %in% data$Species)

RCP45_05_Occ <- RCP45_05_filtered %>% 
  select(
    Species,
    SpeciesID,
    SpecCode,
    CsquareCode,
    CenterLat,
    CenterLong,
    Probability
  )

# 16. Threshold-specific trait matrix — 0.5
traits_05 <- data %>% 
  select(
    Species,habitat,vertical,terrestriality,thermo,feeding,
    diet,max.length,deepwater,body.plan,reproduction
  )

species_05 <- unique(c(
  AqMap_05_Occ$Species,
  RCP45_05_Occ$Species
))

traits_05_filtered <- traits_05 %>%
  filter(Species %in% species_05)

sharks_traits_05 <- traits_05_filtered %>% 
  remove_rownames() %>% 
  column_to_rownames(var = "Species")

# 17. Functional space — 0.5
sp_dist_sharks_05 <- mFD::funct.dist(
  sp_tr = sharks_traits_05,
  tr_cat = sharks_traits_cat,
  metric = "gower",
  scale_euclid = "scale_center",
  ordinal_var = "classic",
  weight_type = "equal",
  stop_if_NA = TRUE
)

fspaces_quality_sharks_05 <- mFD::quality.fspaces(
  sp_dist = sp_dist_sharks_05,
  maxdim_pcoa = 10,
  deviation_weighting = c("absolute", "squared"),
  fdist_scaling = FALSE,
  fdendro = "average"
)

sp_faxes_coord_sharksAqMap_05 <-
  fspaces_quality_sharks_05$details_fspaces$sp_pc_coord

pcoa_05 <- list(
  li = sp_faxes_coord_sharksAqMap_05
)

# 18. Species-grid matrices — 0.5
## Present day
baskets_AqMap_05 <- AqMap_05_Occ %>%
  distinct(CsquareCode, Species) %>%
  mutate(presence = 1) %>%
  pivot_wider(
    names_from = Species,
    values_from = presence,
    values_fill = 0
  ) %>%
  rename(Grid = CsquareCode)

baskets_AqMaps_weights_05 <- baskets_AqMap_05 %>%
  filter(
    rowSums(select(., -Grid)) >= 5
  ) %>% 
  column_to_rownames("Grid")

baskets_sharks_weightsAqMap_05 <- data.matrix(
  baskets_AqMaps_weights_05,
  rownames.force = NA
)

class(baskets_sharks_weightsAqMap_05) <- "numeric"

# Climate change
baskets_RCP45_05 <- RCP45_05_Occ %>%
  distinct(CsquareCode, Species) %>%
  mutate(presence = 1) %>%
  pivot_wider(
    names_from = Species,
    values_from = presence,
    values_fill = 0
  ) %>%
  rename(Grid = CsquareCode)

baskets_RCP45_weights_05 <- baskets_RCP45_05 %>%
  filter(
    rowSums(select(., -Grid)) >= 5
  ) %>% 
  column_to_rownames("Grid")

baskets_sharks_weightsRCP45_05 <- data.matrix(
  baskets_RCP45_weights_05,
  rownames.force = NA
)

class(baskets_sharks_weightsRCP45_05) <- "numeric"

# Climate change + extinctions
RCP45ext_05_Occ <- RCP45_05_Occ %>%
  semi_join(
    Spp_median.CritE %>%
      filter(median >= 77),
    by = "Species"
  )

baskets_RCP45ext_05 <- RCP45ext_05_Occ %>%
  distinct(CsquareCode, Species) %>%
  mutate(presence = 1) %>%
  pivot_wider(
    names_from = Species,
    values_from = presence,
    values_fill = 0
  ) %>%
  rename(Grid = CsquareCode)

baskets_RCP45ext_weights_05 <- baskets_RCP45ext_05 %>%
  filter(
    rowSums(select(., -Grid)) >= 5
  ) %>% 
  column_to_rownames("Grid")

baskets_sharks_weightsRCP45ext_05 <- data.matrix(
  baskets_RCP45ext_weights_05,
  rownames.force = NA
)

class(baskets_sharks_weightsRCP45ext_05) <- "numeric"

# 19. Species lists — 0.5
test1_05 <- apply(
  baskets_sharks_weightsAqMap_05,
  1,
  function(x)
    colnames(baskets_sharks_weightsAqMap_05)[which(x == 1)]
)

test2_05 <- apply(
  baskets_sharks_weightsRCP45_05,
  1,
  function(x)
    colnames(baskets_sharks_weightsRCP45_05)[which(x == 1)]
)

test3_05 <- apply(
  baskets_sharks_weightsRCP45ext_05,
  1,
  function(x)
    colnames(baskets_sharks_weightsRCP45ext_05)[which(x == 1)]
)

# 20. Calculate FV_Sp — 0.5
pb_fvsp1_05 <- progress_bar$new(
  format = "  Calculating FV_Sp present 0.5 [:bar] :percent in :elapsed",
  total = length(test1_05),
  clear = FALSE,
  width = 60
)

pb_fvsp2_05 <- progress_bar$new(
  format = "  Calculating FV_Sp climate 0.5 [:bar] :percent in :elapsed",
  total = length(test2_05),
  clear = FALSE,
  width = 60
)

pb_fvsp3_05 <- progress_bar$new(
  format = "  Calculating FV_Sp climate + extinction 0.5 [:bar] :percent in :elapsed",
  total = length(test3_05),
  clear = FALSE,
  width = 60
)


lala4_05 <- vector("list", length(test1_05))
lala5_05 <- vector("list", length(test2_05))
lala6_05 <- vector("list", length(test3_05))


# Present-day FV_Sp
for (i in seq_along(test1_05)) {
  
  result <- tryCatch(
    get_FV_Sp(
      ax = c(1:4),
      pcoa = pcoa_05,
      Selected_sp = test1_05[[i]]
    ),
    error = function(e) NA
  )
  
  lala4_05[[i]] <- result
  
  pb_fvsp1_05$tick()
}


# Climate change FV_Sp
for (i in seq_along(test2_05)) {
  
  result2 <- tryCatch(
    get_FV_Sp(
      ax = c(1:4),
      pcoa = pcoa_05,
      Selected_sp = test2_05[[i]]
    ),
    error = function(e) NA
  )
  
  lala5_05[[i]] <- result2
  
  pb_fvsp2_05$tick()
}


# Climate change + extinctions FV_Sp
for (i in seq_along(test3_05)) {
  
  result3 <- tryCatch(
    get_FV_Sp(
      ax = c(1:4),
      pcoa = pcoa_05,
      Selected_sp = test3_05[[i]]
    ),
    error = function(e) NA
  )
  
  lala6_05[[i]] <- result3
  
  pb_fvsp3_05$tick()
}

# 21. Extract FRic — 0.5
fric_AqMap_05 <- sapply(
  seq_along(lala4_05),
  function(i) {
    result <- lala4_05[[i]]
    
    if (is.list(result)) {
      x <- result$data$Fric_r
      
      if (length(x) == 0) {
        NA_real_
      } else {
        as.numeric(x)
      }
    } else {
      NA_real_
    }
  }
)

names(fric_AqMap_05) <- names(test1_05)

fric_RCP45_05 <- sapply(
  seq_along(lala5_05),
  function(i) {
    x <- lala5_05[[i]]$data$Fric_r
    
    if (length(x) == 0) {
      NA_real_
    } else {
      as.numeric(x)
    }
  }
)

names(fric_RCP45_05) <- names(test2_05)

fric_RCP45ext_05 <- sapply(
  seq_along(lala6_05),
  function(i) {
    x <- lala6_05[[i]]$data$Fric_r
    
    if (length(x) == 0) {
      NA_real_
    } else {
      as.numeric(x)
    }
  }
)

names(fric_RCP45ext_05) <- names(test3_05)

# 22. Create FRic + SR + coordinate dataframes — 0.5
## Present-day grid information
grid_info_AqMap_05 <- AqMap_05_Occ %>%
  distinct(
    CsquareCode,
    CenterLat,
    CenterLong
  )

# RCP45 grid information
grid_info_RCP45_05 <- RCP45_05_Occ %>%
  distinct(
    CsquareCode,
    CenterLat,
    CenterLong
  )

# RCP45 + extinction grid information
grid_info_RCP45ext_05 <- RCP45ext_05_Occ %>%
  distinct(
    CsquareCode,
    CenterLat,
    CenterLong
  )

# Present day
FD_AqMap_05 <- data.frame(
  Grid = names(fric_AqMap_05),
  fric = as.numeric(fric_AqMap_05),
  row.names = NULL
) %>%
  left_join(
    baskets_AqMaps_weights_05 %>%
      rownames_to_column("Grid") %>%
      mutate(
        sp_richn = rowSums(
          select(., -Grid)
        )
      ) %>%
      select(Grid, sp_richn),
    by = "Grid"
  ) %>%
  left_join(
    grid_info_AqMap_05 %>%
      rename(Grid = CsquareCode),
    by = "Grid"
  )

save(FD_AqMap_05,file="C:/Users/Jack Cooper/Documents/iucn/iucn_sim/Maps/AqMap_05.RData")

# RCP45
FD_RCP45_05 <- data.frame(
  Grid = names(fric_RCP45_05),
  fric = as.numeric(fric_RCP45_05),
  row.names = NULL
) %>%
  left_join(
    baskets_RCP45_weights_05 %>%
      rownames_to_column("Grid") %>%
      mutate(
        sp_richn = rowSums(
          select(., -Grid)
        )
      ) %>%
      select(Grid, sp_richn),
    by = "Grid"
  ) %>%
  left_join(
    grid_info_RCP45_05 %>%
      rename(Grid = CsquareCode),
    by = "Grid"
  )

save(FD_RCP45_05,file="C:/Users/Jack Cooper/Documents/iucn/iucn_sim/Maps/RCP45_05.RData")

# RCP45 + extinction
FD_RCP45ext_05 <- data.frame(
  Grid = names(fric_RCP45ext_05),
  fric = as.numeric(fric_RCP45ext_05),
  row.names = NULL
) %>%
  left_join(
    baskets_RCP45ext_weights_05 %>%
      rownames_to_column("Grid") %>%
      mutate(
        sp_richn = rowSums(
          select(., -Grid)
        )
      ) %>%
      select(Grid, sp_richn),
    by = "Grid"
  ) %>%
  left_join(
    grid_info_RCP45ext_05 %>%
      rename(Grid = CsquareCode),
    by = "Grid"
  )

save(FD_RCP45ext_05,file="C:/Users/Jack Cooper/Documents/iucn/iucn_sim/Maps/RCP45ext_05.RData")

# 23. Calculate shifts — 0.5
FD_shift_RCP45_05 <- FD_AqMap_05 %>%
  inner_join(
    FD_RCP45_05,
    by = "Grid",
    suffix = c("_present", "_rcp45")
  ) %>%
  mutate(
    SR_shift = sp_richn_rcp45 - sp_richn_present,
    FRic_shift = fric_rcp45 - fric_present
  ) %>%
  select(
    Grid,
    CenterLat_present,
    CenterLong_present,
    SR_shift,
    FRic_shift
  ) %>%
  rename(
    CenterLat = CenterLat_present,
    CenterLong = CenterLong_present
  ) %>%
  na.omit()

FD_shift_RCP45ext_05 <- FD_AqMap_05 %>%
  inner_join(
    FD_RCP45ext_05,
    by = "Grid",
    suffix = c("_present", "_rcp45ext")
  ) %>%
  mutate(
    SR_shift = sp_richn_rcp45ext - sp_richn_present,
    FRic_shift = fric_rcp45ext - fric_present
  ) %>%
  select(
    Grid,
    CenterLat_present,
    CenterLong_present,
    SR_shift,
    FRic_shift
  ) %>%
  rename(
    CenterLat = CenterLat_present,
    CenterLong = CenterLong_present
  ) %>%
  na.omit()

# 24. Calculate grid changes — 0.5
sr_shift_RCP45_05 <- FD_shift_RCP45_05 %>%
  mutate(
    SR_change = case_when(
      SR_shift > 0 ~ "Increase",
      SR_shift < 0 ~ "Decrease",
      TRUE ~ "No change"
    )
  ) %>%
  count(SR_change) %>%
  mutate(Proportion = n / sum(n))

fric_shift_RCP45_05 <- FD_shift_RCP45_05 %>%
  mutate(
    FRic_change = case_when(
      FRic_shift > 0 ~ "Increase",
      FRic_shift < 0 ~ "Decrease",
      TRUE ~ "No change"
    )
  ) %>%
  count(FRic_change) %>%
  mutate(Proportion = n / sum(n))

sr_shift_RCP45ext_05 <- FD_shift_RCP45ext_05 %>%
  mutate(
    SR_change = case_when(
      SR_shift > 0 ~ "Increase",
      SR_shift < 0 ~ "Decrease",
      TRUE ~ "No change"
    )
  ) %>%
  count(SR_change) %>%
  mutate(Proportion = n / sum(n))

fric_shift_RCP45ext_05 <- FD_shift_RCP45ext_05 %>%
  mutate(
    FRic_change = case_when(
      FRic_shift > 0 ~ "Increase",
      FRic_shift < 0 ~ "Decrease",
      TRUE ~ "No change"
    )
  ) %>%
  count(FRic_change) %>%
  mutate(Proportion = n / sum(n))

# 25. Compare numbers of cells across thresholds
cat("\n\nNUMBER OF CELLS WITH VALID FRic SHIFTS\n")

cat("\n0.2 RCP45:", nrow(FD_shift_RCP45_02))
cat("\n0.2 RCP45 + extinction:", nrow(FD_shift_RCP45ext_02))

cat("\n0.3 RCP45:", nrow(FD_shift_RCP45))
cat("\n0.3 RCP45 + extinction:", nrow(FD_shift_RCP45ext))

cat("\n0.5 RCP45:", nrow(FD_shift_RCP45_05))
cat("\n0.5 RCP45 + extinction:", nrow(FD_shift_RCP45ext_05))

# 26. Display FRic gain/loss proportions - forms Table S7
fric_shift_RCP45
fric_shift_RCP45ext

fric_shift_RCP45_02
fric_shift_RCP45ext_02

fric_shift_RCP45_05
fric_shift_RCP45ext_05

# 27. Statistical tests on valid cells - produces Table S8
# Function to compare two thresholds using the same grid cells
compare_thresholds <- function(data_02, data_03, data_05) {
  
  # Identify grid cells with valid FRic shifts in ALL three thresholds
  common_cells <- Reduce(
    intersect,
    list(data_02$Grid, data_03$Grid, data_05$Grid)
  )
  
  # Extract FRic shifts for common cells
  comparison <- data.frame(
    Grid = common_cells,
    FRic_02 = data_02$FRic_shift[
      match(common_cells, data_02$Grid)
    ],
    FRic_03 = data_03$FRic_shift[
      match(common_cells, data_03$Grid)
    ],
    FRic_05 = data_05$FRic_shift[
      match(common_cells, data_05$Grid)
    ]
  ) %>%
    drop_na()
  
  # Differences relative to primary threshold (0.3)
  comparison <- comparison %>%
    mutate(
      Difference_02_vs_03 = FRic_02 - FRic_03,
      Difference_05_vs_03 = FRic_05 - FRic_03
    )
  
  # Wilcoxon tests
  test_02 <- wilcox.test(
    comparison$FRic_02,
    comparison$FRic_03,
    paired = TRUE,
    exact = FALSE
  )
  
  test_05 <- wilcox.test(
    comparison$FRic_05,
    comparison$FRic_03,
    paired = TRUE,
    exact = FALSE
  )
  
  # Effect size r = Z / sqrt(N)
  # Calculate Z from the Wilcoxon statistic
  calc_r <- function(test, x, y) {
    
    n <- sum(complete.cases(x, y))
    
    differences <- x - y
    differences <- differences[differences != 0]
    
    n_eff <- length(differences)
    
    W <- as.numeric(test$statistic)
    
    expected_W <- n_eff * (n_eff + 1) / 4
    sd_W <- sqrt(
      n_eff * (n_eff + 1) * (2 * n_eff + 1) / 24
    )
    
    Z <- (W - expected_W) / sd_W
    
    Z / sqrt(n_eff)
  }
  
  # Create results table
  results <- tibble(
    Comparison = c("0.2 vs 0.3", "0.5 vs 0.3"),
    Common_cells = nrow(comparison),
    
    Wilcoxon_V = c(
      as.numeric(test_02$statistic),
      as.numeric(test_05$statistic)
    ),
    
    P_value = c(
      test_02$p.value,
      test_05$p.value
    ),
    
    Median_difference = c(
      median(comparison$Difference_02_vs_03),
      median(comparison$Difference_05_vs_03)
    ),
    
    Mean_difference = c(
      mean(comparison$Difference_02_vs_03),
      mean(comparison$Difference_05_vs_03)
    ),
    
    IQR_difference = c(
      IQR(comparison$Difference_02_vs_03),
      IQR(comparison$Difference_05_vs_03)
    ),
    
    Effect_size_r = c(
      calc_r(
        test_02,
        comparison$FRic_02,
        comparison$FRic_03
      ),
      calc_r(
        test_05,
        comparison$FRic_05,
        comparison$FRic_03
      )
    )
  )
  
  return(results)
}

# RCP 4.5 — climate only
results_RCP45 <- compare_thresholds(
  FD_shift_RCP45_02,
  FD_shift_RCP45,
  FD_shift_RCP45_05
) %>%
  mutate(Scenario = "RCP 4.5")

# RCP 4.5 + extinction
results_RCP45ext <- compare_thresholds(
  FD_shift_RCP45ext_02,
  FD_shift_RCP45ext,
  FD_shift_RCP45ext_05
) %>%
  mutate(Scenario = "RCP 4.5 + extinction")

# Combine results into one table
threshold_sensitivity_results <- bind_rows(
  results_RCP45,
  results_RCP45ext
) %>%
  select(
    Scenario,
    Comparison,
    Common_cells,
    Wilcoxon_V,
    P_value,
    Median_difference,
    Mean_difference,
    IQR_difference,
    Effect_size_r
  )


# Display table
threshold_sensitivity_results
