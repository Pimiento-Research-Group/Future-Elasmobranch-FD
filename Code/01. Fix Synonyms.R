###############################################################################################################
# 01. Fix synonyms
## This R code updates IUCN status of all elasmobranchs and corrects for synonyms
## it produces one Rdata file of the final species list and current IUCN status
###############################################################################################################

# Import packages
library(readxl)
library(tidyverse)
library(janitor)
library(rredlist)
library(purrr)

# API token
Sys.setenv(IUCN_KEY = "1fecbeea639ba430f60510af483c5d4b282e3f097aa98e4613003c9903970df6")
Sys.getenv("IUCN_KEY")
apikey <- Sys.getenv("IUCN_KEY") # needed to access rredlist package

# Load species names downloaded from IUCN
species <- read_xlsx("C:/Users/2022207/Dropbox/Jack's PhD/Chapter 3. Future shark FD/Analyses/IUCN_sim/Species_iucn.xlsx")

# Get IUCN info from downloaded species
species_iucn<-species %>%
  dplyr::select(Scientific_name) %>% 
  mutate(iucn_pull = map(Scientific_name, rl_search, key = apikey))

## clean dataset and gather only iucn statuses
species_iucn_clean <-  species_iucn %>%
  mutate(category = map_chr(iucn_pull, pluck, "result", "category", .default = NA)) %>%
  dplyr::select(Scientific_name, category) %>% 
  bind_cols(species$Scientific_name) %>%
  rename(trait.names=...3)

## Order species alphabetically for checking against trait data later
species_iucn_clean <- species_iucn_clean[order(species_iucn_clean$Scientific_name),]

# identify duplicated names
synonyms <- species_iucn_clean %>%
  group_by(trait.names) %>% 
  filter(n()>1)

# identify synonyms that did not result in duplicated names
which(is.na(match(species_iucn_clean$Scientific_name, species_iucn_clean$trait.names))==T)

# eliminate duplicated and leave status of matching species names
synonyms <-  synonyms[row.names(synonyms) %in% 
                        unique(match(synonyms$trait.names, synonyms$Scientific_name)),]

species_iucn_updated <- species_iucn_clean %>%
  filter(!trait.names %in% synonyms$trait.names) %>%
  bind_rows(synonyms) %>%
  dplyr::select(-Scientific_name)

which(duplicated(species_iucn_updated$trait.names)==T) #no duplicates

# species which names changed in iucn website when we downloaded the data
wrong.names <- species_iucn_updated %>%
  dplyr:: filter(is.na(category))

species_iucn_updated_fixed <- species_iucn_updated %>%
  mutate(category=replace(category, trait.names=="Dipturus australis", "NT")) %>%
  mutate(category=replace(category, trait.names=="Dipturus cerva", "NT")) %>%
  mutate(category=replace(category, trait.names=="Dipturus confusus", "CR")) %>%
  mutate(category=replace(category, trait.names=="Dipturus endeavouri", "NT")) %>%
  mutate(category=replace(category, trait.names=="Dipturus grahami", "LC")) %>%
  mutate(category=replace(category, trait.names=="Dipturus healdi", "LC")) %>%
  mutate(category=replace(category, trait.names=="Zearaja maugeana", "EN")) %>%
  mutate(category=replace(category, trait.names=="Dipturus oculus", "LC")) %>%
  mutate(category=replace(category, trait.names=="Dipturus polyommata", "LC")) %>%
  mutate(category=replace(category, trait.names=="Himantura oxyrhyncha", "EN")) %>%
  mutate(category=replace(category, trait.names=="Myliobatis freminvillei", "VU")) %>%
  mutate(category=replace(category, trait.names=="Narcine lasti", "LC")) %>%
  mutate(category=replace(category, trait.names=="Narcine nelsoni", "LC")) %>%
  mutate(category=replace(category, trait.names=="Narcine ornata", "LC")) %>%
  mutate(category=replace(category, trait.names=="Narcine tasmaniensis", "LC")) %>%
  mutate(category=replace(category, trait.names=="Narcine westraliensis", "LC")) %>%
  mutate(category=replace(category, trait.names=="Psammobatis parvacauda", "LC")) %>%
  mutate(category=replace(category, trait.names=="Raja inornata", "LC")) %>%
  mutate(category=replace(category, trait.names=="Raja stellulata", "LC")) %>%
  mutate(category=replace(category, trait.names=="Taeniura grabata", "NT"))

# add missing species from trait dataset
# load look up table with synonyms
all.names <-read_xlsx("C:/Users/2022207/Dropbox/Jack's PhD/Chapter 3. Future shark FD/Analyses/Lookup_Taxonomy.xlsx", 
                      sheet = "Species") %>%
  dplyr::select(-Family, -Order, -Superorder)

# what species are not already in data?
missing.sp <- all.names %>%
  filter(!Species %in% species_iucn_updated_fixed$trait.names) 

# Identify synonyms
syn.1 =  missing.sp %>%
  dplyr::select(Synonyms.1) %>%
  dplyr::filter(!Synonyms.1=="NA") %>%
  dplyr::rename(Species = Synonyms.1)
syn.2 =  missing.sp %>%
  dplyr::select(Synonyms.2) %>%
  dplyr::filter(!Synonyms.2=="NA")%>%
  dplyr::rename(Species = Synonyms.2)
syn.3 =  missing.sp %>%
  dplyr::select(Synonyms.3) %>%
  dplyr::filter(!Synonyms.3=="NA")%>%
  dplyr::rename(Species = Synonyms.3)
syn.4 =  missing.sp %>%
  dplyr::select(Synonyms.4) %>%
  dplyr::filter(!Synonyms.4=="NA")%>%
  dplyr::rename(Species = Synonyms.4)
syn.5 =  missing.sp %>%
  dplyr::select(Synonyms.5) %>%
  dplyr::filter(!Synonyms.5=="NA")%>%
  dplyr::rename(Species = Synonyms.5)
syn.6 =  missing.sp %>%
  dplyr::select(Synonyms.6) %>%
  dplyr::filter(!Synonyms.6=="NA")%>%
  dplyr::rename(Species = Synonyms.6)
syn.7 =  missing.sp %>%
  dplyr::select(Synonyms.7) %>%
  dplyr::filter(!Synonyms.7=="NA")%>%
  dplyr::rename(Species = Synonyms.7)
syn.8 =  missing.sp %>%
  dplyr::select(Synonyms.8) %>%
  dplyr::filter(!Synonyms.8=="NA")%>%
  dplyr::rename(Species = Synonyms.8)
syn.9 =  missing.sp %>%
  dplyr::select(Synonyms.9) %>%
  dplyr::filter(!Synonyms.9=="NA")%>%
  dplyr::rename(Species = Synonyms.9)
syn.10 =  missing.sp %>%
  dplyr::select(Synonyms.10) %>%
  dplyr::filter(!Synonyms.10=="NA")%>%
  dplyr::rename(Species = Synonyms.10)

# check the iucn category of these species that are missing and their synonyms
missing.sp.syn <- missing.sp %>%
  dplyr::select(Species) %>%
  bind_rows(syn.1,syn.2,syn.3,syn.4,syn.5) %>%
  dplyr::select(Species) %>% # now apply the rl_search function to each species using purrr::map()
  mutate(iucn_pull = map(Species, rl_search, key = apikey)) 

missing.sp.syn_clean <- missing.sp.syn %>%
  mutate(category = map_chr(iucn_pull, pluck, "result", "category", .default = NA)) %>% 
  dplyr::select(Species, category) %>% 
  dplyr::filter(!is.na(category)) %>%
  dplyr::select(category, Species)

intersect(missing.sp.syn_clean$Species, missing.sp$Species)

species_iucn_updated_no.syn <- species_iucn_updated_fixed %>%
  rename(Species = trait.names) %>%
  bind_rows(missing.sp.syn_clean)

# get still missing species
setdiff(all.names$Species, species_iucn_updated_no.syn$Species) 
# these are all not evaluated and are considered junior synonyms of other species

# add not evaluated species
species_iucn_updated_no.syn <- species_iucn_updated_no.syn %>%
  bind_rows(
    tibble(setdiff(all.names$Species,species_iucn_updated_no.syn$Species)) %>%
      setNames("Species") %>%
      mutate(category="NE")
  )

# get their synonyms from iucn
species_iucn_updated_fixed <- species_iucn_updated_fixed %>% 
  rename(Species = trait.names)

missing.syns <- as_tibble(setdiff(all.names$Species, species_iucn_updated_fixed$Species)) %>%
  setNames("Species")%>%
  mutate(iucn_pull = map(Species, rl_synonyms, key = apikey))

missing.syns_clean <- missing.syns %>% 
  mutate(accepted_name = map_chr(iucn_pull, ~ if(length(.x$result$accepted_name) > 0) .x$result$accepted_name[1] else NA,.default = NA)) %>% 
  select(Species, accepted_name) %>%
  filter(!is.na(accepted_name))

# Check if any synonyms now occur in current species list as well as if all accepted names are in list
intersect(missing.syns_clean$Species, species_iucn_updated_fixed$Species)   # False - none in iucn_updated_fixed

all(missing.syns_clean$accepted_name %in% species_iucn_updated_fixed$Species)   # TRUE - no synonyms present

# Remove chimeras from species_list
chim1 <- str_subset(species_iucn_updated_fixed$Species, "Callorhinchus")
chim2 <- str_subset(species_iucn_updated_fixed$Species, "Chimaera")
chim3 <- str_subset(species_iucn_updated_fixed$Species, "Hydrolagus")
chim4 <- str_subset(species_iucn_updated_fixed$Species, "Harriotta")
chim5 <- str_subset(species_iucn_updated_fixed$Species, "Neoharriotta")
chim6 <- str_subset(species_iucn_updated_fixed$Species, "Rhinochimaera")

species_iucn_updated_fixed <- species_iucn_updated_fixed %>%
  filter(!Species %in% chim1 & !Species%in%chim2 & !Species%in%chim3 &
           !Species%in%chim4 & !Species%in%chim5 & !Species%in%chim6)

# Check that your species are all within the list of Dulvy et al. 2021
Dulvy <- read_xlsx("C:/Users/2022207/Dropbox/Jack's PhD/Chapter 3. Future shark FD/Data/Dulvy_2021.xlsx")

missing_species <- setdiff(species_iucn_updated_fixed$Species, Dulvy$`Latin binomial`)

# Display the missing species, if any
if (length(missing_species) > 0) {
  print("Some species are missing:")
  print(missing_species)
} else {
  print("All species are present in Dulvy.")
}

synonym_dataframes <- list(syn.1, syn.2, syn.3, syn.4, syn.5, syn.6, syn.7, syn.8, syn.9, syn.10)

# Check which missing species are present in any of the synonym dataframes
missing_species <- data.frame(Species = missing_species)

missing_species_synonyms <- lapply(synonym_dataframes, function(syn_df) {
  species_present <- missing_species$Species %in% syn_df$Species
  if (any(species_present)) {
    missing_species$Species[species_present]
  } else {
    NULL
  }
})

# Filter out NULL values from the list
missing_species_synonyms <- Filter(Negate(is.null), missing_species_synonyms)

# Display the missing species and their corresponding synonyms
if (length(missing_species_synonyms) > 0) {
  print("Missing species and their synonyms:")
  for (i in seq_along(missing_species_synonyms)) {
    print(paste("Synonym source:", paste0("syn.", i)))
    print(missing_species_synonyms[[i]])
  }
} else {
  print("No missing species found among synonyms.")
}

# Save final data for use in simulations and FD analyses
save(species_iucn_updated_fixed,file = "C:/Users/2022207/Dropbox/Jack's PhD/Chapter 3. Future shark FD/Data/R data/Species_list_IUCN.RData")