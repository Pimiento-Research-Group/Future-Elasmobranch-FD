###############################################################################################################
# 02. IUCN_sim
## This R code updates IUCN status of all elasmobranchs and simulates their future using IUCN_sim
## it uses the Crit E EX mode approach - based on IUCN status extinction probabilities; no GL data
###############################################################################################################

# Import packages
library(readxl)
library(tidyverse)
library(devtools)
library(rredlist)
library(reticulate)
library(iucnsim)
library(tibble)
library(here)

# Load species and iucn data
load(file = "~/Species_list_IUCN.RData")

## Load python functions from IUCN_sim into R
reticulate::source_python(here("~/iucn_sim.py"))

# set reference group and define iucn key
reference_group = "Chondrichthyes"
reference_rank = "class" 
iucn_key="1fecbeea639ba430f60510af483c5d4b282e3f097aa98e4613003c9903970df6"

## Download IUCN history of reference group
outdir = "~/iucn_data" # define where you want the output files to be stored

# get iucn history of reference group, will be written to file
iucn_history_file = get_iucn_history(reference_group=reference_group,
                                     reference_rank=reference_rank,
                                     iucn_key=iucn_key,
                                     outdir=outdir)

# Import IUCN_history 
IUCN_history <- read.table("~\\CHONDRICHTHYES_iucn_history.txt", header = TRUE, sep = "\t")

# Download IUCN history from Dulvy et al. 2024 and convert to history file
D24_hist <- read_xlsx("~/Dulvy_2024.xlsx")

# Isolate species and IUCN history from Dulvy et al. 2024
IUCN_D24 <- D24_hist %>%
  select(species, y1970, y1980, y2005, y2020)

# Standardise IUCN history and 2023 status files by species name
IUCN_D24 <- IUCN_D24 %>%
  arrange(species)

species_iucn_updated_fixed <- species_iucn_updated_fixed %>%
  rename(species = Species) %>%
  arrange(species) %>%
  rename(y2023 = category)

# Form final IUCN history file for analyses
IUCN_D24 <- IUCN_D24 %>%
  left_join(species_iucn_updated_fixed %>% select(species, y2023), by = "species")

IUCN_D24 <- IUCN_D24 %>%
  mutate(
    y2020 = ifelse(
      species %in% c("Carcharhinus obsoletus", "Torpedo suessii", "Urolophus javanicus"),
      "CR",  # Set to CR to avoid error; possibly_extinct taxa will account for this"
      y2020
    )
  )

IUCN_D24_sub <- IUCN_D24 %>%
  filter(!is.na(y2023))

names(IUCN_D24_sub) <- sub("^y", "", names(IUCN_D24))

# Rename to final file
IUCN_D24_final <- IUCN_D24_sub

# Specify the file path for the new file
output_file <- "~\\IUCN_final_Dulvy2024.txt"

# Write the data frame to a tab-separated text file
write.table(IUCN_D24_final, output_file, sep = "\t", quote = FALSE, row.names = FALSE)

# Import IUCN history
IUCN_history_final <- "~\\IUCN_final_Dulvy2024.txt"

## Get IUCN history transitions
counted_status_transition_events = evaluate_iucn_history(IUCN_history_final)

# Get species_list
species_list <- IUCN_D24_final$species

# get most recent status for each taxon in target species list
extant_taxa_current_status = get_most_recent_status_target_species(species_list=species_list,
                                                                   iucn_history_file=IUCN_history_final,
                                                                   iucn_key=iucn_key,
                                                                   outdir=outdir)

# Mark possibly extinct species
possibly_extinct_taxa = get_possibly_extinct_iucn_info(IUCN_history_final,
                                                       outdir=outdir)

# Set current status to 2023
extant_taxa_current_status <- species_iucn_updated_fixed %>%
  select(species, y2023, everything())

# Estimate status transition rates
outdir = "~/transition_rates_D24_CritE"
transition_rates_out = estimate_transition_rates(extant_taxa_current_status,
                                                 IUCN_history_final,
                                                 outdir,
                                                 extinction_probs_mode=0,
                                                 possibly_extinct_list=possibly_extinct_taxa,
                                                 rate_samples=100)

# Simulate future (next 1000 years; 77 years = 2100)
## Extract future status of species
outdir = "~/future_simulations_future_status_D24_CritE"
sim_years = 1000
# Run simulations (supercomputer)
future_sim_output = run_future_sim(transition_rates_out,
                                   outdir,
                                   n_years=sim_years,
                                   n_sim=10000)

# extract the different output items
extinction_times.CritE = future_sim_output[[1]]
future_div_min_max.CritE = future_sim_output[[2]]
status_through_time_trajectories.CritE = future_sim_output[[3]]

# read the resulting pickle file and calculate the probability of each status per year per species
future_status_file = paste0(outdir,'/future_status_array_list.pkl')
future_status_data = summarize_future_status_array(future_status_file)
future_status_probs = future_status_data[[1]]
most_likely_future_status = future_status_data[[2]]

# Estimate extinction rates
outdir = "C:/Users/Jack Cooper/Documents/iucn/iucn_sim/iucn_data/extinction_rates_D24_CritE"
ext_rates = estimate_extinction_rates(extinction_times.CritE,
                                      sim_years,
                                      outdir,
                                      load_from_file=FALSE)

# Produce a dataframe to determine which species are extinct in the future
dat.ext <- read_delim("~/te_all_species.txt")

# Make dat.ext a dataframe & tidy
ext.df <- dat.ext %>% as.data.frame()

ext.df1 <- rbind(colnames(ext.df), ext.df)

ext.df1[1, ] <- lapply(ext.df1[1, ], function(x) ifelse(is.na(as.numeric(x)), NA, as.numeric(x)))

ext.df1[1, 1] <- "Acroteriobatus annulatus"

# Calculate and save median (for spatial analyses)
# Define function to calculate median extinction time, ignoring NA values
calculate_median <- function(x) {
  x <- as.numeric(x[!is.na(x)])  # Remove NAs and ensure values are numeric
  if (length(x) == 0) return(NA)  # Return NA if all values are NA
  return(median(x))  # Calculate and return median
}

# Apply median calculation to each row (species)
ext.df1$median <- sapply(1:nrow(ext.df1), function(i) calculate_median(ext.df1[i, -1]))

# Create a summary table with species name and median extinction time
ext.times.median <- ext.df1 %>%
  select(`Acroteriobatus annulatus`, median) %>%
  mutate(median = as.numeric(median)) %>%
  arrange(median)

# Rename column for clarity
Spp_median.CritE <- ext.times.median %>%
  rename(Species = `Acroteriobatus annulatus`)
save(Spp_median.CritE,file = "~/CritE_median_D24.RData")

# Filter out rows where the median is under 77 - identifies species extinct by 2100 
extinct_spp_2100.CritE <- ext.times.median %>%
  filter(median < 77) %>%
  rename(Species = `Acroteriobatus annulatus`)
save(extinct_spp_2100.CritE,file = "~/CritE_2100_extinct_species_D24.RData")

