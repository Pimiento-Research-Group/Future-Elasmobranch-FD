###############################################################################################################
# 02. IUCN simulations
## This R code updates IUCN status of all elasmobranchs and simulates their future using IUCN_sim
## it produces various outpit files containing extinction simulation results
###############################################################################################################

# Import packages
library(tidyverse)
library(devtools)
library(reticulate)
library(iucnsim)
library(tibble)
library(here)
library(progress)
library(purrr)
library(furrr)

# Load species and iucn data
load(file = "C:/Users/2022207/Dropbox/Jack's PhD/Chapter 3. Future shark FD/Data/R data/Species_list_IUCN.RData")

## Load python functions from IUCN_sim into R
reticulate::install_python(version = '3.11.4')
reticulate::source_python(here("C:/Users/2022207/Dropbox/Jack's PhD/Chapter 3. Future shark FD/IUCN_sim/iucn_sim.py"))

# set reference group and define iucn key
reference_group = "Chondrichthyes"
reference_rank = "class" 
iucn_key="1fecbeea639ba430f60510af483c5d4b282e3f097aa98e4613003c9903970df6"

## Download IUCN history of reference group
outdir = 'data/iucn_sim/iucn_data' # define where you want the output files to be stored

# get iucn history of reference group, will be written to file
iucn_history_file = get_iucn_history(reference_group=reference_group,
                                     reference_rank=reference_rank,
                                     iucn_key=iucn_key,
                                     outdir=outdir)

# Import IUCN_history & check column matches with final species list - add 2023 status to history file
IUCN_history <- read.table("data/iucn_sim/iucn_data\\CHONDRICHTHYES_iucn_history.txt", header = TRUE, sep = "\t")

# Check which species from iucn are missing from history file & filter out
missing_species <- setdiff(species_iucn_updated_fixed$Species, IUCN_history$species)

filtered_species_iucn_updated <- species_iucn_updated_fixed %>%
  filter(Species %in% missing_species)
filtered_species_iucn_updated <- filtered_species_iucn_updated[order(filtered_species_iucn_updated$Species),]
filtered_species_iucn_updated <- filtered_species_iucn_updated %>% rename(X2023 = category)

filtered_species_iucn_history <- species_iucn_updated_fixed %>%
  filter(!Species %in% missing_species)
filtered_species_iucn_history <- filtered_species_iucn_history[order(filtered_species_iucn_history$Species),]
filtered_species_iucn_history <- filtered_species_iucn_history %>% rename(X2023 = category)

# Order species columns and add current status based on filtered species (guarantees status matches)
IUCN_history <- IUCN_history[order(IUCN_history$species),] %>% rename(Species = species)

IUCN_history_2023 <- IUCN_history %>%
  left_join(filtered_species_iucn_history %>% select(Species, X2023), by = c("Species" = "Species"))

IUCN_history_2023$X2023[is.na(IUCN_history_2023$X2023)] <- ""

# Add all 20 other columns from IUCN_history to filtered_species_iucn_updated as NAs & reorder columns
filtered_species_iucn_updated$X2015<-NA;filtered_species_iucn_updated$X2009<-NA;filtered_species_iucn_updated$X2000<-NA;
filtered_species_iucn_updated$X1996<-NA;filtered_species_iucn_updated$X2019<-NA;filtered_species_iucn_updated$X1994<-NA;
filtered_species_iucn_updated$X1990<-NA;filtered_species_iucn_updated$X2018<-NA;filtered_species_iucn_updated$X2007<-NA;
filtered_species_iucn_updated$X2006<-NA;filtered_species_iucn_updated$X2013<-NA;filtered_species_iucn_updated$X2016<-NA;
filtered_species_iucn_updated$X2005<-NA;filtered_species_iucn_updated$X2003<-NA;filtered_species_iucn_updated$X2002<-NA;
filtered_species_iucn_updated$X2004<-NA;filtered_species_iucn_updated$X2017<-NA;filtered_species_iucn_updated$X2011<-NA;
filtered_species_iucn_updated$X2012<-NA;filtered_species_iucn_updated$X2010<-NA

filtered_species_iucn_updated <- filtered_species_iucn_updated %>% 
  select(2, 3:22, 1, everything())

# Merge filtered species with history file by shared columns
IUCN_final <- rbind(IUCN_history_2023, filtered_species_iucn_updated) %>% rename(species = Species)

IUCN_final[is.na(IUCN_final)] <- ""
colnames(IUCN_final)[2:22] <- sub("^X", "", colnames(IUCN_final)[2:22])

# Specify the file path for the new file
output_file <- "C:/Users/2022207/Dropbox/Jack's PhD/Chapter 3. Future shark FD/Analyses/IUCN_sim/data/iucn_sim/iucn_data\\IUCN_final.txt"

# Write the data frame to a tab-separated text file
write.table(IUCN_final, output_file, sep = "\t", quote = FALSE, row.names = FALSE)

# Import IUCN history with 2023 status
IUCN_history_final <- "C:/Users/2022207/Dropbox/Jack's PhD/Chapter 3. Future shark FD/Analyses/IUCN_sim/data/iucn_sim/iucn_data\\IUCN_final.txt"

## Get IUCN history transitions
counted_status_transition_events = evaluate_iucn_history(IUCN_history_final)

# Get species_list
species_list <- species_iucn_updated_fixed$Species

# get most recent status for each taxon in target species list
extant_taxa_current_status = get_most_recent_status_target_species(species_list=species_list,
                                                                   iucn_history_file=IUCN_history_final,
                                                                   iucn_key=iucn_key,
                                                                   outdir=outdir)

# Mark possibly extinct species
possibly_extinct_taxa = get_possibly_extinct_iucn_info(IUCN_history_final,
                                                       outdir=outdir)

# Estimate status transition rates
outdir = "C:/Users/2022207/Dropbox/Jack's PhD/Chapter 3. Future shark FD/Analyses/IUCN_sim/data/iucn_sim/transition_rates"
transition_rates_out = estimate_transition_rates(extant_taxa_current_status,
                                                 IUCN_history_final,
                                                 outdir,
                                                 extinction_probs_mode=1,
                                                 possibly_extinct_list=possibly_extinct_taxa,
                                                 rate_samples=100,
                                                 n_gen = 100000,
                                                 burnin = 1000)

# Simulate future (next 500 years)
## Extract future status of species
outdir = "C:/Users/2022207/Dropbox/Jack's PhD/Chapter 3. Future shark FD/Analyses/IUCN_sim/data/iucn_sim/future_simulations_future_status"
sim_years = 1000
# Run simulations (supercomputer)
future_sim_output = run_future_sim(transition_rates_out,
                                   outdir,
                                   n_years=sim_years,
                                   n_sim=10000,
                                   status_change = 1,
                                   conservation_increase_factor = 1,
                                   threat_increase_factor = 1,
                                   model_unknown_as_lc = 0,
                                   until_n_taxa_extinct = 0,
                                   plot_diversity_trajectory = 1,
                                   plot_status_trajectories = 1,
                                   plot_histograms = 1,
                                   plot_status_piechart = 1,
                                   save_future_status_array = TRUE)

# extract the different output items
extinction_times = future_sim_output[[1]]
save(extinction_times, file = "C:/Users/2022207/Dropbox/Jack's PhD/Chapter 3. Future shark FD/Data/R data/IUCN_sim_extinction_times.RData")
future_div_min_max = future_sim_output[[2]]
status_through_time_trajectories = future_sim_output[[3]]

# read the resulting pickle file and calculate the probability of each status per year per species
future_status_file = paste0(outdir,'/future_status_array_list.pkl')
future_status_data = summarize_future_status_array(future_status_file)
future_status_probs = future_status_data[[1]]
most_likely_future_status = future_status_data[[2]]

# Estimate extinction rates
outdir = "C:/Users/2022207/Dropbox/Jack's PhD/Chapter 3. Future shark FD/Analyses/IUCN_sim/data/iucn_sim/extinction_rates"
ext_rates = estimate_extinction_rates(extinction_times,
                                      sim_years,
                                      outdir,
                                      load_from_file=FALSE)

# Produce a dataframe to determine which species are extinct in the future
dat.ext <- read_delim("C:/Users/2022207/Dropbox/Jack's PhD/Chapter 3. Future shark FD/Analyses/IUCN_sim/data/iucn_sim/future_simulations_future_status/te_all_species.txt")

# Make dat.ext a dataframe & tidy
ext.df <- dat.ext %>% as.data.frame()

ext.df1 <- rbind(colnames(ext.df), ext.df)

ext.df1[1, ] <- lapply(ext.df1[1, ], function(x) ifelse(is.na(as.numeric(x)), NA, as.numeric(x)))

ext.df1[1, 1] <- "Acroteriobatus annulatus"

# Form mode database
ext.df1$iucn <- species_iucn_updated_fixed$category

# Calculate mode extinction time per species across simulations
calculate_mode <- function(x) {
  x <- x[!is.na(x)]  # Remove NAs
  if (length(x) == 0) return(NA)  # Return NA if all values are NA
  ux <- unique(x)
  counts <- tabulate(match(x, ux))
  mode_values <- ux[counts == max(counts)]
  return(mode_values[1])  
}

ext.df1$mode <- sapply(1:nrow(ext.df1), function(i) calculate_mode(ext.df1[i, -1]))

ext.times.mode <- ext.df1 %>%
  select(`Acroteriobatus annulatus`,mode) %>%
  mutate(mode = as.numeric(mode)) %>%
  arrange(mode)
ext.times.mode$mode <- as.numeric(ext.times.mode$mode)

Spp_mode <- ext.times.mode %>% 
  rename(Species = `Acroteriobatus annulatus`)
save(Spp_mode,file = "C:/Users/2022207/Dropbox/Jack's PhD/Chapter 3. Future shark FD/Data/R data/2100_mode.RData")

# Filter out rows where the mode is under 77 - identifies species extinct by 2100
extinct_spp <- ext.times.mode %>%
  filter(mode < 77) %>%
  rename(Species = `Acroteriobatus annulatus`)
save(extinct_spp,file = "C:/Users/2022207/Dropbox/Jack's PhD/Chapter 3. Future shark FD/Data/R data/2100_extinct_species.RData")
