###############################################################################################################
# 06. Figure 1
## This R code plots future functional diversity (FD) based on our simulations
## it produces Figure 1
###############################################################################################################

## Import packages
library(tidyverse)
library(tibble)
library(reshape2)
library(cowplot)
library(ggsci)
library(RColorBrewer)
library(purrr)

## Load species and iucn data
load(file = "C:/Users/2022207/Dropbox/Jack's PhD/Chapter 3. Future shark FD/Data/R data/Species_list_IUCN.RData")

# Load boxplots and null model
load(file = "C:/Users/2022207/Dropbox/Jack's PhD/Chapter 3. Future shark FD/Analyses/Final paper R codes/data/Crit_E/CritE_full_buffer.RData")
load(file = "C:/Users/2022207/Dropbox/Jack's PhD/Chapter 3. Future shark FD/Analyses/Final paper R codes/data/Crit_E/CritE_full_long_buffer.RData")
load(file = "C:/Users/2022207/Dropbox/Jack's PhD/Chapter 3. Future shark FD/Analyses/Final paper R codes/data/Crit_E/Null_long.RData")
load(file = "C:/Users/2022207/Dropbox/Jack's PhD/Chapter 3. Future shark FD/Analyses/Final paper R codes/data/Crit_E/Null_wide.RData")

## Calculate % changes for empirical, buffer and null in species & functional diversity - Crit E EX mode
# Empirical
FDindices_taxon.CritE_long <- melt(FDmetrics_taxonvar.CritE, id.vars= "Scenario")

FDmetrics_Spp.CritE <- FDindices_taxon.CritE_long %>% 
  filter(variable == "sp_richn")
FDmetrics_FRic.CritE <- FDindices_taxon.CritE_long %>% 
  filter(variable == "fric")

FDmetrics_Spp.CritE <- FDmetrics_Spp.CritE %>% 
  mutate(percent_change = ((value - 1100) / 1100) * 100)
FDmetrics_FRic.CritE <- FDmetrics_FRic.CritE %>% 
  mutate(percent_change = ((value - 1) / 1) * 100)

# Calculate mean changes per scenario per metric
mean_percent_change_Spp <- FDmetrics_Spp.CritE %>%
  group_by(Scenario) %>%
  summarise(mean_percent_change = mean(percent_change, na.rm = TRUE))

mean_percent_change_FRic <- FDmetrics_FRic.CritE %>%
  group_by(Scenario) %>%
  summarise(mean_percent_change = mean(percent_change, na.rm = TRUE))

# Null
SPPmetrics_null_Taxonvar.CritE <- FDmetrics_null_long_taxonvar.CritE %>% 
  filter(variable == "sp_richn")
FRicmetrics_null_Taxonvar.CritE <- FDmetrics_null_long_taxonvar.CritE %>% 
  filter(variable == "fric")

SPPmetrics_null_Taxonvar.CritE <- SPPmetrics_null_Taxonvar.CritE %>% 
  mutate(percent_change = ((value - 1100) / 1100) * 100)
FRicmetrics_null_Taxonvar.CritE <- FRicmetrics_null_Taxonvar.CritE %>% 
  mutate(percent_change = ((value - 1) / 1) * 100)

## Plot results
# Filter only 2100 scenario (assuming "Future" = 2100)
FDmetrics_Spp_2100 <- FDmetrics_Spp.CritE %>% filter(Scenario == "Future")
FRicmetrics_empirical_2100 <- FDmetrics_FRic.CritE %>% filter(Scenario == "Future")
FRicmetrics_null_2100 <- FRicmetrics_null_Taxonvar.CritE %>% filter(Scenario == "Future")

# Plot 1: Species Richness (SR) - 2100 only
Spp_CritE_2100 <- ggplot(data = FDmetrics_Spp_2100, aes(x = Scenario, y = percent_change)) +
  geom_boxplot(fill = "#ca4678") +  
  labs(x = "", y = "% change") +
  ggtitle("Species Richness (SR)") +
  theme_bw() +
  theme(
    legend.position = "none",
    axis.text = element_text(size = 8, color = "black"),
    axis.text.x = element_blank(),
    axis.ticks.x = element_blank(),
    axis.title = element_text(size = 10),
    panel.background = element_rect(fill = "white"),
    panel.grid = element_blank(),
    panel.border = element_rect(fill = "transparent"),
    plot.title = element_text(hjust = 0.5, size = 10)
  )
save(Spp_CritE_2100, file = "C:/Users/2022207/Dropbox/Jack's PhD/Chapter 3. Future shark FD/Data/R data/Spp_2100.RData")

# Plot 2: Functional Richness (FRic) - 2100 only
FRic_CritE_2100 <- ggplot(data = FRicmetrics_null_2100, aes(x = Scenario, y = percent_change)) +
  geom_violin(trim = TRUE, fill = 'grey', color = "grey", alpha = 0.25, scale = "width") +
  geom_boxplot(data = FRicmetrics_empirical_2100, aes(x = Scenario, y = percent_change), 
               fill = "#f78c33") +  
  labs(x = "", y = "% change") +
  ggtitle("Functional Richness (FRic)") +
  theme_bw() +
  theme(
    legend.position = "none",
    axis.text = element_text(size = 8, color = "black"),
    axis.text.x = element_blank(),
    axis.ticks.x = element_blank(),
    axis.title = element_text(size = 10),
    panel.background = element_rect(fill = "white"),
    panel.grid = element_blank(),
    panel.border = element_rect(fill = "transparent"),
    plot.title = element_text(hjust = 0.5, size = 10)
  )
save(FRic_CritE_2100, file = "C:/Users/2022207/Dropbox/Jack's PhD/Chapter 3. Future shark FD/Data/R data/FRic_2100.RData")

# Wilcoxon signed-ranked tests
perform_wilcoxon_test <- function(empirical_data, null_data, scenario_name) {
  # Filter empirical and null data for the specific scenario
  empirical_values <- empirical_data %>% 
    filter(Scenario == scenario_name) %>% 
    pull(percent_change)
  
  null_values <- null_data %>% 
    filter(Scenario == scenario_name) %>% 
    pull(percent_change)
  
  # Perform the Wilcoxon signed-rank test
  wilcox.test(empirical_values, null_values, paired = TRUE)
}

# Set scenarios
scenarios <- c("Future")

# Perform tests for FRic
cat("\nWilcoxon test results for FRic:\n")
for (scenario in scenarios) {
  result <- perform_wilcoxon_test(FDmetrics_FRic.CritE, FRicmetrics_null_Taxonvar.CritE, scenario)
  cat(paste("Scenario:", scenario, "\n"))
  print(result)
}
