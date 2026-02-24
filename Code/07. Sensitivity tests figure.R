###############################################################################################################
# 07. Sensitivity tests figure
## This R code plots future FRic changes based on sensitivity tests
## it produces Figure S1 and Table S1
###############################################################################################################

## Import packages
library(tidyverse)
library(tibble)
library(reshape2)
library(cowplot)
library(ggsci)
library(RColorBrewer)
library(purrr)

# Load simulated results
load("~/CritE_full_buffer_minushabitat.RData")
load("~/CritE_full_long_buffer_minushabitat.RData")
load("~/CritE_full_buffer_minusVP.RData")
load("~/CritE_full_long_buffer_minusVP.RData")
load("~/CritE_full_buffer_minusterrest.RData")
load("~/CritE_full_long_buffer_minusterrest.RData")
load("~/CritE_full_buffer_minusthermo.RData")
load("~/CritE_full_long_buffer_minusthermo.RData")
load("~/CritE_full_buffer_minusFM.RData")
load("~/CritE_full_long_buffer_minusFM.RData")
load("~/CritE_full_buffer_minusdiet.RData")
load("~/CritE_full_long_buffer_minusdiet.RData")
load("~/CritE_full_buffer_minusTL.RData")
load("~/CritE_full_long_buffer_minusTL.RData")
load("~/CritE_full_buffer_minusdeepwater.RData")
load("~/CritE_full_long_buffer_minusdeepwater.RData")
load("~/CritE_full_buffer_minusbody.RData")
load("~/CritE_full_long_buffer_minusbody.RData")
load("~/CritE_full_buffer_minusreproduction.RData")
load("~/CritE_full_long_buffer_minusreproduction.RData")

# Calculate % changes for species and FRic
## Minus habitat
FDindices_taxon.habitat_long <- melt(FDmetrics_taxonvar.habitat, id.vars= "Scenario")

FDmetrics_Spp.habitat <- FDindices_taxon.habitat_long %>% 
  filter(variable == "sp_richn")
FDmetrics_FRic.habitat <- FDindices_taxon.habitat_long %>% 
  filter(variable == "fric")

FDmetrics_Spp.habitat <- FDmetrics_Spp.habitat %>% 
  mutate(percent_change = ((value - 1100) / 1100) * 100)
FDmetrics_FRic.habitat <- FDmetrics_FRic.habitat %>% 
  mutate(percent_change = ((value - 1) / 1) * 100)

## Minus vertical
FDindices_taxon.vertical_long <- melt(FDmetrics_taxonvar.vertical, id.vars= "Scenario")

FDmetrics_Spp.vertical <- FDindices_taxon.vertical_long %>% 
  filter(variable == "sp_richn")
FDmetrics_FRic.vertical <- FDindices_taxon.vertical_long %>% 
  filter(variable == "fric")

FDmetrics_Spp.vertical <- FDmetrics_Spp.vertical %>% 
  mutate(percent_change = ((value - 1100) / 1100) * 100)
FDmetrics_FRic.vertical <- FDmetrics_FRic.vertical %>% 
  mutate(percent_change = ((value - 1) / 1) * 100)

## Minus terrest
FDindices_taxon.terrest_long <- melt(FDmetrics_taxonvar.terrest, id.vars= "Scenario")

FDmetrics_Spp.terrest <- FDindices_taxon.terrest_long %>% 
  filter(variable == "sp_richn")
FDmetrics_FRic.terrest <- FDindices_taxon.terrest_long %>% 
  filter(variable == "fric")

FDmetrics_Spp.terrest <- FDmetrics_Spp.terrest %>% 
  mutate(percent_change = ((value - 1100) / 1100) * 100)
FDmetrics_FRic.terrest <- FDmetrics_FRic.terrest %>% 
  mutate(percent_change = ((value - 1) / 1) * 100)

## Minus thermo
FDindices_taxon.thermo_long <- melt(FDmetrics_taxonvar.thermo, id.vars= "Scenario")

FDmetrics_Spp.thermo <- FDindices_taxon.thermo_long %>% 
  filter(variable == "sp_richn")
FDmetrics_FRic.thermo <- FDindices_taxon.thermo_long %>% 
  filter(variable == "fric")

FDmetrics_Spp.thermo <- FDmetrics_Spp.thermo %>% 
  mutate(percent_change = ((value - 1100) / 1100) * 100)
FDmetrics_FRic.thermo <- FDmetrics_FRic.thermo %>% 
  mutate(percent_change = ((value - 1) / 1) * 100)

## Minus feeding
FDindices_taxon.feeding_long <- melt(FDmetrics_taxonvar.feeding, id.vars= "Scenario")

FDmetrics_Spp.feeding <- FDindices_taxon.feeding_long %>% 
  filter(variable == "sp_richn")
FDmetrics_FRic.feeding <- FDindices_taxon.feeding_long %>% 
  filter(variable == "fric")

FDmetrics_Spp.feeding <- FDmetrics_Spp.feeding %>% 
  mutate(percent_change = ((value - 1100) / 1100) * 100)
FDmetrics_FRic.feeding <- FDmetrics_FRic.feeding %>% 
  mutate(percent_change = ((value - 1) / 1) * 100)

## Minus diet
FDindices_taxon.diet_long <- melt(FDmetrics_taxonvar.diet, id.vars= "Scenario")

FDmetrics_Spp.diet <- FDindices_taxon.diet_long %>% 
  filter(variable == "sp_richn")
FDmetrics_FRic.diet <- FDindices_taxon.diet_long %>% 
  filter(variable == "fric")

FDmetrics_Spp.diet <- FDmetrics_Spp.diet %>% 
  mutate(percent_change = ((value - 1100) / 1100) * 100)
FDmetrics_FRic.diet <- FDmetrics_FRic.diet %>% 
  mutate(percent_change = ((value - 1) / 1) * 100)

## Minus maxTL
FDindices_taxon.maxTL_long <- melt(FDmetrics_taxonvar.TL, id.vars= "Scenario")

FDmetrics_Spp.maxTL <- FDindices_taxon.maxTL_long %>% 
  filter(variable == "sp_richn")
FDmetrics_FRic.maxTL <- FDindices_taxon.maxTL_long %>% 
  filter(variable == "fric")

FDmetrics_Spp.maxTL <- FDmetrics_Spp.maxTL %>% 
  mutate(percent_change = ((value - 1100) / 1100) * 100)
FDmetrics_FRic.maxTL <- FDmetrics_FRic.maxTL %>% 
  mutate(percent_change = ((value - 1) / 1) * 100)

## Minus deepwater
FDindices_taxon.deepwater_long <- melt(FDmetrics_taxonvar.deepwater, id.vars= "Scenario")

FDmetrics_Spp.deepwater <- FDindices_taxon.deepwater_long %>% 
  filter(variable == "sp_richn")
FDmetrics_FRic.deepwater <- FDindices_taxon.deepwater_long %>% 
  filter(variable == "fric")

FDmetrics_Spp.deepwater <- FDmetrics_Spp.deepwater %>% 
  mutate(percent_change = ((value - 1100) / 1100) * 100)
FDmetrics_FRic.deepwater <- FDmetrics_FRic.deepwater %>% 
  mutate(percent_change = ((value - 1) / 1) * 100)

## Minus bodyplan
FDindices_taxon.bodyplan_long <- melt(FDmetrics_taxonvar.bodyplan, id.vars= "Scenario")

FDmetrics_Spp.bodyplan <- FDindices_taxon.bodyplan_long %>% 
  filter(variable == "sp_richn")
FDmetrics_FRic.bodyplan <- FDindices_taxon.bodyplan_long %>% 
  filter(variable == "fric")

FDmetrics_Spp.bodyplan <- FDmetrics_Spp.bodyplan %>% 
  mutate(percent_change = ((value - 1100) / 1100) * 100)
FDmetrics_FRic.bodyplan <- FDmetrics_FRic.bodyplan %>% 
  mutate(percent_change = ((value - 1) / 1) * 100)

## Minus reproduction
FDindices_taxon.reproduction_long <- melt(FDmetrics_taxonvar.reproduction, id.vars= "Scenario")

FDmetrics_Spp.reproduction <- FDindices_taxon.reproduction_long %>% 
  filter(variable == "sp_richn")
FDmetrics_FRic.reproduction <- FDindices_taxon.reproduction_long %>% 
  filter(variable == "fric")

FDmetrics_Spp.reproduction <- FDmetrics_Spp.reproduction %>% 
  mutate(percent_change = ((value - 1100) / 1100) * 100)
FDmetrics_FRic.reproduction <- FDmetrics_FRic.reproduction %>% 
  mutate(percent_change = ((value - 1) / 1) * 100)

# Calculate mean changes for species and FRic per scenario
## Minus habitat
mean_percent_change_Spp.habitat <- FDmetrics_Spp.habitat %>%
  group_by(Scenario) %>%
  summarise(mean_percent_change = mean(percent_change, na.rm = TRUE))

mean_percent_change_FRic.habitat <- FDmetrics_FRic.habitat %>%
  group_by(Scenario) %>%
  summarise(mean_percent_change = mean(percent_change, na.rm = TRUE))

## Minus vertical
mean_percent_change_Spp.vertical <- FDmetrics_Spp.vertical %>%
  group_by(Scenario) %>%
  summarise(mean_percent_change = mean(percent_change, na.rm = TRUE))

mean_percent_change_FRic.vertical <- FDmetrics_FRic.vertical %>%
  group_by(Scenario) %>%
  summarise(mean_percent_change = mean(percent_change, na.rm = TRUE))

## Minus terrest
mean_percent_change_Spp.terrest <- FDmetrics_Spp.terrest %>%
  group_by(Scenario) %>%
  summarise(mean_percent_change = mean(percent_change, na.rm = TRUE))

mean_percent_change_FRic.terrest <- FDmetrics_FRic.terrest %>%
  group_by(Scenario) %>%
  summarise(mean_percent_change = mean(percent_change, na.rm = TRUE))

## Minus thermo
mean_percent_change_Spp.thermo <- FDmetrics_Spp.thermo %>%
  group_by(Scenario) %>%
  summarise(mean_percent_change = mean(percent_change, na.rm = TRUE))

mean_percent_change_FRic.thermo <- FDmetrics_FRic.thermo %>%
  group_by(Scenario) %>%
  summarise(mean_percent_change = mean(percent_change, na.rm = TRUE))

## Minus feeding
mean_percent_change_Spp.feeding <- FDmetrics_Spp.feeding %>%
  group_by(Scenario) %>%
  summarise(mean_percent_change = mean(percent_change, na.rm = TRUE))

mean_percent_change_FRic.feeding <- FDmetrics_FRic.feeding %>%
  group_by(Scenario) %>%
  summarise(mean_percent_change = mean(percent_change, na.rm = TRUE))

## Minus diet
mean_percent_change_Spp.diet <- FDmetrics_Spp.diet %>%
  group_by(Scenario) %>%
  summarise(mean_percent_change = mean(percent_change, na.rm = TRUE))

mean_percent_change_FRic.diet <- FDmetrics_FRic.diet %>%
  group_by(Scenario) %>%
  summarise(mean_percent_change = mean(percent_change, na.rm = TRUE))

## Minus maxTL
mean_percent_change_Spp.maxTL <- FDmetrics_Spp.maxTL %>%
  group_by(Scenario) %>%
  summarise(mean_percent_change = mean(percent_change, na.rm = TRUE))

mean_percent_change_FRic.maxTL <- FDmetrics_FRic.maxTL %>%
  group_by(Scenario) %>%
  summarise(mean_percent_change = mean(percent_change, na.rm = TRUE))

## Minus deepwater
mean_percent_change_Spp.deepwater <- FDmetrics_Spp.deepwater %>%
  group_by(Scenario) %>%
  summarise(mean_percent_change = mean(percent_change, na.rm = TRUE))

mean_percent_change_FRic.deepwater <- FDmetrics_FRic.deepwater %>%
  group_by(Scenario) %>%
  summarise(mean_percent_change = mean(percent_change, na.rm = TRUE))

## Minus bodyplan
mean_percent_change_Spp.bodyplan <- FDmetrics_Spp.bodyplan %>%
  group_by(Scenario) %>%
  summarise(mean_percent_change = mean(percent_change, na.rm = TRUE))

mean_percent_change_FRic.bodyplan <- FDmetrics_FRic.bodyplan %>%
  group_by(Scenario) %>%
  summarise(mean_percent_change = mean(percent_change, na.rm = TRUE))

## Minus reproduction
mean_percent_change_Spp.reproduction <- FDmetrics_Spp.reproduction %>%
  group_by(Scenario) %>%
  summarise(mean_percent_change = mean(percent_change, na.rm = TRUE))

mean_percent_change_FRic.reproduction <- FDmetrics_FRic.reproduction %>%
  group_by(Scenario) %>%
  summarise(mean_percent_change = mean(percent_change, na.rm = TRUE))

# Plot results
## Filter only 2100 scenario (assuming "Future" = 2100)
### Species richness
FDmetrics_Spp_2100.habitat <- FDmetrics_Spp.habitat %>% filter(Scenario == "Future")
FDmetrics_Spp_2100.vertical <- FDmetrics_Spp.vertical %>% filter(Scenario == "Future")
FDmetrics_Spp_2100.terrest <- FDmetrics_Spp.terrest %>% filter(Scenario == "Future")
FDmetrics_Spp_2100.thermo <- FDmetrics_Spp.thermo %>% filter(Scenario == "Future")
FDmetrics_Spp_2100.feeding <- FDmetrics_Spp.feeding %>% filter(Scenario == "Future")
FDmetrics_Spp_2100.diet <- FDmetrics_Spp.diet %>% filter(Scenario == "Future")
FDmetrics_Spp_2100.maxTL <- FDmetrics_Spp.maxTL %>% filter(Scenario == "Future")
FDmetrics_Spp_2100.deepwater <- FDmetrics_Spp.deepwater %>% filter(Scenario == "Future")
FDmetrics_Spp_2100.bodyplan <- FDmetrics_Spp.bodyplan %>% filter(Scenario == "Future")
FDmetrics_Spp_2100.reproduction <- FDmetrics_Spp.reproduction %>% filter(Scenario == "Future")

FDmetrics_Spp_2100_all <- bind_rows(
  FDmetrics_Spp_2100.habitat      %>% mutate(trait_removed = "Habitat"),
  FDmetrics_Spp_2100.vertical     %>% mutate(trait_removed = "Vertical"),
  FDmetrics_Spp_2100.terrest      %>% mutate(trait_removed = "Terrestrial"),
  FDmetrics_Spp_2100.thermo       %>% mutate(trait_removed = "Thermal"),
  FDmetrics_Spp_2100.feeding      %>% mutate(trait_removed = "Feeding"),
  FDmetrics_Spp_2100.diet         %>% mutate(trait_removed = "Diet"),
  FDmetrics_Spp_2100.maxTL        %>% mutate(trait_removed = "Max TL"),
  FDmetrics_Spp_2100.deepwater    %>% mutate(trait_removed = "Deepwater"),
  FDmetrics_Spp_2100.bodyplan     %>% mutate(trait_removed = "Body plan"),
  FDmetrics_Spp_2100.reproduction %>% mutate(trait_removed = "Reproduction")
)

### FRic
FDmetrics_FRic_2100.habitat <- FDmetrics_FRic.habitat %>% filter(Scenario == "Future")
FDmetrics_FRic_2100.vertical <- FDmetrics_FRic.vertical %>% filter(Scenario == "Future")
FDmetrics_FRic_2100.terrest <- FDmetrics_FRic.terrest %>% filter(Scenario == "Future")
FDmetrics_FRic_2100.thermo <- FDmetrics_FRic.thermo %>% filter(Scenario == "Future")
FDmetrics_FRic_2100.feeding <- FDmetrics_FRic.feeding %>% filter(Scenario == "Future")
FDmetrics_FRic_2100.diet <- FDmetrics_FRic.diet %>% filter(Scenario == "Future")
FDmetrics_FRic_2100.maxTL <- FDmetrics_FRic.maxTL %>% filter(Scenario == "Future")
FDmetrics_FRic_2100.deepwater <- FDmetrics_FRic.deepwater %>% filter(Scenario == "Future")
FDmetrics_FRic_2100.bodyplan <- FDmetrics_FRic.bodyplan %>% filter(Scenario == "Future")
FDmetrics_FRic_2100.reproduction <- FDmetrics_FRic.reproduction %>% filter(Scenario == "Future")

FDmetrics_FRic_2100_all <- bind_rows(
  FDmetrics_FRic_2100.habitat      %>% mutate(trait_removed = "Habitat"),
  FDmetrics_FRic_2100.vertical     %>% mutate(trait_removed = "Vertical"),
  FDmetrics_FRic_2100.terrest      %>% mutate(trait_removed = "Terrestriality"),
  FDmetrics_FRic_2100.thermo       %>% mutate(trait_removed = "Thermo"),
  FDmetrics_FRic_2100.feeding      %>% mutate(trait_removed = "Feeding"),
  FDmetrics_FRic_2100.diet         %>% mutate(trait_removed = "Diet"),
  FDmetrics_FRic_2100.maxTL        %>% mutate(trait_removed = "Max TL"),
  FDmetrics_FRic_2100.deepwater    %>% mutate(trait_removed = "Deepwater"),
  FDmetrics_FRic_2100.bodyplan     %>% mutate(trait_removed = "Body plan"),
  FDmetrics_FRic_2100.reproduction %>% mutate(trait_removed = "Reproduction")
) %>%
  mutate(
    trait_removed = factor(
      trait_removed,
      levels = c(
        "Habitat", "Vertical", "Terrestriality", "Thermo", "Feeding", "Diet",
        "Max TL", "Deepwater", "Body plan", "Reproduction"))
  )

## Plot species richness changes - should match
Spp_sensitivity <- ggplot(
  FDmetrics_Spp_2100_all,
  aes(x = trait_removed, y = percent_change)
) +
  geom_boxplot(fill = "#ca4678", 
               width = 0.7, outlier.size = 0.8) +
  labs(
    x = "Trait removed",
    y = "% change in species richness",
    title = "Species richness (SR) sensitivity tests"
  ) +
  theme_bw() +
  theme(
    axis.text.x = element_text(angle = 45, hjust = 1, size = 8),
    axis.text.y = element_text(size = 8),
    axis.title = element_text(size = 10),
    panel.grid = element_blank(),
    panel.border = element_rect(fill = NA),
    plot.title = element_text(hjust = 0.5, size = 10)
  )

## Plot FRic changes
FRic_sensitivity <- ggplot(
  FDmetrics_FRic_2100_all,
  aes(x = trait_removed, y = percent_change, fill = trait_removed)
) +
  geom_boxplot(
    width = 0.7,
    outlier.size = 0.8
  ) +
  scale_fill_manual(
    values = c(
      "Habitat"         = "#FED439FF",
      "Vertical"        = "#709AE1FF",
      "Terrestriality"  = "#8A9197FF",
      "Thermo"          = "#D2AF81FF",
      "Feeding"         = "#FD7446FF",
      "Diet"            = "#D5E4A2FF",
      "Max TL"          = "#197EC0FF",
      "Deepwater"       = "#F05C3BFF",
      "Body plan"       = "#46732EFF",
      "Reproduction"    = "#71D0F5FF"
    )
  ) +
  labs(
    x = "Trait removed",
    y = "FRic % change",
    title = ""
  ) +
  theme_bw() +
  theme(
    axis.text.x = element_text(angle = 45, hjust = 1, size = 8),
    axis.text.y = element_text(size = 8),
    axis.title  = element_text(size = 10),
    panel.grid  = element_blank(),
    panel.border = element_rect(fill = NA),
    plot.title  = element_text(hjust = 0.5, size = 10)
  ) +
  guides(fill = "none")

