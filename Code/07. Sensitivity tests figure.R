###############################################################################################################
# 07. Sensitivity tests figure
## This R code plots future FRic changes based on sensitivity tests
## it produces Figure S3 and Table S3
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
load("C:/Users/Jack Cooper/Documents/iucn/iucn_sim/Robustness tests/Minus_habitat/CritE_full_buffer_minushabitat.RData")
load("C:/Users/Jack Cooper/Documents/iucn/iucn_sim/Robustness tests/Minus_habitat/CritE_full_long_buffer_minushabitat.RData")
load("C:/Users/Jack Cooper/Documents/iucn/iucn_sim/Robustness tests/Minus_vertical/CritE_full_buffer_minusVP.RData")
load("C:/Users/Jack Cooper/Documents/iucn/iucn_sim/Robustness tests/Minus_vertical/CritE_full_long_buffer_minusVP.RData")
load("C:/Users/Jack Cooper/Documents/iucn/iucn_sim/Robustness tests/Minus_terrest/CritE_full_buffer_minusterrest.RData")
load("C:/Users/Jack Cooper/Documents/iucn/iucn_sim/Robustness tests/Minus_terrest/CritE_full_long_buffer_minusterrest.RData")
load("C:/Users/Jack Cooper/Documents/iucn/iucn_sim/Robustness tests/Minus_thermo/CritE_full_buffer_minusthermo.RData")
load("C:/Users/Jack Cooper/Documents/iucn/iucn_sim/Robustness tests/Minus_thermo/CritE_full_long_buffer_minusthermo.RData")
load("C:/Users/Jack Cooper/Documents/iucn/iucn_sim/Robustness tests/Minus_FM/CritE_full_buffer_minusFM.RData")
load("C:/Users/Jack Cooper/Documents/iucn/iucn_sim/Robustness tests/Minus_FM/CritE_full_long_buffer_minusFM.RData")
load("C:/Users/Jack Cooper/Documents/iucn/iucn_sim/Robustness tests/Minus_diet/CritE_full_buffer_minusdiet.RData")
load("C:/Users/Jack Cooper/Documents/iucn/iucn_sim/Robustness tests/Minus_diet/CritE_full_long_buffer_minusdiet.RData")
load("C:/Users/Jack Cooper/Documents/iucn/iucn_sim/Robustness tests/Minus_maxTL/CritE_full_buffer_minusTL.RData")
load("C:/Users/Jack Cooper/Documents/iucn/iucn_sim/Robustness tests/Minus_maxTL/CritE_full_long_buffer_minusTL.RData")
load("C:/Users/Jack Cooper/Documents/iucn/iucn_sim/Robustness tests/Minus_deepwater/CritE_full_buffer_minusdeepwater.RData")
load("C:/Users/Jack Cooper/Documents/iucn/iucn_sim/Robustness tests/Minus_deepwater/CritE_full_long_buffer_minusdeepwater.RData")
load("C:/Users/Jack Cooper/Documents/iucn/iucn_sim/Robustness tests/Minus_bodyplan/CritE_full_buffer_minusbody.RData")
load("C:/Users/Jack Cooper/Documents/iucn/iucn_sim/Robustness tests/Minus_bodyplan/CritE_full_long_buffer_minusbody.RData")
load("C:/Users/Jack Cooper/Documents/iucn/iucn_sim/Robustness tests/Minus_reproduction/CritE_full_buffer_minusreproduction.RData")
load("C:/Users/Jack Cooper/Documents/iucn/iucn_sim/Robustness tests/Minus_reproduction/CritE_full_long_buffer_minusreproduction.RData")

# Calculate % changes for species and FRic
## Minus habitat
FDindices_taxon.habitat_long <- reshape2::melt(FDmetrics_taxonvar.habitat, id.vars= "Scenario")

FDmetrics_Spp.habitat <- FDindices_taxon.habitat_long %>% 
  filter(variable == "sp_richn")
FDmetrics_FRic.habitat <- FDindices_taxon.habitat_long %>% 
  filter(variable == "fric")

FDmetrics_Spp.habitat <- FDmetrics_Spp.habitat %>% 
  mutate(percent_change = ((value - 1100) / 1100) * 100)
FDmetrics_FRic.habitat <- FDmetrics_FRic.habitat %>% 
  mutate(percent_change = ((value - 1) / 1) * 100)

## Minus vertical
FDindices_taxon.vertical_long <- reshape2::melt(FDmetrics_taxonvar.vertical, id.vars= "Scenario")

FDmetrics_Spp.vertical <- FDindices_taxon.vertical_long %>% 
  filter(variable == "sp_richn")
FDmetrics_FRic.vertical <- FDindices_taxon.vertical_long %>% 
  filter(variable == "fric")

FDmetrics_Spp.vertical <- FDmetrics_Spp.vertical %>% 
  mutate(percent_change = ((value - 1100) / 1100) * 100)
FDmetrics_FRic.vertical <- FDmetrics_FRic.vertical %>% 
  mutate(percent_change = ((value - 1) / 1) * 100)

## Minus terrest
FDindices_taxon.terrest_long <- reshape2::melt(FDmetrics_taxonvar.terrest, id.vars= "Scenario")

FDmetrics_Spp.terrest <- FDindices_taxon.terrest_long %>% 
  filter(variable == "sp_richn")
FDmetrics_FRic.terrest <- FDindices_taxon.terrest_long %>% 
  filter(variable == "fric")

FDmetrics_Spp.terrest <- FDmetrics_Spp.terrest %>% 
  mutate(percent_change = ((value - 1100) / 1100) * 100)
FDmetrics_FRic.terrest <- FDmetrics_FRic.terrest %>% 
  mutate(percent_change = ((value - 1) / 1) * 100)

## Minus thermo
FDindices_taxon.thermo_long <- reshape2::melt(FDmetrics_taxonvar.thermo, id.vars= "Scenario")

FDmetrics_Spp.thermo <- FDindices_taxon.thermo_long %>% 
  filter(variable == "sp_richn")
FDmetrics_FRic.thermo <- FDindices_taxon.thermo_long %>% 
  filter(variable == "fric")

FDmetrics_Spp.thermo <- FDmetrics_Spp.thermo %>% 
  mutate(percent_change = ((value - 1100) / 1100) * 100)
FDmetrics_FRic.thermo <- FDmetrics_FRic.thermo %>% 
  mutate(percent_change = ((value - 1) / 1) * 100)

## Minus feeding
FDindices_taxon.feeding_long <- reshape2::melt(FDmetrics_taxonvar.feeding, id.vars= "Scenario")

FDmetrics_Spp.feeding <- FDindices_taxon.feeding_long %>% 
  filter(variable == "sp_richn")
FDmetrics_FRic.feeding <- FDindices_taxon.feeding_long %>% 
  filter(variable == "fric")

FDmetrics_Spp.feeding <- FDmetrics_Spp.feeding %>% 
  mutate(percent_change = ((value - 1100) / 1100) * 100)
FDmetrics_FRic.feeding <- FDmetrics_FRic.feeding %>% 
  mutate(percent_change = ((value - 1) / 1) * 100)

## Minus diet
FDindices_taxon.diet_long <- reshape2::melt(FDmetrics_taxonvar.diet, id.vars= "Scenario")

FDmetrics_Spp.diet <- FDindices_taxon.diet_long %>% 
  filter(variable == "sp_richn")
FDmetrics_FRic.diet <- FDindices_taxon.diet_long %>% 
  filter(variable == "fric")

FDmetrics_Spp.diet <- FDmetrics_Spp.diet %>% 
  mutate(percent_change = ((value - 1100) / 1100) * 100)
FDmetrics_FRic.diet <- FDmetrics_FRic.diet %>% 
  mutate(percent_change = ((value - 1) / 1) * 100)

## Minus maxTL
FDindices_taxon.maxTL_long <- reshape2::melt(FDmetrics_taxonvar.TL, id.vars= "Scenario")

FDmetrics_Spp.maxTL <- FDindices_taxon.maxTL_long %>% 
  filter(variable == "sp_richn")
FDmetrics_FRic.maxTL <- FDindices_taxon.maxTL_long %>% 
  filter(variable == "fric")

FDmetrics_Spp.maxTL <- FDmetrics_Spp.maxTL %>% 
  mutate(percent_change = ((value - 1100) / 1100) * 100)
FDmetrics_FRic.maxTL <- FDmetrics_FRic.maxTL %>% 
  mutate(percent_change = ((value - 1) / 1) * 100)

## Minus deepwater
FDindices_taxon.deepwater_long <- reshape2::melt(FDmetrics_taxonvar.deepwater, id.vars= "Scenario")

FDmetrics_Spp.deepwater <- FDindices_taxon.deepwater_long %>% 
  filter(variable == "sp_richn")
FDmetrics_FRic.deepwater <- FDindices_taxon.deepwater_long %>% 
  filter(variable == "fric")

FDmetrics_Spp.deepwater <- FDmetrics_Spp.deepwater %>% 
  mutate(percent_change = ((value - 1100) / 1100) * 100)
FDmetrics_FRic.deepwater <- FDmetrics_FRic.deepwater %>% 
  mutate(percent_change = ((value - 1) / 1) * 100)

## Minus bodyplan
FDindices_taxon.bodyplan_long <- reshape2::melt(FDmetrics_taxonvar.bodyplan, id.vars= "Scenario")

FDmetrics_Spp.bodyplan <- FDindices_taxon.bodyplan_long %>% 
  filter(variable == "sp_richn")
FDmetrics_FRic.bodyplan <- FDindices_taxon.bodyplan_long %>% 
  filter(variable == "fric")

FDmetrics_Spp.bodyplan <- FDmetrics_Spp.bodyplan %>% 
  mutate(percent_change = ((value - 1100) / 1100) * 100)
FDmetrics_FRic.bodyplan <- FDmetrics_FRic.bodyplan %>% 
  mutate(percent_change = ((value - 1) / 1) * 100)

## Minus reproduction
FDindices_taxon.reproduction_long <- reshape2::melt(FDmetrics_taxonvar.reproduction, id.vars= "Scenario")

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

## Plot species richness changes - a test to show all 10 should match
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

## Plot FRic changes - produces Figure S3
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

# Secondary analysis: trait removal tests against primary analysis
# Load primary analysis
load(file = "C:/Users/Jack Cooper/Documents/iucn/iucn_sim/Crit E analyses/CritE_full_buffer.RData")
load(file = "C:/Users/Jack Cooper/Documents/iucn/iucn_sim/Crit E analyses/CritE_full_long_buffer.RData")

FDindices_taxon.CritE_long <- reshape2::melt(FDmetrics_taxonvar.CritE, id.vars= "Scenario")

FDmetrics_Spp.CritE <- FDindices_taxon.CritE_long %>% 
  filter(variable == "sp_richn")
FDmetrics_FRic.CritE <- FDindices_taxon.CritE_long %>% 
  filter(variable == "fric")

FDmetrics_Spp_2100 <- FDmetrics_Spp.CritE %>% filter(Scenario == "Future")
FRicmetrics_empirical_2100 <- FDmetrics_FRic.CritE %>% filter(Scenario == "Future")

FRicmetrics_empirical_2100 <- FRicmetrics_empirical_2100 %>%
  mutate(percent_change = ((value - 1) / 1) * 100)

# Compare primary analysis to all 10 senstivity tests - pairwise wilcox
FRic_primary <- FRicmetrics_empirical_2100 %>%
  mutate(sim_id = row_number()) %>%
  select(sim_id, percent_change) %>%
  rename(All_traits = percent_change)

FRic_sensitivity_data <- FDmetrics_FRic_2100_all %>%
  group_by(trait_removed) %>%
  mutate(sim_id = row_number()) %>%
  ungroup() %>%
  select(sim_id, trait_removed, percent_change) %>%
  pivot_wider(
    names_from = trait_removed,
    values_from = percent_change
  )

nrow(FRic_primary)
nrow(FRic_sensitivity_data)

FRic_comparison <- FRic_sensitivity_data %>%
  inner_join(FRic_primary, by = "sim_id")

# Identify the ten trait-exclusion scenarios
trait_names <- setdiff(
  names(FRic_comparison),
  c("sim_id", "All_traits")
)

# Run paired comparisons
FRic_sensitivity_results <- map_dfr(trait_names, function(trait) {
  
  # Extract paired values
  x <- FRic_comparison$All_traits
  y <- FRic_comparison[[trait]]
  
  # Remove missing values
  valid <- is.finite(x) & is.finite(y)
  x <- x[valid]
  y <- y[valid]
  
  # Calculate paired differences: primary minus sensitivity
  diffs <- x - y
  diffs_nonzero <- diffs[diffs != 0]
  N <- length(diffs_nonzero)
  
  # Wilcoxon signed-rank test
  wt <- wilcox.test(
    x, y,
    paired = TRUE,
    alternative = "two.sided",
    exact = FALSE,
    correct = TRUE
  )
  
  V <- as.numeric(wt$statistic)
  
  # Standardised Z statistic (tie-adjusted)
  ranks <- rank(abs(diffs_nonzero), ties.method = "average")
  mean_V <- sum(ranks) / 2
  sd_V <- sqrt(sum(ranks^2) / 4)
  
  Z <- (V - mean_V - sign(V - mean_V) * 0.5) / sd_V
  
  # Effect size
  r <- Z / sqrt(N)
  
  # Results table
  tibble(
    Trait_removed = trait,
    Common_pairs = length(diffs),
    Effective_n = N,
    Mean_paired_difference = mean(diffs),
    Median_paired_difference = median(diffs),
    IQR = IQR(diffs),
    Wilcoxon_V = V,
    Z_statistic = Z,
    Effect_size_r = r,
    P_value = wt$p.value
  )
}) %>%
  mutate(
    P_adjusted = p.adjust(P_value, method = "holm")
  )

# Display results
print(FRic_sensitivity_results, n = Inf)