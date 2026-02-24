###############################################################################################################
# 08. Overfitting tests
## This R code plots future functional diversity (FD) based on fewer simulations
## it produces Figure S2 and Table S2
###############################################################################################################

## Import packages
library(tidyverse)
library(tibble)
library(reshape2)
library(cowplot)
library(ggsci)
library(RColorBrewer)
library(purrr)

# Load boxplots and null model
load(file = "~/CritE_full_buffer.RData")
load(file = "~/CritE_full_long_buffer.RData")
load(file = "~/Null_wide_CritE.RData")
load(file = "~/Null_long_CritE.RData")

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

# Subset the first 1000 
# Empirical
FDmetrics_Spp.1000 <- FDmetrics_Spp.CritE %>%
  slice(1:2000)
FDmetrics_FRic.1000 <- FDmetrics_FRic.CritE %>%
  slice(1:2000)

# Calculate mean changes per metric in this subset
mean_percent_change_Spp.1000 <- FDmetrics_Spp.1000 %>%
  group_by(Scenario) %>%
  summarise(mean_percent_change = mean(percent_change, na.rm = TRUE))

mean_percent_change_FRic.1000 <- FDmetrics_FRic.1000 %>%
  group_by(Scenario) %>%
  summarise(mean_percent_change = mean(percent_change, na.rm = TRUE))

# Null
SPPmetrics_null_Taxonvar.1000 <- FDmetrics_null_long_taxonvar.CritE %>% 
  filter(variable == "sp_richn") %>% 
  slice(1:2000)
FRicmetrics_null_Taxonvar.1000 <- FDmetrics_null_long_taxonvar.CritE %>% 
  filter(variable == "fric") %>% 
  slice(1:2000)

SPPmetrics_null_Taxonvar.1000 <- SPPmetrics_null_Taxonvar.1000 %>% 
  mutate(percent_change = ((value - 1100) / 1100) * 100)
FRicmetrics_null_Taxonvar.1000 <- FRicmetrics_null_Taxonvar.1000 %>% 
  mutate(percent_change = ((value - 1) / 1) * 100)

## Plot results
# Filter only future
FDmetrics_Spp_1000 <- FDmetrics_Spp.1000 %>% filter(Scenario == "Future")
FRicmetrics_empirical_1000 <- FDmetrics_FRic.1000 %>% filter(Scenario == "Future")
FRicmetrics_null_1000 <- FRicmetrics_null_Taxonvar.1000 %>% filter(Scenario == "Future")

# Plot 1: Species Richness (SR) - 2100 only
Spp_1000 <- ggplot(data = FDmetrics_Spp_1000, aes(x = Scenario, y = percent_change)) +
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

# Plot 2: Functional Richness (FRic) - 2100 only
FRic_1000 <- ggplot(data = FRicmetrics_null_1000, aes(x = Scenario, y = percent_change)) +
  geom_violin(trim = TRUE, fill = 'grey', color = "grey", alpha = 0.25, scale = "width") +
  geom_boxplot(data = FRicmetrics_empirical_1000, aes(x = Scenario, y = percent_change), 
               fill = "#f78c33") +  
  labs(x = "", y = "FRic % change") +
  ggtitle("") +
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

# Wilcoxon signed-ranked tests
perform_wilcoxon_extended <- function(empirical_data, null_data, scenario_name, 
                                      alternative = "two.sided",
                                      correct = TRUE) {
  
  # Extract paired values
  empirical_values <- empirical_data %>%
    filter(Scenario == scenario_name) %>%
    pull(percent_change)
  
  null_values <- null_data %>%
    filter(Scenario == scenario_name) %>%
    pull(percent_change)
  
  # Ensure equal length
  stopifnot(length(empirical_values) == length(null_values))
  
  # Compute paired differences
  diffs <- empirical_values - null_values
  
  # Remove zero differences
  diffs_nonzero <- diffs[diffs != 0]
  N <- length(diffs_nonzero)
  
  # Run Wilcoxon signed-rank test
  wt <- wilcox.test(empirical_values,
                    null_values,
                    paired = TRUE,
                    alternative = alternative,
                    exact = FALSE,
                    correct = correct)
  
  V <- as.numeric(wt$statistic)
  
  # Expected mean and SD under null
  mean_V <- N * (N + 1) / 4
  sd_V <- sqrt(N * (N + 1) * (2 * N + 1) / 24)
  
  # Continuity correction
  cc <- if (correct) 0.5 else 0
  
  if (alternative == "two.sided") {
    Z <- (V - mean_V - sign(V - mean_V) * cc) / sd_V
  } else if (alternative == "greater") {
    Z <- (V - mean_V - cc) / sd_V
  } else {
    Z <- (V - mean_V + cc) / sd_V
  }
  
  # Effect size
  r <- Z / sqrt(N)
  
  # Median paired difference
  median_diff <- median(diffs)
  
  # Create results table
  results_table <- tibble(
    Scenario = scenario_name,
    Wilcoxon_V = V,
    Sample_size_nonzero_pairs = N,
    Test_tail = alternative,
    P_value = wt$p.value,
    P_value_type = "Asymptotic (normal approximation)",
    Z_statistic = Z,
    Effect_size_r = r,
    Median_paired_difference = median_diff
  )
  
  # Convert to character BEFORE pivoting to avoid type conflict
  results_table %>%
    mutate(across(everything(), as.character)) %>%
    pivot_longer(
      cols = everything(),
      names_to = "Statistic",
      values_to = "result"
    )
}

# Tabulate outputs
cat("\nExtended Wilcoxon results for FRic:\n")

FRic_results_1000 <- perform_wilcoxon_extended(
  FDmetrics_FRic.1000,
  FRicmetrics_null_Taxonvar.1000,
  scenario_name = "Future",
  alternative = "two.sided"
)

print(FRic_results_1000)

### Test 2: 2500 simulations
# Subset  
# Empirical
FDmetrics_Spp.2500 <- FDmetrics_Spp.CritE %>%
  slice(1:5000)
FDmetrics_FRic.2500 <- FDmetrics_FRic.CritE %>%
  slice(1:5000)

# Calculate mean changes per metric in this subset
mean_percent_change_Spp_2500 <- FDmetrics_Spp.2500 %>%
  group_by(Scenario) %>%
  summarise(mean_percent_change = mean(percent_change, na.rm = TRUE))

mean_percent_change_FRic_2500 <- FDmetrics_FRic.2500 %>%
  group_by(Scenario) %>%
  summarise(mean_percent_change = mean(percent_change, na.rm = TRUE))

# Null
SPPmetrics_null_Taxonvar.2500 <- FDmetrics_null_long_taxonvar.CritE %>% 
  filter(variable == "sp_richn") %>% 
  slice(1:5000)
FRicmetrics_null_Taxonvar.2500 <- FDmetrics_null_long_taxonvar.CritE %>% 
  filter(variable == "fric") %>% 
  slice(1:5000)

SPPmetrics_null_Taxonvar.2500 <- SPPmetrics_null_Taxonvar.2500 %>% 
  mutate(percent_change = ((value - 1100) / 1100) * 100)
FRicmetrics_null_Taxonvar.2500 <- FRicmetrics_null_Taxonvar.2500 %>% 
  mutate(percent_change = ((value - 1) / 1) * 100)

## Plot results
# Filter only future
FDmetrics_Spp_2500 <- FDmetrics_Spp.2500 %>% filter(Scenario == "Future")
FRicmetrics_empirical_2500 <- FDmetrics_FRic.2500 %>% filter(Scenario == "Future")
FRicmetrics_null_2500 <- FRicmetrics_null_Taxonvar.2500 %>% filter(Scenario == "Future")

# Plot 1: Species Richness (SR) - 2100 only
Spp_2500 <- ggplot(data = FDmetrics_Spp_2500, aes(x = Scenario, y = percent_change)) +
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

# Plot 2: Functional Richness (FRic) - 2100 only
FRic_2500 <- ggplot(data = FRicmetrics_null_2500, aes(x = Scenario, y = percent_change)) +
  geom_violin(trim = TRUE, fill = 'grey', color = "grey", alpha = 0.25, scale = "width") +
  geom_boxplot(data = FRicmetrics_empirical_2500, aes(x = Scenario, y = percent_change), 
               fill = "#f78c33") +  
  labs(x = "", y = "FRic % change") +
  ggtitle("") +
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

# Tabulate outputs
cat("\nExtended Wilcoxon results for FRic:\n")

FRic_results_2500 <- perform_wilcoxon_extended(
  FDmetrics_FRic.2500,
  FRicmetrics_null_Taxonvar.2500,
  scenario_name = "Future",
  alternative = "two.sided"
)

print(FRic_results_2500)

# Test 3: 5000 simulations
# Subset  
# Empirical
FDmetrics_Spp.5000 <- FDmetrics_Spp.CritE %>%
  slice(1:10000)
FDmetrics_FRic.5000 <- FDmetrics_FRic.CritE %>%
  slice(1:10000)

# Calculate mean changes per metric in this subset
mean_percent_change_Spp_5000 <- FDmetrics_Spp.5000 %>%
  group_by(Scenario) %>%
  summarise(mean_percent_change = mean(percent_change, na.rm = TRUE))

mean_percent_change_FRic_5000 <- FDmetrics_FRic.5000 %>%
  group_by(Scenario) %>%
  summarise(mean_percent_change = mean(percent_change, na.rm = TRUE))

# Null
SPPmetrics_null_Taxonvar.5000 <- FDmetrics_null_long_taxonvar.CritE %>% 
  filter(variable == "sp_richn") %>% 
  slice(1:10000)
FRicmetrics_null_Taxonvar.5000 <- FDmetrics_null_long_taxonvar.CritE %>% 
  filter(variable == "fric") %>% 
  slice(1:10000)

SPPmetrics_null_Taxonvar.5000 <- SPPmetrics_null_Taxonvar.5000 %>% 
  mutate(percent_change = ((value - 1100) / 1100) * 100)
FRicmetrics_null_Taxonvar.5000 <- FRicmetrics_null_Taxonvar.5000 %>% 
  mutate(percent_change = ((value - 1) / 1) * 100)

## Plot results
# Filter only future
FDmetrics_Spp_5000 <- FDmetrics_Spp.5000 %>% filter(Scenario == "Future")
FRicmetrics_empirical_5000 <- FDmetrics_FRic.5000 %>% filter(Scenario == "Future")
FRicmetrics_null_5000 <- FRicmetrics_null_Taxonvar.5000 %>% filter(Scenario == "Future")

# Plot 1: Species Richness (SR) - 2100 only
Spp_5000 <- ggplot(data = FDmetrics_Spp_5000, aes(x = Scenario, y = percent_change)) +
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

# Plot 2: Functional Richness (FRic) - 2100 only
FRic_5000 <- ggplot(data = FRicmetrics_null_5000, aes(x = Scenario, y = percent_change)) +
  geom_violin(trim = TRUE, fill = 'grey', color = "grey", alpha = 0.25, scale = "width") +
  geom_boxplot(data = FRicmetrics_empirical_5000, aes(x = Scenario, y = percent_change), 
               fill = "#f78c33") +  
  labs(x = "", y = "FRic % change") +
  ggtitle("") +
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

# Tabulate outputs
cat("\nExtended Wilcoxon results for FRic:\n")

FRic_results_5000 <- perform_wilcoxon_extended(
  FDmetrics_FRic.5000,
  FRicmetrics_null_Taxonvar.5000,
  scenario_name = "Future",
  alternative = "two.sided"
)

print(FRic_results_5000)

# Form database of all tests together
FRic_empirical_all <- bind_rows(
  FRicmetrics_empirical_1000 %>% mutate(Simulations = "1000"),
  FRicmetrics_empirical_2500 %>% mutate(Simulations = "2500"),
  FRicmetrics_empirical_5000 %>% mutate(Simulations = "5000")
)

FRic_null_all <- bind_rows(
  FRicmetrics_null_1000 %>% mutate(Simulations = "1000"),
  FRicmetrics_null_2500 %>% mutate(Simulations = "2500"),
  FRicmetrics_null_5000 %>% mutate(Simulations = "5000")
)

FRic_empirical_all$Simulations <- factor(
  FRic_empirical_all$Simulations,
  levels = c("1000", "2500", "5000")
)

FRic_null_all$Simulations <- factor(
  FRic_null_all$Simulations,
  levels = c("1000", "2500", "5000")
)

# Plot all tests together - produces Figure S2
FRic_combined <- ggplot() +
  
  # Null violins (fixed grey)
  geom_violin(
    data = FRic_null_all,
    aes(x = Simulations, y = percent_change),
    trim = TRUE,
    fill = "grey",
    color = "grey",
    alpha = 0.25,
    scale = "width"
  ) +
  
  # Empirical boxplots (mapped fill)
  geom_boxplot(
    data = FRic_empirical_all,
    aes(x = Simulations, y = percent_change, fill = Simulations),
    width = 0.2
  ) +
  
  scale_fill_manual(
    values = c(
      "1000" = "#CC0C00FF",
      "2500" = "#FFCD00FF",
      "5000" = "#00B5E2FF"
    )
  ) +
  
  labs(
    x = "Number of simulations",
    y = "FRic % change"
  ) +
  
  theme_bw() +
  theme(
    legend.position = "none",
    axis.text = element_text(size = 8, color = "black"),
    axis.title = element_text(size = 10),
    panel.background = element_rect(fill = "white"),
    panel.grid = element_blank(),
    panel.border = element_rect(fill = "transparent")
  )

