###############################################################################################################
# 06. Figure 1
## This R code plots future functional diversity (FD) based on our simulations
## it produces Figure 1A & D
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

## Calculate % changes in species & functional diversity - Crit E EX mode
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
save(Spp_CritE_2100, file = "~/Spp_2100.RData")

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
save(FRic_CritE_2100, file = "~/FRic_2100.RData")

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

FRic_results <- perform_wilcoxon_extended(
  FDmetrics_FRic.CritE,
  FRicmetrics_null_Taxonvar.CritE,
  scenario_name = "Future",
  alternative = "two.sided"
)


print(FRic_results)
