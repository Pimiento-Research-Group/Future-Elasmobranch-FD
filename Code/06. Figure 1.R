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
FDmetrics_FSpe.CritE <- FDindices_taxon.CritE_long %>% 
  filter(variable == "fsp")
FDmetrics_FUn.CritE <- FDindices_taxon.CritE_long %>% 
  filter(variable == "fun")

FDmetrics_Spp.CritE <- FDmetrics_Spp.CritE %>% 
  mutate(percent_change = ((value - 1100) / 1100) * 100)
FDmetrics_FRic.CritE <- FDmetrics_FRic.CritE %>% 
  mutate(percent_change = ((value - 1) / 1) * 100)
FDmetrics_FSpe.CritE <- FDmetrics_FSpe.CritE %>% 
  mutate(percent_change = ((value - 0.2000877) / 0.2000877) * 100)
FDmetrics_FUn.CritE <- FDmetrics_FUn.CritE %>% 
  mutate(percent_change = ((value - 0.007140415) / 0.007140415) * 100)

# Calculate mean changes per scenario per metric
mean_percent_change_Spp <- FDmetrics_Spp.CritE %>%
  group_by(Scenario) %>%
  summarise(mean_percent_change = mean(percent_change, na.rm = TRUE))

mean_percent_change_FRic <- FDmetrics_FRic.CritE %>%
  group_by(Scenario) %>%
  summarise(mean_percent_change = mean(percent_change, na.rm = TRUE))

mean_percent_change_FSpe <- FDmetrics_FSpe.CritE %>%
  group_by(Scenario) %>%
  summarise(mean_percent_change = mean(percent_change, na.rm = TRUE))

mean_percent_change_FUn <- FDmetrics_FUn.CritE %>%
  group_by(Scenario) %>%
  summarise(mean_percent_change = mean(percent_change, na.rm = TRUE))

# Null
SPPmetrics_null_Taxonvar.CritE <- FDmetrics_null_long_taxonvar.CritE %>% 
  filter(variable == "sp_richn")
FRicmetrics_null_Taxonvar.CritE <- FDmetrics_null_long_taxonvar.CritE %>% 
  filter(variable == "fric")
FSpemetrics_null_Taxonvar.CritE <- FDmetrics_null_long_taxonvar.CritE %>% 
  filter(variable == "fsp")
FUnmetrics_null_Taxonvar.CritE <- FDmetrics_null_long_taxonvar.CritE %>% 
  filter(variable == "fun")

SPPmetrics_null_Taxonvar.CritE <- SPPmetrics_null_Taxonvar.CritE %>% 
  mutate(percent_change = ((value - 1100) / 1100) * 100)
FRicmetrics_null_Taxonvar.CritE <- FRicmetrics_null_Taxonvar.CritE %>% 
  mutate(percent_change = ((value - 1) / 1) * 100)
FSpemetrics_null_Taxonvar.CritE <- FSpemetrics_null_Taxonvar.CritE %>% 
  mutate(percent_change = ((value - 0.2000877) / 0.2000877) * 100)
FUnmetrics_null_Taxonvar.CritE <- FUnmetrics_null_Taxonvar.CritE %>% 
  mutate(percent_change = ((value - 0.007140415) / 0.007140415) * 100)

## Plot Crit E EX mode results
# Species richness
Spp_CritE <- ggplot() +
  scale_x_discrete(limits=c("Future", "AT")) +  
  geom_boxplot(data = subset(FDmetrics_Spp.CritE, Scenario!="Present"), 
               aes(x = Scenario, y = percent_change), fill = "#d95f02") +  
  scale_color_manual(values=c("Future"="#000000", "AT"="#000000")) +  
  scale_fill_manual(values=c("Future"="#000000", "AT"="#000000")) +  
  labs(x = "", y = "Species % change") +
  ggtitle("") +
  theme_minimal() +
  theme_bw() +
  theme(legend.position = "none",
        axis.text = element_text(size= 8, color= "black"),
        axis.title= element_text(size= 10), 
        panel.background= element_rect(fill= "white"),
        panel.grid = element_blank(), 
        panel.border = element_rect(fill= "transparent"),
        axis.text.x = element_blank(),
        plot.title = element_text(hjust = 0.5))

# Functional Richness
FRic_CritE <- ggplot(data=subset(FRicmetrics_null_Taxonvar.CritE, Scenario!="Present"), 
                     aes(x= Scenario, y= percent_change)) +
  geom_violin(trim=TRUE, fill='grey', color="grey", alpha = 0.25, scale = "width") +
  scale_x_discrete(limits=c("Future", "AT")) +  
  geom_boxplot(data = subset(FDmetrics_FRic.CritE, Scenario!="Present"), 
               aes(x = Scenario, y = percent_change), fill = "#7570b3") +  
  scale_color_manual(values=c("Future"="#000000", "AT"="#000000")) + 
  scale_fill_manual(values=c("Future"="#000000", "AT"="#000000")) +  
  labs(x = "", y = "FRic % change") +
  ggtitle("") +
  theme_minimal() +
  theme_bw() +
  theme(legend.position = "none",
        axis.text = element_text(size= 8, color= "black"),
        axis.title= element_text(size= 10), 
        panel.background= element_rect(fill= "white"),
        panel.grid = element_blank(), 
        panel.border = element_rect(fill= "transparent"),
        axis.text.x = element_blank(),
        plot.title = element_text(hjust = 0.5))

# Functional Uniqueness 
FUn_CritE <- ggplot(data=subset(FUnmetrics_null_Taxonvar.CritE, Scenario!="Present"), 
                    aes(x= Scenario, y= percent_change)) +
  geom_violin(trim=TRUE, fill='grey', color="grey", alpha = 0.25, scale = "width") +
  scale_x_discrete(limits=c("Future", "AT"), labels=c("2100", "IUCN AT")) +  
  geom_boxplot(data = subset(FDmetrics_FUn.CritE, Scenario!="Present"), 
               aes(x = Scenario, y = percent_change), fill = "#e7298a") +  
  scale_color_manual(values=c("Future"="#000000", "AT"="#000000")) +
  scale_fill_manual(values=c("Future"="#000000", "AT"="#000000")) +  
  labs(x = "Future Scenario", y = "FUn % change") +
  ggtitle("") +
  theme_minimal() +
  theme_bw() +
  theme(legend.position = "none") +
  theme(axis.text = element_text(size= 8, color= "black"),
        axis.title= element_text(size= 10), 
        panel.background= element_rect(fill= "white")) +
  theme(panel.grid = element_blank(), 
        panel.border = element_rect(fill= "transparent"))

# Functional Specialisation 
FSp_CritE <- ggplot(data=subset(FSpemetrics_null_Taxonvar.CritE, Scenario!="Present"), 
                    aes(x= Scenario, y= percent_change)) +
  geom_violin(trim=TRUE, fill='grey', color="grey", alpha = 0.25, scale = "width") +
  scale_x_discrete(limits=c("Future", "AT"), labels=c("2100", "IUCN AT")) +  
  geom_boxplot(data = subset(FDmetrics_FSpe.CritE, Scenario!="Present"), 
               aes(x = Scenario, y = percent_change), fill = "#e6ab02") +  
  scale_color_manual(values=c("Future"="#000000", "AT"="#000000")) +  
  scale_fill_manual(values=c("Future"="#000000", "AT"="#000000")) +  
  labs(x = "Future Scenario", y = "FSp % change") +
  theme_minimal() +
  theme_bw() +
  theme(legend.position = "none") +
  theme(axis.text = element_text(size= 8, color= "black"),
        axis.title= element_text(size= 10), 
        panel.background= element_rect(fill= "white")) +
  theme(panel.grid = element_blank(), 
        panel.border = element_rect(fill= "transparent"))

# Make Crit E EX mode column
Fig.CritE <- plot_grid(Spp_CritE,
                       FRic_CritE,
                       FUn_CritE,
                       FSp_CritE,
                       labels= c("(a)","(b)","(c)","(d)"), 
                       label_size = 10,align = "hv", label_fontface = "bold",  nrow=2)

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
scenarios <- c("Future", "AT")

# Perform tests for FRic
cat("\nWilcoxon test results for FRic:\n")
for (scenario in scenarios) {
  result <- perform_wilcoxon_test(FDmetrics_FRic.CritE, FRicmetrics_null_Taxonvar.CritE, scenario)
  cat(paste("Scenario:", scenario, "\n"))
  print(result)
}

# Perform tests for FUn
cat("\nWilcoxon test results for FUn:\n")
for (scenario in scenarios) {
  result <- perform_wilcoxon_test(FDmetrics_FUn.CritE, FUnmetrics_null_Taxonvar.CritE, scenario)
  cat(paste("Scenario:", scenario, "\n"))
  print(result)
}

# Perform tests for FSp
cat("\nWilcoxon test results for FSp:\n")
for (scenario in scenarios) {
  result <- perform_wilcoxon_test(FDmetrics_FSpe.CritE, FSpemetrics_null_Taxonvar.CritE, scenario)
  cat(paste("Scenario:", scenario, "\n"))
  print(result)
}