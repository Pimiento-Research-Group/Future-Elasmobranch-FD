###############################################################################################################
# 07. Future FD plotting
## This R code plots future functional diversity (FD) based on simulated IUCN-based extinctions
## It incorporates code loaded from the supercomputer, but the results of code 6 can also be loaded
## it produces Figure 1
###############################################################################################################

# Import packages
library(tidyverse)
library(tibble)
library(reshape2)
library(cowplot)
library(ggsci)
library(RColorBrewer)
library(reshape2)
library(purrr)

## No. of extinctions across iucn simulations
# Load iucnsim data for Figure 1
dat.ext <- read_delim("C:/Users/2022207/Dropbox/Jack's PhD/Chapter 3. Future shark FD/Analyses/IUCN_sim/data/iucn_sim/future_simulations_future_status/te_all_species.txt")

# Make dat.ext a dataframe & tidy
ext.df <- dat.ext %>% as.data.frame()

ext.df1 <- rbind(colnames(ext.df), ext.df)

ext.df1[1, ] <- lapply(ext.df1[1, ], function(x) ifelse(is.na(as.numeric(x)), NA, as.numeric(x)))

ext.df1[1, 1] <- "Acroteriobatus annulatus"

# Calculate mode extinction time per species across simulations
calculate_mode <- function(x) {
  x <- x[!is.na(x)]  # Remove NAs
  if (length(x) == 0) return(NA)  # Return NA if all values are NA
  ux <- unique(x)
  counts <- tabulate(match(x, ux))
  mode_values <- ux[counts == max(counts)]
  return(mode_values[1])  # Return the first mode value
}

ext.df1$mode <- sapply(1:nrow(ext.df1), function(i) calculate_mode(ext.df1[i, -1]))

# Calculate no. extinctions per scenario for each simulation and make boxplot (out of long version?)
ext <- ext.df1
# Define the thresholds
thresholds <- c(0, 77, 100, 200, 300, 400, 500)
threshold_names <- paste0("Ext_", c("Present", "2100", "Y100", "Y200", "Y300", "Y400", "Y500"))

# Create a function to count values below a threshold
count_below_threshold <- function(df, threshold) {
  apply(df, 2, function(x) sum(as.numeric(x) < threshold, na.rm = TRUE))
}

# Create new rows for each threshold and append to the dataframe
for (i in 1:length(thresholds)) {
  new_row <- count_below_threshold(ext[-1, ], thresholds[i])  # Exclude the first row for counting
  new_row <- c(threshold_names[i], new_row)  # Add the row name 
  ext <- rbind(ext, new_row)  # Append the new row to the dataframe
}

# Filter to form dataframe
colnames(ext)[1] <- "Scenario"
filtered_ext <- ext %>% filter(Scenario %in% threshold_names)

# Long form for plotting
extFig_1 <- melt(filtered_ext, id.vars= "Scenario")
extFig_1$value <- as.numeric(extFig_1$value)
extFig_1$percent_change <- extFig_1$value / 1100 * 100

extFig <- tail(extFig_1, -7)

# Fig 1
Fig_1a <- ggplot(data = extFig, aes(x = Scenario, y = percent_change)) +
  scale_x_discrete(limits=c("Ext_Present", "Ext_2100", "Ext_Y100", "Ext_Y200", "Ext_Y300", "Ext_Y400", "Ext_Y500"),
                   labels=c("Present", "2100", "100 years", "200 years", "300 years", "400 years", "500 years")) +
  geom_violin(trim=TRUE, fill='grey', color="grey", alpha = 0.25, scale = "width") +
  labs(x = "", y = "% extinctions") +
  geom_point(aes(x = 2, y = 4), color = "#1f78b4", fill = "#1f78b4", shape = 21, size = 5) + #mode
  geom_point(aes(x = 3, y = 5.27), color = "#6a3d9a", fill = "#6a3d9a", shape = 21, size = 5) +
  geom_point(aes(x = 4, y = 11.63), color = "#b19cd9", fill = "#b19cd9", shape = 21, size = 5) +
  geom_point(aes(x = 5, y = 19.81), color = "#ffa07a", fill = "#ffa07a", shape = 21, size = 5) + #Colour = alpha("#683646", 0.8)
  geom_point(aes(x = 6, y = 26.82), color = "#ff7f00", fill = "#ff7f00", shape = 21, size = 5) +  #"#683646"
  geom_point(aes(x = 7, y = 32.45), color = "#e31a1c", fill = "#e31a1c", shape = 21, size = 5) + #"#AA1638"
  ggtitle("")+
  theme_minimal() +
  theme_bw() +
  theme(legend.position = "none") +
  theme(axis.text = element_text(size= 6.5, color= "black"),
        axis.title= element_text(size= 8), panel.background= element_rect(fill= "white"))+
  theme(panel.grid = element_blank(), 
        panel.border = element_rect(fill= "transparent")) +
  theme(axis.text.x = element_blank())+
  theme(plot.title = element_text(hjust = 0.5)) + 
  scale_y_reverse()

# Load future FD results from code 4
load(file = "C:/Users/2022207/Dropbox/Jack's PhD/Chapter 3. Future shark FD/Data/R data/Future_FD.RData")

# Replicate 1000 times to match with null distribution
Res_FDmetrics_TaxonVar <- map_dfr(1:1000, ~ FDindices_taxon)

# Convert to long
FDmetrics_long_TaxonVar<- melt(Res_FDmetrics_TaxonVar, id.vars= "Scenario")
save(FDmetrics_long_TaxonVar,file = "C:/Users/2022207/Dropbox/Jack's PhD/Chapter 3. Future shark FD/Data/R data/Future_Taxon_variation_long_metrics.RData")

# Load null supercomputer results
load(file = "C:/Users/2022207/Dropbox/Jack's PhD/Chapter 3. Future shark FD/Data/R data/ModeImps1_500_Future_null_long_metrics.RData")
FDmetrics_null_long_taxonvar_1_500 <- FDmetrics_null_long_taxonvar
load(file = "C:/Users/2022207/Dropbox/Jack's PhD/Chapter 3. Future shark FD/Data/R data/ModeImps501_1000_Future_null_long_metrics.RData")
FDmetrics_null_long_taxonvar_501_1000 <- FDmetrics_null_long_taxonvar

# Join
FDmetrics_null_long_taxonvar <- rbind(FDmetrics_null_long_taxonvar_1_500, FDmetrics_null_long_taxonvar_501_1000)
save(FDmetrics_null_long_taxonvar,file = "C:/Users/2022207/Dropbox/Jack's PhD/Chapter 3. Future shark FD/Data/R data/Merged_mode_Future_null_long_metrics.RData")

# Load wide supercomputer results (null) and calculate medians
load(file = "C:/Users/2022207/Dropbox/Jack's PhD/Chapter 3. Future shark FD/Data/R data/ModeImps1_500_Future_null_metrics.RData")
FDmetrics_null_wide_TaxonVar_1_500 <- Null_FDmetrics_taxonvar
load(file = "C:/Users/2022207/Dropbox/Jack's PhD/Chapter 3. Future shark FD/Data/R data/ModeImps501_1000_Future_null_metrics.RData")
FDmetrics_null_wide_TaxonVar_501_1000 <- Null_FDmetrics_taxonvar
Null_FDmetrics_taxonvar <- rbind(FDmetrics_null_wide_TaxonVar_1_500, FDmetrics_null_wide_TaxonVar_501_1000)
save(Null_FDmetrics_taxonvar,file = "C:/Users/2022207/Dropbox/Jack's PhD/Chapter 3. Future shark FD/Data/R data/Merged_mode_Future_null_metrics.RData")

# Calculate medians and % changes
Future_Null_var <- Null_FDmetrics_taxonvar %>% 
  group_by(Scenario) %>%
  summarise(Sp_mean = mean(sp_richn),
            Sp_med = median(sp_richn),
            Sp_sd = sd(sp_richn),
            FRic_mean = mean(fric),
            FRic_med = median(fric),
            FRic_sd = sd(fric),
            Fori_mean = mean(fori),
            Fori_med = median(fori),
            Fori_sd = sd(fori),
            Fspe_mean = mean(fspe),
            Fspe_med = median(fspe),
            Fspe_sd = sd(fspe),
            Fun_mean = mean(fun),
            Fun_med = median(fun),
            Fun_sd = sd(fun)) %>% 
  mutate(Sp_decrease = ifelse(Scenario != "Present", (Sp_med - Sp_med[Scenario == "Present"]) / Sp_med[Scenario == "Present"] * 100, 0),
         FRic_decrease = ifelse(Scenario != "Present", (FRic_med - FRic_med[Scenario == "Present"]) / FRic_med[Scenario == "Present"] * 100, 0),
         Fori_decrease = ifelse(Scenario != "Present", (Fori_med - Fori_med[Scenario == "Present"]) / Fori_med[Scenario == "Present"] * 100, 0),
         Fspe_decrease = ifelse(Scenario != "Present", (Fspe_med - Fspe_med[Scenario == "Present"]) / Fspe_med[Scenario == "Present"] * 100, 0),
         Fun_decrease = ifelse(Scenario != "Present", (Fun_med - Fun_med[Scenario == "Present"]) / Fun_med[Scenario == "Present"] * 100, 0)
  )

Future_Null_var$Scenario <- ordered(Future_Null_var$Scenario, levels=c("Present","Future","Y100","Y200","Y300","Y400","Y500"))

## Filter by metric for both empirical & null models
SPPmetrics_null_Taxonvar <- FDmetrics_null_long_taxonvar %>% 
  filter(variable == "sp_richn")
FRicmetrics_null_Taxonvar <- FDmetrics_null_long_taxonvar %>% 
  filter(variable == "fric")
FOrimetrics_null_Taxonvar <- FDmetrics_null_long_taxonvar %>% 
  filter(variable == "fori")
FSpemetrics_null_Taxonvar <- FDmetrics_null_long_taxonvar %>% 
  filter(variable == "fspe")
FUnmetrics_null_Taxonvar <- FDmetrics_null_long_taxonvar %>% 
  filter(variable == "fun")

FDmetrics_Spp <- FDmetrics_long_TaxonVar %>% 
  filter(variable == "sp_richn")
FDmetrics_FRic <- FDmetrics_long_TaxonVar %>% 
  filter(variable == "fric")
FDmetrics_FOri <- FDmetrics_long_TaxonVar %>% 
  filter(variable == "fori")
FDmetrics_FSpe <- FDmetrics_long_TaxonVar %>% 
  filter(variable == "fspe")
FDmetrics_FUn <- FDmetrics_long_TaxonVar %>% 
  filter(variable == "fun")

# Calculate % change for plots
## Empirical
FDmetrics_Spp <- FDmetrics_Spp %>% 
  mutate(percent_change = ((value - 1100) / 1100) * 100)
FDmetrics_FRic <- FDmetrics_FRic %>% 
  mutate(percent_change = ((value - 1) / 1) * 100)
FDmetrics_FOri <- FDmetrics_FOri %>% 
  mutate(percent_change = ((value - 0.01766838) / 0.01766838) * 100)
FDmetrics_FSpe <- FDmetrics_FSpe %>% 
  mutate(percent_change = ((value - 0.2746819) / 0.2746819) * 100)
FDmetrics_FUn <- FDmetrics_FUn %>% 
  mutate(percent_change = ((value - 0.004141561) / 0.004141561) * 100)
## Null
SPPmetrics_null_Taxonvar <- SPPmetrics_null_Taxonvar %>% 
  mutate(percent_change = ((value - 1100) / 1100) * 100)
FRicmetrics_null_Taxonvar <- FRicmetrics_null_Taxonvar %>% 
  mutate(percent_change = ((value - 1) / 1) * 100)
FOrimetrics_null_Taxonvar <- FOrimetrics_null_Taxonvar %>% 
  mutate(percent_change = ((value - 0.01766838) / 0.01766838) * 100)
FSpemetrics_null_Taxonvar <- FSpemetrics_null_Taxonvar %>% 
  mutate(percent_change = ((value - 0.2746819) / 0.2746819) * 100)
FUnmetrics_null_Taxonvar <- FUnmetrics_null_Taxonvar %>% 
  mutate(percent_change = ((value - 0.004141561) / 0.004141561) * 100)

# For plotting single dot
FDmetrics_FRic_subset <- FDmetrics_FRic[1:7, ]
FDmetrics_FUn_subset <- FDmetrics_FUn[1:7, ]
FDmetrics_FSpe_subset <- FDmetrics_FSpe[1:7, ]

# Plots
## Functional richness
FRic_null_variation <- ggplot(data=FRicmetrics_null_Taxonvar, aes(x= Scenario, y= percent_change)) +
  geom_violin(trim=TRUE, fill='grey', color="grey", alpha = 0.25, scale = "width") +
  scale_x_discrete(limits=c("Present", "Future", "Y100", "Y200", "Y300", "Y400", "Y500")) +
  geom_point(data = FDmetrics_FRic_subset, 
             aes(x = Scenario, y = percent_change,
                 color = Scenario, fill = Scenario),
             size = 5, shape = 21) +
  scale_color_manual(values=c("Present"="#add8e6", "Future"="#1f78b4", "Y100"="#6a3d9a", 
                              "Y200"="#b19cd9", "Y300"="#ffa07a", "Y400"="#ff7f00", "Y500"="#e31a1c")) +
  scale_fill_manual(values=c("Present"="#add8e6", "Future"="#1f78b4", "Y100"="#6a3d9a", 
                             "Y200"="#b19cd9", "Y300"="#ffa07a", "Y400"="#ff7f00", "Y500"="#e31a1c")) +
  labs(x = "", y = "FRic % change") +
  ggtitle("") +
  theme_minimal() +
  theme_bw() +
  theme(legend.position = "none",
        axis.text = element_text(size= 6.5, color= "black"),
        axis.title= element_text(size= 8), 
        panel.background= element_rect(fill= "white"),
        panel.grid = element_blank(), 
        panel.border = element_rect(fill= "transparent"),
        axis.text.x = element_blank(),
        plot.title = element_text(hjust = 0.5))

# Functional Uniqueness
FUn_null_variation <- ggplot(data=FUnmetrics_null_Taxonvar, aes(x= Scenario, y= percent_change)) +
  geom_violin(trim=TRUE, fill='grey', color="grey", alpha = 0.25, scale = "width") +
  scale_x_discrete(limits=c("Present", "Future", "Y100", "Y200", "Y300", "Y400", "Y500"),
                   labels=c("Present", "2100", "100 years", "200 years", "300 years", "400 years", "500 years")) +
  geom_point(data = FDmetrics_FUn_subset, 
             aes(x = Scenario, y = percent_change, color = Scenario, fill = Scenario), 
             size = 5, shape = 21) +
  scale_color_manual(values=c("Present"="#add8e6", "Future"="#1f78b4", "Y100"="#6a3d9a", 
                              "Y200"="#b19cd9", "Y300"="#ffa07a", "Y400"="#ff7f00", "Y500"="#e31a1c")) +
  scale_fill_manual(values=c("Present"="#add8e6", "Future"="#1f78b4", "Y100"="#6a3d9a", 
                             "Y200"="#b19cd9", "Y300"="#ffa07a", "Y400"="#ff7f00", "Y500"="#e31a1c")) +
  labs(x = "", y = "FUn % change") +
  ggtitle("") +
  theme_minimal() +
  theme_bw() +
  theme(legend.position = "none",
        axis.text = element_text(size= 6.5, color= "black"),
        axis.title= element_text(size= 8), 
        panel.background= element_rect(fill= "white"),
        panel.grid = element_blank(), 
        panel.border = element_rect(fill= "transparent"),
        axis.text.x = element_blank(),
        plot.title = element_text(hjust = 0.5))

# Functional Specialisation
FSp_null_variation <- ggplot(data=FSpemetrics_null_Taxonvar, aes(x= Scenario, y= percent_change)) +
  geom_violin(trim=TRUE, fill='grey', color="grey", alpha = 0.25, scale = "width") +
  scale_x_discrete(limits=c("Present", "Future", "Y100", "Y200", "Y300", "Y400", "Y500"),
                   labels=c("Present", "2100", "100 years", "200 years", "300 years", "400 years", "500 years")) +
  geom_point(data = FDmetrics_FSpe_subset, 
             aes(x = Scenario, y = percent_change, color = Scenario, fill = Scenario), 
             size = 5, shape = 21) +
  scale_color_manual(values=c("Present"="#add8e6", "Future"="#1f78b4", "Y100"="#6a3d9a", 
                              "Y200"="#b19cd9", "Y300"="#ffa07a", "Y400"="#ff7f00", "Y500"="#e31a1c")) +
  scale_fill_manual(values=c("Present"="#add8e6", "Future"="#1f78b4", "Y100"="#6a3d9a", 
                             "Y200"="#b19cd9", "Y300"="#ffa07a", "Y400"="#ff7f00", "Y500"="#e31a1c")) +
  labs(x = "Scenario", y = "FSp % change") +
  theme_minimal() +
  theme_bw() +
  theme(legend.position = "none") +
  theme(axis.text = element_text(size= 6.5, color= "black"),
        axis.title= element_text(size= 8), 
        panel.background= element_rect(fill= "white")) +
  theme(panel.grid = element_blank(), 
        panel.border = element_rect(fill= "transparent"))

# Final plot - plots Figure 1
Fig <- plot_grid(Fig_1a,
                 FRic_null_variation,
                 FUn_null_variation,
                 FSp_null_variation,
                 labels= c("(a)","(b)","(c)","(d)"), 
                 label_size = 10,align = "hv", label_fontface = "bold",  nrow=4)