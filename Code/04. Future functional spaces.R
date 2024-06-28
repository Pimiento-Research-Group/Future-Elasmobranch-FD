###############################################################################################################
# 04. Future functional spaces
## This R code forms future functional spaces based on modal extinction times from IUCN simulations
## it produces produces Figure S5
###############################################################################################################

# Import packages
library(tidyverse)
library(tibble)
library(scales)
library(cowplot)
library(doBy)
library(reshape2)
library(RColorBrewer)
library(mFD)
library(ggsci)
library(data.table)
library(funrar)
library(viridis)

# Load data
load(file="C:/Users/2022207/Dropbox/Jack's PhD/Chapter 3. Future shark FD/Data/R data/Data.RData")

# Load mode IUCN_sim results
load(file = "C:/Users/2022207/Dropbox/Jack's PhD/Chapter 3. Future shark FD/Data/R data/2100_mode.RData")

# Make future functional space
# Occurrence data
Occ_Future <- data %>% 
  select(Species)

Occ_Future$Present <- 1

# Set up future scenarios
## 2100
Occ_Future <- Occ_Future %>%
  left_join(Spp_mode, by = "Species") %>%
  mutate(Future = ifelse(is.na(mode), 1, ifelse(mode <= 77, 0, 1))) %>%
  select(Species,Present,Future)

## 100 years
Occ_Future <- Occ_Future %>%
  left_join(Spp_mode, by = "Species") %>%
  mutate(Y100 = ifelse(is.na(mode), 1, ifelse(mode <= 100, 0, 1))) %>%
  select(Species,Present,Future,Y100)

## 200 years
Occ_Future <- Occ_Future %>%
  left_join(Spp_mode, by = "Species") %>%
  mutate(Y200 = ifelse(is.na(mode), 1, ifelse(mode <= 200, 0, 1))) %>%
  select(Species,Present,Future,Y100,Y200)

## 300 years
Occ_Future <- Occ_Future %>%
  left_join(Spp_mode, by = "Species") %>%
  mutate(Y300 = ifelse(is.na(mode), 1, ifelse(mode <= 300, 0, 1))) %>%
  select(Species,Present,Future,Y100,Y200,Y300)

## 400 years
Occ_Future <- Occ_Future %>%
  left_join(Spp_mode, by = "Species") %>%
  mutate(Y400 = ifelse(is.na(mode), 1, ifelse(mode <= 400, 0, 1))) %>%
  select(Species,Present,Future,Y100,Y200,Y300,Y400)

## 500 years
Occ_Future <- Occ_Future %>%
  left_join(Spp_mode, by = "Species") %>%
  mutate(Y500 = ifelse(is.na(mode), 1, ifelse(mode <= 500, 0, 1))) %>%
  select(Species,Present,Future,Y100,Y200,Y300,Y400,Y500)

baskets <- Occ_Future %>% 
  remove_rownames %>% 
  column_to_rownames(var="Species")

# Transpose to form final occurrence matrix
baskets_sharks_weights <- t(baskets)

# Species-trait matrix
traits <- data %>% 
  select(Species,habitat,vertical,terrestriality,thermo,feeding,diet,max.length)

sharks_traits <- traits %>% 
  remove_rownames %>% 
  column_to_rownames(var="Species")

# Form trait category matrix
sharks_traits_cat <- tibble(trait_name = c("habitat","vertical","terrestriality","thermo","feeding","diet","max.length"),
                            trait_type = c("N", "N", "O", "O", "O", "N", "Q"))

# Summarise dataset
sharks_traits_summ <- mFD::sp.tr.summary(
  tr_cat     = sharks_traits_cat,   
  sp_tr      = sharks_traits, 
  stop_if_NA = TRUE)

# Reveal number of states for each trait
sharks_traits_summ$tr_summary_list

# make assemblage data into a matrix & ensure all matrix values are numeric
baskets_sharks_weights <- data.matrix(baskets_sharks_weights, rownames.force = NA)
class(baskets_sharks_weights) <- "numeric"

#summarise dataset
asb_sp_sharks_summ <- mFD::asb.sp.summary(asb_sp_w = baskets_sharks_weights)

# retrieve species occurrences for all assemblages:
asb_sp_sharks_occ <- asb_sp_sharks_summ$asb_sp_occ

# Construct trait distance matrix using species
sp_dist_sharks <- mFD::funct.dist(
  sp_tr         = sharks_traits,
  tr_cat        = sharks_traits_cat,
  metric        = "gower",
  scale_euclid  = "scale_center",
  ordinal_var   = "classic",
  weight_type   = "equal",
  stop_if_NA    = TRUE)

# Assess quality
fspaces_quality_sharks <- mFD::quality.fspaces(
  sp_dist             = sp_dist_sharks,
  maxdim_pcoa         = 10,
  deviation_weighting = c("absolute","squared"),
  fdist_scaling       = FALSE,
  fdendro             = "average")

# Return coordinates of each axis
sp_faxes_coord_sharks <- fspaces_quality_sharks$"details_fspaces"$"sp_pc_coord"

# Calculate FD metrics
alpha_fd_indices_sharks <- mFD::alpha.fd.multidim(
  sp_faxes_coord   = sp_faxes_coord_sharks[ , c("PC1", "PC2", "PC3", "PC4")],
  asb_sp_w         = baskets_sharks_weights,
  ind_vect         = c("fric", "fori", "fspe"),
  scaling          = TRUE,
  check_input      = TRUE,
  details_returned = TRUE)

fd_ind_values_sharks <- alpha_fd_indices_sharks$"functional_diversity_indices"

## Form dataframe
FDindices_taxon <- as.data.frame(fd_ind_values_sharks) %>% 
  tibble::rownames_to_column("Scenario")
FDindices_taxon$Scenario <- ordered(FDindices_taxon$Scenario, levels=c("Present","Future","Y100","Y200","Y300","Y400","Y500"))

# Calculate functional uniqueness per species
Funct_Uni <- uniqueness(baskets_sharks_weights, as.matrix(sp_dist_sharks)) %>% 
  as.data.frame()
Funct_Uni$Present <- baskets$Present
Funct_Uni$Future <- baskets$Future
Funct_Uni$Y100 <- baskets$Y100
Funct_Uni$Y200 <- baskets$Y200
Funct_Uni$Y300 <- baskets$Y300
Funct_Uni$Y400 <- baskets$Y400
Funct_Uni$Y500 <- baskets$Y500

# Calculate functional uniqueness per scenario
Present_FUn <- mean(Funct_Uni$Ui)
Future_FUn_df <- Funct_Uni %>% 
  filter(Future == "1") %>% 
  select(Ui)
Future_FUn <- mean(Future_FUn_df$Ui)

Y100_FUn_df <- Funct_Uni %>% 
  filter(Y100 == "1") %>% 
  select(Ui)
Y100_FUn <- mean(Y100_FUn_df$Ui)

Y200_FUn_df <- Funct_Uni %>% 
  filter(Y200 == "1") %>% 
  select(Ui)
Y200_FUn <- mean(Y200_FUn_df$Ui)

Y300_FUn_df <- Funct_Uni %>% 
  filter(Y300 == "1") %>% 
  select(Ui)
Y300_FUn <- mean(Y300_FUn_df$Ui)

Y400_FUn_df <- Funct_Uni %>% 
  filter(Y400 == "1") %>% 
  select(Ui)
Y400_FUn <- mean(Y400_FUn_df$Ui)

Y500_FUn_df <- Funct_Uni %>% 
  filter(Y500 == "1") %>% 
  select(Ui)
Y500_FUn <- mean(Y500_FUn_df$Ui)

# Add uniqueness to final FD results data frame
Fun <- data.frame(
  Fun = c(Present_FUn,Future_FUn,Y100_FUn,Y200_FUn,Y300_FUn,Y400_FUn,Y500_FUn),
  row.names = c("Present","Future","Y100","Y200","Y300","Y400","Y500")
)

FDindices_taxon$fun <- Fun$Fun

# Save future FD results for plotting
save(FDindices_taxon, file = "C:/Users/2022207/Dropbox/Jack's PhD/Chapter 3. Future shark FD/Data/R data/Future_FD.RData")

# Form future functional spaces
limits_PC1 <- c(-0.4, 0.4)
limits_PC2 <- c(-0.3, 0.5)
limits_PC3 <- c(-0.3, 0.4)
limits_PC4 <- c(-0.6, 0.3)

## Plot
Present_space <- mFD::alpha.multidim.plot(
  output_alpha_fd_multidim = alpha_fd_indices_sharks,
  plot_asb_nm              = "Present",
  ind_nm                   = c("fric", "fori", "fspe"),
  faxes                    = NULL,
  faxes_nm                 = NULL,
  range_faxes              = c(NA, NA),
  color_bg                 = NA,
  shape_sp                 = c(pool = 24, asb1 = 21),
  size_sp                  = c(pool = 1.5, asb1 = 1),
  color_sp                 = c(pool = "black", asb1 = "#F8766D"),
  color_vert               = c(pool = "black", asb1 = "#F8766D"),
  fill_sp                  = c(pool = "black", asb1 = "#F8766D"), 
  fill_vert                = c(pool = "black", asb1 = "#F8766D"),
  color_ch                 = c(pool = 'grey50', asb1 = "#F8766D"),  # Set dark grey outline for convex hull
  fill_ch                  = c(pool = "white", asb1 = "#F8766D"),
  alpha_ch                 = c(pool = 1, asb1 = 0.1),
  shape_centroid_fdis      = c(asb1 = 22),
  shape_centroid_fdiv      = c(asb1 = 22),
  shape_centroid_fspe      = 23,
  color_centroid_fspe      = "black",
  size_sp_nm               = 3, 
  color_sp_nm              = "black",
  plot_sp_nm               = NULL,
  fontface_sp_nm           = "plain",
  save_file                = FALSE,
  check_input              = TRUE)

Present_a <- Present_space$fric$PC1_PC2 +
  labs(x = "PCoA1", y = "PCoA2")+
  coord_fixed(ratio = 1, xlim = limits_PC1, ylim = limits_PC2)

Present_b <- Present_space$fric$PC3_PC4 +
  labs(x = "PCoA3", y = "PCoA4")+
  coord_fixed(ratio = 1, xlim = limits_PC3, ylim = limits_PC4)

Future_space <- mFD::alpha.multidim.plot(
  output_alpha_fd_multidim = alpha_fd_indices_sharks,
  plot_asb_nm              = "Future",
  ind_nm                   = c("fric", "fori", "fspe"),
  faxes                    = NULL,
  faxes_nm                 = NULL,
  range_faxes              = c(NA, NA),
  color_bg                 = NA,
  shape_sp                 = c(pool = 24, asb1 = 21),
  size_sp                  = c(pool = 1.5, asb1 = 1),
  color_sp                 = c(pool = "black", asb1 = "#F8766D"),
  color_vert               = c(pool = "black", asb1 = "#F8766D"),
  fill_sp                  = c(pool = "black", asb1 = "#F8766D"),  
  fill_vert                = c(pool = "black", asb1 = "#F8766D"),
  color_ch                 = c(pool = 'grey50', asb1 = "#F8766D"),  # Set dark grey outline for convex hull
  fill_ch                  = c(pool = "white", asb1 = "#F8766D"),
  alpha_ch                 = c(pool = 1, asb1 = 0.1),
  shape_centroid_fdis      = c(asb1 = 22),
  shape_centroid_fdiv      = c(asb1 = 22),
  shape_centroid_fspe      = 23,
  color_centroid_fspe      = "black",
  size_sp_nm               = 3, 
  color_sp_nm              = "black",
  plot_sp_nm               = NULL,
  fontface_sp_nm           = "plain",
  save_file                = FALSE,
  check_input              = TRUE)

Future_2100a <- Future_space$fric$PC1_PC2 +
  labs(x = "PCoA1", y = "PCoA2")+
  coord_fixed(ratio = 1, xlim = limits_PC1, ylim = limits_PC2)

Future_2100b <- Future_space$fric$PC3_PC4 +
  labs(x = "PCoA3", y = "PCoA4")+
  coord_fixed(ratio = 1, xlim = limits_PC3, ylim = limits_PC4)

Y100_space <- mFD::alpha.multidim.plot(
  output_alpha_fd_multidim = alpha_fd_indices_sharks,
  plot_asb_nm              = "Y100",
  ind_nm                   = c("fric", "fori", "fspe"),
  faxes                    = NULL,
  faxes_nm                 = NULL,
  range_faxes              = c(NA, NA),
  color_bg                 = NA,
  shape_sp                 = c(pool = 24, asb1 = 21),
  size_sp                  = c(pool = 1.5, asb1 = 1),
  color_sp                 = c(pool = "black", asb1 = "#F8766D"),
  color_vert               = c(pool = "black", asb1 = "#F8766D"),
  fill_sp                  = c(pool = "black", asb1 = "#F8766D"),  
  fill_vert                = c(pool = "black", asb1 = "#F8766D"),
  color_ch                 = c(pool = 'grey50', asb1 = "#F8766D"),  # Set dark grey outline for convex hull
  fill_ch                  = c(pool = "white", asb1 = "#F8766D"),
  alpha_ch                 = c(pool = 1, asb1 = 0.1),
  shape_centroid_fdis      = c(asb1 = 22),
  shape_centroid_fdiv      = c(asb1 = 22),
  shape_centroid_fspe      = 23,
  color_centroid_fspe      = "black",
  size_sp_nm               = 3, 
  color_sp_nm              = "black",
  plot_sp_nm               = NULL,
  fontface_sp_nm           = "plain",
  save_file                = FALSE,
  check_input              = TRUE)

Future_100a <- Y100_space$fric$PC1_PC2 +
  labs(x = "PCoA1", y = "PCoA2")+
  coord_fixed(ratio = 1, xlim = limits_PC1, ylim = limits_PC2)

Future_100b <- Y100_space$fric$PC3_PC4 +
  labs(x = "PCoA3", y = "PCoA4")+
  coord_fixed(ratio = 1, xlim = limits_PC3, ylim = limits_PC4)

Y200_space <- mFD::alpha.multidim.plot(
  output_alpha_fd_multidim = alpha_fd_indices_sharks,
  plot_asb_nm              = "Y200",
  ind_nm                   = c("fric", "fori", "fspe"),
  faxes                    = NULL,
  faxes_nm                 = NULL,
  range_faxes              = c(NA, NA),
  color_bg                 = NA,
  shape_sp                 = c(pool = 24, asb1 = 21),
  size_sp                  = c(pool = 1.5, asb1 = 1),
  color_sp                 = c(pool = "black", asb1 = "#F8766D"),
  color_vert               = c(pool = "black", asb1 = "#F8766D"),
  fill_sp                  = c(pool = "black", asb1 = "#F8766D"),  
  fill_vert                = c(pool = "black", asb1 = "#F8766D"),
  color_ch                 = c(pool = 'grey50', asb1 = "#F8766D"),  # Set dark grey outline for convex hull
  fill_ch                  = c(pool = "white", asb1 = "#F8766D"),
  alpha_ch                 = c(pool = 1, asb1 = 0.1),
  shape_centroid_fdis      = c(asb1 = 22),
  shape_centroid_fdiv      = c(asb1 = 22),
  shape_centroid_fspe      = 23,
  color_centroid_fspe      = "black",
  size_sp_nm               = 3, 
  color_sp_nm              = "black",
  plot_sp_nm               = NULL,
  fontface_sp_nm           = "plain",
  save_file                = FALSE,
  check_input              = TRUE)

## Plot future space
Future_200a <- Y200_space$fric$PC1_PC2 +
  labs(x = "PCoA1", y = "PCoA2")+
  coord_fixed(ratio = 1, xlim = limits_PC1, ylim = limits_PC2)

Future_200b <- Y200_space$fric$PC3_PC4 +
  labs(x = "PCoA3", y = "PCoA4")+
  coord_fixed(ratio = 1, xlim = limits_PC3, ylim = limits_PC4)

Y300_space <- mFD::alpha.multidim.plot(
  output_alpha_fd_multidim = alpha_fd_indices_sharks,
  plot_asb_nm              = "Y300",
  ind_nm                   = c("fric", "fori", "fspe"),
  faxes                    = NULL,
  faxes_nm                 = NULL,
  range_faxes              = c(NA, NA),
  color_bg                 = NA,
  shape_sp                 = c(pool = 24, asb1 = 21),
  size_sp                  = c(pool = 1.5, asb1 = 1),
  color_sp                 = c(pool = "black", asb1 = "#F8766D"),
  color_vert               = c(pool = "black", asb1 = "#F8766D"),
  fill_sp                  = c(pool = "black", asb1 = "#F8766D"),  
  fill_vert                = c(pool = "black", asb1 = "#F8766D"),
  color_ch                 = c(pool = 'grey50', asb1 = "#F8766D"),  # Set dark grey outline for convex hull
  fill_ch                  = c(pool = "white", asb1 = "#F8766D"),
  alpha_ch                 = c(pool = 1, asb1 = 0.1),
  shape_centroid_fdis      = c(asb1 = 22),
  shape_centroid_fdiv      = c(asb1 = 22),
  shape_centroid_fspe      = 23,
  color_centroid_fspe      = "black",
  size_sp_nm               = 3, 
  color_sp_nm              = "black",
  plot_sp_nm               = NULL,
  fontface_sp_nm           = "plain",
  save_file                = FALSE,
  check_input              = TRUE)

Future_300a <- Y300_space$fric$PC1_PC2 +
  labs(x = "PCoA1", y = "PCoA2")+
  coord_fixed(ratio = 1, xlim = limits_PC1, ylim = limits_PC2)

Future_300b <- Y300_space$fric$PC3_PC4 +
  labs(x = "PCoA3", y = "PCoA4")+
  coord_fixed(ratio = 1, xlim = limits_PC3, ylim = limits_PC4)

Y400_space <- mFD::alpha.multidim.plot(
  output_alpha_fd_multidim = alpha_fd_indices_sharks,
  plot_asb_nm              = "Y400",
  ind_nm                   = c("fric", "fori", "fspe"),
  faxes                    = NULL,
  faxes_nm                 = NULL,
  range_faxes              = c(NA, NA),
  color_bg                 = NA,
  shape_sp                 = c(pool = 24, asb1 = 21),
  size_sp                  = c(pool = 1.5, asb1 = 1),
  color_sp                 = c(pool = "black", asb1 = "#F8766D"),
  color_vert               = c(pool = "black", asb1 = "#F8766D"),
  fill_sp                  = c(pool = "black", asb1 = "#F8766D"),  
  fill_vert                = c(pool = "black", asb1 = "#F8766D"),
  color_ch                 = c(pool = 'grey50', asb1 = "#F8766D"),  # Set dark grey outline for convex hull
  fill_ch                  = c(pool = "white", asb1 = "#F8766D"),
  alpha_ch                 = c(pool = 1, asb1 = 0.1),
  shape_centroid_fdis      = c(asb1 = 22),
  shape_centroid_fdiv      = c(asb1 = 22),
  shape_centroid_fspe      = 23,
  color_centroid_fspe      = "black",
  size_sp_nm               = 3, 
  color_sp_nm              = "black",
  plot_sp_nm               = NULL,
  fontface_sp_nm           = "plain",
  save_file                = FALSE,
  check_input              = TRUE)

Future_400a <- Y400_space$fric$PC1_PC2 +
  labs(x = "PCoA1", y = "PCoA2")+
  coord_fixed(ratio = 1, xlim = limits_PC1, ylim = limits_PC2)

Future_400b <- Y400_space$fric$PC3_PC4 +
  labs(x = "PCoA3", y = "PCoA4")+
  coord_fixed(ratio = 1, xlim = limits_PC3, ylim = limits_PC4)

Y500_space <- mFD::alpha.multidim.plot(
  output_alpha_fd_multidim = alpha_fd_indices_sharks,
  plot_asb_nm              = "Y500",
  ind_nm                   = c("fric", "fori", "fspe"),
  faxes                    = NULL,
  faxes_nm                 = NULL,
  range_faxes              = c(NA, NA),
  color_bg                 = NA,
  shape_sp                 = c(pool = 24, asb1 = 21),
  size_sp                  = c(pool = 1.5, asb1 = 1),
  color_sp                 = c(pool = "black", asb1 = "#F8766D"),
  color_vert               = c(pool = "black", asb1 = "#F8766D"),
  fill_sp                  = c(pool = "black", asb1 = "#F8766D"), 
  fill_vert                = c(pool = "black", asb1 = "#F8766D"),
  color_ch                 = c(pool = 'grey50', asb1 = "#F8766D"),  # Set dark grey outline for convex hull
  fill_ch                  = c(pool = "white", asb1 = "#F8766D"),
  alpha_ch                 = c(pool = 1, asb1 = 0.1),
  shape_centroid_fdis      = c(asb1 = 22),
  shape_centroid_fdiv      = c(asb1 = 22),
  shape_centroid_fspe      = 23,
  color_centroid_fspe      = "black",
  size_sp_nm               = 3, 
  color_sp_nm              = "black",
  plot_sp_nm               = NULL,
  fontface_sp_nm           = "plain",
  save_file                = FALSE,
  check_input              = TRUE)

Future_500a <- Y500_space$fric$PC1_PC2 +
  labs(x = "PCoA1", y = "PCoA2")+
  coord_fixed(ratio = 1, xlim = limits_PC1, ylim = limits_PC2)

Future_500b <- Y500_space$fric$PC3_PC4 +
  labs(x = "PCoA3", y = "PCoA4")+
  coord_fixed(ratio = 1, xlim = limits_PC3, ylim = limits_PC4)

# Plot all future spaces
Sup_Fig <- plot_grid(Present_a,Present_b,
                     Future_2100a,Future_2100b,
                     Future_100a,Future_100b,
                     Future_200a,Future_200b,
                     Future_300a,Future_300b,
                     Future_400a,Future_400b,
                     Future_500a,Future_500b,
                   labels = c("(a)","(b)","(c)","(d)","(e)","(f)","(g)","(h)","(i)","(j)","(k)","(l)","(m)","(n)"),
                   label_size = 12,align = "hv", label_fontface = "bold", hjust = -0.15,  nrow=7)

# Create a plot to extract the legend
legend_plot <- ggplot() +
  geom_point(aes(x = 1, y = 1, shape = "extant"), color = "#F8766D", size = 3) +
  geom_point(aes(x = 1, y = 1, shape = "extinct"), color = "black", size = 3) +
  scale_shape_manual(name = NULL, values = c(extant = 16, extinct = 17), labels = c("extant", "extinct")) +
  guides(shape = guide_legend(override.aes = list(color = c("#F8766D", "black")))) +
  theme_void() +
  theme(legend.position = "top", legend.direction = "horizontal", legend.title = element_blank())

# Extract the legend
legend <- get_legend(legend_plot)

# Produce final figure - produces Figure S5
final_sup_fig <- plot_grid(legend, Sup_Fig, ncol = 1, rel_heights = c(0.1, 1))
