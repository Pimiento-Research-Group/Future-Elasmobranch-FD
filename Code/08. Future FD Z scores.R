###############################################################################################################
# 08. Future FD Z-scores
## This R code calculates if changes in functional diversity are greater than expected by random chance
## it produces Table S4
###############################################################################################################

# Import packages
library(tidyverse)
library(tibble)
library(scales)
library(doBy)
library(reshape2)
library(data.table)
library(writexl)

# Load empirical results and null distribution
load(file = "C:/Users/2022207/Dropbox/Jack's PhD/Chapter 3. Future shark FD/Data/R data/Future_Taxon_variation_long_metrics.RData")

load(file = "C:/Users/2022207/Dropbox/Jack's PhD/Chapter 3. Future shark FD/Data/R data/Future_null_metrics.RData")
load(file = "C:/Users/2022207/Dropbox/Jack's PhD/Chapter 3. Future shark FD/Data/R data/Future_null_long_metrics.RData")

# Test for normality - to determine if Z scores are based on mean or median
FRicmetrics_null_Taxonvar <- FDmetrics_null_long_taxonvar %>% 
  filter(variable == "fric")
FOrimetrics_null_Taxonvar <- FDmetrics_null_long_taxonvar %>% 
  filter(variable == "fori")
FSpemetrics_null_Taxonvar <- FDmetrics_null_long_taxonvar %>% 
  filter(variable == "fspe")
FUnmetrics_null_Taxonvar <- FDmetrics_null_long_taxonvar %>% 
  filter(variable == "fun")

shapiro.test(FRicmetrics_null_Taxonvar$value[1:3500])
shapiro.test(FRicmetrics_null_Taxonvar$value[3501:7000])
shapiro.test(FSpemetrics_null_Taxonvar$value[1:3500])
shapiro.test(FSpemetrics_null_Taxonvar$value[3501:7000])
shapiro.test(FUnmetrics_null_Taxonvar$value[1:3500])
shapiro.test(FUnmetrics_null_Taxonvar$value[3501:7000])   # None are normally distributed; so use median

FDmetrics_FRic <- FDmetrics_long_TaxonVar %>% 
  filter(variable == "fric")
FDmetrics_FOri <- FDmetrics_long_TaxonVar %>% 
  filter(variable == "fori")
FDmetrics_FSpe <- FDmetrics_long_TaxonVar %>% 
  filter(variable == "fspe")
FDmetrics_FUn <- FDmetrics_long_TaxonVar %>% 
  filter(variable == "fun")

shapiro.test(FDmetrics_FRic$value[1:3500])
shapiro.test(FDmetrics_FRic$value[3501:7000])
shapiro.test(FDmetrics_FSpe$value[1:3500])
shapiro.test(FDmetrics_FSpe$value[3501:7000])
shapiro.test(FDmetrics_FUn$value[1:3500]) 
shapiro.test(FDmetrics_FUn$value[3501:7000])              # Non-normal; Use median; all points are the same under mode

# Z-score calculations
# Calculate Z scores between empirical analyses and null model per scenario
# Form empty dataframe
emp.null<- as.data.frame(matrix(data= NA,nrow= 7, ncol= 10, dimnames= list(c("Present","Future","Y100","Y200","Y300","Y400","Y500"),
                                                                           c("Sp_richn","fric","fori","fspe","fun",
                                                                             "Sp_richn_Z","fric_Z","fori_Z","fspe_Z","fun_Z"))))
# Separate by epoch for each analysis
## Empirical 
FDmetrics_taxonvar_P<- Res_FDmetrics_TaxonVar[grep("Present", Res_FDmetrics_TaxonVar$Scenario), (2:6)]
FDmetrics_taxonvar_F<- Res_FDmetrics_TaxonVar[grep("Future", Res_FDmetrics_TaxonVar$Scenario), (2:6)]
FDmetrics_taxonvar_100<- Res_FDmetrics_TaxonVar[grep("Y100", Res_FDmetrics_TaxonVar$Scenario), (2:6)]
FDmetrics_taxonvar_200<- Res_FDmetrics_TaxonVar[grep("Y200", Res_FDmetrics_TaxonVar$Scenario), (2:6)]
FDmetrics_taxonvar_300<- Res_FDmetrics_TaxonVar[grep("Y300", Res_FDmetrics_TaxonVar$Scenario), (2:6)]
FDmetrics_taxonvar_400<- Res_FDmetrics_TaxonVar[grep("Y400", Res_FDmetrics_TaxonVar$Scenario), (2:6)]
FDmetrics_taxonvar_500<- Res_FDmetrics_TaxonVar[grep("Y500", Res_FDmetrics_TaxonVar$Scenario), (2:6)]

## Null model
Null_FDmetrics_taxonvar_P<- Null_FDmetrics_taxonvar[grep("Present", Null_FDmetrics_taxonvar$Scenario), (2:6)]
Null_FDmetrics_taxonvar_F<- Null_FDmetrics_taxonvar[grep("Future", Null_FDmetrics_taxonvar$Scenario), (2:6)]
Null_FDmetrics_taxonvar_100<- Null_FDmetrics_taxonvar[grep("Y100", Null_FDmetrics_taxonvar$Scenario), (2:6)]
Null_FDmetrics_taxonvar_200<- Null_FDmetrics_taxonvar[grep("Y200", Null_FDmetrics_taxonvar$Scenario), (2:6)]
Null_FDmetrics_taxonvar_300<- Null_FDmetrics_taxonvar[grep("Y300", Null_FDmetrics_taxonvar$Scenario), (2:6)]
Null_FDmetrics_taxonvar_400<- Null_FDmetrics_taxonvar[grep("Y400", Null_FDmetrics_taxonvar$Scenario), (2:6)]
Null_FDmetrics_taxonvar_500<- Null_FDmetrics_taxonvar[grep("Y500", Null_FDmetrics_taxonvar$Scenario), (2:6)]

# Make empirical dataframes - means or medians of empirical & null models
FDemp <- Res_FDmetrics_TaxonVar %>% 
  group_by(Scenario) %>%
  summarise(Sp_med = median(sp_richn),
            FRic_med = median(fric),
            Fori_med = median(fori),
            Fspe_med = median(fspe),
            FUn_med = median(fun))
FDemp <- FDemp %>% 
  column_to_rownames(var = "Scenario")

FDnull <- Null_FDmetrics_taxonvar %>% 
  group_by(Scenario) %>%
  summarise(Sp_med = median(sp_richn),
            FRic_med = median(fric),
            Fori_med = median(fori),
            Fspe_med = median(fspe),
            FUn_med = median(fun))
FDnull <- FDnull %>% 
  column_to_rownames(var = "Scenario")

# Calculate empirical differences; use FDemp as scenarios are rownames rather than a column
for(e in 1:5) { 
  for (g in 1:7){
    emp.null[g,e]<- as.numeric(FDemp[g,e])-as.numeric(FDnull[g,e])
    
  }
}


# write function to calculate slopes for null model values
slopes_fun<- function(data1, data2){
  output<- matrix(data= NA, nrow= dim(data1)[1], ncol= dim(data1)[2])
  for (k in 1:5) {
    output[,k]<- data2[,k]- data1[,k]
  }
  return(output)
}

P_null_Slopes<- slopes_fun(FDmetrics_taxonvar_P, Null_FDmetrics_taxonvar_P)
F_null_Slopes<- slopes_fun(FDmetrics_taxonvar_F, Null_FDmetrics_taxonvar_F)
F_null_Slopes_100<- slopes_fun(FDmetrics_taxonvar_100, Null_FDmetrics_taxonvar_100)
F_null_Slopes_200<- slopes_fun(FDmetrics_taxonvar_200, Null_FDmetrics_taxonvar_200)
F_null_Slopes_300<- slopes_fun(FDmetrics_taxonvar_300, Null_FDmetrics_taxonvar_300)
F_null_Slopes_400<- slopes_fun(FDmetrics_taxonvar_400, Null_FDmetrics_taxonvar_400)
F_null_Slopes_500<- slopes_fun(FDmetrics_taxonvar_500, Null_FDmetrics_taxonvar_500)

# Calculate Z scores
for (i in 1:5){
  emp.null[1,i+5]<- (emp.null[1,i]- median(P_null_Slopes[,i], na.rm = TRUE))/sd(P_null_Slopes[,i], na.rm= TRUE)}

for (i in 1:5){
  emp.null[2,i+5]<- (emp.null[2,i]- median(F_null_Slopes[,i],  na.rm = TRUE))/sd(F_null_Slopes[,i],  na.rm = TRUE)}

for (i in 1:5){
  emp.null[3,i+5]<- (emp.null[3,i]- median(F_null_Slopes_100[,i],  na.rm = TRUE))/sd(F_null_Slopes_100[,i],  na.rm = TRUE)}

for (i in 1:5){
  emp.null[4,i+5]<- (emp.null[4,i]- median(F_null_Slopes_200[,i],  na.rm = TRUE))/sd(F_null_Slopes_200[,i],  na.rm = TRUE)}

for (i in 1:5){
  emp.null[5,i+5]<- (emp.null[5,i]- median(F_null_Slopes_300[,i],  na.rm = TRUE))/sd(F_null_Slopes_300[,i],  na.rm = TRUE)}

for (i in 1:5){
  emp.null[6,i+5]<- (emp.null[6,i]- median(F_null_Slopes_400[,i],  na.rm = TRUE))/sd(F_null_Slopes_400[,i],  na.rm = TRUE)}

for (i in 1:5){
  emp.null[7,i+5]<- (emp.null[7,i]- median(F_null_Slopes_500[,i],  na.rm = TRUE))/sd(F_null_Slopes_500[,i],  na.rm = TRUE)}


emp.null

# Extract results as excel worksheet - produces Table SX
save(emp.null, file = "C:/Users/2022207/Dropbox/Jack's PhD/Chapter 3. Future shark FD/Data/R data/Mode_Future Z-scores.Rdata")
write_xlsx(emp.null, "C:/Users/2022207/Dropbox/Jack's PhD/Chapter 3. Future shark FD/Data/R data/Mode_Future Z-scores.xlsx")
