# R code for the analysis of: Future impacts of overfishing and climate change on elasmobranch functional diversity

## Authors
Jack A. Cooper, Gregor H. Mathes, John N. Griffin, Fabien Leprieur, Kristin Kaschner, Kathleen Kesner-Reyes & Catalina Pimiento

## Introduction
This folder contains R scripts for all analyses exploring the future functional diversity of sharks and rays. 

The following packages are required for:
- Data manipulation & handling: ```readxl```, ```tidyverse```, ```janitor```, ```rredlist```, ```reticulate```,```tibble```, ```here```, ```doBy```, ```reshape2```, and ```data.table```.
- Analyses: ```purrr```, ```devtools```, ```iucnsim```, ```progress```, ```mFD```, ```foreach```, ```parallel```, ```doParallel```, and ```geometry```.
- Visualisation: ```scales```, ```cowplot```, ```RColorBrewer```, ```ggsci```, ```viridis```, ```ggplotify```, and ```wesanderson```.

These packages can be downloaded and installed using the following commands:
``` {r}
install.packages("readxl")
install.packages("tidyverse")
install.packages("janitor")
install.packages("rredlist")
install.packages("purrr")
install.packages("devtools")
install.packages("reticulate")
install.packages("iucnsim")
install.packages("tibble")
install.packages("here")
install.packages("progress")
install.packages("scales")
install.packages("cowplot")
install.packages("doBy")
install.packages("reshape2")
install.packages("RColorBrewer")
install.packages("mFD")
install.packages("ggsci")
install.packages("viridis")
install.packages("data.table")
install.packages("foreach")
install.packages("parallel")
install.packages("doParallel")
install.packages("ggplotify")
install.packages("geometry")
install.packages("wesanderson")
```
## Folders
All code and data are placed in the following folders:

Code: containing all R scripts

Data: containing Rdata and xlsx inputs and outputs from the code. Rdata files can be loaded directly into replications of our analyses using the load() command

## Codes
A total of 11 R scripts were produced containing all analyses. The scripts are as follows:

01. **Fix Synonyms**: This code updates IUCN status of all elasmobranchs and corrects for synonyms
02. **IUCN_sim**: This code simulates the future IUCN status, and therefore extinction times, of all elasmobranchs using the iucn_sim program
03. **Empirical analyses**: This code calculates functional diversity in our extinction scenarios based on our simulations
04. **Null analyses**: This code calculates null functional diversity in our extinction scenarios based on our simulations
05. **Median functional space**: This code forms the functional spaces of our extinction scenarios - producing Figure S1 and S2
06. **Figure 1**: This code visualises changes in functional diversity in our extinction scenarios following our simulations - producing Figure 1
07. **AquaMaps functional space**: This code calculates the functional diversity occupation of AquaMaps - producing Figure S3
08. **AquaMaps FD analyses**: This code analyses global functional diversity on the spatial scale from AquaMaps
09. **Map figures**: This code visualises all map-based analyses - producing the map components of Figure 1 and Figures 2 and S4
10. **AquaMaps supplementary FD analyses**: This code analyses global functional diversity on alternative climate change scenarios
11. **Supplementary map figures**: This code visualises supplementary map-based analyses - produces Figures S5-S7

