# R code for the analysis of: Overfishing extinguishes climate change refugia for elasmobranch functional diversity

## Authors
Jack A. Cooper, Gregor H. Mathes, John N. Griffin, Kristin Kaschner, Kathleen Kesner-Reyes, Fabien Leprieur & Catalina Pimiento

## Introduction
This folder contains R scripts for all analyses exploring the future functional diversity of sharks and rays. 

The following packages are required for:
- Data manipulation & handling: ```readxl```, ```tidyverse```, ```janitor```, ```rredlist```, ```reticulate```,```tibble```, ```here```, ```doBy```, ```reshape2```, and ```data.table```.
- Analyses: ```purrr```, ```devtools```, ```iucnsim```, ```progress```, ```mFD```, ```foreach```, ```parallel```, ```doParallel```, and ```geometry```.
- Visualisation: ```scales```, ```cowplot```, ```RColorBrewer```, ```ggsci```, ```viridis```, ```ggplotify```, ```ggthemes```, and ```wesanderson```.

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
install.packages("ggthemes")
install.packages("geometry")
install.packages("wesanderson")
```
## Folders
All code and data are placed in the following folders:

Code: containing all R scripts

Data: containing Rdata and xlsx inputs and outputs from the code. Rdata files can be loaded directly into replications of our analyses using the load() command

## Codes
A total of 13 R scripts were produced containing all analyses. A sub-folder of codes for individual trait-removal sensitivity tests is also included. The scripts in the main folder are as follows:

01. **Fix Synonyms**: This code updates IUCN status of all elasmobranchs and corrects for synonyms
02. **IUCN_sim**: This code simulates the future IUCN status, and therefore extinction times, of all elasmobranchs using the iucn_sim program
03. **Empirical analyses**: This code calculates functional diversity in our extinction scenarios based on our simulations
04. **Null analyses**: This code calculates null functional diversity in our extinction scenarios based on our simulations
05. **Median functional spaces**: This code forms the functional spaces of our extinction scenario - producing Figure S3 and S4
06. **Figure 1**: This code visualises changes in functional diversity following our simulations - producing Figure 1A, D
07. **Sensitivity tests figure**: This code visualises our sensitivity tests for trait selection - producing Figure S1
08. **Overfitting tests**: This code visualises changes in functional diversity following our simulations - producing Figure S2
09. **AquaMaps functional space**: This code calculates the functional diversity occupation of AquaMaps - producing Figure S5
10. **AquaMaps FD analyses**: This code analyses global functional diversity on the spatial scale from AquaMaps
11. **Map figures**: This code visualises all map-based analyses - producing Figure 1, 2 and S6
12. **AquaMaps supplementary FD analyses**: This code analyses global functional diversity on alternative climate change scenarios
13. **Supplementary map figures**: This code visualises supplementary map-based analyses - produces Figures S7-S9

