# R code for the analysis of: Temporal and spatial decline in functional diversity of sharks and rays under simulated extinctions

## Authors
Jack A. Cooper, Gregor H. Mathes, John N. Griffin, Fabien Leprieur, Kristin Kaschner, Kathleen Kesner-Reyes & Catalina Pimiento

## Introduction
This folder contains R scripts for all analyses exploring the future functional diversity of sharks and rays. 

The following packages are required for:
- Data manipulation & handling: ```readxl```, ```tidyverse```, ```janitor```, ```rredlist```, ```reticulate```,```tibble```, ```here```, ```doBy```, ```reshape2```, and ```data.table```.
- Analyses: ```purrr```, ```devtools```, ```iucnsim```, ```progress```, ```mFD```, ```foreach```, ```parallel```, ```doParallel```, ```geometry```, ```ade4```, and ```lme4```.
- Visualisation: ```scales```, ```cowplot```, ```RColorBrewer```, ```ggsci```, ```viridis```, ```ggplotify```, ```ggthemes```, ```wesanderson```, ```colourspace```, and ```visreg```.

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
install.packages("ade4")
install.packages("geometry")
install.packages("ggthemes")
install.packages("wesanderson")
install.packages("colorspace")
install.packages("visreg")
install.packages("lme4")
```
## Folders
All code and data are placed in the following folders:

Code: containing all R scripts

Data: containing Rdata and xlsx inputs and outputs from the code. Rdata files can be loaded directly into replications of our analyses using the load() command

## Codes
A total of 20 R scripts were produced containing all analyses. The scripts are as follows:

01. **Fix Synonyms**: This code updates IUCN status of all elasmobranchs and corrects for synonyms
02. **IUCN_sim**: This code simulates the future IUCN status, and therefore extinction times, of all elasmobranchs using the iucn_sim program
03. **Empirical analyses**: This code calculates functional diversity in our extinction scenarios based on our simulations
04. **Null analyses**: This code calculates null functional diversity in our extinction scenarios based on our simulations
05. **Median functional space**: This code forms the functional spaces of our extinction scenarios - producing Figure S1 and S3
06. **Figure 1**: This code visualises changes in functional diversity in our extinction scenarios following our simulations - producing Figure 1
07. **Median GLMMs**: This code uses GLMMs to assess if and when the most distinct species are expected to be lost in our scenarios - producing Figure S2 and Table 1
08. **AquaMaps functional space**: This code calculates the functional diversity occupation of AquaMaps - producing Figure S4
09. **AquaMaps present-day analyses**: This code analyses global functional diversity in the present
10. **AquaMaps RCP 2.6 analyses**: This code analyses global functional diversity under the RCP 2.6 scenario
11. **AquaMaps RCP 4.5 analyses**: This code analyses global functional diversity under the RCP 4.5 scenario
12. **AquaMaps RCP 8.5 analyses**: This code analyses global functional diversity under the RCP 8.5 scenario
13. **AquaMaps present-day maps**: This code plots the global distribution of present-day functional diversity - producing Figure S5
14. **RCP 2.6 maps**: This code plots the functional diversity shift under RCP 2.6 climate change
15. **RCP 2.6 extinction maps**: This code plots the functional diversity shift under RCP 2.6 climate change while accounting for median extinction times recorded by IUCN_sim
16. **RCP 4.5 maps**: This code plots the functional diversity shift under RCP 4.5 climate change
17. **RCP 4.5 extinction maps**: This code plots the functional diversity shift under RCP 4.5 climate change while accounting for median extinction times recorded by IUCN_sim
18. **RCP 8.5 maps**: This code plots the functional diversity shift under RCP 8.5 climate change
19. **RCP 8.5 extinction maps**: This code plots the functional diversity shift under RCP 8.5 climate change while accounting for median extinction times recorded by IUCN_sim
20. **AquaMaps shift maps**: This code plots all shift maps from AquaMaps - producing Figures 2, S6 and S7
