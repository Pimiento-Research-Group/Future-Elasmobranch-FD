# R code for the analysis of: The future of elasmobranch functional diversity under a changing climate and simulated extinctions

## Authors
Jack A. Cooper, Gregor H. Mathes & Catalina Pimiento

## Introduction
This folder contains R scripts for all analyses exploring the future functional diversity of sharks and rays in 2100 and beyond. 

The following packages are required for:
- Data manipulation & handling: ```readxl```, ```tidyverse```, ```janitor```, ```rredlist```, ```tibble```, ```here```, ```doBy```, ```reshape2```, ```data.table```, and ```writexl```.
- Analyses: ```purrr```, ```devtools```, ```iucnsim```, ```progress```, ```furrr```, ```mFD```, ```foreach```, ```parallel```, ```doParallel```, ```funrar```, ```ade4```, and ```lme4```.
- Visualisation: ```scales```, ```cowplot```, ```RColorBrewer```, ```ggsci```, ```ggplotify```, ```ggthemes```, ```wesanderson```, ```colourspace```, ```visreg```, ```plotly```, ```gridExtra``` and ```ggpubr```.

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
install.packages("furrr")
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
install.packages("funrar")
install.packages("ggplotify")
install.packages("ade4")
install.packages("writexl")
install.packages("ggthemes")
install.packages("wesanderson")
install.packages("colorspace")
install.packages("visreg")
install.packages("lme4")
install.packages("plotly")
install.packages("gridExtra")
install.packages("ggpubr")
```

## Codes
A total of 26 R scripts were produced containing all analyses. The scripts are as follows:

01. **Fix Synonyms**: This code updates IUCN status of all elasmobranchs and corrects for synonyms
02. **IUCN simulations_final**: This code simulates the future IUCN status, and therefore extinction times, of all elasmobranchs using the iucn_sim program
03. **Functional space**: This code creates the functional space & assesses functional extent of AquaMaps
04. **Future functional spaces**: This code forms future functional spaces based on modal extinction times from IUCN simulations
05. **Future FUn-FSp scores**: This code forms future functional spaces based on functional uniqueness (FUn) and functional specialisation (FSp) scores
06. **Null model**: This code calculates null distribution of future functional diversity analyses; adapted from code run in the supercomputer
07. **Future FD plotting**: This code plots future functional diversity (FD) based on simulated IUCN-based extinctions - producing Figure 1 and Table 1
08. **Future FD Z-scores**: This code calculates if changes in functional diversity are greater than expected by random chance
09. **AquaMaps present-day analyses**: This code analyses global functional diversity in the present
10. **AquaMaps RCP 2.6 analyses**: This code analyses global functional diversity under the RCP 2.6 scenario
11. **AquaMaps RCP 4.5 analyses**: This code analyses global functional diversity under the RCP 4.5 scenario
12. **AquaMaps RCP 8.5 analyses**: This code analyses global functional diversity under the RCP 8.5 scenario
13. **AquaMaps present-day maps**: This code plots the global distribution of present-day functional diversity
14. **RCP 2.6 maps**: This code plots the functional diversity shift under RCP 2.6 climate change
15. **RCP 2.6 extinction maps**: This code plots the functional diversity shift under RCP 2.6 climate change while accounting for mode extinction times recorded by IUCN_sim
16. **RCP 4.5 maps**: This code plots the functional diversity shift under RCP 4.5 climate change
17. **RCP 4.5 extinction maps**: This code plots the functional diversity shift under RCP 4.5 climate change while accounting for mode extinction times recorded by IUCN_sim
18. **RCP 8.5 maps**: This code plots the functional diversity shift under RCP 8.5 climate change
19. **RCP 8.5 extinction maps**: This code plots the functional diversity shift under RCP 8.5 climate change while accounting for mode extinction times recorded by IUCN_sim
20. **AquaMaps shift maps**: This code plots all shift maps from AquaMaps - producing Figures 2 and 3
21. **Trait and environmental GLMs**: This code uses GLMs to analyses relationships between traits & extinction; and environmental parameters & occurrence in AquaMaps
22. **Environmental drivers of RCP 2.6 FD**: This code analyses the relationship between environmental AquaMaps data and FD in 2100 under RCP 2.6 climate change and extinctions
23. **Environmental drivers of RCP 4.5 FD**: This code analyses the relationship between environmental AquaMaps data and FD in 2100 under RCP 4.5 climate change and extinctions
24. **Environmental drivers of RCP 8.5 FD**: This code analyses the relationship between environmental AquaMaps data and FD in 2100 under RCP 8.5 climate change and extinctions
25. **Extinction GLMs**: This code assesses the relationship between ranked FUSE species and simulated extinction - producing Figure 4
26. **FUSE ranking extinctions**: This code identifies the simulated species extinctions expected to have the most adverse effects on functional diversity in 2100
