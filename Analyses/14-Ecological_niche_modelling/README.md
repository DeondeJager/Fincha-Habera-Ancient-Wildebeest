# Code for ecological niche models
Here we provide the code to run the models, but not scripts to prepare the occurrence and climate input files, which are quite large (~3 Gb compressed, ~10 Gb uncompressed). 

These are instead available, with a detailed README as Supplementary Data S4 on Zenodo: https://doi.org/10.5281/zenodo.18865979.

## Scripts
- `WildebeestBlueENM.R`: Run the blue wildebeest ENMs. Author: Nicholas A. Freymueller
- `WildebeestBlackENM.R`: Run the black wildebeest ENMs. Author: Nicholas A. Freymueller
- `NAF_Functions.R`: Functions required by the two ENM scripts. Author: Nicholas A. Freymueller
- `filterGBIF.R`: Filter the initial occurrences downloaded from GBIF. Author: Deon de Jager
  - Processes files `Blue_wildebeest_0034647-240626123714530_GBIF.tsv` and `Black_wildebeest_0034563-240626123714530_GBIF.tsv`.
