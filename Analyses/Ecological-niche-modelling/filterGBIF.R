# Filter GBIF occurrences of black and blue wildebeest
# Author: Deon de Jager

## Based on this post: https://data-blog.gbif.org/post/gbif-filtering-guide/ by John Waller.

#-- GBIF download citations --#
# Black wildebeest (Connochaetes gnou): GBIF.org (24 July 2024) GBIF Occurrence Download  https://doi.org/10.15468/dl.wm5x3f
# Blue wildebeest (Connochaetes taurinus): GBIF.org (24 July 2024) GBIF Occurrence Download  https://doi.org/10.15468/dl.bahh78

#-- Load packages
library(rgbif)
library(tidyverse)
#install.packages("CoordinateCleaner")
library(CoordinateCleaner)

#-- Set working directory
#setwd("path/to/occurrenceData")

#-- Import occurrence data
## Note: Some basic filtering was done on GBIF and the occurences downloaded as a simple CSV, which is actually a TSV (tab-separated)
## GBIF filters: 
### Black wildebeest: "Species = Connochaetes gnou (Zimmermann, 1780)" AND "Continent = AFRICA" AND "Has coordinate = true" AND "Has geospatial issue = false" AND "Occurrence status = present"
### Blue wildebeest: "Species = Connochaetes taurinus (Burchell, 1823)" AND "Continent = AFRICA" AND "Has coordinate = true" AND "Has geospatial issue = false" AND "Occurrence status = present"

# Read TSV
Ct_occ <- read_tsv("Blue_wildebeest_0034647-240626123714530_GBIF.tsv", col_names = TRUE, na = "")
Cg_occ <- read_tsv("Black_wildebeest_0034563-240626123714530_GBIF.tsv", col_names = TRUE, na = "")
## Namibia's country code is NA, which is by default treated as missing data, thus have to specify that only empty entries are missing data.

# Set lowercase column names to work with CoordinateCleaner
Cg_occ <- Cg_occ %>%
  setNames(tolower(names(.)))
Ct_occ <- Ct_occ %>%
  setNames(tolower(names(.)))

#-- Filter occurrences
## Black wildebeest
Cg_occ_filter <- Cg_occ %>% 
  filter(countrycode %in% c("ZA", "SZ", "LS")) %>% # keep only records from South Africa, eSwatini, and Lesotho
  filter(establishmentmeans != "introduced" | is.na(establishmentmeans)) %>% # remove records from introduced populations (if recorded), but keep missing values
  filter(!basisofrecord %in% c("FOSSIL_SPECIMEN", "MATERIAL_CITATION", "LIVING_SPECIMEN")) %>% # remove fossils and living specimens (often in zoos)
  filter(coordinateprecision < 0.01 | is.na(coordinateprecision)) %>% # keeps missing values
  filter(coordinateuncertaintyinmeters < 10000 | is.na(coordinateuncertaintyinmeters)) %>% # keeps missing values
  filter(!coordinateuncertaintyinmeters %in% c(301,3036,999,9999)) %>% # common default values that are inaccurate
  filter(!decimallatitude == 0 | !decimallongitude == 0) %>% # remove points along prime meridian or equator
  cc_cen(lon = "decimallongitude", lat = "decimallatitude", buffer = 2000) %>% # remove country centroids within 2km 
  cc_cap(lon = "decimallongitude", lat = "decimallatitude", buffer = 2000) %>% # remove capitals centroids within 2km
  cc_inst(lon = "decimallongitude", lat = "decimallatitude", buffer = 2000) %>% # remove zoo and herbaria within 2km 
  #cc_sea(lon = "decimallongitude", lat = "decimallatitude", value = "flagged") %>% # Flag records from ocean. Check output for accuracy (tends to remove those close to ocean). Remove manually if in ocean.
  #distinct(decimallongitude,decimallatitude,specieskey,datasetkey, .keep_all = TRUE) %>% # Not done here, as spatial thinning is done later.
  glimpse() # look at results of pipeline

# Write csv
write_csv(Cg_occ_filter, file = "Black_wildebeest_0034563-240626123714530_GBIF_filterR.csv", na = "")

## Blue wildebeest
Ct_occ_filter <- Ct_occ %>% 
  filter(establishmentmeans != "introduced" | is.na(establishmentmeans)) %>% # remove records from introduced populations (if recorded), but keep missing values
  filter(collectioncode != "MammalsAngola" | is.na(collectioncode)) %>% # these records are a gridded dataset. Thus remove, but keep missing values. 
  filter(!basisofrecord %in% c("FOSSIL_SPECIMEN", "MATERIAL_CITATION", "LIVING_SPECIMEN", "OCCURRENCE", "MACHINE_OBSERVATION", "OBSERVATION")) %>% # remove fossils, living specimens (often in zoos), the invalid basis "OCCURRENCE", and the vague "OBSERVATION". This leaves "HUMAN_OBSERVATION", "PRESERVED_SPECIMEN", and "MATERIAL_SAMPLE" (DNA-based).
  filter(coordinateprecision < 0.01 | is.na(coordinateprecision)) %>% # keeps missing values
  filter(coordinateuncertaintyinmeters < 10000 | is.na(coordinateuncertaintyinmeters)) %>% # keeps missing values
  filter(!coordinateuncertaintyinmeters %in% c(301,3036,999,9999)) %>% # common default values that are inaccurate
  filter(!decimallatitude == 0 | !decimallongitude == 0) %>% # remove points along prime meridian or equator
  cc_cen(lon = "decimallongitude", lat = "decimallatitude", buffer = 2000) %>% # remove country centroids within 2km 
  cc_cap(lon = "decimallongitude", lat = "decimallatitude", buffer = 2000) %>% # remove capitals centroids within 2km
  cc_inst(lon = "decimallongitude", lat = "decimallatitude", buffer = 2000) %>% # remove zoo and herbaria within 2km 
  #cc_sea(lon = "decimallongitude", lat = "decimallatitude", value = "flagged") #%>% # Flag records from ocean. Check output for accuracy (tends to remove those close to ocean). Remove manually if in ocean.
  #distinct(decimallongitude,decimallatitude,specieskey,datasetkey, .keep_all = TRUE) %>% # Not done here, as spatial thinning is done later.
  glimpse() # look at results of pipeline

# Write csv
write_csv(Ct_occ_filter, file = "Blue_wildebeest_0034647-240626123714530_GBIF_filterR.csv", na = "")

# Resulting occurrences were then manually filtered based on several criteria, but predominantly whether they likely represent an introduced population outside the natural range of the species. See manuscript supplementary text for more detail.