# Plot elevation profile

# Set working directory
setwd("C:/Users/pzx702/Documents/MSCAFellowship2021/Manuscripts/2022_Wildebeest_AncientEastAfrica/Map")

# Load packages
library(tidyverse)

# Import data
profile <- read_csv("Wildebeest_n132_metadata_ElevationProfile_Filtered.csv", col_names = T)
samples <- read_csv("Wildebeest_n132_metadata_ElevationProfile_GenomeSamples.csv", col_names = T)
sites <- read_csv("Wildebeest_Localities_metadata_Elevation.csv", col_names = T)
## Lock in order of sites for ggplot
sites$Locality <- factor(sites$Locality, levels = sites$Locality)

# Colour palette
myPalette <- c("#f2933b", "#cace3d", "#8dc5a6", "#cc95ff", "#ff96ca", "#8585c4", "#5c5c5c")

# Plot elevation profile
ggplot() +
  geom_line(data = profile, aes(distance/1000, elevation), colour = "darkgrey", linewidth = 0.5) +
  geom_point(data = samples, aes(distance/1000, elevation, fill = Subspecies_common), pch = 21, size = 5, stroke = 1.5) +
  scale_fill_manual(values = myPalette,
                      breaks = c("Bale wildebeest", "Western white-bearded", "Eastern white-bearded", "Nyassa", "Cookson", "Brindled", "Black"),
                      name = "Taxon") +
  #scale_x_continuous(name = "Distance (km)",
  #                   labels = scales::label_comma()) +
  scale_x_reverse(name = "North to South path through sampling localities",
                  breaks = NULL) +
  scale_y_continuous(name = "Elevation (m a.s.l.)",
                     labels = scales::label_comma()) +
  theme_bw(base_size = 18) +
  theme(legend.position = "none")
ggsave("Elevation_profile.png", dpi = 300, width = 290, height = 100, units = "mm")
ggsave("Elevation_profile.pdf", width = 290, height = 100, units = "mm")
ggsave("Elevation_profile.svg", width = 290, height = 100, units = "mm")

# Plot elevation of sites only
ggplot(sites) +
  geom_point(aes(Locality, Elevation, fill = Subspecies_common), pch = 21, size = 5, stroke = 0.7) +
  scale_fill_manual(values = myPalette,
                    breaks = c("Bale wildebeest", "Western white-bearded", "Eastern white-bearded", "Nyassa", "Cookson", "Brindled", "Black"),
                    name = "Taxon") +
  scale_y_continuous(name = "Elevation (m a.s.l.)",
                     labels = scales::label_comma()) +
  scale_x_discrete(name = "Localities (north to south)",
                   labels = NULL,
                   breaks = NULL) +
  theme_bw(base_size = 24) +
  theme(legend.position = "none", 
        axis.ticks.x = element_blank())
ggsave("Genome_Samples_Localities_Elevation.png", dpi = 300)
ggsave("Genome_Samples_Localities_Elevation_narrow.png", dpi = 300, width = 290, height = 100, units = "mm")
ggsave("Genome_Samples_Localities_Elevation_narrow.svg", width = 290, height = 100, units = "mm")

# Plot elevation of sites by subspecies
## Set order of subspecies for plotting
sites$Subspecies_common <- factor(sites$Subspecies_common, levels = c("Bale wildebeest", "Western white-bearded", "Eastern white-bearded", "Nyassa", "Cookson", "Brindled", "Black"))
ggplot(sites) +
  geom_point(aes(Subspecies_common, Elevation, fill = Subspecies_common), pch = 21, size = 5, stroke = 0.7) +
  scale_fill_manual(values = myPalette,
                    breaks = c("Bale wildebeest", "Western white-bearded", "Eastern white-bearded", "Nyassa", "Cookson", "Brindled", "Black"),
                    name = "Taxon") +
  scale_y_continuous(name = "Elevation (m a.s.l.)",
                     labels = scales::label_comma()) +
  scale_x_discrete(name = "Localities by taxon",
                   labels = NULL,
                   breaks = NULL) +
  theme_bw(base_size = 24) +
  theme(legend.position = "none", 
        axis.ticks.x = element_blank())
ggsave("Genome_Samples_LocalitiesByTaxon_Elevation.png", dpi = 300)
ggsave("Genome_Samples_LocalitiesByTaxon_Elevation_narrow.png", dpi = 300, width = 290/2, height = 100, units = "mm")
ggsave("Genome_Samples_LocalitiesByTaxon_Elevation_narrow.svg", width = 290/2, height = 100, units = "mm")
