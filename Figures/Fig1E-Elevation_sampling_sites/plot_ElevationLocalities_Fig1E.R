# Plot elevation of sample localities
## See the README.md file in this folder for a description of how the input data were prepared.
## Produces Fig. 1E

# Set working directory
setwd("./")

# Load packages
library(tidyverse)

# Import data
sites <- read_csv("Wildebeest_Localities_metadata_Elevation.csv", col_names = T)

## Lock in order of sites for ggplot
#sites$Locality <- factor(sites$Locality, levels = sites$Locality)

# Colour palette
myPalette <- c("#f2933b", "#cace3d", "#8dc5a6", "#cc95ff", "#ff96ca", "#8585c4", "#5c5c5c")

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
ggsave("Genome_Samples_LocalitiesByTaxon_Elevation.svg")
