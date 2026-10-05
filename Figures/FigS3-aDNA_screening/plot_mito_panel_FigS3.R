# Plot screening data mapped to mito_panel for species ID
## Produces fig. S3

# Set working directory
setwd("./")

## Load packages
library(tidyverse)
library(reshape2)
library(writexl)
library(ggforce)
library(ggh4x) # to duplicate discrete axis in ggplot (currently only natively available for continuous scales)
library(grid)

# Import data
screening <- read_csv("FinchaHabera_mito_panel_summary_long_TableS5.csv", col_names = T)

# Create plotting order (Bovidae ordered according to Bibi 2013 tree from top to bottom, then the rest with approx. increasing distance from Bovidae)
screening$reference <- factor(screening$reference, levels = c(NA, "ElandCommon", "KuduGreater", "KuduLesser", "NyalaMountain", "BuffaloCape", "Impala", "RoanAntelope",
                                                              "Hartebeest", "WildebeestBlue", "Aoudad", "Goat", "IbexNubian", "IbexWalia", "ReedbuckMt", "Waterbuck",
                                                              "Klipspringer", "DuikerCommon", "DikdikSalts", "Gerenuk", "GazelleGrants", "Giraffe", "Hippopotamus", 
                                                              "Bushpig", "Warthog", "ZebraPlains", "EthiopianWolf","Leopard", "Hyena", "Human", "BaboonHamadryas", 
                                                              "BaboonOlive", "RootRatGiant", "Elephant"))
screening$SampleID <- factor(screening$SampleID, levels = c("Fossil936", "B1ancient", "B2ancient", "LBancient", "IBancient"))
screening$Family <- factor(screening$Family, levels = c(NA, "Bovidae", "Giraffidae", "Hippopotamidae", "Suidae", "Equidae", "Canidae", "Felidae", "Hyaenidae",
                                                        "Hominidae", "Cercopithecidae", "Spalacidae", "Elephantidae"))


# Plot ----
## Plotting number of unique (i.e. duplicates removed) reads mapping to each reference mitogenome
### All samples, facet wrap
ggplot(filter(screening, variable == "hits_unique")) +
  geom_col(aes(x = reference, y = value, fill = Family)) +
  geom_text(aes(x = reference, y = value, label = value), size = 2.5, hjust = -0.2, angle = 90) +
  labs(x = "Reference mitogenome", y = "No. of unique reads mapped") +
  coord_cartesian(clip = "off") +
  facet_wrap(~SampleID + CGG, scales = "free_y", ncol = 2, nrow = 5) +
  theme_minimal(base_size = 12) +
  theme(strip.text = element_text(face = "bold", size = 10),
        legend.position = "top",
        axis.text.x = element_text(size = 6)) +
  guides(fill = guide_legend(nrow = 2)) +
  scale_x_discrete(guide = guide_axis(angle = 45))
ggsave("panel_FigS5.png", width = 205, height = 292, units = "mm", dpi = 300)
ggsave("panel_FigS5.pdf", width = 205, height = 292, units = "mm")
