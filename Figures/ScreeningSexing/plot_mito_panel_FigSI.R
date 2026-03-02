# Plot screening data mapped to mito_panel for species ID

# Set working directory
#setwd("")

## Load packages
library(tidyverse)
library(reshape2)
library(writexl)
library(ggforce)
library(ggh4x) # to duplicate discrete axis in ggplot (currently only natively available for continuous scales)
library(grid)

# Import data
screening <- read_csv("FinchaHabera_mito_panel_summary_long_TableSI.csv", col_names = T)

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
  #  geom_vline(xintercept = 9.5, linetype = "dashed") +
  #  geom_vline(xintercept = 7.5, linetype = "dashed") +
  coord_cartesian(clip = "off") +
  facet_wrap(~SampleID + CGG, scales = "free_y", ncol = 2, nrow = 5) +
  theme_minimal(base_size = 12) +
  theme(strip.text = element_text(face = "bold", size = 10),
        legend.position = "top",
        axis.text.x = element_text(size = 6)) +
  guides(fill = guide_legend(nrow = 2)) +
  scale_x_discrete(guide = guide_axis(angle = 45))
  ggtheme_classic2()
ggsave("panel_FigSI.png", width = 205, height = 292, units = "mm", dpi = 300)
ggsave("panel_FigSI.pdf", width = 205, height = 292, units = "mm")

### Fossil936 only
ggplot(filter(screening, variable == "hits_unique", SampleID == "Fossil936")) +
  geom_col(aes(x = reference, y = value, fill = Family)) +
  scale_fill_discrete(guide = guide_legend(reverse = TRUE)) +
  #geom_text(aes(x = reference, y = value+150, label = value), size = 5) +
  labs(x = "Reference mitogenome", y = "No. of unique reads mapped") +
  theme_minimal(base_size = 14) +
  scale_y_continuous(limits = c(0, 3000),
                     labels = scales::label_comma()) +
  coord_flip()
ggsave("Fossil936.png", width = 7.5, height = 5, units = "in", dpi = 300)
ggsave("Fossil936.pdf", width = 7.5, height = 5, units = "in")
