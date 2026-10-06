# Plot emu pca output (pca based on pseudohap files)
## NOte: The pseudohap files exclude transitions and were filtered for maf=0.05
## Produces basis for Fig. S8C

# Set working directory
setwd("./")

# Load packages
library(tidyverse)
library(janitor)
library(scales)
library(ggpubr)

# Load data ----
metadata <- read_csv("../metadata.csv", col_names = TRUE)
## n82_fossil - fossil + modern (excl. SRR27955277_trim_Mut), blue and black wildebeest ----
eigenvec_n82_fossil <- read_table("n82_fossil_Goat_autosomes_noHart.emu.eigvecs", col_names = T)
eigenvec_n82_fossil <- eigenvec_n82_fossil %>%
  clean_names() %>%
  select(!number_fid) # remove first column as it is redundant with the second column in this case
## Join datasets
PCs_n82_fossil <- left_join(eigenvec_n82_fossil, metadata, by = c("iid"="SRR"))
PCs_n82_fossil <- rename(PCs_n82_fossil, SRR = "iid") # Rename first column

## n72_fossil - fossil + modern blue (excl. SRR27955277_trim_Mut & black wildebeest) ----
eigenvec_n72_fossil <- read_table("n72_fossil_Goat_autosomes_noBlack_noHart.emu.eigvecs", col_names = T)
eigenvec_n72_fossil <- eigenvec_n72_fossil %>%
  clean_names() %>%
  select(!number_fid) # remove first column as it is redundant with the second column in this case
## Join datasets
PCs_n72_fossil <- left_join(eigenvec_n72_fossil, metadata, by = c("iid"="SRR"))
PCs_n72_fossil <- rename(PCs_n72_fossil, SRR = "iid") # Rename first column


# Colour palette ----
myPalette <- c("#f2933b", "#cace3d", "#8dc5a6", "#cc95ff", "#ff96ca", "#8585c4", "#5c5c5c")

# Plot ----

## n82_fossil - fossil + modern (excl. SRR27955277_trim_Mut), blue and black wildebeest ----
# Note: The Fossil936 triangle was later filled in Inkscape to be seen more clearly, as seen in Fig. 2C and Fig. S10.
# Wildebeest silhouettes also added in Inkscape.

## PC1 v PC2
### For Fig. S8
p_n82_fossil_1_2 <- ggplot(PCs_n82_fossil) +
  geom_point(aes(pc1, pc2, colour = Colour_Subspecies), pch = 2, size = 10) +
  scale_colour_manual(values = myPalette,
                      breaks = c("#f2933b", "#cace3d", "#8dc5a6", "#cc95ff", "#ff96ca", "#8585c4", "#5c5c5c"),
                      labels =  c("Bale wildebeest", "Western white-bearded", "Eastern white-bearded", "Nyassa", "Cookson", "Brindled", "Black"),
                      name = "Taxon") +
  theme_bw(base_size = 18) +
  labs(x = "PC1 (31.4%)", y = "PC2 (18.5%)") +
  theme(legend.position = "none")
p_n82_fossil_1_2$layers[[1]]$aes_params$stroke <- 1.5 # Change stroke width of the points
p_n82_fossil_1_2
ggsave(plot = p_n82_fossil_1_2, "emu_pca_n82_fossil_PC1_2_FigS6.png", dpi = 300, width = 5.57, height = 3.72, units = "in")
ggsave(plot = p_n82_fossil_1_2, "emu_pca_n82_fossil_PC1_2_FigS6.pdf", width = 5.57, height = 3.72, units = "in")


## PC1 v PC3 - not included in manuscript
#p_n82_fossil_1_3 <- ggplot(PCs_n82_fossil) +
#  geom_point(aes(pc1, pc3, colour = Colour_Subspecies), pch = 2, size = 10) +
#  scale_colour_manual(values = myPalette,
#                      breaks = c("#f2933b", "#cace3d", "#8dc5a6", "#cc95ff", "#ff96ca", "#8585c4", "#5c5c5c"),
#                      labels =  c("Bale wildebeest", "Western white-bearded", "Eastern white-bearded", "Nyassa", "Cookson", "Brindled", "Black"),
#                      name = "Taxon") +
#  theme_bw(base_size = 16) +
#  labs(x = "PC1 (31.4%)", y = "PC3 (12.9%)") +
#  theme(legend.position = "none")
#p_n82_fossil_1_3$layers[[1]]$aes_params$stroke <- 1.5 # Change stroke width of the points
#p_n82_fossil_1_3
#ggsave(plot = p_n82_fossil_1_3, "emu_pca_n82_fossil_PC1_3.png", dpi = 300)

## PC2 v PC3 - not included in manuscript
#p_n82_fossil_2_3 <- ggplot(PCs_n82_fossil) +
#  geom_point(aes(pc2, pc3, colour = Colour_Subspecies), pch = 2, size = 10) +
#  scale_colour_manual(values = myPalette,
#                      breaks = c("#f2933b", "#cace3d", "#8dc5a6", "#cc95ff", "#ff96ca", "#8585c4", "#5c5c5c"),
#                      labels =  c("Bale wildebeest", "Western white-bearded", "Eastern white-bearded", "Nyassa", "Cookson", "Brindled", "Black"),
#                      name = "Taxon") +
#  theme_bw(base_size = 16) +
#  labs(x = "PC2 (18.5%)", y = "PC3 (12.9%)") +
#  theme(legend.position = "none")
#p_n82_fossil_2_3$layers[[1]]$aes_params$stroke <- 1.5 # Change stroke width of the points
#p_n82_fossil_2_3
#ggsave(plot = p_n82_fossil_2_3, "emu_pca_n82_fossil_PC2_3.png", dpi = 300)

## n72_fossil - fossil + modern blue (excl. SRR27955277_trim_Mut & black wildebeest) ----

## PC1 v PC2
### For Fig. S8
p_n72_fossil_1_2 <- ggplot(PCs_n72_fossil) +
  geom_point(aes(pc1, pc2, colour = Colour_Subspecies), pch = 2, size = 10) +
  scale_colour_manual(values = myPalette,
                      breaks = c("#f2933b", "#cace3d", "#8dc5a6", "#cc95ff", "#ff96ca", "#8585c4"),
                      labels =  c("Bale wildebeest", "Western white-bearded", "Eastern white-bearded", "Nyassa", "Cookson", "Brindled"),
                      name = "Taxon") +
  theme_bw(base_size = 18) +
  labs(x = "PC1 (20.3%)", y = "PC2 (13.5%)") +
  theme(legend.position = "none")
p_n72_fossil_1_2$layers[[1]]$aes_params$stroke <- 1.5 # Change stroke width of the points
p_n72_fossil_1_2
ggsave(plot = p_n72_fossil_1_2, "emu_pca_n72_fossil_PC1_2_FigS6.png", dpi = 300, width = 5.57, height = 3.72, units = "in")
ggsave(plot = p_n72_fossil_1_2, "emu_pca_n72_fossil_PC1_2_FigS6.pdf", width = 5.57, height = 3.72, units = "in")

## PC1 v PC3 - not included in manuscript
#p_n72_fossil_1_3 <- ggplot(PCs_n72_fossil) +
#  geom_point(aes(pc1, pc3, colour = Colour_Subspecies), pch = 2, size = 10) +
#  scale_colour_manual(values = myPalette,
#                      breaks = c("#f2933b", "#cace3d", "#8dc5a6", "#cc95ff", "#ff96ca", "#8585c4"),
#                      labels =  c("Bale wildebeest", "Western white-bearded", "Eastern white-bearded", "Nyassa", "Cookson", "Brindled"),
#                      name = "Taxon") +
#  theme_bw(base_size = 16) +
#  labs(x = "PC1 (20.3%)", y = "PC3 (7.3%)") +
#  theme(legend.position = "none")
#p_n72_fossil_1_3$layers[[1]]$aes_params$stroke <- 1.5 # Change stroke width of the points
#p_n72_fossil_1_3
#ggsave(plot = p_n72_fossil_1_3, "emu_pca_n72_fossil_PC1_3.png", dpi = 300)

## PC2 v PC3 - not included in manuscript
#p_n72_fossil_2_3 <- ggplot(PCs_n72_fossil) +
#  geom_point(aes(pc2, pc3, colour = Colour_Subspecies), pch = 2, size = 10) +
#  scale_colour_manual(values = myPalette,
#                      breaks = c("#f2933b", "#cace3d", "#8dc5a6", "#cc95ff", "#ff96ca", "#8585c4"),
#                      labels =  c("Bale wildebeest", "Western white-bearded", "Eastern white-bearded", "Nyassa", "Cookson", "Brindled"),
#                      name = "Taxon") +
#  theme_bw(base_size = 16) +
#  labs(x = "PC2 (13.5%)", y = "PC3 (7.3%)") +
#  theme(legend.position = "none")
#p_n72_fossil_2_3$layers[[1]]$aes_params$stroke <- 1.5 # Change stroke width of the points
#p_n72_fossil_2_3
#ggsave(plot = p_n72_fossil_2_3, "emu_pca_n72_fossil_PC2_3.png", dpi = 300)

