# Plot emu pca output (pca based on pseudohap files)
## Note: The pseudohap files exclude transitions and were filtered for maf=0.05
## Produces basis for Fig. 2C and Fig. S10

# Set working directory
setwd("./")

# Load packages
library(tidyverse)
library(janitor)
library(scales)
library(ggpubr)

# Load data ----
metadata <- read_csv("../metadata.csv", col_names = TRUE)

## Note: The percentage variance explained for each PC can be found in the *.eigenval files.

## n82_fossil - fossil + modern blue and black wildebeest ----
eigenvec_n82_fossil <- read_table("n82_fossil_BWD_HiC_autosomes_noHart.emu.eigvecs", col_names = T)
eigenvec_n82_fossil <- eigenvec_n82_fossil %>%
  clean_names() %>%
  select(!number_fid) # remove first column as it is redundant with the second column in this case
## Join datasets
PCs_n82_fossil <- left_join(eigenvec_n82_fossil, metadata, by = c("iid"="SRR"))
PCs_n82_fossil <- rename(PCs_n82_fossil, SRR = "iid") # Rename first column

## n72_fossil - fossil + modern blue (excl. SRR27955277_damaged & black wildebeest) ----
eigenvec_n72_fossil <- read_table("n72_fossil_BWD_HiC_autosomes_noBlack_noHart.emu.eigvecs", col_names = T)
eigenvec_n72_fossil <- eigenvec_n72_fossil %>%
  clean_names() %>%
  select(!number_fid) # remove first column as it is redundant with the second column in this case
## Join datasets
PCs_n72_fossil <- left_join(eigenvec_n72_fossil, metadata, by = c("iid"="SRR"))
PCs_n72_fossil <- rename(PCs_n72_fossil, SRR = "iid") # Rename first column

## n81_modern - SRR27955277_damaged (excl. Fossil936), blue and black wildebeest (no Hart, as with all PCA datasets) ----
eigenvec_n81_modern <- read_table("n81_modern_BWD_HiC_autosomes_noHart.emu.eigvecs", col_names = T)
eigenvec_n81_modern <- eigenvec_n81_modern %>%
  clean_names() %>%
  select(!number_fid) # remove first column as it is redundant with the second column in this case
## Join datasets
PCs_n81_modern <- left_join(eigenvec_n81_modern, metadata, by = c("iid"="SRR"))
PCs_n81_modern <- rename(PCs_n81_modern, SRR = "iid") # Rename first column

## n81_modern_noBlack - SRR27955277_damaged and blue wildebeest (excl. Fossil936 and black wildebeest) ----
eigenvec_n81_modern_noBlack <- read_table("n81_modern_BWD_HiC_autosomes_noBlack_noHart.emu.eigvecs", col_names = T)
eigenvec_n81_modern_noBlack <- eigenvec_n81_modern_noBlack %>%
  clean_names() %>%
  select(!number_fid) # remove first column as it is redundant with the second column in this case
## Join datasets
PCs_n81_modern_noBlack <- left_join(eigenvec_n81_modern_noBlack, metadata, by = c("iid"="SRR"))
PCs_n81_modern_noBlack <- rename(PCs_n81_modern_noBlack, SRR = "iid") # Rename first column

# Colour palette ----
myPalette <- c("#f2933b", "#cace3d", "#8dc5a6", "#cc95ff", "#ff96ca", "#8585c4", "#5c5c5c") # For Fossil936 datasets
myPalette1 <- c("#780000", "#cace3d", "#8dc5a6", "#cc95ff", "#ff96ca", "#8585c4", "#5c5c5c") # For 277_damaged datasets


# Plot ----

## n82_fossil - fossil + modern (excl. SRR27955277_trim_Mut), blue and black wildebeest ----
# Note: The Fossil936 triangle was later filled in Inkscape to be seen more clearly, as seen in Fig. 2C and Fig. S10.
# Wildebeest silhouettes also added in Inkscape.

## PC1 v PC2
### For Fig. 2
p_fig2.1 <- ggplot(PCs_n82_fossil) +
  geom_point(aes(pc1, pc2, colour = Colour_Subspecies), pch = 2, size = 10) +
  scale_colour_manual(values = myPalette,
                      breaks = c("#f2933b", "#cace3d", "#8dc5a6", "#cc95ff", "#ff96ca", "#8585c4", "#5c5c5c"),
                      labels =  c("Bale wildebeest", "Western white-bearded", "Eastern white-bearded", "Nyassa", "Cookson", "Brindled", "Black"),
                      name = "Taxon") +
  scale_y_continuous(limits = c(-0.2, 0.15)) +
  scale_x_continuous(limits = c(-0.11, 0.33)) +
  theme_bw(base_size = 18) +
  labs(x = "PC1 (32.5%)", y = "PC2 (18.8%)") +
  theme(legend.position = "right")
p_fig2.1$layers[[1]]$aes_params$stroke <- 1.5 # Change stroke width of the points
p_fig2.1
ggsave(plot = p_fig2.1, "emu_pca_n82_fossil_PC1_2_fig2.1.png", dpi = 300)
ggsave(plot = p_fig2.1, "emu_pca_n82_fossil_PC1_2_fig2.1.pdf")

### Get legend from p_fig2.1 (change theme(legend.position = "bottom") as needed ("none" to remove))
p_fig2.1.legend <- get_legend(p_fig2.1)
p_fig2.1.legend <- as_ggplot(p_fig2.1.legend)
p_fig2.1.legend
ggsave(plot = p_fig2.1.legend, "emu_pca_n82_fossil_legend.png", width = 200, height = 30, units = "mm", dpi = 300)
ggsave(plot = p_fig2.1.legend, "emu_pca_n82_fossil_legend.pdf", width = 200, height = 30, units = "mm")
ggsave(plot = p_fig2.1.legend, "emu_pca_n82_fossil_legend_vert.pdf", width = 75, height = 90, units = "mm")

## PC1 v PC3
p_n82_fossil_1_3 <- ggplot(PCs_n82_fossil) +
  geom_point(aes(pc1, pc3, colour = Colour_Subspecies), pch = 2, size = 10) +
  scale_colour_manual(values = myPalette,
                      breaks = c("#f2933b", "#cace3d", "#8dc5a6", "#cc95ff", "#ff96ca", "#8585c4", "#5c5c5c"),
                      labels =  c("Bale wildebeest", "Western white-bearded", "Eastern white-bearded", "Nyassa", "Cookson", "Brindled", "Black"),
                      name = "Taxon") +
  scale_x_continuous(name = "PC1 (32.5%)", limits = c(-0.1, 0.31)) +
  scale_y_continuous(name = "PC3 (13.0%)", limits = c(-0.1, 0.35)) +
  theme_bw(base_size = 16) +
  theme(legend.position = "none")
p_n82_fossil_1_3$layers[[1]]$aes_params$stroke <- 1.5 # Change stroke width of the points
p_n82_fossil_1_3
ggsave(plot = p_n82_fossil_1_3, "emu_pca_n82_fossil_PC1_3.png", dpi = 300)
ggsave(plot = p_n82_fossil_1_3, "emu_pca_n82_fossil_PC1_3_noLeg.pdf", width = 5.57, height = 3.72, units = "in")

## PC2 v PC3
p_n82_fossil_2_3 <- ggplot(PCs_n82_fossil) +
  geom_point(aes(pc2, pc3, colour = Colour_Subspecies), pch = 2, size = 10) +
  scale_colour_manual(values = myPalette,
                      breaks = c("#f2933b", "#cace3d", "#8dc5a6", "#cc95ff", "#ff96ca", "#8585c4", "#5c5c5c"),
                      labels =  c("Bale wildebeest", "Western white-bearded", "Eastern white-bearded", "Nyassa", "Cookson", "Brindled", "Black"),
                      name = "Taxon") +
  scale_x_continuous(name = "PC2 (18.8%)", limits = c(-0.17, 0.135)) +
  scale_y_continuous(name = "PC3 (13.0%)", limits = c(-0.1, 0.35)) +
  theme_bw(base_size = 16) +
  theme(legend.position = "none")
p_n82_fossil_2_3$layers[[1]]$aes_params$stroke <- 1.5 # Change stroke width of the points
p_n82_fossil_2_3
ggsave(plot = p_n82_fossil_2_3, "emu_pca_n82_fossil_PC2_3.png", dpi = 300)
ggsave(plot = p_n82_fossil_2_3, "emu_pca_n82_fossil_PC2_3_noLeg.pdf", width = 5.57, height = 3.72, units = "in")

## n72_fossil - fossil + modern blue (excl. SRR27955277_trim_Mut & black wildebeest) ----
## PC1 v PC2
### For Fig. 2
p_fig2.2 <- ggplot(PCs_n72_fossil) +
  geom_point(aes(pc1, pc2, colour = Colour_Subspecies), pch = 2, size = 10) +
  scale_colour_manual(values = myPalette,
                      breaks = c("#f2933b", "#cace3d", "#8dc5a6", "#cc95ff", "#ff96ca", "#8585c4"),
                      labels =  c("Bale wildebeest", "Western white-bearded", "Eastern white-bearded", "Nyassa", "Cookson", "Brindled"),
                      name = "Taxon") +
  scale_y_continuous(limits = c(-0.305, 0.12)) +
  scale_x_continuous(limits = c(-0.11, 0.2)) +
  theme_bw(base_size = 18) +
  labs(x = "PC1 (20.2%)", y = "PC2 (13.6%)") +
  theme(legend.position = "none")
p_fig2.2$layers[[1]]$aes_params$stroke <- 1.5 # Change stroke width of the points
p_fig2.2
ggsave(plot = p_fig2.2, "emu_pca_n72_fossil_PC1_2_fig2.2.png", dpi = 300)
ggsave(plot = p_fig2.2, "emu_pca_n72_fossil_PC1_2_fig2.2.pdf", dpi = 300)

## PC1 v PC3
p_n72_fossil_1_3 <- ggplot(PCs_n72_fossil) +
  geom_point(aes(pc1, pc3, colour = Colour_Subspecies), pch = 2, size = 10) +
  scale_colour_manual(values = myPalette,
                      breaks = c("#f2933b", "#cace3d", "#8dc5a6", "#cc95ff", "#ff96ca", "#8585c4"),
                      labels =  c("Bale wildebeest", "Western white-bearded", "Eastern white-bearded", "Nyassa", "Cookson", "Brindled"),
                      name = "Taxon") +
  scale_x_continuous(name = "PC1 (20.2%)", limits = c(-0.115, 0.2)) +
  scale_y_continuous(name = "PC3 (7.1%)", limits = c(-0.16, 0.32)) +
  theme_bw(base_size = 16) +
  theme(legend.position = "none")
p_n72_fossil_1_3$layers[[1]]$aes_params$stroke <- 1.5 # Change stroke width of the points
p_n72_fossil_1_3
ggsave(plot = p_n72_fossil_1_3, "emu_pca_n72_fossil_PC1_3.png", dpi = 300)
ggsave(plot = p_n72_fossil_1_3, "emu_pca_n72_fossil_PC1_3_noLeg.pdf", width = 5.57, height = 3.72, units = "in")

## PC2 v PC3
p_n72_fossil_2_3 <- ggplot(PCs_n72_fossil) +
  geom_point(aes(pc2, pc3, colour = Colour_Subspecies), pch = 2, size = 10) +
  scale_colour_manual(values = myPalette,
                      breaks = c("#f2933b", "#cace3d", "#8dc5a6", "#cc95ff", "#ff96ca", "#8585c4"),
                      labels =  c("Bale wildebeest", "Western white-bearded", "Eastern white-bearded", "Nyassa", "Cookson", "Brindled"),
                      name = "Taxon") +
  scale_x_continuous(name = "PC2 (13.6%)", limits = c(-0.305, 0.12)) +
  scale_y_continuous(name = "PC3 (7.1%)", limits = c(-0.16, 0.315)) +
  theme_bw(base_size = 16) +
  theme(legend.position = "none")
p_n72_fossil_2_3$layers[[1]]$aes_params$stroke <- 1.5 # Change stroke width of the points
p_n72_fossil_2_3
ggsave(plot = p_n72_fossil_2_3, "emu_pca_n72_fossil_PC2_3.png", dpi = 300)
ggsave(plot = p_n72_fossil_2_3, "emu_pca_n72_fossil_PC2_3_noLeg.pdf", width = 5.57, height = 3.72, units = "in")



# Modern dataset ----
## n81_modern - SRR27955277_damaged (excl. Fossil936), blue and black wildebeest (no Hart, as with all PCA datasets) ----
# Note: The outline colour of the 277_damaged triangle was later changed in Inkscape and it was filled, to match what is seen in Fig. S10.

## PC1 v PC2
p_n81_modern_1_2 <- ggplot(PCs_n81_modern) + 
  geom_point(aes(pc1, pc2, colour = Colour_Subspecies), pch = 24, size = 10) +
  scale_colour_manual(values = myPalette1,
                      breaks = c("#780000", "#cace3d", "#8dc5a6", "#cc95ff", "#ff96ca", "#8585c4", "#5c5c5c"),
                      labels =  c("277_damaged", "Western white-bearded", "Eastern white-bearded", "Nyassa", "Cookson", "Brindled", "Black"),
                      name = "Taxon") +
  theme_bw(base_size = 16) +
  labs(x = "PC1 (35.2%)", y = "PC2 (20.4%)") +
  theme(legend.position = "none")
p_n81_modern_1_2$layers[[1]]$aes_params$stroke <- 1.5 # Change stroke width of the points
p_n81_modern_1_2
ggsave(plot = p_n81_modern_1_2, "emu_pca_n81_modern_PC1_2.png", width = 5.57, height = 3.72, units = "in", dpi = 300)
ggsave(plot = p_n81_modern_1_2, "emu_pca_n81_modern_PC1_2.svg", width = 5.57, height = 3.72, units = "in")

### Get legend from p_n81_modern_1_2 (change theme(legend.position = "bottom") as needed ("none" to remove))
p_n81_modern_1_2.legend <- get_legend(p_n81_modern_1_2)
p_n81_modern_1_2.legend <- as_ggplot(p_n81_modern_1_2.legend)
p_n81_modern_1_2.legend
ggsave(plot = p_n81_modern_1_2.legend, "emu_pca_n81_modern_legend.svg", width = 200, height = 30, units = "mm")

## PC1 v PC3
p_n81_modern_1_3 <- ggplot(PCs_n81_modern) + 
  geom_point(aes(pc1, pc3, colour = Colour_Subspecies), pch = 24, size = 10) +
  scale_colour_manual(values = myPalette1,
                      breaks = c("#780000", "#cace3d", "#8dc5a6", "#cc95ff", "#ff96ca", "#8585c4", "#5c5c5c"),
                      labels =  c("277_damaged", "Western white-bearded", "Eastern white-bearded", "Nyassa", "Cookson", "Brindled", "Black"),
                      name = "Taxon") +
  theme_bw(base_size = 16) +
  labs(x = "PC1 (35.2%)", y = "PC3 (14.4%)") +
  theme(legend.position = "none")
p_n81_modern_1_3$layers[[1]]$aes_params$stroke <- 1.5 # Change stroke width of the points
p_n81_modern_1_3
ggsave(plot = p_n81_modern_1_3, "emu_pca_n81_modern_PC1_3.png", width = 5.57, height = 3.72, units = "in", dpi = 300)
ggsave(plot = p_n81_modern_1_3, "emu_pca_n81_modern_PC1_3.svg", width = 5.57, height = 3.72, units = "in")

## PC2 v PC3
p_n81_modern_2_3 <- ggplot(PCs_n81_modern) + 
  geom_point(aes(pc2, pc3, colour = Colour_Subspecies), pch = 24, size = 10) +
  scale_colour_manual(values = myPalette1,
                      breaks = c("#780000", "#cace3d", "#8dc5a6", "#cc95ff", "#ff96ca", "#8585c4", "#5c5c5c"),
                      labels =  c("277_damaged", "Western white-bearded", "Eastern white-bearded", "Nyassa", "Cookson", "Brindled", "Black"),
                      name = "Taxon") +
  theme_bw(base_size = 16) +
  labs(x = "PC2 (20.4%)", y = "PC3 (14.4%)") +
  theme(legend.position = "none")
p_n81_modern_2_3$layers[[1]]$aes_params$stroke <- 1.5 # Change stroke width of the points
p_n81_modern_2_3
ggsave(plot = p_n81_modern_2_3, "emu_pca_n81_modern_PC2_3.png", width = 5.57, height = 3.72, units = "in", dpi = 300)
ggsave(plot = p_n81_modern_2_3, "emu_pca_n81_modern_PC2_3.svg", width = 5.57, height = 3.72, units = "in")

## n81_modern_noBlack - SRR27955277_damaged and blue wildebeest (excl. Fossil936 and black wildebeest) ----
## PC1 v PC2
p_n81_modern_noBlack_1_2 <- ggplot(PCs_n81_modern_noBlack) + 
  geom_point(aes(pc1, pc2, colour = Colour_Subspecies), pch = 24, size = 10) +
  scale_colour_manual(values = myPalette1,
                      breaks = c("#780000", "#cace3d", "#8dc5a6", "#cc95ff", "#ff96ca", "#8585c4"),
                      labels =  c("277_damaged", "Western white-bearded", "Eastern white-bearded", "Nyassa", "Cookson", "Brindled"),
                      name = "Taxon") +
  theme_bw(base_size = 16) +
  labs(x = "PC1 (21.8%)", y = "PC2 (15.1%)") +
  theme(legend.position = "none")
p_n81_modern_noBlack_1_2$layers[[1]]$aes_params$stroke <- 1.5 # Change stroke width of the points
p_n81_modern_noBlack_1_2
ggsave(plot = p_n81_modern_noBlack_1_2, "emu_pca_n81_modern_noBlack_PC1_2.png", width = 5.57, height = 3.72, units = "in", dpi = 300)
ggsave(plot = p_n81_modern_noBlack_1_2, "emu_pca_n81_modern_noBlack_PC1_2.svg", width = 5.57, height = 3.72, units = "in")

## PC1 v PC3
p_n81_modern_noBlack_1_3 <- ggplot(PCs_n81_modern_noBlack) + 
  geom_point(aes(pc1, pc3, colour = Colour_Subspecies), pch = 24, size = 10) +
  scale_colour_manual(values = myPalette1,
                      breaks = c("#780000", "#cace3d", "#8dc5a6", "#cc95ff", "#ff96ca", "#8585c4"),
                      labels =  c("277_damaged", "Western white-bearded", "Eastern white-bearded", "Nyassa", "Cookson", "Brindled"),
                      name = "Taxon") +
  theme_bw(base_size = 16) +
  labs(x = "PC1 (21.8%)", y = "PC3 (7.8%)") +
  theme(legend.position = "none")
p_n81_modern_noBlack_1_3$layers[[1]]$aes_params$stroke <- 1.5 # Change stroke width of the points
p_n81_modern_noBlack_1_3
ggsave(plot = p_n81_modern_noBlack_1_3, "emu_pca_n81_modern_noBlack_PC1_3.png", width = 5.57, height = 3.72, units = "in", dpi = 300)
ggsave(plot = p_n81_modern_noBlack_1_3, "emu_pca_n81_modern_noBlack_PC1_3.svg", width = 5.57, height = 3.72, units = "in")

## PC2 v PC3
p_n81_modern_noBlack_2_3 <- ggplot(PCs_n81_modern_noBlack) + 
  geom_point(aes(pc2, pc3, colour = Colour_Subspecies), pch = 24, size = 10) +
  scale_colour_manual(values = myPalette1,
                      breaks = c("#780000", "#cace3d", "#8dc5a6", "#cc95ff", "#ff96ca", "#8585c4"),
                      labels =  c("277_damaged", "Western white-bearded", "Eastern white-bearded", "Nyassa", "Cookson", "Brindled"),
                      name = "Taxon") +
  theme_bw(base_size = 16) +
  labs(x = "PC2 (15.1%)", y = "PC3 (7.8%)") +
  theme(legend.position = "none")
p_n81_modern_noBlack_2_3$layers[[1]]$aes_params$stroke <- 1.5 # Change stroke width of the points
p_n81_modern_noBlack_2_3
ggsave(plot = p_n81_modern_noBlack_2_3, "emu_pca_n81_modern_noBlack_PC2_3.png", width = 5.57, height = 3.72, units = "in", dpi = 300)
ggsave(plot = p_n81_modern_noBlack_2_3, "emu_pca_n81_modern_noBlack_PC2_3.svg", width = 5.57, height = 3.72, units = "in")
