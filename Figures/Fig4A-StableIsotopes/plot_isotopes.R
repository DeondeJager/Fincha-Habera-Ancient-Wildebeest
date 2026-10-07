# Plot stable isotope data

# Produces Fig. 4A.

# Set working directory
setwd("./")

# Load packages
library(tidyverse)
library(ggpubr)
library(ggtext)

# Load data
## Read-in the saved filtered data with corrected d13C values:
iso_Alc <- read_csv("TableS10_StableIsotopes_dataset.csv", col_names = T, skip = 1, n_max = 129)

## Get some info/stats
max(unique(iso_Alc$Age_upper)) # oldest = 104,000 years
min(unique(iso_Alc$Age_lower)) # youngest = 14,000 years
iso_Alc %>% 
  group_by(Taxon_group) %>%
  summarise(n())
# A tibble: 3 x 2
#Taxon_group           `n()`
#<chr>                 <int>
#  1 Alcelaphini             117
#2 Alcelaphini_Bale          2
#3 Cephalophini             10

# Plot
## Plot with d13C_diet values and perc_C4 axis
p13C <- ggplot(iso_Alc, aes(delta13C_diet, Taxon_group)) +
  annotate(geom = "rect", xmin = -15.1, xmax = -11.7, ymin = 0, ymax = 3.5, fill = "#343a40", colour = NA, alpha = 0.5) + 
  annotate(geom = "rect", xmin = -22.1, xmax = -15.1, ymin = 0, ymax = 3.5, fill = "#6c757d", colour = NA, alpha = 0.5) +
  annotate(geom = "rect", xmin = -26.6, xmax = -22.1, ymin = 0, ymax = 3.5, fill = "#ced4da", colour = NA, alpha = 0.5) +
  annotate(geom = "segment", x = -11.7+1.2, xend = -11.7+1.2, y = 0, yend = 3.5, linetype = "dashed") +
  annotate(geom = "segment", x = -26.6-2.2, xend = -26.6-2.2, y = 0, yend = 3.5, linetype = "dashed") +
  annotate(geom = "text", x = -13.4, y = 3.35, label = "C[4]", parse = TRUE, size = 8) +
  annotate(geom = "text", x = -18.6, y = 3.35, label = "C[3]/C[4]", parse = TRUE, size = 8) +
  annotate(geom = "text", x = -25, y = 3.35, label = "C[3]", parse = TRUE, size = 8) +
  annotate(geom = "text", x = -8.75, y = 3, label = "42 ka", angle = 45, size = 6) +
  annotate(geom = "text", x = -8.75, y = 2, label = "104-14 ka", angle = 45, size = 6) +
  annotate(geom = "text", x = -8.75, y = 1, label = "78-14 ka", angle = 45, size = 6) +
  geom_boxplot(outlier.shape = NA) +
  scale_y_discrete(limits = c("Cephalophini", "Alcelaphini", "Alcelaphini_Bale"),
                   labels = c("Cephalophini\n(Browsers)", "Alcelaphini\n(Grazers)", "Fossil936\n(Grazer)"),
                   name = NULL) +
  geom_jitter(aes(fill = Elevation_m_DEM), shape = 21, size = 5, alpha = 0.65) +
  scale_fill_viridis_c(option = "C", direction = 1, name = "Elevation\n(m a.s.l.)") +
  theme_pubr(base_size = 16) +
  xlab(paste0("\u03b4", "<sup>13</sup>C<sub>diet</sub>")) +
  theme(legend.position = "right", legend.direction = "vertical") +
  theme(axis.title.x = element_markdown()) +
  scale_x_continuous(breaks = round(c(-26.6, -23.62, -20.64, -17.66, -14.68, -11.7), 1), 
                     sec.axis = dup_axis(name = paste0("Percent C", "<sub>4</sub> biomass in diet"),
                                         labels = c(0, 20, 40, 60, 80, 100))) +
  coord_cartesian(clip = "off")
#  labs(subtitle = "Late Pleistocene fauna - East Africa.\nAlcelaphini = Grazers, Cephalophini = Browsers")
p13C
ggsave("d13C_diet_EA_noRR.png", p13C, width = 3028, height = 2356, units = "px", dpi = 300)
ggsave("d13C_diet_EA_noRR.svg", p13C, width = 3028, height = 2356, units = "px")
