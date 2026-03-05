# Plot stable isotope data

# Set working directory
#setwd("")

# Load packages
library(tidyverse)
library(ggpubr)
library(ggtext)

# Load data
## Read-in the saved filtered data with corrected d13C values:
iso_Alc <- read_csv("iso_Alc.csv", col_names = T)

## Get some info/stats
max(unique(iso_Alc$Age_upper)) # oldest = 104,000 years
min(unique(iso_Alc$Age_lower)) # youngest = 14,000 years
iso_Alc %>% 
  group_by(Taxon_group) %>%
  summarise(n())
# A tibble: 5 x 2
#Taxon_group           `n()`
#<chr>                 <int>
#  1 Alcelaphini             117
#2 Alcelaphini_Bale          2
#3 Cephalophini             10
#4 Tachyoryctini_Bale       17
#5 Tachyoryctini_Karungu    17

# Plot
## Plot with corrected d13C values & diet categories following Robinson (2022), who in turn followed Cerling et al. (2015)
### Corrected means the collagen-derived measurements have been adjusted to be comparable with the carbonate-derived measurements in the reference dataset
### Boxplot with ggplot
ggplot(iso_Alc, aes(delta13C_corrected, Taxon_group)) +
  annotate(geom = "rect", xmin = -1, xmax = 4, ymin = 0, ymax = 5.5, fill = "#343a40", colour = NA, alpha = 0.5) +
  annotate(geom = "rect", xmin = -8, xmax = -1, ymin = 0, ymax = 5.5, fill = "#6c757d", colour = NA, alpha = 0.5) +
  annotate(geom = "rect", xmin = -14, xmax = -8, ymin = 0, ymax = 5.5, fill = "#ced4da", colour = NA, alpha = 0.5) +
  annotate(geom = "text", x = 1.5, y = 5.35, label = "C[4]", parse = TRUE, size = 8) +
  annotate(geom = "text", x = -4.5, y = 5.14, label = "C[3]/C[4]", parse = TRUE, size = 8, angle = 90) +
  annotate(geom = "text", x = -11, y = 5.35, label = "C[3]", parse = TRUE, size = 8) +
  annotate(geom = "text", x = 6.75, y = 5, label = "42 ka", angle = 45, size = 5) +
  annotate(geom = "text", x = 6.75, y = 4, label = "104-14 ka", angle = 45, size = 5) +
  annotate(geom = "text", x = 6.75, y = 3, label = "78-14 ka", angle = 45, size = 5) +
  annotate(geom = "text", x = 6.75, y = 2, label = "47-31 ka", angle = 45, size = 5) +
  annotate(geom = "text", x = 6.75, y = 1, label = "94-45 ka", angle = 45, size = 5) +
  geom_boxplot(outlier.shape = NA) +
  scale_y_discrete(limits = c("Tachyoryctini_Karungu", "Tachyoryctini_Bale", "Cephalophini", "Alcelaphini", "Alcelaphini_Bale"),
                   labels = c("East African root-rat\n(Karungu, Kenya)", "Giant root-rat\n(Fincha Habera)", "Cephalophinae\n(Browsers)", "Alcelaphinae\n(Grazers)", "Fossil936\n(Fincha Habera)"),
                   name = NULL) +
  geom_jitter(aes(fill = Elevation_m_DEM), shape = 21, size = 3, alpha = 0.65) +
  scale_fill_viridis_c(option = "C", direction = 1, name = "Elevation\n(m)") +
  theme_pubr(base_size = 16) +
  #scale_x_continuous(breaks = waiver(), n.breaks = 8) +
  xlab(paste0("\u03b4", "<sup>13</sup>C")) +
  theme(legend.position = "right", legend.direction = "vertical") +
  #  theme(axis.text.y = element_text(face = c("italic", "italic", "plain", "plain", "plain"))) +
  theme(axis.title.x = element_markdown()) +
  scale_x_continuous(n.breaks = 6) +
  coord_cartesian(clip = "off")
#  labs(subtitle = "Late Pleistocene fauna - East Africa.\nAlcelaphini = Grazers, Cephalophini = Browsers")
ggsave("d13C_corrected_EA.png", width = 3028, height = 2356, units = "px", dpi = 300)
ggsave("d13C_corrected_EA.pdf", width = 3028, height = 2356, units = "px")
