# Plot East African fossil wildebeest occurrences through time

# Produces Fig. 3E

# Set working directory
setwd("./")

# Load packages
library(tidyverse)
library(janitor)

# Load data
fossils = read_csv("Faith2019_wildebeest_combined.csv", col_names = T) # This is essentially Table S17
fossils = fossils %>% 
  clean_names()

# Plot
## Shaded points
### Define age categories
fossils <- fossils %>%
  mutate(age_cat = case_when(mean_age >= 1.7 ~ "cat5",
                             mean_age < 1.7 & mean_age >= 1.35 ~ "cat4",
                             mean_age < 1.35 & mean_age >= 0.102 ~ "cat3",
                             mean_age < 0.102 & mean_age >= 0.0117 ~ "cat2",
                             mean_age < 0.0117 ~ "cat1"),
         .after = mean_age)
### Define colour palette for point shading
myPalette <- c("#67000d", "#d32020", "#fb7050", "#fcbea5", "#fff5f0")

### Plot
ggplot(fossils) +
  geom_point(aes(mean_age, country, fill = age_cat), shape = 21, size = 4) +
  scale_fill_discrete(palette = rev(myPalette)) +
  theme_bw(base_size = 16) +
  theme(panel.grid.major.x = element_blank(),
        panel.grid.minor.x = element_blank(),
        legend.position = "none") +
  scale_x_continuous(n.breaks = 10, name = "Mean age (Mya)") +
  scale_y_discrete(name = "Country", limits = rev)
ggsave("EA_Fossils_AllWildebeest_Faith2019_colpoints.png", width = 278.268, height = 50, units = "mm", dpi = 300)
ggsave("EA_Fossils_AllWildebeest_Faith2019_colpoints.svg", width = 278.268, height = 50, units = "mm")

