# Plot SeXY results
## Ratios (Autosomes:X) of bootstrap replicates

# Set working directory
#setwd("")

## Load packages
library(tidyverse)

## Import data
Fossilratios <- read_csv("SexCoverage_Wildebeest.csv", col_names = T)

## Make boxplot
ggplot(Fossilratios) +
  annotate("text", label = "Undetermined", x = 0.58, y = 0.75, size = 6) +
  annotate("text", label = "Male", x = 0.58, y = 0.6, size = 6) +
  annotate("text", label = "Female", x = 0.58, y = 0.9, size = 6) +
  annotate("rect", xmin = -Inf, xmax = Inf, ymin = 0.7, ymax = 0.8, alpha = 0.5, fill = "grey") +
  annotate("rect", xmin = -Inf, xmax = Inf, ymin = 0.5, ymax = 0.7, alpha = 0.5, fill = "lightblue") +
  annotate("rect", xmin = -Inf, xmax = Inf, ymin = 0.8, ymax = 1, alpha = 0.5, fill = "pink") +
  geom_boxplot(aes(Species, Ratio)) +
  geom_jitter(aes(Species, Ratio), colour = "black", size = 4, alpha =  0.5) +
  scale_y_continuous(limits = c(0,1)) +
  theme_classic(base_size = 24) +
  scale_x_discrete(labels = rep(NULL, 10)) +
  labs(x = "Fossil936", y = "Coverage ratio\n(X chr:Autosomes)") #+
#  theme(axis.title.x = element_blank())
ggsave("Fossil_wildebeest_sex_boxplot.png", dpi = 300, width = 7.5, height = 5)
ggsave("Fossil_wildebeest_sex_boxplot.pdf", width = 7.5, height = 5)
 