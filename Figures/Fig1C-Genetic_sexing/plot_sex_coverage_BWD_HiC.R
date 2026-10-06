# Plot SeXY results
## Ratios (Autosomes:X) of bootstrap replicates
## Produces Fig. 1C

# Set working directory
setwd("./")

## Load packages
library(tidyverse)

## Import data
Fossilratios <- read_csv("SexCoverage_Fossil936_BWD_HiC.csv", col_names = T)

## Make boxplot
ggplot(Fossilratios) +
  annotate("text", label = "Undetermined", x = 0.58, y = 0.75, size = 6) +
  annotate("text", label = "Male", x = 0.58, y = 0.6, size = 6) +
  annotate("text", label = "Female", x = 0.58, y = 0.9, size = 6) +
  annotate("rect", xmin = -Inf, xmax = Inf, ymin = 0.7, ymax = 0.8, alpha = 0.5, fill = "grey") +
  annotate("rect", xmin = -Inf, xmax = Inf, ymin = 0.5, ymax = 0.7, alpha = 0.5, fill = "lightblue") +
  annotate("rect", xmin = -Inf, xmax = Inf, ymin = 0.8, ymax = 1, alpha = 0.5, fill = "pink") +
  geom_boxplot(aes(SampleID, Ratio)) +
  geom_jitter(aes(SampleID, Ratio), fill = "#f2933b", colour = "#232323", size = 4, pch = 21) +
  scale_y_continuous(limits = c(0,1)) +
  theme_classic(base_size = 24) +
  theme(text = element_text(family = "arial")) +
  scale_x_discrete(labels = rep(NULL, 10)) +
  scale_y_continuous(limits = c(0.49, 1)) +
  labs(x = "Fossil936", y = "Coverage ratio\n(X chr:Autosomes)") #+
#  theme(axis.title.x = element_blank())
ggsave("Fossil936_BWD_HiC_sex_boxplot.png", dpi = 300, width = 7.5, height = 5)
ggsave("Fossil936_BWD_HiC_sex_boxplot.svg", width = 7.5, height = 5)
 