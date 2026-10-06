# Test heterozygosity parameters with 277, 277_3X, and 277_damaged, and HiC_scaffold_1 

# Load packages
library(tidyverse)

# Set working directory
setwd("./")

# Import data
## I just manually prepped the below file from the .ml files
SRR277_het <- read_csv("SRR27955277_tests.csv", col_names = TRUE, col_types = "cccfcn")
## Re-order for plotting
SRR277_het$Short_name <- factor(SRR277_het$Short_name, 
                                levels = c("277_baq0", "277_baq1", "277_baq2", "277_3X_baq0", "277_3X_baq0_minD3", "277_3X_baq1", "277_3X_baq1_minD3", "277_3X_baq2", 
                                           "277_3X_baq2_minD3", "277_damaged_baq0", "277_damaged_baq0_minD3", "277_damaged_baq1", "277_damaged_baq1_minD3", "277_damaged_baq2", 
                                           "277_damaged_baq2_minD3", "277_damaged_baq2_C50", "277_damaged_baq2_trim4", "277_damaged_baq2_trim10", "277_damaged_baq2_fold"))

# Plot
## Quick plot
ggplot(SRR277_het) +
  geom_boxplot(aes(Short_name, Prop_het, colour = baq)) +
  scale_x_discrete(guide = guide_axis(angle = 45), name = "Sample and parameters") +
  scale_y_continuous(name = "Heterozygosity") +
  geom_vline(xintercept = 9.5, linetype = "dashed") +
  annotate("text", x = 14.5, y = 0.0062, label = "Ancient damage added") +
  annotate("text", x = 4.5, y = 0.0062, label = "No damage") #+
  #labs(subtitle = "Heterozygosity estimates of SRR27955277 (HiC_scaffold_1, 10 Mb windows)") +
#  theme(plot.margin = margin(0,0.1,0,1, "cm"))
ggsave(filename = "het_tests_all.png", plot = last_plot(), width = 200, height = 150, units = "mm", dpi = 300)
