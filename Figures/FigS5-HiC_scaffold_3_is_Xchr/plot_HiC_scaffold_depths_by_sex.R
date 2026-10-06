# Plot known male and female coverage ratios per HiC_scaffold
## Produces fig. S5

# Set working directory
setwd("./")

# Load packages
library(tidyverse)

# Import data
female <- read_table("SRR27955224_female_HiC_scaffolds_coverage.txt", col_names = c("depth", "scaffold"))
female <- female %>%
  mutate(sample = "SRR27955224") %>%
  mutate(sex = "Female")
male <- read_table("SSRR27955227_male_HiC_scaffolds_coverage.txt", col_names = c("depth", "scaffold"))
male <- male %>%
  mutate(sample = "SRR27955227") %>%
  mutate(sex = "Male")
## Join datasets
depth_data <- bind_rows(female, male)
depth_data <- depth_data %>%
  mutate(category = if_else(scaffold == "HiC_scaffold_3", "chrX", "autosome"))
depth_data$scaffold <- factor(depth_data$scaffold, levels = unique(depth_data$scaffold))

# Plot
ggplot(depth_data) +
  geom_boxplot(aes(category, depth, colour = sex)) +
  scale_y_continuous(name = "Mean depth of coverage",
                     limits = c(12,24),
                     breaks = c(12, 14, 16, 18, 20, 22, 24)) +
  scale_x_discrete(name = "HiC scaffold",
                   labels = c("HiC_scaffold_1,2,4-29", "HiC_scaffold_3")) +
  scale_colour_discrete(name = "Sample",
                        labels = c(("SRR27955224\n(Female)"), "SRR27955227\n(Male)")) +
  theme_classic(base_size = 18) +
  theme(legend.position = "top")
ggsave("HiC_scaffold_depths_by_sex.png", dpi = 300)
ggsave("HiC_scaffold_depths_by_sex.pdf", dpi = 300)

# Get means
mean_depths <- depth_data %>%
  group_by(category, sex) %>%
  summarise(mean = mean(depth)) %>%
  ungroup()
write_csv(mean_depths, "mean_depths.csv")
