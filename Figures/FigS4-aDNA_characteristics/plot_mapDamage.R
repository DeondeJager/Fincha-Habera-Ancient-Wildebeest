# Plot mapDamage misincorporation and read length distribution patterns for Fossil936 libraries & simulated aDNA sample (277_damaged)
## Produces fig. S4

####---- Set working directory ----####
setwd("./")

####---- Load packages ----####
library(tidyverse)
library(writexl)
library(ggpubr)
library(cowplot)
library(patchwork)

####---- Get and organise misincorporation data ----####
#### Load in previously processed data ####
## If done before, then just load the saved csv files for plotting, otherwise run steps under "Process misincorporation data from mapDamage" and "Process read length distribution data from mapDamage" below.
mapDamage <- read_csv("mapDamage_misincorporation_summary_R.csv", col_names = T)
rl <- read_csv("mapDamage_lgdistribution_summary_R.csv", col_names = T)

#### Process misincorporation data from mapDamage ####
# Code by DdJ and ChatGPT
## The code loops through each sample and uses its "misincorporation.txt" file to calculate the G>A and C>T frequencies per position,
## and then stores this in a dataframe with the SeqID, which can then be used for plotting.

# 1. Get list of all misincorporation files to process (in current working directory)
file_list <- list.files(".", pattern = "_misincorporation.txt$")

# 2. Initialize an empty list to store each file's processed data
results_list <- list()

# 3. Loop through each file
for (file in file_list) {
  
  SeqID <- str_remove(file, "_misincorporation.txt$")

  # Find the line number where "Chr" first appears in *misincorporation.txt file
  first_chr_line <- grep("^Chr", readLines(file), ignore.case = FALSE)[1]
  # Read the TSV starting from that line
  # subtract 1 because skip counts lines *before* reading
  data <- read_tsv(file, skip = first_chr_line - 1, col_names = TRUE)
  
  # Perform calculations on each row
  processed_data <- data %>%
    #rowwise() %>%
    mutate(GtoA = `G>A`/G) %>% 
    mutate(CtoT = `C>T`/C) %>%
    #ungroup() %>%
    select(End, Std, Pos, GtoA, CtoT) %>%
    mutate(SeqID = SeqID) # Add a columns with the SeqID

  # Add the processed data to the results list
  results_list[[length(results_list) + 1]] <- processed_data
}

# 4. Combine all the processed data into a single tibble by row-binding
mapDamage <- bind_rows(results_list)
# View the final combined tibble
print(mapDamage)

# 5. Average the frequencies per position across the two strands (+ and -)
tmp1 <- mapDamage %>% 
  group_by(End, Pos, SeqID) %>% 
  summarise(GtoA_mean = mean(GtoA)) %>%
  arrange(SeqID)
tmp2 <- mapDamage %>% 
  group_by(End, Pos, SeqID) %>% 
  summarise(CtoT_mean = mean(CtoT)) %>%
  arrange(SeqID)
# Overwrite previous dataframe with new one with means across positions
mapDamage <- inner_join(tmp1, tmp2)
print(mapDamage)

# 6. Save as file
write_csv(mapDamage, "mapDamage_misincorporation_summary_R.csv")

#### Process read length distribution data from mapDamage ####
# 1. Get list of all read length files to process (in current working directory)
file_list <- list.files(".", pattern = "_lgdistribution.txt$")

# 2. Initialize an empty list to store each file's processed data
results_list <- list()

# 3. Loop through each file
for (file in file_list) {
  
  SeqID <- str_remove(file, "_lgdistribution.txt$")
  
  # Find the line number where "Std" first appears in *lgdistribution.txt file
  first_chr_line <- grep("^Std", readLines(file), ignore.case = FALSE)[1]
  # Read the TSV starting from that line
  # subtract 1 because skip counts lines *before* reading
  data <- read_tsv(file, skip = first_chr_line - 1, col_names = TRUE)
  
  # Combine + and - strand occurrence counts
  processed_data <- data %>%
    group_by(Length) %>%
    summarise(Occurences = sum(Occurences)) %>%
    mutate(SeqID = SeqID)
  
  # Add the processed data to the results list
  results_list[[length(results_list) + 1]] <- processed_data
}

# 4. Combine all the processed data into a single tibble by row-binding
rl <- bind_rows(results_list)
# View the final combined tibble
print(rl)

# 5. Save as file
write_csv(rl, "mapDamage_lgdistribution_summary_R.csv")



####---- Plot damage patterns and read length dist together ----####
## Get max value to set y-axis limits
filter(mapDamage, Pos <= 25) %>% ungroup() %>% summarise(max(CtoT_mean, na.rm = T)) # 0.398

## Explicitly set order of libraries for plotting order
### Misincorporation
mapDamage$SeqID <- factor(mapDamage$SeqID,
                          levels = c("a9361SE", "277_damaged", "a9361_HiSeq", "a9361_NovaSeq", "a9362_HiSeq", "a9362_NovaSeq", "SCR116_Screening", "SCR116_DeepSeq", "SCR117_Screening", "SCR117_DeepSeq"))
# Verify it worked
levels(mapDamage$SeqID)

### Read length
rl$SeqID <- factor(rl$SeqID,
                   levels = c("a9361SE", "277_damaged", "a9361_HiSeq", "a9361_NovaSeq", "a9362_HiSeq", "a9362_NovaSeq", "SCR116_Screening", "SCR116_DeepSeq", "SCR117_Screening", "SCR117_DeepSeq"))
# Verify it worked
levels(rl$SeqID)

## Make individual plots and then add them together with cowplot
# A. Make an empty list to store each library's plots in
results_list <- list()
# B. Loop through each SeqID in the dataframe and make plots - note we use "levels" instead of "unique" in the for loop, as unique() does not respect the order of levels
for (lib_name in levels(mapDamage$SeqID)) {
  # Make plots
  ## 5p end (CtoT)
  p1 <- ggplot(filter(mapDamage, Pos <= 25) %>% filter(SeqID == lib_name) %>% filter(End == "5p")) +
    theme_pubr() +
    geom_line(aes(Pos, CtoT_mean), colour = "red", linewidth = 1) +
    geom_line(aes(Pos, GtoA_mean), colour = "blue", linewidth = 1) +
    scale_y_continuous(limits = c(0,0.4), breaks = c(0.0, 0.1, 0.2, 0.3, 0.4), name = NULL) +
    scale_x_continuous(breaks = seq(from = 0, to = 25, by = 5), name = NULL)
  ## 3p end (GtoA)
  p2 <- ggplot(filter(mapDamage, Pos <= 25) %>% filter(SeqID == lib_name) %>% filter(End == "3p")) +
    theme_pubr() +
    geom_line(aes(Pos, CtoT_mean), colour = "red", linewidth = 1) +
    geom_line(aes(Pos, GtoA_mean), colour = "blue", linewidth = 1) +
    scale_x_reverse(breaks = seq(from = 0, to = 25, by = 5), labels = c("0", "-5", "-10", "-15", "-20", "-25"), name = NULL) +
    scale_y_continuous(limits = c(0,0.4), breaks = c(0.0, 0.1, 0.2, 0.3, 0.4), name = NULL, position = "right", labels = NULL)
  ## Join together in a side-by-side plot with joint title
  p3 <- cowplot::plot_grid(p1, p2, ncol = 2, align = "hv")
  p4 <- ggdraw(p3) + 
    draw_label(lib_name, fontface = "bold", size = 14, x = 0.5, y = 0.9)
  ## Plot length
  p_rl <- ggplot(filter(rl, SeqID == lib_name)) +
    theme_pubr() +
    geom_col(aes(Length, Occurences)) +
    scale_x_continuous(breaks = seq(from = 30, to = 100, by = 10),
                       labels = c("30", "", "50", "", "70", "", "90", ""),
                       name = "Read length (bp)") +
    coord_cartesian(xlim = c(30, 100)) +
    scale_y_continuous(labels = scales::label_number(scale_cut = scales::cut_short_scale()),
                       n.breaks = 3,
                       name = "Count")
  # Join misincorporation and read length plots
  p4_rl <- ggdraw() +
    draw_plot(p4) +
    draw_plot(p_rl, width = 0.5, height = 0.5, hjust = -0.5, vjust = -0.5)
  
  # Add all the plots to the results list
  results_list[[length(results_list) + 1]] <- p4_rl
}
results_list

# C. Combined plots into a grid
p5 <- cowplot::plot_grid(plotlist = results_list, ncol = 2, nrow = 5, align = "hv")
p5

# D. Get legend to use as common legend - extract it from another plot
## Transform dataframe into long format for End (5p, 3p)
mapDamage_long <- pivot_longer(mapDamage, cols = c("GtoA_mean", "CtoT_mean"), names_to = "Damage", values_to = "Frequency")
mapDamage_long$End <- factor(mapDamage_long$End, levels = c("5p", "3p"))
p_legend <- ggplot() +
  theme_pubr() +
  geom_line(data = filter(mapDamage_long, Pos <= 10) %>% filter(SeqID == "277_damaged"), 
            mapping = aes(Pos, Frequency, colour = Damage), linewidth = 1) +
  scale_colour_manual(values = c("red", "blue"), labels = c("T", "A"), name = "Nucleotide:") +
  theme(legend.title = element_text(size = 16), legend.text = element_text(size = 14),
        legend.margin = margin(t = 0, r = 0, b = 0, l = 0)) +
  facet_wrap(~End)
p_legend
leg <- ggpubr::get_legend(p_legend)
leg

# E. Arrange all plots in a grid with legend
p6 <- ggarrange(p5, ncol = 1, nrow = 1, labels = NULL, align = "hv",
                common.legend = TRUE, legend.grob = leg, legend = "top")
p6

# F. Add x and y titles (common for all plots)
p6_ann <- annotate_figure(p6, left = text_grob("Frequency", size = 16, rot = 90), bottom = text_grob("Read position", size = 16))
p6_ann

# G. Save final plot
ggsave("plot_mapDamage_lgdist_FigS4.png", p6_ann, width = 205, height = 292, units = "mm", dpi = 300)
ggsave("plot_mapDamage_lgdist_FigS4.pdf", p6_ann, width = 205, height = 292, units = "mm")
