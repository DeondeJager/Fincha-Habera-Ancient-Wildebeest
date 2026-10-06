# Calculate and plot heterozygosity of wildebeest samples
## Produces Fig. 2F

# Load packages
library(tidyverse)
library(ggpubr)
library(scales)

################################################################################
# Plot heterozygosity - see "Data prep" below for how to generate the input
################################################################################
setwd("./")

# Load in data
mean_het_manual <- read_csv("mean_het_forFossil_manual.csv", col_names = TRUE)
## This contains the het estimates across the genome for Fossil936 combined with the per-sample heterozygosity estimates of the modern samples, generated in the "Data prep" section at the end of this file.

## Re-order subspecies from north to south
mean_het_manual$Subspecies_common <- factor(mean_het_manual$Subspecies_common, 
                                            levels = c("Fossil936", "Western white-bearded", "Eastern white-bearded", "Nyassa", "Cookson", "Brindled", "Black"))

### Plot depth datasets separately
myPalette <- c("#f2933b", "#cace3d", "#8dc5a6", "#cc95ff", "#ff96ca", "#8585c4", "#5c5c5c")
high_depth <-  mean_het_manual %>%
  filter(SRR != "SRR27955277_trim_Mut") %>% 
  ggplot() +
    geom_jitter(aes(Subspecies_common, mean_het, colour = Subspecies_common), pch = 19, size = 2, alpha = 0.5) +
    geom_boxplot(aes(Subspecies_common, mean_het, fill = Subspecies_common, colour = Subspecies_common), alpha = 0.5) +
    #annotate("text", x = 2, y = 0.00065, label ="****", size = 8) +
    #annotate("text", x = 3, y = 0.00056, label ="****", size = 8) +
    #annotate("text", x = 4, y = 0.00048, label ="***", size = 8) +
    #annotate("text", x = 5, y = 0.00049, label ="**", size = 8) +
    #annotate("text", x = 6, y = 0.00068, label ="***", size = 8) +
    #annotate("text", x = 7, y = 0.00056, label ="****", size = 8) +
    scale_fill_discrete(type = myPalette) +
    scale_colour_discrete(type = myPalette) +
    theme_bw(base_size = 18) +
    theme(legend.position = "none") +
    scale_y_continuous(limits = c(0,7e-04), labels = scales::label_number()) +
    scale_x_discrete(guide = guide_axis(angle = 45),
                     labels = c("Fossil936", "Western\nwhite-bearded", "Eastern\nwhite-bearded", "Nyassa", "Cookson", "Brindled", "Black")) +
    labs(x = "Taxon", y = "Genome-wide heterozygosity")
high_depth
ggsave("plot_het_BWD_HiC_highDepth.png", high_depth, dpi = 300, height = 5, width = 5, units = "in")
ggsave("plot_het_BWD_HiC_highDepth.svg", high_depth, height = 5, width = 5, units = "in")


# Statistical analyses
## Wilcox tests
Black <- mean_het_manual %>%
  filter(Subspecies_common == "Black") %>%
  select(mean_het)
Brindled <- mean_het_manual %>%
  filter(Subspecies_common == "Brindled") %>%
  select(mean_het)
Cookson <- mean_het_manual %>%
  filter(Subspecies_common == "Cookson") %>%
  select(mean_het)
Nyassa <- mean_het_manual %>%
  filter(Subspecies_common == "Nyassa") %>%
  select(mean_het)
Ewb <- mean_het_manual %>%
  filter(Subspecies_common == "Eastern white-bearded") %>%
  select(mean_het)
Wwb <- mean_het_manual %>%
  filter(Subspecies_common == "Western white-bearded") %>%
  select(mean_het)
Fossil936 <- mean_het_manual %>%
  filter(Subspecies_common == "Fossil936") %>%
  select(mean_het)

wilcox.test(Fossil936$mean_het, Black$mean_het, paired = FALSE)
wilcox.test(Fossil936$mean_het, Brindled$mean_het, paired = FALSE)
wilcox.test(Fossil936$mean_het, Cookson$mean_het, paired = FALSE)
wilcox.test(Fossil936$mean_het, Nyassa$mean_het, paired = FALSE)
wilcox.test(Fossil936$mean_het, Ewb$mean_het, paired = FALSE)
wilcox.test(Fossil936$mean_het, Wwb$mean_het, paired = FALSE)


################################################################################
# Data prep - only need to do this once for the modern samples (high depth)
################################################################################
# Set working directory
## For high-depth analysis
setwd("./Liu2024_het/")

# Import data
## Read file names into list
het_files <- list.files(pattern = ".ml", full.names = FALSE)
## Load sample metadata
metadata <- read_csv("metadata.csv", col_names = TRUE) # High-depth

## Get sample names by making an empty vector, then removing "_BWD_HiC_noTrans.10Mb.ml" 
## from the filename and store this as the sample name
samples <- vector("character", length(het_files))
for (i in 1:length(het_files)) {
  samples[i] <- sub("_BWD_HiC_noTrans.10Mb.ml*", "", het_files[i])
}

# Calculate heterozygosity for each modern sample
## Loop through files and calculate proportion of heterozygous sites (Prop_het) for each window 
## and then the mean Prop_het across all windows, for each individual
list_mean <- list()
for (i in 1:length(het_files)){
  
  # Read data
  df <- read_table(het_files[i], , col_names = c("Hom1", "Het", "Hom2"))
  
  # Calculate heterozygosity for each row and add to data frame
  df <- df %>% mutate(Prop_het = Het/(Hom1+Het+Hom2))
  
  # Save mean in data frame
  list_mean[[i]] <- data.frame(SRR = samples[i],
                               mean_het = mean(df$Prop_het))
}
## Will give lots of warnings because there is a fourth empty column in each ml file that we ignore.
## Bind all results together and save as data frame
mean_het_df <- do.call(rbind, list_mean)

# Combine mean_het into data frame with subspecies info
mean_het_df <- left_join(mean_het_df, metadata, by = "SRR") %>%
  select("SRR", "mean_het", "Subspecies_common")

# Quick plot
ggplot(mean_het_df) +
  geom_boxplot(aes(Subspecies_common, mean_het, fill = Subspecies_common)) +
  geom_jitter(aes(Subspecies_common, mean_het, fill = Subspecies_common), pch = 21, size = 2)

## Plot the fossil wildebeest boxplot as representing het across windows in genome, as opposed to mean of subspecies
## Save mean het dataframe as csv to edit manually and then re-import to R
write_csv(mean_het_df, "mean_het_forFossil.csv")
## Manual editing was adding the Prop_het sites per window for Fossil936 to the CSV instead of using a single point mean het for it.
## This file is then saved as mean_het_forFossil_manual.csv and is used to make the final plots (see first section of script)

