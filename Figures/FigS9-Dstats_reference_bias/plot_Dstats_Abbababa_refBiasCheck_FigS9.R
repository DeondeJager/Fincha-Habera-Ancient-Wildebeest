# Plot D-stat results from reference bias check abbababa (individual-based Dstats) in ANGSD
## Produces Fig. S9.
## Note that SRR27955277_trim_Mut == 277_damaged in the code below and input files

# Set working directory
setwd("./")

# Load packages
library(tidyverse)
library(janitor)
library(ggpubr)

# Load data ----
## Goat
### noMaf, baq1 (1 chromosome)
n10_Goat.refBiasCheck.rmTrans.noMaf.baq1 <- read_table("n10_Goat.refBiasCheck.rmTrans.noMaf.baq1.jackknife.txt", col_names = T)
n10_Goat.refBiasCheck.rmTrans.noMaf.baq1 <- n10_Goat.refBiasCheck.rmTrans.noMaf.baq1 %>%
  remove_empty(which = "cols") # Remove empty column "X10"

## BWD_HiC
### noMaf, baq1 (1 chromosome)
n10_BWD_HiC.refBiasCheck.rmTrans.noMaf.baq1 <- read_table("n10_BWD_HiC.refBiasCheck.rmTrans.noMaf.baq1.jackknife.txt", col_names = T)
n10_BWD_HiC.refBiasCheck.rmTrans.noMaf.baq1 <- n10_BWD_HiC.refBiasCheck.rmTrans.noMaf.baq1 %>%
  remove_empty(which = "cols") # Remove empty column "X10"

# Plot Dstat by subspecies
## Goat ----
#--- noMaf, baq1 (1 chromosome) ---#
## Filter dataset for damaged modern in H1 and undamaged in H2 position
n10_Goat.refBiasCheck.rmTrans.noMaf.baq1.1 <- n10_Goat.refBiasCheck.rmTrans.noMaf.baq1 %>%
  filter(H1 == "SRR27955277_trim_Mut" & H2 == "SRR27955277")
## Add (sub)species
n10_Goat.refBiasCheck.rmTrans.noMaf.baq1.1 <- n10_Goat.refBiasCheck.rmTrans.noMaf.baq1.1 %>%
  mutate(taxon = case_match(H3, 
                            "SRR27955210" ~ "Black", 
                            "SRR27955211" ~ "E_Monduli", 
                            "SRR27955222" ~ "B_Etosha", 
                            "SRR27955228" ~ "E_Nairobi", 
                            "SRR27955239" ~ "W_Serengeti", 
                            "SRR27955244" ~ "Nyassa", 
                            "SRR27955347" ~ "Cookson"))
## Collapse into subspecies names (the E_*, W_*, B_* names were for previous analyses)
n10_Goat.refBiasCheck.rmTrans.noMaf.baq1.1 <- n10_Goat.refBiasCheck.rmTrans.noMaf.baq1.1 %>%
  filter(taxon != "E_Nairobi") %>%
  mutate(taxon = str_replace_all(taxon, c("E_Monduli" = "Eastern white-bearded", "W_Serengeti" = "Western white-bearded", "B_Etosha" = "Brindled")))

## Make h3 an ordered factor, for plotting
n10_Goat.refBiasCheck.rmTrans.noMaf.baq1.1$taxon <- factor(n10_Goat.refBiasCheck.rmTrans.noMaf.baq1.1$taxon, 
                                                        levels = c("Black", "Brindled", "Cookson", "Nyassa", "Western white-bearded", "Eastern white-bearded"))
Dplot_goat <- ggplot(n10_Goat.refBiasCheck.rmTrans.noMaf.baq1.1) +
  geom_vline(xintercept = 0, linetype = "dashed") +
  geom_point(aes(x = Dstat, y = taxon, colour = Z), size = 6) +
  geom_linerange(aes(x = Dstat, y = taxon, xmin = Dstat-SE, xmax = Dstat+SE)) +
  scale_colour_gradient2(limits = c(-2, 30), mid = "darkgrey") +
  #  scale_colour_manual(values = myPalette1,
#                      breaks = c("E_Amboseli", "E_Nairobi", "E_Monduli", "Nyassa", "Cookson", "B_Etosha", "Black"),
#                      name = "H2") + # Define colours 
#  guides(colour = "none") + # Remove legend
  theme_bw(base_size = 18) +
  scale_x_continuous(limits = c(-1,1), name ="D (O,H3;277,277_damaged)") +
  scale_y_discrete(name = "H3", labels = c("Black", "Brindled", "Cookson", "Nyassa", "Western\nwhite-bearded", "Eastern\nwhite-bearded")) +
  #  labs(subtitle = "Damaged modern at H1, Modern at H2, Ref = Goat") +
  theme(legend.position = "none")
Dplot_goat
ggsave("Goat/n10_Goat.refBiasCheck.rmTrans.noMaf.baq1_H1H2.png", dpi = 300)
  

## BWD_HiC ----
#--- noMaf, baq1 (1 chromosome) ---#
## Filter dataset for damaged modern in H1 and undamaged in H2 position
n10_BWD_HiC.refBiasCheck.rmTrans.noMaf.baq1.1 <- n10_BWD_HiC.refBiasCheck.rmTrans.noMaf.baq1 %>%
  filter(H1 == "SRR27955277_trim_Mut" & H2 == "SRR27955277")
## Add (sub)species
n10_BWD_HiC.refBiasCheck.rmTrans.noMaf.baq1.1 <- n10_BWD_HiC.refBiasCheck.rmTrans.noMaf.baq1.1 %>%
  mutate(taxon = case_match(H3, 
                            "SRR27955210" ~ "Black", 
                            "SRR27955211" ~ "E_Monduli", 
                            "SRR27955222" ~ "B_Etosha", 
                            "SRR27955228" ~ "E_Nairobi", 
                            "SRR27955239" ~ "W_Serengeti", 
                            "SRR27955244" ~ "Nyassa", 
                            "SRR27955347" ~ "Cookson"))
## Collapse into subspecies names (the E_*, W_*, B_* names were for previous analyses)
n10_BWD_HiC.refBiasCheck.rmTrans.noMaf.baq1.1 <- n10_BWD_HiC.refBiasCheck.rmTrans.noMaf.baq1.1 %>%
  filter(taxon != "E_Nairobi") %>%
  mutate(taxon = str_replace_all(taxon, c("E_Monduli" = "Eastern white-bearded", "W_Serengeti" = "Western white-bearded", "B_Etosha" = "Brindled")))

## Make h3 an ordered factor, for plotting
n10_BWD_HiC.refBiasCheck.rmTrans.noMaf.baq1.1$taxon <- factor(n10_BWD_HiC.refBiasCheck.rmTrans.noMaf.baq1.1$taxon, 
                                                           levels = c("Black", "Brindled", "Cookson", "Nyassa", "Western white-bearded", "Eastern white-bearded"))

Dplot_BWD_HiC <- ggplot(n10_BWD_HiC.refBiasCheck.rmTrans.noMaf.baq1.1) +
  geom_vline(xintercept = 0, linetype = "dashed") +
  geom_point(aes(x = Dstat, y = taxon, colour = Z), size = 6) +
  geom_linerange(aes(x = Dstat, y = taxon, xmin = Dstat-SE, xmax = Dstat+SE)) +
  scale_colour_gradient2(limits = c(-2, 30), mid = "darkgrey") +
  #  scale_colour_manual(values = myPalette1,
  #                      breaks = c("E_Amboseli", "E_Nairobi", "E_Monduli", "Nyassa", "Cookson", "B_Etosha", "Black"),
  #                      name = "H2") + # Define colours 
  #  guides(colour = "none") + # Remove legend
  theme_bw(base_size = 18) +
  scale_x_continuous(limits = c(-1,1), name ="D (O,H3;277,277_damaged)") +
  scale_y_discrete(name = "H3", labels = c("Black", "Brindled", "Cookson", "Nyassa", "Western\nwhite-bearded", "Eastern\nwhite-bearded")) +
#  labs(subtitle = "Damaged modern at H1, Modern at H2, Ref = BWD_HiC") +
  theme(legend.position = "none")
Dplot_BWD_HiC
ggsave("BWD_HiC/n10_BWD_HiC.refBiasCheck.rmTrans.noMaf.baq1_H1H2.png", dpi = 300)
Dplot_legend <- get_legend(Dplot_BWD_HiC) # Change theme(legend.position = "none") to theme(legend.position = "bottom") in previous code section and re-run for this to work

## Combine Dstat plots
Dstats_plots <- ggarrange(Dplot_goat, Dplot_BWD_HiC, ncol = 1, nrow = 2, legend.grob = Dplot_legend, legend = "bottom",
                          labels = c("C", "D"), font.label = list(size = 18))
Dstats_plots
ggsave("n10_Dstats.refBiasCheck.rmTrans.noMaf.baq1_H1H2.png", Dstats_plots, dpi = 300)


# Plot Dstat between damaged and undamaged modern
## Goat
### Filter for rows containing only undamaged modern (no damaged modern involved)
undamaged <- n10_Goat.refBiasCheck.rmTrans.noMaf.baq1 %>%
  filter(H1 != "SRR27955277_trim_Mut" & H2 != "SRR27955277_trim_Mut" & H3 != "SRR27955277_trim_Mut") %>%
  filter(H1 == "SRR27955277" | H2 == "SRR27955277" | H3 == "SRR27955277") %>%
  filter(H1 != "SRR27955228" & H2 != "SRR27955228" & H3 != "SRR27955228") # Removed SRR27955228 (E_Nairobi), to have the same data as the plots above
### Filter for rows containing only damaged modern (no undamaged modern involved)
damaged <- n10_Goat.refBiasCheck.rmTrans.noMaf.baq1 %>%
  filter(H1 != "SRR27955277" & H2 != "SRR27955277" & H3 != "SRR27955277") %>%
  filter(H1 == "SRR27955277_trim_Mut" | H2 == "SRR27955277_trim_Mut" | H3 == "SRR27955277_trim_Mut") %>%
  filter(H1 != "SRR27955228" & H2 != "SRR27955228" & H3 != "SRR27955228") # Removed SRR27955228 (E_Nairobi), to have the same data as the plots above

### Prepend text to the headers of each dataframe
names(undamaged) <- paste("Original", names(undamaged), sep = "_")
names(damaged) <- paste("Simulated", names(damaged), sep = "_")
### Combine the data
combined_data <- cbind(undamaged, damaged)
#write_csv(combined_data, "goat_combined_Dstats.csv")
### Plot
goat_plot <- ggplot(combined_data, aes(Original_Dstat, Simulated_Dstat, fill = Simulated_Dstat - Original_Dstat)) +
  geom_point(shape = 21, size = 3) +
#  geom_text(aes(label = Original_H3), vjust = -0.6, size = 3) +
  scale_x_continuous(limits = c(-1,1), name = "277 Dstat") +
  scale_y_continuous(limits = c(-1,1), name = "277_damaged Dstat") +
  scale_fill_gradient2(low = "blue", mid = "gray70", high = "red",
                       midpoint = 0, name = expression(Delta~Dstat~"(Dmg - Undmg)"),
                       limits = c(-0.05, 0.3)) +
  geom_abline(intercept = 0, slope = 1, linetype = "dashed") +
  theme_bw(base_size = 18) +
  theme(legend.position = "none") #+
  #labs(subtitle = "Goat reference, chr1")
goat_plot
ggsave("n10_Goat.refBiasCheck.rmTrans.noMaf.baq1.png", goat_plot, dpi = 300)

## BWD_HiC
### Filter for rows containing only undamaged modern (no damaged modern involved)
undamaged <- n10_BWD_HiC.refBiasCheck.rmTrans.noMaf.baq1 %>%
  filter(H1 != "SRR27955277_trim_Mut" & H2 != "SRR27955277_trim_Mut" & H3 != "SRR27955277_trim_Mut") %>%
  filter(H1 == "SRR27955277" | H2 == "SRR27955277" | H3 == "SRR27955277") %>%
  filter(H1 != "SRR27955228" & H2 != "SRR27955228" & H3 != "SRR27955228") # Removed SRR27955228 (E_Nairobi), to have the same data as the plots above
### Filter for rows containing only damaged modern (no undamaged modern involved)
damaged <- n10_BWD_HiC.refBiasCheck.rmTrans.noMaf.baq1 %>%
  filter(H1 != "SRR27955277" & H2 != "SRR27955277" & H3 != "SRR27955277") %>%
  filter(H1 == "SRR27955277_trim_Mut" | H2 == "SRR27955277_trim_Mut" | H3 == "SRR27955277_trim_Mut") %>%
  filter(H1 != "SRR27955228" & H2 != "SRR27955228" & H3 != "SRR27955228") # Removed SRR27955228 (E_Nairobi), to have the same data as the plots above
### Prepend text to the headers of each dataframe
names(undamaged) <- paste("Original", names(undamaged), sep = "_")
names(damaged) <- paste("Simulated", names(damaged), sep = "_")
### Combine the data
combined_data <- cbind(undamaged, damaged)
#write_csv(combined_data, "BWD_HiC_combined_Dstats.csv")
### Plot
BWD_HiC_plot <- ggplot(combined_data, aes(Original_Dstat, Simulated_Dstat, fill = Simulated_Dstat - Original_Dstat)) +
  geom_point(shape = 21, size = 3) +
  #geom_text(aes(label = Original_H3), vjust = -0.6, size = 3) +
  scale_x_continuous(limits = c(-1,1), name = "277 Dstat") +
  scale_y_continuous(limits = c(-1,1), name = "277_damaged Dstat") +
  scale_fill_gradient2(low = "blue", mid = "gray70", high = "red",
                       midpoint = 0, name = expression(Delta~Dstat~"(Dmg - Undmg)"),
                       limits = c(-0.05, 0.3)) +
  geom_abline(intercept = 0, slope = 1, linetype = "dashed") +
  theme_bw(base_size = 18) +
  theme(legend.position = "none") #+
  #labs(subtitle = "BWD_HiC reference, HiC_scaffold_1")
BWD_HiC_plot
ggsave("n10_BWD_HiC.refBiasCheck.rmTrans.noMaf.baq1.png", BWD_HiC_plot, dpi = 300)
common_legend <- get_legend(BWD_HiC_plot) # Change theme(legend.position = "none") to theme(legend.position = "bottom") in previous code section and re-run for this to work

## Combine compare plots
compare_plots <- ggarrange(goat_plot, BWD_HiC_plot, ncol = 1, nrow = 2, legend.grob = common_legend, legend = "bottom",
                           labels = c("A", "B"), font.label = list(size = 18))
compare_plots
ggsave("n10_compare.refBiasCheck.rmTrans.noMaf.baq1.png", compare_plots, dpi = 300)

# Combine all plots
ggarrange(compare_plots, Dstats_plots, ncol = 2, nrow = 1)
ggsave("FigureSXV_RefBias.png", dpi = 300, bg = "white")
ggsave("FigureSXV_RefBias.svg", width = 10.2, height = 7.92, bg = "white")
