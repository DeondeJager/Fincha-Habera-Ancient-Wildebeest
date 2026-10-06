# qpGraph analysis with ADMIXTOOLS2 - subspecies
## Tutorial: https://uqrmaie1.github.io/admixtools/articles/admixtools.html
## find_graphs parameters are from: Maier, R. et al. (2023). On the limits of fitting complex models of population history to f-statistics. eLife, 12, e85492. https://doi.org/10.7554/eLife.85492 

# Set working directory
setwd("/projects/lorenzen/people/KUID/wildebeest/admixture_graph/BWD_HiC/admixtools2/subspecies")

# Load packages
#install.packages("devtools") # if "devtools" is not installed already
#devtools::install_github("uqrmaie1/admixtools")
library(admixtools)
library(magrittr)
library(tidyverse)
library(gridExtra)

# 1. Calculate f2 statistics & load results into R ----
prefix = "n82_fossil_BWD_HiC_autosomes.packedancestrymap"
f2_dir = "f2_results"
## Calculate f2 - only has to be done once, then can be loaded from file with f2_from_precomp()
#extract_f2(prefix, f2_dir,
#           adjust_pseudohaploid = TRUE,
#           auto_only = FALSE,
#           overwrite = TRUE)
## Load precomputed f2 results
f2_blocks = f2_from_precomp(f2_dir)

# 2. Perform fully automated graph exploration - to be compared afterwards ----
## The code below performs 100 independent find_graphs runs (the for loop) for a given number of admixture events and:
## i. stores the winning model (lowest LL) of each independent run in an object, and
## ii. plots the winning models of each independent run together on a single page with their scores, so we can evaluate how many independent runs converged on the same winning model.

## NEXT: The best model for each number of admix events must then be compared to all others using:
## a. Out-of-sample scores (to fairly compare models of different complexitites), and 
## b. Bootstrap resampling (of SNPs) for each graph fit to determine whether the different scores from a. are significantly different from each other.
## See: https://uqrmaie1.github.io/admixtools/articles/graphs.html#comparing-the-fits-of-different-graphs

## For 6 admix events ----
opt_results_6 <- list() # Initialize empty list to store results
winner_6 <- list() # Initialize list for winners
plots_6 <- list() # Initialize list for winners plots

for (i in 1:100) {
  opt_results_6[[i]] <- find_graphs(f2_blocks, 
                                    numadmix = 6, 
                                    outpop = "Hartebeest",
                                    mutfuns = namedList(spr_leaves, spr_all, swap_leaves, move_admixedge_once, flipadmix_random, place_root_random, mutate_n),
                                    numgraphs = 10,
                                    stop_gen = 10000,
                                    stop_gen2 = 30,
                                    opt_worst_residual = FALSE,
                                    reject_f4z = 0,
                                    diag = 1e-04,
                                    lsqmode = FALSE,
                                    verbose = TRUE)
  
  winner_6[[i]] <- opt_results_6[[i]] %>% slice_min(score, with_ties = FALSE)
  
  print(paste("Iteration", i, "winner score:", round(winner_6[[i]]$score[[1]], 3)), quote = FALSE)
  
  # Save each plot with title showing the score
  plots_6[[i]] <- plot_graph(winner_6[[i]]$edges[[1]],
                             title = paste("Score:", round(winner_6[[i]]$score[[1]], 3)))
}
## Combine plots into grid
combined_plots_gridArrange <- marrangeGrob(grobs = plots_6, ncol = 2, nrow = 5) # Use marrangeGrob for multiple page pdf
## Save grid as pdf for easy comparison across replicates
ggsave("winner_graphs_grid_6.pdf", combined_plots_gridArrange, width = 12, height = 18) # multiple page pdf
## Save models, winners and plots as Rdata files
save(opt_results_6, file = "opt_results_6.Rdata")
save(winner_6, file = "winner_6.Rdata")
save(plots_6, file = "plots_6.Rdata")

