# Get scores and do tests of graphs after running find_graphs with ADMIXTOOLS2
## Got help from ChatGPT with this code
## Tutorial: https://uqrmaie1.github.io/admixtools/articles/admixtools.html
## find_graphs parameters are from: Maier, R. et al. (2023). On the limits of fitting complex models of population history to f-statistics. eLife, 12, e85492. https://doi.org/10.7554/eLife.85492 

# This script analyses results and produces plots for investigation. 
# For the final plots included in the manuscript, see the script "admixtools2_qpGraph_BWD_HiC_plots_final_Fig3_FigS14.R" in the "Figures" subfolder of the GitHub repo: https://github.com/DeondeJager/Fincha-Habera-Ancient-Wildebeest/tree/main/Figures

# Set working directory
setwd("./")

# Load packages
#install.packages("devtools") # if "devtools" is not installed already
#devtools::install_github("uqrmaie1/admixtools")
library(admixtools)
library(magrittr)
library(tidyverse)
library(gridExtra)
library(igraph)
#install.packages("shiny", dependencies = TRUE)
#library(shiny)

# Set high number of digits for printing precision, as some likelihood scores appear identical for the first few decimals, but are actually different.
options(digits = 16)

# 1. Calculate/load f2 statistics & load results into R ----
prefix = "n82_fossil_BWD_HiC_autosomes.packedancestrymap"
f2_dir = "f2_results"
#extract_f2(prefix, f2_dir,
#           adjust_pseudohaploid = TRUE,
#           auto_only = FALSE,
#           overwrite = TRUE)
## Load results
f2_blocks = f2_from_precomp(f2_dir)

# 2. Compare models within complexity levels (number of admixture events) ----
## Load results of model exploration
load(file = 'winner_0.Rdata')
load(file = 'winner_1.Rdata')
load(file = 'winner_2.Rdata')
load(file = 'winner_3.Rdata')
load(file = 'winner_4.Rdata')
load(file = 'winner_5.Rdata')
load(file = 'winner_6.Rdata')

# 3. Deduplicate models based on identical graph structures, using their hash values ----
deduplicate_graph_list <- function(winner_list) {
  # Step 1: Compute hashes
  hashes <- map_chr(winner_list, ~ graph_hash(.x$graph[[1]]))
  
  # Step 2: Find indices of first occurrence of each unique hash
  unique_indices <- match(unique(hashes), hashes)
  
  # Step 3: Keep only the unique graphs
  winner_list[unique_indices]
}
# Check the deduplication worked
winner_0_unique <- deduplicate_graph_list(winner_0)
length(winner_0)        # 100 original graphs
length(winner_0_unique) # Fewer (only unique ones)

winner_1_unique <- deduplicate_graph_list(winner_1)
length(winner_1)        # 100 original graphs
length(winner_1_unique) # Fewer (only unique ones)

winner_2_unique <- deduplicate_graph_list(winner_2)
length(winner_2)        # 100 original graphs
length(winner_2_unique) # Fewer (only unique ones)

winner_3_unique <- deduplicate_graph_list(winner_3)
length(winner_3)        # 100 original graphs
length(winner_3_unique) # Fewer (only unique ones)

winner_4_unique <- deduplicate_graph_list(winner_4)
length(winner_4)        # 100 original graphs
length(winner_4_unique) # Fewer (only unique ones)

winner_5_unique <- deduplicate_graph_list(winner_5)
length(winner_5)        # 100 original graphs
length(winner_5_unique) # Fewer (only unique ones)

winner_6_unique <- deduplicate_graph_list(winner_6)
length(winner_6)        # 100 original graphs
length(winner_6_unique) # Fewer (only unique ones)

# 4. Rank unique models by score and sort them by rank using the sort_and_rank_by_score function defined below ----
## sort_and_rank_by_score function - takes the dataframe list (e.g. winner_0) as input
sort_and_rank_by_score <- function(df_list) {
  # Step 1: Add original index to each tibble
  df_list_with_index <- lapply(seq_along(df_list), function(i) {
    x <- df_list[[i]]
    x$orig_index <- i
    return(x)
  })
  
  # Step 2: Extract score and orig_index
  scores <- sapply(df_list_with_index, function(x) x$score)
  orig_indices <- sapply(df_list_with_index, function(x) x$orig_index)
  
  # Step 3: Determine sort order
  ordering <- order(scores, orig_indices)
  
  # Step 4: Sort and assign rank
  ranked_list <- lapply(seq_along(ordering), function(rank) {
    x <- df_list_with_index[[ordering[rank]]]
    x$rank <- rank
    return(x)
  })
  
  return(ranked_list)
}

## 0 admix events 
winner_0_ranked <- sort_and_rank_by_score(winner_0_unique)
#### Plot models
plots_0 <- list() # Initialize list for models plots
for (i in 1:length(winner_0_ranked)){
  plots_0[[i]] <- plot_graph(winner_0_ranked[[i]]$edges[[1]],
                             title = paste("0 admix events\n",
                                           "Score:", winner_0_ranked[[i]]$score[[1]], 
                                           "Original index:", winner_0_ranked[[i]]$orig_index[[1]], 
                                           "Rank:", winner_0_ranked[[i]]$rank[[1]]))
}
plots_0
plots_arrange <- marrangeGrob(grobs = plots_0, ncol = 2, nrow = 3) 
## Save grid as pdf for easy comparison across replicates
ggsave("winner_0_ranked.pdf", plots_arrange, width = 12, height = 18) 
#### Store scores in a tibble and save as file
winner_0_ranked_summary_df <- winner_0_ranked %>%
  map_df(~ tibble(
    score = format(.x$score, digits = 16, scientific = FALSE),
    orig_index = round(.x$orig_index),
    rank = round(.x$rank)
  ))
write_csv(winner_0_ranked_summary_df, "winner_0_ranked_summary.csv")

### 1 admix events 
winner_1_ranked <- sort_and_rank_by_score(winner_1_unique)
#### Plot models
plots_1 <- list() # Initialize list for models plots
for (i in 1:length(winner_1_ranked)){
  plots_1[[i]] <- plot_graph(winner_1_ranked[[i]]$edges[[1]],
                             title = paste("1 admix events\n",
                                           "Score:", winner_1_ranked[[i]]$score[[1]], 
                                           "Original index:", winner_1_ranked[[i]]$orig_index[[1]], 
                                           "Rank:", winner_1_ranked[[i]]$rank[[1]]))
}
plots_1
plots_arrange <- marrangeGrob(grobs = plots_1, ncol = 2, nrow = 3) 
## Save grid as pdf for easy comparison across replicates
ggsave("winner_1_ranked.pdf", plots_arrange, width = 12, height = 18) 
#### Store scores in a tibble and save as file
winner_1_ranked_summary_df <- winner_1_ranked %>%
  map_df(~ tibble(
    score = format(.x$score, digits = 16, scientific = FALSE),
    orig_index = round(.x$orig_index),
    rank = round(.x$rank)
  ))
write_csv(winner_1_ranked_summary_df, "winner_1_ranked_summary.csv")

### 2 admix events 
winner_2_ranked <- sort_and_rank_by_score(winner_2_unique)
#### Plot models
plots_2 <- list() # Initialize list for models plots
for (i in 1:length(winner_2_ranked)){
  plots_2[[i]] <- plot_graph(winner_2_ranked[[i]]$edges[[1]],
                             title = paste("2 admix events\n",
                                           "Score:", winner_2_ranked[[i]]$score[[1]], 
                                           "Original index:", winner_2_ranked[[i]]$orig_index[[1]], 
                                           "Rank:", winner_2_ranked[[i]]$rank[[1]]))
}
plots_2
plots_arrange <- marrangeGrob(grobs = plots_2, ncol = 2, nrow = 3) 
## Save grid as pdf for easy comparison across replicates
ggsave("winner_2_ranked.pdf", plots_arrange, width = 12, height = 18) 
#### Store scores in a tibble and save as file
winner_2_ranked_summary_df <- winner_2_ranked %>%
  map_df(~ tibble(
    score = format(.x$score, digits = 16, scientific = FALSE),
    orig_index = round(.x$orig_index),
    rank = round(.x$rank)
  ))
write_csv(winner_2_ranked_summary_df, "winner_2_ranked_summary.csv")

### 3 admix events 
winner_3_ranked <- sort_and_rank_by_score(winner_3_unique)
#### Plot models
plots_3 <- list() # Initialize list for models plots
for (i in 1:length(winner_3_ranked)){
  plots_3[[i]] <- plot_graph(winner_3_ranked[[i]]$edges[[1]],
                             title = paste("3 admix events\n",
                                           "Score:", winner_3_ranked[[i]]$score[[1]], 
                                           "Original index:", winner_3_ranked[[i]]$orig_index[[1]], 
                                           "Rank:", winner_3_ranked[[i]]$rank[[1]]))
}
plots_3
plots_arrange <- marrangeGrob(grobs = plots_3, ncol = 2, nrow = 3) 
## Save grid as pdf for easy comparison across replicates
ggsave("winner_3_ranked.pdf", plots_arrange, width = 12, height = 18) 
#### Store scores in a tibble and save as file
winner_3_ranked_summary_df <- winner_3_ranked %>%
  map_df(~ tibble(
    score = format(.x$score, digits = 16, scientific = FALSE),
    orig_index = round(.x$orig_index),
    rank = round(.x$rank)
  ))
write_csv(winner_3_ranked_summary_df, "winner_3_ranked_summary.csv")

### 4 admix events 
winner_4_ranked <- sort_and_rank_by_score(winner_4_unique)
#### Plot models
plots_4 <- list() # Initialize list for models plots
for (i in 1:length(winner_4_ranked)){
  plots_4[[i]] <- plot_graph(winner_4_ranked[[i]]$edges[[1]],
                             title = paste("4 admix events\n",
                                           "Score:", winner_4_ranked[[i]]$score[[1]], 
                                           "Original index:", winner_4_ranked[[i]]$orig_index[[1]], 
                                           "Rank:", winner_4_ranked[[i]]$rank[[1]]))
}
plots_4
plots_arrange <- marrangeGrob(grobs = plots_4, ncol = 2, nrow = 3) 
## Save grid as pdf for easy comparison across replicates
ggsave("winner_4_ranked.pdf", plots_arrange, width = 12, height = 18) 
#### Store scores in a tibble and save as file
winner_4_ranked_summary_df <- winner_4_ranked %>%
  map_df(~ tibble(
    score = format(.x$score, digits = 16, scientific = FALSE),
    orig_index = round(.x$orig_index),
    rank = round(.x$rank)
  ))
write_csv(winner_4_ranked_summary_df, "winner_4_ranked_summary.csv")

### 5 admix events 
winner_5_ranked <- sort_and_rank_by_score(winner_5_unique)
#### Plot models
plots_5 <- list() # Initialize list for models plots
for (i in 1:length(winner_5_ranked)){
  plots_5[[i]] <- plot_graph(winner_5_ranked[[i]]$edges[[1]],
                             title = paste("5 admix events\n",
                                           "Score:", winner_5_ranked[[i]]$score[[1]], 
                                           "Original index:", winner_5_ranked[[i]]$orig_index[[1]], 
                                           "Rank:", winner_5_ranked[[i]]$rank[[1]]))
}
plots_5
plots_arrange <- marrangeGrob(grobs = plots_5, ncol = 2, nrow = 3) 
## Save grid as pdf for easy comparison across replicates
ggsave("winner_5_ranked.pdf", plots_arrange, width = 12, height = 18) 
#### Store scores in a tibble and save as file
winner_5_ranked_summary_df <- winner_5_ranked %>%
  map_df(~ tibble(
    score = format(.x$score, digits = 16, scientific = FALSE),
    orig_index = round(.x$orig_index),
    rank = round(.x$rank)
  ))
write_csv(winner_5_ranked_summary_df, "winner_5_ranked_summary.csv")

### 6 admix events 
winner_6_ranked <- sort_and_rank_by_score(winner_6_unique)
### Plot models
plots_6 <- list() # Initialize list for models plots
for (i in 1:length(winner_6_ranked)){
  plots_6[[i]] <- plot_graph(winner_6_ranked[[i]]$edges[[1]],
                             title = paste("6 admix events\n",
                                           "Score:", winner_6_ranked[[i]]$score[[1]], 
                                           "Original index:", winner_6_ranked[[i]]$orig_index[[1]], 
                                           "Rank:", winner_6_ranked[[i]]$rank[[1]]))
}
plots_6
plots_arrange <- marrangeGrob(grobs = plots_6, ncol = 2, nrow = 3) 
## Save grid as pdf for easy comparison across replicates
ggsave("winner_6_ranked.pdf", plots_arrange, width = 12, height = 18) 
#### Store scores in a tibble and save as file
winner_6_ranked_summary_df <- winner_6_ranked %>%
  map_df(~ tibble(
    score = format(.x$score, digits = 16, scientific = FALSE),
    orig_index = round(.x$orig_index),
    rank = round(.x$rank)
  ))
write_csv(winner_6_ranked_summary_df, "winner_6_ranked_summary.csv")


# 5. Do statistical comparisons of best/all models for each complexity (# admix events) ----
## Sample f2_blocks
nblocks = dim(f2_blocks)[3]
train = sample(1:nblocks, round(nblocks/2))

# 5.1. Get out-of-sample scores for all unique models
## List of names of all winner lists
winner_names <- paste0("winner_", 0:6, "_ranked")

# Initialize empty list to store all OOS results
oos_results <- list()

# Loop through each winner list and evaluate OOS scores
for (name in winner_names) {
  # Get the actual list object
  winner_list <- get(name)
  
  # Run qpgraph on each model in the list
  oos_list <- map(winner_list, ~ qpgraph(
    data = f2_blocks[, , train],
    graph = .x$edges[[1]],
    f2_blocks_test = f2_blocks[, , -train]
  ))
  
  # Store under dynamic name 
  assign(
    paste0("oos_", sub("winner_", "", name)),
    oos_list
  )
  
  # Combine OOS results with metadata (orig_index, rank)
  oos_df <- map2_dfr(oos_list, winner_list, ~ tibble(
    model = name,
    orig_index = .y$orig_index,
    rank = .y$rank,
    orig_score = .y$score,
    oos_score = .x$score,
    oos_score_test = .x$score_test
  ))
  
  # Save into a cumulative list
  oos_results[[name]] <- oos_df
}

# Combine all model results into one final dataframe
final_oos_scores <- bind_rows(oos_results)
write_tsv(final_oos_scores, "final_oos_scores.tsv")

## 5.2. Calculate worst residuals
# Loop through each winner group (0 to 6)
for (i in 0:6) {
  
  winner_list <- get(paste0("winner_", i, "_ranked"))
  n_models <- length(winner_list)
  
  # Run qpgraph for each model, saving the result (with return_fstats = TRUE)
  qpgraph_results <- map(winner_list, ~ qpgraph(
    data = f2_blocks,
    graph = .x$edges[[1]],
    numstart = 1000,
    return_fstats = TRUE
  ))
  
  # Save the results as an RData file
  save(qpgraph_results, file = paste0("winner_", i, "_ranked_qpGraphs.Rdata"))
  
  # Extract worst residuals
  worst_residuals <- tibble(
    model = paste0("W", i, "_M", seq_len(n_models)),
    worst_residual = map_dbl(qpgraph_results, ~ .x$worst_residual)
  )
  
  # Save worst residuals to CSV
  write_csv(worst_residuals, paste0("winner_", i, "_worst_residuals.csv"))
}

#--- Step 5.3 takes very long and was submitted as a 14-day job on the cluster, which still did not finish and so had to be run for each admix level separately by changing "6" to "5" to "4", etc on lines 364-365 below.

# 5.3. Bootstrap-resampled graph fits to statistically test whether different models within each winner group are significantly different or not
# Group model IDs and models by winner category (e.g. W0, W1, ..., W6)
#winner_ranked_lists <- map(0:6, ~ get(paste0("winner_", ., "_ranked")))
#names(winner_ranked_lists) <- paste0("W", 0:6)
winner_ranked_lists <- map(6, ~ get(paste0("winner_", ., "_ranked"))) # 6 admix events only
names(winner_ranked_lists) <- paste0("W", 6) # 6 admix events only

# Flatten all models into one list
all_models <- flatten(winner_ranked_lists)
# Generate model IDs like W0_M1, W0_M2, ..., W6_Mn depending on number of models per winner_X
model_ids <- unlist(
  map2(
    names(winner_ranked_lists), 
    map_int(winner_ranked_lists, length), 
    ~ paste0(.x, "_M", seq_len(.y))
  )
)

names(all_models) <- model_ids
# Loop through each winner category
for (winner_name in names(winner_ranked_lists)) {
  # Filter model IDs for this winner only
  models_in_category <- model_ids[startsWith(model_ids, winner_name)]
  
  # Generate all pairs within this category
  model_pairs <- combn(models_in_category, 2, simplify = FALSE)
  
  # Compare pairs and extract p_emp values
  pemp_results <- map_dfr(model_pairs, function(pair) {
    m1 <- all_models[[pair[1]]]$edges[[1]]
    m2 <- all_models[[pair[2]]]$edges[[1]]
    
    fits <- qpgraph_resample_multi(f2_blocks, list(m1, m2), nboot = 100)
    fits_compared <- compare_fits(fits[[1]]$score_test, fits[[2]]$score_test)
    
    tibble(
      model1 = pair[1],
      model2 = pair[2],
      p_emp = fits_compared$p_emp
    )
  })
  
  # Make symmetric matrix for plotting
  pemp_plot_data <- bind_rows(
    pemp_results,
    pemp_results %>% rename(model1 = model2, model2 = model1)
  )
  
  diag_entries <- tibble(
    model1 = models_in_category,
    model2 = models_in_category,
    p_emp = NA_real_
  )
  
  pemp_full <- bind_rows(pemp_plot_data, diag_entries)
  
  # Save long and wide CSVs for this winner category
  write_csv(pemp_full, paste0(winner_name, "_pemp_matrix_long.csv"))
  
  pemp_wide <- pemp_full %>% pivot_wider(names_from = model2, values_from = p_emp)
  write_csv(pemp_wide, paste0(winner_name, "_pemp_matrix_square.csv"))
  
  # Plot heatmap
  plot <- ggplot(pemp_full, aes(x = model1, y = model2, fill = p_emp)) +
    geom_tile(color = "white") +
    scale_fill_viridis_c(option = "plasma", direction = -1, na.value = "grey90", name = "p_emp") +
    theme_minimal() +
    coord_fixed() +
    theme(
      axis.text.x = element_text(angle = 90, vjust = 0.5, hjust = 1),
      axis.title = element_blank()
    ) +
    geom_text(aes(label = round(p_emp, 2)), size = 2)
  
  ggsave(paste0(winner_name, "_pemp_matrix_plot.pdf"), plot = plot,
         width = 180, height = 180, units = "mm")
}

