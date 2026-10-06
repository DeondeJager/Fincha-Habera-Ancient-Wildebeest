# Plot admixture models for manuscript

# Produces Fig. 3A-B (at the very end of the file) and plots for Fig. S14.

# Set working directory
setwd("./")

# Load packages
library(admixtools)
library(magrittr)
library(tidyverse)
library(gridExtra)
library(igraph)

# Load data
model_stats <- read_csv("final_oos_scores_worst_residuals.csv", col_names = T)

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

### 1 admix events 
winner_1_ranked <- sort_and_rank_by_score(winner_1_unique)

### 2 admix events 
winner_2_ranked <- sort_and_rank_by_score(winner_2_unique)

### 3 admix events 
winner_3_ranked <- sort_and_rank_by_score(winner_3_unique)

### 4 admix events 
winner_4_ranked <- sort_and_rank_by_score(winner_4_unique)

### 5 admix events 
winner_5_ranked <- sort_and_rank_by_score(winner_5_unique)

### 6 admix events 
winner_6_ranked <- sort_and_rank_by_score(winner_6_unique)


# Plot all----
# Group model IDs and models by winner category (e.g. W0, W1, ..., W6) 
winner_ranked_lists <- map(0:6, ~ get(paste0("winner_", ., "_ranked")))
names(winner_ranked_lists) <- paste0("W", 0:6)
# Flatten all models into one list
all_models <- flatten(winner_ranked_lists)
# Generate model IDs like W0_M1, W0_M2, ..., W6_Mn depending on number of models per winner_X
model_ids <- unlist(map2(names(winner_ranked_lists),
                         map_int(winner_ranked_lists, length), ~ paste0(.x, "_M", seq_len(.y))))
names(all_models) <- model_ids

# Plot best model per complexity level, and all other models that are not significantly different from it (based on pemp results)
# 0 admix events
nodiff_models_0 <- all_models[c("W0_M1", "W0_M2")]
nodiff_models_0_stats <- model_stats %>% 
  filter(model %in% c("W0_M1", "W0_M2"))
plots_0 <- list()
for (i in 1:length(nodiff_models_0)){
  plots_0[[i]] <- plot_graph(nodiff_models_0[[i]]$edges[[1]],
                             title = paste(nodiff_models_0_stats$model[[i]],"\n",
                                           "Score:", round(nodiff_models_0_stats$orig_score[[i]], 2), 
                                           "Out-of-sample score:", round(nodiff_models_0_stats$oos_score[[i]], 2), 
                                           "Worst-residual:", round(nodiff_models_0_stats$worst_residual[[i]], 2)))
}
plots_0
## Save as pdf
plots_arrange <- marrangeGrob(grobs = plots_0, ncol = 2, nrow = 3) 
ggsave("final_nodiff_models_0.pdf", plots_arrange, width = 12, height = 18)

# 1 admix events
nodiff_models_1 <- all_models[c("W1_M1", "W1_M2", "W1_M3", "W1_M4")]
nodiff_models_1_stats <- model_stats %>% 
  filter(model %in% c("W1_M1", "W1_M2", "W1_M3", "W1_M4"))
plots_1 <- list()
for (i in 1:length(nodiff_models_1)){
  plots_1[[i]] <- plot_graph(nodiff_models_1[[i]]$edges[[1]],
                             title = paste(nodiff_models_1_stats$model[[i]],"\n",
                                           "Score:", round(nodiff_models_1_stats$orig_score[[i]], 2), 
                                           "Out-of-sample score:", round(nodiff_models_1_stats$oos_score[[i]], 2), 
                                           "Worst-residual:", round(nodiff_models_1_stats$worst_residual[[i]], 2)))
}
plots_1
## Save as pdf
plots_arrange <- marrangeGrob(grobs = plots_1, ncol = 2, nrow = 3) 
ggsave("final_nodiff_models_1.pdf", plots_arrange, width = 12, height = 18)

# 2 admix events
nodiff_models_2 <- all_models[c("W2_M1", "W2_M2", "W2_M3")]
nodiff_models_2_stats <- model_stats %>% 
  filter(model %in% c("W2_M1", "W2_M2", "W2_M3"))
plots_2 <- list()
for (i in 1:length(nodiff_models_2)){
  plots_2[[i]] <- plot_graph(nodiff_models_2[[i]]$edges[[1]],
                             title = paste(nodiff_models_2_stats$model[[i]],"\n",
                                           "Score:", round(nodiff_models_2_stats$orig_score[[i]], 2), 
                                           "Out-of-sample score:", round(nodiff_models_2_stats$oos_score[[i]], 2), 
                                           "Worst-residual:", round(nodiff_models_2_stats$worst_residual[[i]], 2)))
}
plots_2
## Save as pdf
plots_arrange <- marrangeGrob(grobs = plots_2, ncol = 2, nrow = 3) 
ggsave("final_nodiff_models_2.pdf", plots_arrange, width = 12, height = 18)

# 3 admix events
nodiff_models_3 <- all_models[c("W3_M1", "W3_M2", "W3_M3")]
nodiff_models_3_stats <- model_stats %>% 
  filter(model %in% c("W3_M1", "W3_M2", "W3_M3"))
plots_3 <- list()
for (i in 1:length(nodiff_models_3)){
  plots_3[[i]] <- plot_graph(nodiff_models_3[[i]]$edges[[1]],
                             title = paste(nodiff_models_3_stats$model[[i]],"\n",
                                           "Score:", round(nodiff_models_3_stats$orig_score[[i]], 2), 
                                           "Out-of-sample score:", round(nodiff_models_3_stats$oos_score[[i]], 2), 
                                           "Worst-residual:", round(nodiff_models_3_stats$worst_residual[[i]], 2)))
}
plots_3
## Save as pdf
plots_arrange <- marrangeGrob(grobs = plots_3, ncol = 2, nrow = 3) 
ggsave("final_nodiff_models_3.pdf", plots_arrange, width = 12, height = 18)

# 4 admix events
nodiff_models_4 <- all_models[c("W4_M1", "W4_M2", "W4_M3", "W4_M4", "W4_M5", "W4_M6", "W4_M7", "W4_M8", "W4_M9", "W4_M10", "W4_M11", "W4_M12",
                                "W4_M13", "W4_M14", "W4_M15", "W4_M16", "W4_M17", "W4_M19", "W4_M20", "W4_M22")]
nodiff_models_4_stats <- model_stats %>% 
  filter(model %in% c("W4_M1", "W4_M2", "W4_M3", "W4_M4", "W4_M5", "W4_M6", "W4_M7", "W4_M8", "W4_M9", "W4_M10", "W4_M11", "W4_M12",
                      "W4_M13", "W4_M14", "W4_M15", "W4_M16", "W4_M17", "W4_M19", "W4_M20", "W4_M22"))
plots_4 <- list()
for (i in 1:length(nodiff_models_4)){
  plots_4[[i]] <- plot_graph(nodiff_models_4[[i]]$edges[[1]],
                             title = paste(nodiff_models_4_stats$model[[i]],"\n",
                                           "Score:", round(nodiff_models_4_stats$orig_score[[i]], 2), 
                                           "Out-of-sample score:", round(nodiff_models_4_stats$oos_score[[i]], 2), 
                                           "Worst-residual:", round(nodiff_models_4_stats$worst_residual[[i]], 2)))
}
plots_4
## Save as pdf
plots_arrange <- marrangeGrob(grobs = plots_4, ncol = 2, nrow = 3) 
ggsave("final_nodiff_models_4.pdf", plots_arrange, width = 12, height = 18)

# 5 admix events
nodiff_models_5 <- all_models[c("W5_M1", "W5_M2", "W5_M3", "W5_M4", "W5_M5", "W5_M6", "W5_M7", "W5_M8", "W5_M9", "W5_M10", "W5_M11", "W5_M12",
                                "W5_M13", "W5_M14", "W5_M15", "W5_M16", "W5_M17", "W5_M19", "W5_M20", "W5_M21", "W5_M23", "W5_M26", "W5_M28", "W5_M31",
                                "W5_M33", "W5_M34", "W5_M35", "W5_M36", "W5_M37", "W5_M38", "W5_M42", "W5_M47", "W5_M48", "W5_M49")]
nodiff_models_5_stats <- model_stats %>% 
  filter(model %in% c("W5_M1", "W5_M2", "W5_M3", "W5_M4", "W5_M5", "W5_M6", "W5_M7", "W5_M8", "W5_M9", "W5_M10", "W5_M11", "W5_M12",
                      "W5_M13", "W5_M14", "W5_M15", "W5_M16", "W5_M17", "W5_M19", "W5_M20", "W5_M21", "W5_M23", "W5_M26", "W5_M28", "W5_M31",
                      "W5_M33", "W5_M34", "W5_M35", "W5_M36", "W5_M37", "W5_M38", "W5_M42", "W5_M47", "W5_M48", "W5_M49"))
plots_5 <- list()
for (i in 1:length(nodiff_models_5)){
  plots_5[[i]] <- plot_graph(nodiff_models_5[[i]]$edges[[1]],
                             title = paste(nodiff_models_5_stats$model[[i]],"\n",
                                           "Score:", round(nodiff_models_5_stats$orig_score[[i]], 2), 
                                           "Out-of-sample score:", round(nodiff_models_5_stats$oos_score[[i]], 2), 
                                           "Worst-residual:", round(nodiff_models_5_stats$worst_residual[[i]], 2)))
}
plots_5
## Save as pdf
plots_arrange <- marrangeGrob(grobs = plots_5, ncol = 2, nrow = 3) 
ggsave("final_nodiff_models_5.pdf", plots_arrange, width = 12, height = 18)

# 6 admix events
nodiff_models_6 <- all_models[c("W6_M1", "W6_M2", "W6_M3", "W6_M4", "W6_M5", "W6_M6", "W6_M7", "W6_M8", "W6_M9", "W6_M10", "W6_M11", "W6_M12",
                                "W6_M13", "W6_M14", "W6_M15", "W6_M16", "W6_M17", "W6_M18", "W6_M19", "W6_M20", "W6_M21", "W6_M22", "W6_M23", "W6_M24", 
                                "W6_M25", "W6_M26", "W6_M27", "W6_M28", "W6_M29", "W6_M30", "W6_M31", "W6_M32", "W6_M33", "W6_M34", "W6_M35", "W6_M36", 
                                "W6_M37", "W6_M38", "W6_M39", "W6_M40", "W6_M41", "W6_M42", "W6_M43", "W6_M44", "W6_M45", "W6_M46", "W6_M47", "W6_M49",
                                "W6_M50", "W6_M51", "W6_M52", "W6_M54", "W6_M55", "W6_M56", "W6_M59", "W6_M60", "W6_M61", "W6_M62", "W6_M63", "W6_M64",
                                "W6_M65", "W6_M68", "W6_M73", "W6_M79", "W6_M83", "W6_M84", "W6_M85", "W6_M92")]
nodiff_models_6_stats <- model_stats %>% 
  filter(model %in% c("W6_M1", "W6_M2", "W6_M3", "W6_M4", "W6_M5", "W6_M6", "W6_M7", "W6_M8", "W6_M9", "W6_M10", "W6_M11", "W6_M12",
                      "W6_M13", "W6_M14", "W6_M15", "W6_M16", "W6_M17", "W6_M18", "W6_M19", "W6_M20", "W6_M21", "W6_M22", "W6_M23", "W6_M24", 
                      "W6_M25", "W6_M26", "W6_M27", "W6_M28", "W6_M29", "W6_M30", "W6_M31", "W6_M32", "W6_M33", "W6_M34", "W6_M35", "W6_M36", 
                      "W6_M37", "W6_M38", "W6_M39", "W6_M40", "W6_M41", "W6_M42", "W6_M43", "W6_M44", "W6_M45", "W6_M46", "W6_M47", "W6_M49",
                      "W6_M50", "W6_M51", "W6_M52", "W6_M54", "W6_M55", "W6_M56", "W6_M59", "W6_M60", "W6_M61", "W6_M62", "W6_M63", "W6_M64",
                      "W6_M65", "W6_M68", "W6_M73", "W6_M79", "W6_M83", "W6_M84", "W6_M85", "W6_M92"))
plots_6 <- list()
for (i in 1:length(nodiff_models_6)){
  plots_6[[i]] <- plot_graph(nodiff_models_6[[i]]$edges[[1]],
                             title = paste(nodiff_models_6_stats$model[[i]],"\n",
                                           "Score:", round(nodiff_models_6_stats$orig_score[[i]], 2), 
                                           "Out-of-sample score:", round(nodiff_models_6_stats$oos_score[[i]], 2), 
                                           "Worst-residual:", round(nodiff_models_6_stats$worst_residual[[i]], 2)))
}
plots_6
## Save as pdf
plots_arrange <- marrangeGrob(grobs = plots_6, ncol = 2, nrow = 3) 
ggsave("final_nodiff_models_6.pdf", plots_arrange, width = 12, height = 18)


# Plot 3 admix events for Figure 3 ----
p_W3_M1 <- plot_graph(winner_3_ranked[[1]]$edges[[1]])
p_W3_M1
ggsave("p_W3_M1.svg", p_W3_M1, width = 90, height = 65, units = "mm")
p_W3_M2 <- plot_graph(winner_3_ranked[[2]]$edges[[1]])
p_W3_M2
ggsave("p_W3_M2.svg", p_W3_M2, width = 90, height = 65, units = "mm")
## These were futher edited for visual clarity in Inkscape.