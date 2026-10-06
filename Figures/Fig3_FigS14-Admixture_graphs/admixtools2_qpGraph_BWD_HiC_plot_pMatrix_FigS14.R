# Maunually re-plot pairwise matrices produced by admixtools2_qpGraph_BWD_HiC_signif.R

# Set working directory
setwd("C:/Users/pzx702/Documents/MSCAFellowship2021/Manuscripts/2022_Wildebeest_AncientEastAfrica/Results/admixture_graph/BWD_HiC/admixtools2/subspecies")

# Load packages
library(tidyverse)
library(ggpubr)

# Load data
W0 <- read_csv("W0_pemp_matrix_long.csv")
W1 <- read_csv("W1_pemp_matrix_long.csv")
W2 <- read_csv("W2_pemp_matrix_long.csv")
W3 <- read_csv("W3_pemp_matrix_long.csv")
W4 <- read_csv("W4_pemp_matrix_long.csv")
W5 <- read_csv("W5_pemp_matrix_long.csv")
W6 <- read_csv("W6_pemp_matrix_long.csv")

# Order models numerically for plotting
## Extract unique model names in correct numeric order
order_models <- function(x) {
  ordered_models <- x %>%
  mutate(model_num = as.numeric(str_extract(model1, "\\d+"))) %>%
  arrange(model_num) %>%
  pull(model1) %>%
  unique()

  # Apply to both axes
  x$model1 <- factor(x$model1, levels = ordered_models)
  x$model2 <- factor(x$model2, levels = ordered_models)
  
  return(x)
} 

W0 <- order_models(W0)
W1 <- order_models(W1)
W2 <- order_models(W2)
W3 <- order_models(W3)
W4 <- order_models(W4)
W5 <- order_models(W5)
W6 <- order_models(W6)
  
# Plot
## W0
p_W0 <- ggplot(W0, aes(x = model1, y = model2, fill = if_else(p_emp<0.05, "TRUE", "FALSE", missing = "NA"))) +
  geom_tile(color = "white", linewidth = 1) +
  #  scale_fill_viridis_c(option = "plasma", direction = -1, na.value = "grey90", name = "p_emp") +
  scale_fill_manual(values = c("TRUE" = "gold", "FALSE" = "purple3", "NA" = "grey90"),
                    name = "p-value",
                    breaks = c("TRUE", "FALSE", "NA"),
                    labels = c("<0.05", ">0.05", "NA")) +
  theme_minimal(base_size = 16) +
  coord_fixed() +
  labs(subtitle = "0 admix events") +
  theme(axis.text.x = element_text(angle = 90, vjust = 0.5, hjust = 1),
        axis.title = element_blank()) #+
#  geom_text(aes(label = round(p_emp, 2)), size = 2)
p_W0
ggsave("W0_pemp_matrix_plot_manual.pdf", p_W0, width = 180, height = 180, units = "mm")
ggsave("W0_pemp_matrix_plot_manual.png", p_W0, width = 180, height = 180, units = "mm", dpi = 300)

## W1
p_W1 <- ggplot(W1, aes(x = model1, y = model2, fill = if_else(p_emp<0.05, "TRUE", "FALSE", missing = "NA"))) +
  geom_tile(color = "white", linewidth = 1) +
  #  scale_fill_viridis_c(option = "plasma", direction = -1, na.value = "grey90", name = "p_emp") +
  scale_fill_manual(values = c("TRUE" = "gold", "FALSE" = "purple3", "NA" = "grey90"),
                    name = "p-value",
                    breaks = c("TRUE", "FALSE", "NA"),
                    labels = c("<0.05", ">0.05", "NA")) +
  theme_minimal(base_size = 16) +
  coord_fixed() +
  labs(subtitle = "1 admix event") +
  theme(axis.text.x = element_text(angle = 90, vjust = 0.5, hjust = 1),
        axis.title = element_blank()) #+
#  geom_text(aes(label = round(p_emp, 2)), size = 2)
p_W1
ggsave("W1_pemp_matrix_plot_manual.pdf", p_W1, width = 180, height = 180, units = "mm")
ggsave("W1_pemp_matrix_plot_manual.png", p_W1, width = 180, height = 180, units = "mm", dpi = 300)


## W2
p_W2 <- ggplot(W2, aes(x = model1, y = model2, fill = if_else(p_emp<0.05, "TRUE", "FALSE", missing = "NA"))) +
  geom_tile(color = "white", linewidth = 1) +
#  scale_fill_viridis_c(option = "plasma", direction = -1, na.value = "grey90", name = "p_emp") +
  scale_fill_manual(values = c("TRUE" = "gold", "FALSE" = "purple3", "NA" = "grey90"),
                    name = "p-value",
                    breaks = c("TRUE", "FALSE", "NA"),
                    labels = c("<0.05", ">0.05", "NA")) +
  theme_minimal(base_size = 16) +
  coord_fixed() +
  labs(subtitle = "2 admix events") +
  theme(axis.text.x = element_text(angle = 90, vjust = 0.5, hjust = 1),
        axis.title = element_blank()) #+
#  geom_text(aes(label = round(p_emp, 2)), size = 2)
p_W2
ggsave("W2_pemp_matrix_plot_manual.pdf", p_W2, width = 180, height = 180, units = "mm")
ggsave("W2_pemp_matrix_plot_manual.png", p_W2, width = 180, height = 180, units = "mm", dpi = 300)

## W3
p_W3 <- ggplot(W3, aes(x = model1, y = model2, fill = if_else(p_emp<0.05, "TRUE", "FALSE", missing = "NA"))) +
  geom_tile(color = "white", linewidth = 0.5) +
  #  scale_fill_viridis_c(option = "plasma", direction = -1, na.value = "grey90", name = "p_emp") +
  scale_fill_manual(values = c("TRUE" = "gold", "FALSE" = "purple3", "NA" = "grey90"),
                    name = "p-value",
                    breaks = c("TRUE", "FALSE", "NA"),
                    labels = c("<0.05", ">0.05", "NA")) +
  theme_minimal(base_size = 16) +
  coord_fixed() +
  labs(subtitle = "3 admix events") +
  theme(axis.text.x = element_text(angle = 90, vjust = 0.5, hjust = 1),
        axis.title = element_blank(),
        axis.text = element_text(size = 9)) #+
#  geom_text(aes(label = round(p_emp, 2)), size = 2)
p_W3
ggsave("W3_pemp_matrix_plot_manual.pdf", p_W3, width = 180, height = 180, units = "mm")
ggsave("W3_pemp_matrix_plot_manual.png", p_W3, width = 180, height = 180, units = "mm", dpi = 300)

## W4
p_W4 <- ggplot(W4, aes(x = model1, y = model2, fill = if_else(p_emp<0.05, "TRUE", "FALSE", missing = "NA"))) +
  geom_tile(color = "white", linewidth = 0.25) +
  #  scale_fill_viridis_c(option = "plasma", direction = -1, na.value = "grey90", name = "p_emp") +
  scale_fill_manual(values = c("TRUE" = "gold", "FALSE" = "purple3", "NA" = "grey90"),
                    name = "p-value",
                    breaks = c("TRUE", "FALSE", "NA"),
                    labels = c("<0.05", ">0.05", "NA")) +
  theme_minimal(base_size = 16) +
  coord_fixed() +
  labs(subtitle = "4 admix events") +
  theme(axis.text.x = element_text(angle = 90, vjust = 0.5, hjust = 1),
        axis.title = element_blank(),
        axis.text = element_text(size = 6)) #+
#  geom_text(aes(label = round(p_emp, 2)), size = 2)
p_W4
ggsave("W4_pemp_matrix_plot_manual.pdf", p_W4, width = 190, height = 190, units = "mm")
ggsave("W4_pemp_matrix_plot_manual.png", p_W4, width = 190, height = 190, units = "mm", dpi = 300)

## W5
p_W5 <- ggplot(W5, aes(x = model1, y = model2, fill = if_else(p_emp<0.05, "TRUE", "FALSE", missing = "NA"))) +
  geom_tile(color = "white", linewidth = 0.25) +
#  scale_fill_viridis_c(option = "plasma", direction = -1, na.value = "grey90", name = "p_emp") +
  scale_fill_manual(values = c("TRUE" = "gold", "FALSE" = "purple3", "NA" = "grey90"),
                    name = "p-value",
                    breaks = c("TRUE", "FALSE", "NA"),
                    labels = c("<0.05", ">0.05", "NA")) +
  theme_minimal(base_size = 16) +
  coord_fixed() +
  labs(subtitle = "5 admix events") +
  theme(axis.text.x = element_text(angle = 90, vjust = 0.5, hjust = 1),
        axis.title = element_blank(),
        axis.text = element_text(size = 5)) #+
#  geom_text(aes(label = round(p_emp, 2)), size = 2)
p_W5
ggsave("W5_pemp_matrix_plot_manual.pdf", p_W5, width = 190, height = 190, units = "mm")
ggsave("W5_pemp_matrix_plot_manual.png", p_W5, width = 190, height = 190, units = "mm", dpi = 300)

## W6
p_W6 <- ggplot(W6, aes(x = model1, y = model2, fill = if_else(p_emp<0.05, "TRUE", "FALSE", missing = "NA"))) +
  geom_tile(color = "white", linewidth = 0.25) +
  #  scale_fill_viridis_c(option = "plasma", direction = -1, na.value = "grey90", name = "p_emp") +
  scale_fill_manual(values = c("TRUE" = "gold", "FALSE" = "purple3", "NA" = "grey90"),
                    name = "p-value",
                    breaks = c("TRUE", "FALSE", "NA"),
                    labels = c("<0.05", ">0.05", "NA")) +
  theme_minimal(base_size = 16) +
  coord_fixed() +
  labs(subtitle = "6 admix events") +
  theme(axis.text.x = element_text(angle = 90, vjust = 0.5, hjust = 1),
        axis.title = element_blank(),
        axis.text = element_text(size = 5)) #+
#  geom_text(aes(label = round(p_emp, 2)), size = 2)
p_W6
ggsave("W6_pemp_matrix_plot_manual.pdf", p_W6, width = 190, height = 190, units = "mm")
ggsave("W6_pemp_matrix_plot_manual.png", p_W6, width = 190, height = 190, units = "mm", dpi = 300)
common_legend <- get_legend(p_W6)

## Combine
plots <- list(p_W0, p_W1, p_W2, p_W3, p_W4, p_W5, p_W6)
p_W <- ggarrange(plotlist = plots, ncol = 2, nrow = 4, legend.grob = common_legend, legend.position = "top")
ggsave("pemp_plots_W0toW6.pdf", width = 205, height = 290, units = "mm")
