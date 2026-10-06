# Calculate f4 stats in admixtools2
## See: https://uqrmaie1.github.io/admixtools/reference/qpdstat.html
## And: https://uqrmaie1.github.io/admixtools/articles/fstats.html
## And: https://compvar-workshop.readthedocs.io/en/latest/contents/03_f3stats/f3stats.html#f4-statistics

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
library(scales)
library(openxlsx)
#install.packages("openxlsx")
#install.packages("shiny", dependencies = TRUE)
#library(shiny)

# Set high number of digits for printing precision, as some likelihood scores appear identical for the first few decimals, but are actually different.
options(digits = 4)

# 1. Calculate/load f2 statistics & load results into R ----
prefix = "n82_fossil_BWD_HiC_autosomes.packedancestrymap"
f2_dir = "f2_results"
#extract_f2(prefix, f2_dir,
#           adjust_pseudohaploid = TRUE,
#           auto_only = FALSE,
#           overwrite = TRUE)
## Load results
f2_blocks = f2_from_precomp(f2_dir)

# 2. Test for effects of potential baises in the data ----
# See: https://uqrmaie1.github.io/admixtools/articles/fstats.html#biases

## 2.1. Use the same SNPs for every f4-statistic (maxmiss = 0, the default) - the most conservative option and the one used for the admixture graph analysis
### Calculate f4 stats from precomputed f2 stats
### Using f2
f4_Hart_Black_rest_Fossil936 <- qpdstat(f2_blocks,
        pop1 = "Hartebeest",
        pop2 = "Black",
        pop3 = c("Wwb", "Ewb", "Cookson", "Nyassa", "Brindled"),
        pop4 = "Fossil936",
        afprod = FALSE)
write_csv(f4_Hart_Black_rest_Fossil936, "f4_Hart_Black_rest_Fossil936.csv")
#### Plot
# Order factors
f4_Hart_Black_rest_Fossil936$pop3 <- factor(f4_Hart_Black_rest_Fossil936$pop3, 
                                            levels = c("Brindled", "Cookson", "Nyassa", "Wwb", "Ewb"))
# Define palette
myPalette = c("#8585c4", "#ff96ca", "#cc95ff", "#cace3d", "#8dc5a6")
# Plot
ggplot(f4_Hart_Black_rest_Fossil936) +
  geom_vline(xintercept = 0, linetype = "dashed") +
  geom_point(aes(x = est, y = pop3, colour = pop3), size = 5) +
  geom_linerange(aes(x = est, y = pop3, xmin = est-se, xmax = est+se)) +
  geom_label(aes(x = -0.002, y = pop3, label = paste0("Z = ", round(z, digits = 2))), fontface = "bold", size = 4, label.size = 0) +
  scale_colour_manual(values = myPalette,
                      breaks = c("Brindled", "Cookson", "Nyassa", "Wwb", "Ewb"),
                      name = "H3") + # Define colours 
  guides(colour = "none") + # Remove legend
  theme_bw(base_size = 18) +
#  scale_x_continuous(limits = c(-0.001, 0.0075)) +
  #scale_x_continuous(limits = c(-0.05, 0.05)) +
  labs(x = "f4(Hartebeest, Black; H3, Fossil936)", y = "H3")
ggsave("f4_Hart_Black_rest_Fossil936.png", dpi = 300, width = 6.74, height = 4.89, units = "in")

## 2.2. Use different SNPs for each f4-statistic, using genotype files directly (allsnps = TRUE)
f4_Hart_Black_rest_Fossil936_geno <- qpdstat(prefix,
        pop1 = "Hartebeest",
        pop2 = "Black",
        pop3 = c("Wwb", "Ewb", "Cookson", "Nyassa", "Brindled"),
        pop4 = "Fossil936",
        allsnps = TRUE)
## Slight change: Results in the D-stat between Black and Wwb to be non-significant (Z = -2.65), whereas in 2.1. it was significant (Z = -3.19)
write_csv(f4_Hart_Black_rest_Fossil936_geno, "f4_Hart_Black_rest_Fossil936_geno.csv")
#### Plot
# Order factors
f4_Hart_Black_rest_Fossil936_geno$pop3 <- factor(f4_Hart_Black_rest_Fossil936_geno$pop3, 
                                                 levels = c("Brindled", "Cookson", "Nyassa", "Wwb", "Ewb"))
# Plot
ggplot(f4_Hart_Black_rest_Fossil936_geno) +
  geom_vline(xintercept = 0, linetype = "dashed") +
  geom_point(aes(x = est, y = pop3, colour = pop3), size = 5) +
  geom_linerange(aes(x = est, y = pop3, xmin = est-se, xmax = est+se)) +
  geom_label(aes(x = -0.002, y = pop3, label = paste0("Z = ", round(z, digits = 2))), fontface = "bold", size = 4, label.size = 0) +
  scale_colour_manual(values = myPalette,
                      breaks = c("Brindled", "Cookson", "Nyassa", "Wwb", "Ewb"),
                      name = "H3") + # Define colours 
  guides(colour = "none") + # Remove legend
  theme_bw(base_size = 18) +
  #  scale_x_continuous(limits = c(-0.001, 0.0075)) +
  #scale_x_continuous(limits = c(-0.05, 0.05)) +
  labs(x = "f4(Hartebeest, Black; H3, Fossil936)", y = "H3")
ggsave("f4_Hart_Black_rest_Fossil936_geno.png", dpi = 300, width = 6.74, height = 4.89, units = "in")

## 2.3. Use different SNPs for each f2-statistic (maxmiss > 0)
# Requires recalculating f2 statistics
prefix = "n82_fossil_BWD_HiC_autosomes.packedancestrymap"
f2_dir_miss = "f2_results_miss"
#extract_f2(prefix, f2_dir_miss,
#           adjust_pseudohaploid = TRUE,
#           auto_only = FALSE,
#           overwrite = TRUE,
#           maxmiss = 0.2) # Max missingness of 20% - just randomly chosen
## Load results
f2_blocks_miss = f2_from_precomp(f2_dir_miss)
afprod_blocks = f2_from_precomp(f2_dir_miss, afprod = TRUE) # Use afprod in cases where maxmiss > 0.

# do f4
f4_Hart_Black_rest_Fossil936_miss0.2 <- qpdstat(f2_blocks_miss,
                                        pop1 = "Hartebeest",
                                        pop2 = "Black",
                                        pop3 = c("Wwb", "Ewb", "Cookson", "Nyassa", "Brindled"),
                                        pop4 = "Fossil936",
                                        afprod = FALSE)
write_csv(f4_Hart_Black_rest_Fossil936_miss0.2, "f4_Hart_Black_rest_Fossil936_miss0.2.csv")
#### Plot
# Order factors
f4_Hart_Black_rest_Fossil936_miss0.2$pop3 <- factor(f4_Hart_Black_rest_Fossil936_miss0.2$pop3, 
                                            levels = c("Brindled", "Cookson", "Nyassa", "Wwb", "Ewb"))
# Plot
ggplot(f4_Hart_Black_rest_Fossil936_miss0.2) +
  geom_vline(xintercept = 0, linetype = "dashed") +
  geom_point(aes(x = est, y = pop3, colour = pop3), size = 5) +
  geom_linerange(aes(x = est, y = pop3, xmin = est-se, xmax = est+se)) +
  geom_label(aes(x = -0.002, y = pop3, label = paste0("Z = ", round(z, digits = 2))), fontface = "bold", size = 4, label.size = 0) +
  scale_colour_manual(values = myPalette,
                      breaks = c("Brindled", "Cookson", "Nyassa", "Wwb", "Ewb"),
                      name = "H3") + # Define colours 
  guides(colour = "none") + # Remove legend
  theme_bw(base_size = 18) +
  #  scale_x_continuous(limits = c(-0.001, 0.0075)) +
  #scale_x_continuous(limits = c(-0.05, 0.05)) +
  labs(x = "f4(Hartebeest, Black; H3, Fossil936)", y = "H3")
ggsave("f4_Hart_Black_rest_Fossil936_miss0.2.png", dpi = 300, width = 6.74, height = 4.89, units = "in")

### 2.3.1. Using allele frequency products when missing data were allowed (afprod = TRUE)
f4_Hart_Black_rest_Fossil936_afprod <- qpdstat(afprod_blocks,
                                               pop1 = "Hartebeest",
                                               pop2 = "Black",
                                               pop3 = c("Wwb", "Ewb", "Cookson", "Nyassa", "Brindled"),
                                               pop4 = "Fossil936",
                                               afprod = TRUE)
write_csv(f4_Hart_Black_rest_Fossil936_afprod, "f4_Hart_Black_rest_Fossil936_afprod.csv")
#### Plot
# Order factors
f4_Hart_Black_rest_Fossil936_afprod$pop3 <- factor(f4_Hart_Black_rest_Fossil936_afprod$pop3, 
                                                   levels = c("Brindled", "Cookson", "Nyassa", "Wwb", "Ewb"))
# Plot
ggplot(f4_Hart_Black_rest_Fossil936_afprod) +
  geom_vline(xintercept = 0, linetype = "dashed") +
  geom_point(aes(x = est, y = pop3, colour = pop3), size = 5) +
  geom_linerange(aes(x = est, y = pop3, xmin = est-se, xmax = est+se)) +
  geom_label(aes(x = -0.002, y = pop3, label = paste0("Z = ", round(z, digits = 2))), fontface = "bold", size = 4, label.size = 0) +
  scale_colour_manual(values = myPalette,
                      breaks = c("Brindled", "Cookson", "Nyassa", "Wwb", "Ewb"),
                      name = "H3") + # Define colours 
  guides(colour = "none") + # Remove legend
  theme_bw(base_size = 18) +
  #  scale_x_continuous(limits = c(-0.001, 0.0075)) +
  #scale_x_continuous(limits = c(-0.05, 0.05)) +
  labs(x = "f4(Hartebeest, Black; H3, Fossil936)", y = "H3")
ggsave("f4_Hart_Black_rest_Fossil936_afprod.png", dpi = 300, width = 6.74, height = 4.89, units = "in")
### Both 3.2. and 3.2.1. have slightly significant false positives, indicating missing data is having an effect and it is best to use the maxmiss = 0 option.

# 3. f4 calculations to test gene flow scenarios (based on maxmiss = 0 data) ----

## 3.1. Testing allele sharing between Fossil936 and Wwb (Western white-bearded)
f4_Hart_Fossil936_rest_Wwb <- qpdstat(f2_blocks,
        pop1 = "Hartebeest",
        pop2 = "Fossil936",
        pop3 = c("Ewb", "Cookson", "Nyassa", "Brindled", "Black"),
        pop4 = "Wwb")
write_csv(f4_Hart_Fossil936_rest_Wwb, "f4_Hart_Fossil936_rest_Wwb.csv")
### Plot
myPalette2 = c("#8dc5a6", "#ff96ca", "#cc95ff", "#8585c4", "#5c5c5c")
# Order factors
f4_Hart_Fossil936_rest_Wwb$pop3 <- factor(f4_Hart_Fossil936_rest_Wwb$pop3, 
                                                   levels = c("Ewb", "Cookson", "Nyassa", "Brindled", "Black"))
ggplot(f4_Hart_Fossil936_rest_Wwb) +
  geom_vline(xintercept = 0, linetype = "dashed") +
  geom_point(aes(x = est, y = pop3, colour = pop3), size = 5) +
  geom_linerange(aes(x = est, y = pop3, xmin = est-se, xmax = est+se)) +
  geom_label(aes(x = 0.04, y = pop3, label = paste0("Z = ", round(z, digits = 2))), fontface = "bold", size = 4, label.size = 0) +
  scale_colour_manual(values = myPalette2,
                      breaks = c("Ewb", "Cookson", "Nyassa", "Brindled", "Black"),
                      name = "H3") + # Define colours 
  guides(colour = "none") + # Remove legend
  theme_bw(base_size = 18) +
  #  scale_x_continuous(limits = c(-0.001, 0.0075)) +
  #scale_x_continuous(limits = c(-0.05, 0.05)) +
  labs(x = "f4(Hartebeest, Fossil936; H3, Wwb)", y = "H3")
ggsave("f4_Hart_Fossil936_rest_Wwb.png", dpi = 300, width = 6.74, height = 4.89, units = "in")
### RESULT: All values positive, indicating shared alleles (gene flow) between Fossil936 and Wwb, regardless of which taxon is at H3

## 3.2. Testing to see what happens if we switch Wwb and Ewb
f4_Hart_Fossil936_rest_Ewb <- qpdstat(f2_blocks,
                                      pop1 = "Hartebeest",
                                      pop2 = "Fossil936",
                                      pop3 = c("Wwb", "Cookson", "Nyassa", "Brindled", "Black"),
                                      pop4 = "Ewb")
write_csv(f4_Hart_Fossil936_rest_Ewb, "f4_Hart_Fossil936_rest_Ewb.csv")
### Plot
f4_Hart_Fossil936_rest_Ewb <- read_csv("f4_Hart_Fossil936_rest_Ewb.csv", col_names = T)
myPalette3 = c("#cace3d", "#ff96ca", "#cc95ff", "#8585c4", "#5c5c5c")
# Order factors
f4_Hart_Fossil936_rest_Ewb$pop3 <- factor(f4_Hart_Fossil936_rest_Ewb$pop3, 
                                          levels = c("Wwb", "Cookson", "Nyassa", "Brindled", "Black"))
ggplot(f4_Hart_Fossil936_rest_Ewb) +
  geom_vline(xintercept = 0, linetype = "dashed") +
  geom_point(aes(x = est, y = pop3, colour = pop3), size = 10) +
  geom_linerange(aes(x = est, y = pop3, xmin = est-se, xmax = est+se)) +
  geom_label(aes(x = 0.04, y = pop3, label = paste0("Z = ", round(z, digits = 2))), fontface = "bold", size = 4, label.size = 0) +
  scale_colour_manual(values = myPalette3,
                      breaks = c("Wwb", "Cookson", "Nyassa", "Brindled", "Black"),
                      name = "H3") + # Define colours 
  guides(colour = "none") + # Remove legend
  theme_bw(base_size = 18) +
  #  scale_x_continuous(limits = c(-0.001, 0.0075)) +
  #scale_x_continuous(limits = c(-0.05, 0.05)) +
  labs(x = "f4 (Hartebeest, Fossil936; H3, Eastern white-bearded)", y = "H3")
ggsave("f4_Hart_Fossil936_rest_Ewb.png", dpi = 300, width = 6.74, height = 4.89, units = "in")
ggsave("f4_Hart_Fossil936_rest_Ewb.svg", width = 6.74, height = 4.89, units = "in")
ggsave("f4_Hart_Fossil936_rest_Ewb_flat.png", dpi = 300, width = 12, height = 4.89, units = "in")
ggsave("f4_Hart_Fossil936_rest_Ewb_flat.svg", width = 12, height = 4.89, units = "in")


### RESULT: f4 is now negative between Fossil936 and Wwb, showing that indeed there is more allele sharing between Fossil936 and Wwb than between Fossil936 & Ewb

## 3.3. Testing allele sharing between Wwb and Fossil936 and Cookson (i.e. is Wwb a ~50/50 lineage between the two? Yes: f4 = 0, No: f4 != 0)
f4_Hart_Wwb_Fossil936_Cookson <- qpdstat(f2_blocks,
                                      pop1 = "Hartebeest",
                                      pop2 = "Wwb",
                                      pop3 = "Fossil936",
                                      pop4 = "Cookson")
write_csv(f4_Hart_Wwb_Fossil936_Cookson, "f4_Hart_Wwb_Fossil936_Cookson.csv")
# f4 = -0.001356037, z = -7.490136. Therefore, Wwb shares an excess of alleles with Fossil936, as compared to Cookson, and thus this result does not support the 50/50 hypothesis

# 4. f3 statistics to test the two different scenarios from qpGraph (exemplified by W3_M1 and W3_M2) ----
## 4.1. Testing W3_M1 scenario: Wwb is admixed ~50/50 between ancestor of Cookson/Nyassa and ancestor of Fossil936
qp3pop(f2_blocks,
   pop1 = "Brindled",
   pop2 = "Black",
   pop3 = c("Cookson", "Nyassa", "Ewb", "Wwb", "Fossil936"))

## 4.2. Testing W3_M2 scenario: Fossil936 is admixed ~70/30 between Wwb ancestor of all blue wildebeest
f3(f2_blocks,
   pop1 = "Fossil936",
   pop2 = "Wwb",
   pop3 = "Brindled")

# 5. qpAdm to test different scenarios from qpGraph (exemplified by W3_M1 and W3_M2) ----
## 5.1. Testing W3_M1 scenario: Wwb is admixed ~50/50 between ancestor of Cookson/Nyassa and ancestor of Fossil936 - exclude Ewb, as it is in a clade with Wwb.
### It is also not appropriate to include Nyassa in this scenario as a right (reference) pop, since it shares a more common recent ancestor with Cookson (a source/left pop) 
### than the proposed split/gene flow event of the target pop (Wwb). This violates the assumptions of qpAdm and the model will be identified as implausible (scen1 below).
### Therefore, have to include hartebeest as a reference pop in order to have more ref pops than source pops in the analysis (scen1.1 below). 
### See: Harney et al. (2021) supplementary info pg. 32 (qpAdm User Guide)
scen1 <-  qpadm(f2_blocks,
      left = c("Fossil936", "Cookson"),
      right = c("Brindled", "Nyassa", "Black"),
      target = "Wwb")
scen1
# p-value in scen1$rankdrop for the full mode (top row) is 7.42e-249 (i.e. <0.05),
# so null hypothesis is rejected (i.e. Wwb is not an admixed pop between Fossil936 and Cookson). However, this violates the assumptions of qpAdm! Therefore, use hartebeest instead of Nyassa
scen1.1 <-  qpadm(f2_blocks,
                left = c("Fossil936", "Cookson"),
                right = c("Brindled", "Black", "Hartebeest"),
                target = "Wwb")
scen1.1
# p-value in scen1$rankdrop for the full mode (top row) is 0.71 (i.e. >0.05),
# so null hypothesis is not rejected (i.e. Wwb can be an admixed pop between Fossil936 and Cookson, with admix props of 47.9% from Fossil936 and 52.1% from Cookson).
# This is the correct model to use for scen1 (scen1.1) and is used for the qpadm_rotate analysis below.
# Note, when Ewb in included, the admix proportions (scen1$weights) are also outside the biologically relevant range of 0-1 (see: https://doi.org/10.1093/genetics/iyaa045), potentially because this violates the assumptions of qpAdm.
# See: https://github.com/uqrmaie1/admixtools/issues/78#issuecomment-2608287154
scen1.miss <- qpadm(prefix,
                    left = c("Fossil936", "Cookson"),
                    right = c("Brindled", "Black", "Hartebeest"),
                    target = "Wwb",
                    allsnps = TRUE)
scen1.miss
#i "allsnps = TRUE" uses different SNPs for each f4-statistic
#Number of SNPs used for each f4-statistic:
#  pop1      pop2     pop3       pop4        n
#1  Wwb   Cookson Brindled      Black 10085930
#2  Wwb   Cookson Brindled Hartebeest  9902274
#3  Wwb Fossil936 Brindled      Black 10071752
#4  Wwb Fossil936 Brindled Hartebeest  9888961
# Null hypothesis also not rejected when using allsnps = TRUE (to include more data (SNPs)), with p = 0.177 and admix props of 48.1% (Fossil936) and 51.9% (Cookson).

## 5.2. Testing W3_M2 scenario: Fossil936 is admixed ~70/30 between Wwb ancestor of all blue wildebeest  - exclude Ewb, as it is in a clade with Wwb
scen2 <-  qpadm(f2_blocks,
                left = c("Black", "Wwb"),
                right = c("Brindled", "Nyassa", "Cookson"),
                target = "Fossil936")
scen2
# p-value in scen1$rankdrop for the full mode (top row) is 0.723 (i.e. >0.05), so null hypothesis is not rejected (i.e. Fossil936 being an admixed pop between Wwb and an ancestral lineage (represented by black wildebeest) cannot be rejected)
# See: https://github.com/uqrmaie1/admixtools/issues/78#issuecomment-2608287154
# Additionally, the admixture proportions are: Black (ancestral) = 0.207, and Wwb = 0.793.
# Note: When Ewb is included, the result remains consistent, with only a slight change in the p-value (0.540) and admix proportions (0.202, 0.798, for Black and Wwb, respectively)
# Consequently, W3_M2 is preferred over W3_M1
scen2.miss <- qpadm(prefix,
                    left = c("Black", "Wwb"),
                    right = c("Nyassa", "Cookson", "Brindled"),
                    target = "Fossil936",
                    allsnps = TRUE)
scen2.miss
#i "allsnps = TRUE" uses different SNPs for each f4-statistic
#Number of SNPs used for each f4-statistic:
#  pop1  pop2   pop3     pop4        n
#1 Fossil936 Black Nyassa Brindled 10071053
#2 Fossil936 Black Nyassa  Cookson 10069616
#3 Fossil936   Wwb Nyassa Brindled 10075091
#4 Fossil936   Wwb Nyassa  Cookson 10073645
# Same result as above, with p-value = 0.651 and admix props of 0.212 (black) and 0.788 (Wwb)

### 5.3. Sanity check - Brindled is admixed between black (7%) and ancestor of Cookson (93%), but excluding Fossil936 (so only modern pops)
scen3 <- qpadm(f2_blocks,
               left = c("Cookson", "Black"),
               right = c("Nyassa", "Ewb", "Wwb"),
               target = "Brindled")
scen3
# Null hypothesis is not rejected: i.e. Brindled is admixed (p = 0.416), but proportions are a bit off compared to the qpGraph models: 0.325 (Cookson) and 0.675 (Black),
# but this scenario is actually difficult to model, because it feels like we don't have ideal reference (right) populations that are outside blue and black wildebeest, so 
# all the reference populations share some drift with Brindled and Cookson (I might be way off here though), and because gene flow from Black into Brindled liekly occurred after
# split between Brindled and Cookson, so admixture proportions from Black into Brindled are elevated due to this model violation (see Harney et al. (2021) pg. 11-12)
# Harney at al. (2021): https://doi.org/10.1093/genetics/iyaa045
scen3.miss <- qpadm(prefix,
               left = c("Cookson", "Black"),
               right = c("Nyassa", "Ewb", "Wwb"),
               target = "Brindled",
               allsnps = TRUE)
scen3.miss
#i "allsnps = TRUE" uses different SNPs for each f4-statistic
#Number of SNPs used for each f4-statistic:
#  pop1    pop2   pop3 pop4        n
#1 Brindled   Black Nyassa  Ewb 10086854
#2 Brindled   Black Nyassa  Wwb 10086854
#3 Brindled Cookson Nyassa  Ewb 10089407
#4 Brindled Cookson Nyassa  Wwb 10089407
# Null hypothesis is not rejected (so Brindled modelled as admix between Cookson and Black is not rejected (p = 0.0908), with admix props of 0.326 for Cookson and 0.674 for Black)

### 5.4. Sanity check - Cookson is admixed between ancestral Brindled (8%) and ancestral Nyassa (92%)
scen4 <- qpadm(f2_blocks,
               left = c("Brindled", "Nyassa"),
               right = c("Black", "Wwb", "Ewb"),
               target = "Cookson")
scen4
# Null hypothesis is rejected (p = 2.65e-62), but admix proportions are more in line with what is expected: 0.121 (Brindled) and 0.879 (Nyassa)
scen4.miss <- qpadm(prefix,
               left = c("Brindled", "Nyassa"),
               right = c("Black", "Wwb", "Ewb"),
               target = "Cookson",
               allsnps = TRUE)
scen4.miss
# Same result as scen4: p = 1.44e-53, admix props: 0.117 (Brindled), 0.883 (Nyassa)

# 5.5. Use qpadm_rotate() to test many combinations of populations ----
## 5.5.1. Testing all models for Wwb as an admixed pop
scen1.all <- qpadm_rotate(f2_blocks, 
             leftright = c("Fossil936", "Cookson", "Brindled", "Black"), 
             target = "Wwb", 
             rightfix = "Hartebeest")
## Save as excel file (need to convert first two columns to character strings)
scen1.all %>%
  mutate(across(c("left", "right"), as.character)) %>%
  writexl::write_xlsx("qpAdm_scen1.all.xlsx", format_headers = FALSE)
  
scen1.all.full <- qpadm_rotate(f2_blocks, 
                          leftright = c("Fossil936", "Cookson", "Brindled", "Black"), 
                          target = "Wwb", 
                          rightfix = "Hartebeest",
                          full_results = TRUE)
# One plausible and non-significant p-value (0.719) model: Fossil936 and Cookson ancestors contribute ~48% and 52% to Wwb. Null hypothesis is not rejected. Same as scen1.1 above

## 5.5.2. Testing all models for Fossil936 as an admixed pop
scen2.all <- qpadm_rotate(f2_blocks, 
                          leftright = c("Wwb", "Nyassa", "Brindled", "Black"), 
                          target = "Fossil936", 
                          rightfix = "Cookson")
## Save as excel file (need to convert first two columns to character strings)
scen2.all %>%
  mutate(across(c("left", "right"), as.character)) %>%
  writexl::write_xlsx("qpAdm_scen2.all.xlsx", format_headers = FALSE)

scen2.all.full <- qpadm_rotate(f2_blocks, 
                               leftright = c("Wwb", "Nyassa", "Brindled", "Black"), 
                               target = "Fossil936", 
                               rightfix = "Cookson",
                               full_results = TRUE)
# Only one plausible and non-significant result: Fossil936 is admixed between ancestral Wwb and ancestral Blue (represented by Black), same as single result above (5.2)
# p-value = 0.7233027, admix props: 0.7932316 (Wwb) and 0.2067684 
# Notably: Wwb and Brindled as sources is implausible (admix props outside of 0-1), which perhaps can be interpreted as aligning with W3_M2 showing gene flow from and ancestral blue,
# pop that predates the split of all blue subspecies (therefore the qpAdm model is only supported when Black is the source, because it is the only source pop outside Blue, and thus is represents ancestral Blue in this scenario)

## 5.5.3. Testing all models for Brindled as an admixed pop
scen3.all <- qpadm_rotate(f2_blocks, 
                          leftright = c("Cookson", "Nyassa", "Ewb", "Black"), 
                          target = "Brindled", 
                          rightfix = "Wwb")
## Save as excel file (need to convert first two columns to character strings)
scen3.all %>%
  mutate(across(c("left", "right"), as.character)) %>%
  writexl::write_xlsx("qpAdm_scen3.all.xlsx", format_headers = FALSE)

scen3.all.full <- qpadm_rotate(f2_blocks, 
                          leftright = c("Cookson", "Nyassa", "Ewb", "Black"), 
                          target = "Brindled", 
                          rightfix = "Wwb",
                          full_results = TRUE)
# Only one plausible and non-significant result: Cookson-Black as sources (p = 0.4256971), with admix props of 0.3247033 (Cookson) and 0.6752967 (Black)

## 5.5.4. Testing all models for Cookson as an admixed pop
scen4.all <- qpadm_rotate(f2_blocks, 
                          leftright = c("Brindled", "Nyassa", "Wwb", "Ewb"), 
                          target = "Cookson", 
                          rightfix = "Black")
## Save as excel file (need to convert first two columns to character strings)
scen4.all %>%
  mutate(across(c("left", "right"), as.character)) %>%
  writexl::write_xlsx("qpAdm_scen4.all.xlsx", format_headers = FALSE)

scen4.all.full <- qpadm_rotate(f2_blocks, 
                               leftright = c("Brindled", "Nyassa", "Wwb", "Ewb"), 
                               target = "Cookson", 
                               rightfix = "Black",
                               full_results = TRUE)
# Same as scen4 above: Brindled-Nyassa is plausible but rejected (p-value = 2.645632e-62), admix props of 0.1209549 (Brindled) and 0.8790451 (Nyassa)

# Save as scen*.all.full as excel files - one sheet per model
## Define function
save_nested_scenarios_excel <- function(nested_tbl, file) {
  # Validate input
  if (!all(c("weights", "f4", "rankdrop", "popdrop") %in% names(nested_tbl))) {
    stop("Input tibble must contain list-columns: weights, f4, rankdrop, and popdrop.")
  }
  
  wb <- createWorkbook()
  
  for (i in seq_len(nrow(nested_tbl))) {
    # Create safe sheet name (max 31 chars)
    sheet_name <- paste0("model_", i)
    addWorksheet(wb, sheet_name)
    
    row_pos <- 1
    
    # Helper function to write one section
    write_section <- function(title, data) {
      nonlocal_row_pos <- row_pos  # capture for this section
      writeData(wb, sheet_name, paste0("Table: ", title), startRow = nonlocal_row_pos)
      nonlocal_row_pos <- nonlocal_row_pos + 1
      writeData(wb, sheet_name, data, startRow = nonlocal_row_pos)
      return(nonlocal_row_pos + nrow(data) + 2)
    }
    
    # Write all 4 tables
    row_pos <- write_section("Weights",  nested_tbl$weights[[i]])
    row_pos <- write_section("f4",       nested_tbl$f4[[i]])
    row_pos <- write_section("Rankdrop", nested_tbl$rankdrop[[i]])
    row_pos <- write_section("Popdrop",  nested_tbl$popdrop[[i]])
    
    # Auto-adjust columns
    setColWidths(wb, sheet = sheet_name, cols = 1:20, widths = "auto")
  }
  
  # Save file
  saveWorkbook(wb, file, overwrite = TRUE)
  message("File saved to: ", normalizePath(file))
}
## Save scens
save_nested_scenarios_excel(scen1.all.full, "qpAdm_scen1.all.full.xlsx")
save_nested_scenarios_excel(scen2.all.full, "qpAdm_scen2.all.full.xlsx")
save_nested_scenarios_excel(scen3.all.full, "qpAdm_scen3.all.full.xlsx")
save_nested_scenarios_excel(scen4.all.full, "qpAdm_scen4.all.full.xlsx")
