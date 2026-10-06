# Compare het for 277, 277_3X, and 277_damaged across the whole genome

# Set working directory
setwd("./")

# Load packages
library(tidyverse)
library(ggpubr)

# Load data
het <- read_csv("SRR27955277_wholeGenomeComparisons.csv", col_names = TRUE)

# Test whether data are normally distributed (if yes, do t-test, if no, do wilcoxon test)
x <- filter(het, Sample == "SRR27955277") # Sample size = 167
shapiro.test(x$Het) # p-value = 0.2411; data are normally distributed
mean(x$Het) # 0.0005090513

x <- filter(het, Sample == "SRR27955277_3X") # Sample size = 155
shapiro.test(x$Het) # p-value = 0.03956 # data are not normally distributed
mean(x$Het) # 0.0002901853

x <- filter(het, Sample == "SRR27955277_damaged") # Sample size = 126
shapiro.test(x$Het) # p-value = 0.07442 # Data are normally distributed...only just
mean(x$Het) # 0.0005275009

## Thus, use Wilcoxon test, as not all data are normally distributed

# Plot with stats
myComparisons <- list(c("SRR27955277_damaged", "SRR27955277"), c("SRR27955277_damaged", "SRR27955277_3X"), c("SRR27955277", "SRR27955277_3X"))
ggboxplot(het, x = "Sample", y = "Het") +
  stat_compare_means(method = "wilcox.test", comparisons = myComparisons) +
  scale_x_discrete(labels = c("277" , "277_3X", "277_damaged"),
                   name = "Sample") +
  scale_y_continuous(name = "Genome-wide heterozygosity")
ggsave("SRR27955277_HetCompare_FigS7.png", dpi = 300)

# Get Wilcoxon test statistics
s277 <- filter(het, Sample == "SRR27955277")
s277_3X <- filter(het, Sample == "SRR27955277_3X")
s277_damaged <- filter(het, Sample == "SRR27955277_damaged")

## 277 vs 277_3X
wilcox.test(s277$Het, s277_3X$Het)
#Wilcoxon rank sum test with continuity correction
#
#data:  s277$Het and s277_3X$Het
#W = 25297, p-value < 2.2e-16
#alternative hypothesis: true location shift is not equal to 0

## 277 vs 277_damaged
wilcox.test(s277$Het, s277_damaged$Het)
#Wilcoxon rank sum test with continuity correction
#
#data:  s277$Het and s277_damaged$Het
#W = 8763, p-value = 0.01437
#alternative hypothesis: true location shift is not equal to 0

## 277_damaged vs 277_3X
wilcox.test(s277_damaged$Het, s277_3X$Het)
#Wilcoxon rank sum test with continuity correction

#data:  s277_damaged$Het and s277_3X$Het
#W = 19427, p-value < 2.2e-16
#alternative hypothesis: true location shift is not equal to 0
