# Check convergence of MCMC runs with convenience R package
## Mainly used to check that mulitple independent runs have converged to the same posterior using the KS statistic
## See: https://revbayes.github.io/tutorials/convergence/

## Load packages
library(convenience)

## Set working directory
#setwd("")

## Check convergence
### Note: Only using log files (excl. tree files), as convenience has trouble with the coupledMCMC (MC3) tree files
checkRunI_MC3 <- checkConvergence(format = "beast",
                                  list_files = c("RunI1_MC3_convenience.log", "RunI2_MC3_convenience.log", "RunI3_MC3_convenience.log"))
checkRunI_MC3

## Plot Kolmogorov-Smirnov statistic for continuous parameters
pdf("KS.RunI_MC3.pdf", width = 11.1, height = 11.3)
plotKS(checkRunI_MC3)
dev.off()
### All parameters are below the threshold. Thus all runs have converged to the same posterior.
