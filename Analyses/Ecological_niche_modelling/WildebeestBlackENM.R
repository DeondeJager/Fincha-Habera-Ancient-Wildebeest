

##### Packages -------------------------------------------------------------------------------------------------------------

# loading in the requisite packages
packs <- c('data.table', 'DescTools', 'dplyr', 'enmSdmX', 'ggplot2', 'matrixStats',
           'mgcv', 'predicts', 'rJava', 'sf', 'splitstackshape', 'stringr')

# loading in the packages
lapply(packs, library, character.only=T)

# giving the package versions
packs.df <- as.data.frame(matrix(NA, nrow=length(packs), ncol=2))
colnames(packs.df) <- c('pkg.name', 'pkg.version')
for(i in 1:length(packs)){
  packs.df[i,1] <- packs[i]
  packs.df[i,2] <- as.character(packageVersion(packs[i]))
}
packs.df

# test if the maxent.jar file is in the predicts package folder
list.files( system.file('java', package = 'predicts'), pattern='jar' )
# [1] "dismo.jar"  "maxent.jar"



##### Loading in the Data --------------------------------------------------------------------------------------------------

# where data in environmental space is
setwd("~/Library/CloudStorage/Box-Box/ADL_POSTDOC/_Projects/_Wildebeest/Occ_Data/ForModeling")
list.files()
#  [1] "back_black_modern.csv"            "back_blue_modern.csv"            
#  [3] "back_fossil_mu.csv"               "back_fossil.csv"                 
#  [5] "BlackVariableCorrelation.pdf"     "BlueVariableCorrelation.pdf"     
#  [7] "occ_black_modern.csv"             "occ_blue_modern.csv"             
#  [9] "occ_fossil_mu.csv"                "occ_fossil.csv"                  
# [11] "site_finchahabera_mu_sd.csv"      "site_finchahabera_mu.csv"        
# [13] "site_finchahabera.csv"            "WildebeestEnvironmentalSpace.png"


### loading in data
# the fossil data has been merged - split into constituent species using the `ID` column
occ_fossil <- fread('occ_fossil_mu.csv')
back_fossil <- fread('back_fossil_mu.csv')

##  black wildebeest
#   background data
back_black_modern <- fread('back_black_modern.csv')
back_black_fossil <- back_fossil[str_detect(back_fossil$ID, pattern='Cg_'),] # Cg indicates Connochaetes gnou
nrow(back_black_modern); nrow(back_black_fossil)
# [1] 5897
# [1] 81042
#   occurrence data
occ_black_modern <- fread('occ_black_modern.csv')
occ_black_fossil <- occ_fossil[str_detect(occ_fossil$ID, pattern='Cg_'),] # Cg indicates Connochaetes gnou
nrow(occ_black_modern); nrow(occ_black_fossil)
# [1] 73
# [1] 90
## fincha habera data
fh_mu <- fread('site_finchahabera_mu.csv')
fh_mu_sd <- fread('site_finchahabera_mu_sd.csv')


### there are way more backgrounds cells for the fossils than there are for the modern
# select 10% of the background grid cells to so the ratio of FossilGridCells:ModernGridCells more closely aligns with FossilOccs:ModernOccs
back_black_fossil_list <- split(back_black_fossil, f=back_black_fossil$ID)
length(back_black_fossil_list)
# [1] 90
for(i in 1:90){
  set.seed(20260202)
  a <- back_black_fossil_list[[i]]
  a <- a[sample(1:nrow(a), floor(nrow(a)/10), replace = FALSE),]
  back_black_fossil_list[[i]] <- a
}
back_black_fossil <- do.call('rbind', back_black_fossil_list)
nrow(back_black_fossil)
# [1] 8077

### merging occ data and environmental data together
# black wildebeest 
back_black <- rbind(back_black_modern, back_black_fossil)
occ_black <- rbind(occ_black_modern, occ_black_fossil)


# color for fincha habera
fh_col <- '#E5984E'

PerformanceAnalytics::chart.Correlation(occ_black[,9:18])



##### Setting up the models ------------------------------------------------------------------------------------------------



# sourcing the R code to run all the functions
setwd("~/Library/CloudStorage/Box-Box/ADL_POSTDOC/_Projects/_Wildebeest/Rcode")
source("NAF_Functions.R")

# working directory for models
setwd("~/Library/CloudStorage/Box-Box/ADL_POSTDOC/_Projects/_Wildebeest/WildebeestModels/BlackModels")

# create null df - empty dataframe that null data will go into later
blackW_null <- create_null_df(taxon_name = 'blackW', 
                              timeBin = 'Quat', 
                              extent = 'Africa'
                              )

# create summary df
blackW_summary <- create_summary_df(taxon_name = 'blackW', 
                                    timeBin = 'Quat', 
                                    extent = 'Africa',
                                    fold = 'train',
                                    f_class = c('LQP', 'Q'),
                                    beta_values = c(0.025, 0.05, 0.1, 0.25, 0.5, 0.75, 1, 1.5, 2, 2.5, 3, 3.5, 4, 4.5, 5)
                                    )

make_maxent_folders(blackW_summary)

# run null model
blackW_null  <- null_aic(null_df = blackW_null,
                         occs = occ_black,
                         background = back_black,
                         first_occ_col = 9)
#   taxon timeBin extent  fold features betas   n k      ln_L      AIC    AICc Boyce   lambdas
# 1 blackW    Quat Africa train       NA    NA 163 1 -1853.793 3709.586 3709.61    NA Intercept


##### Optimizing Model Parameters ------------------------------------------------------------------------------------------

colnames(occ_black)
#  [1] "taxon"    "x"        "y"        "ID"       "prob"     "timeBin"  "extent"   "cellID"   "bio1"     "bio4"    
# [11] "bio10"    "bio11"    "bio12"    "bio13"    "bio14"    "bio15"    "bio16"    "bio17"    "presence" "fold" 

# find most optimized hyperparameter settings
setwd("~/Library/CloudStorage/Box-Box/ADL_POSTDOC/_Projects/_Wildebeest/WildebeestModels/BlackModels")
blackW_optim <- optimize_maxent_likelihood(summary_df = blackW_summary,   # name of the output summary file
                                           occs = occ_black, # species occurrences
                                           background = back_black, # background data
                                           presenceCol = 19, # `presence` column number
                                           predicCols = c(10:11,17:18), # column numbers of predictor variables
                                           Boyce = 'cr' # smooth boyce spline method (Liu et al. 2024 Ecography)
                                           )
# this takes ~30 seconds total
blackW_optim |> dplyr::arrange(AICc) |> View()

# the model with LQP feature classes and beta multiplier of 0.1 seems to be the best based on AICc
# Boyce index looks really nice too
# bio17 (column 18) was displaying concave-down behavior, which is non-sensical

# re-running without bio17

# create summary df
blackW_summary <- create_summary_df(taxon_name = 'blackW', 
                                    timeBin = 'Quat', 
                                    extent = 'Africa',
                                    fold = 'train',
                                    f_class = c('LQP', 'Q'),
                                    beta_values = c(0.025, 0.05, 0.1, 0.25, 0.5, 0.75, 1, 1.5, 2, 2.5, 3, 3.5, 4, 4.5, 5)
)
# find most optimized hyperparameter settings
setwd("~/Library/CloudStorage/Box-Box/ADL_POSTDOC/_Projects/_Wildebeest/WildebeestModels/BlackModels")
blackW_optim <- optimize_maxent_likelihood(summary_df = blackW_summary,   # name of the output summary file
                                           occs = occ_black, # species occurrences
                                           background = back_black, # background data
                                           presenceCol = 19, # `presence` column number
                                           predicCols = c(10:11,17), # column numbers of predictor variables
                                           Boyce = 'cr' # smooth boyce spline method (Liu et al. 2024 Ecography)
)
# this takes ~20 seconds total
blackW_optim |> dplyr::arrange(AICc) |> View()

# the model with LQP feature classes and beta multiplier of 0.75 seems to be the best
# bio10 has > 80 % contribution - try with bio10 and bio11


# create summary df
blackW_summary <- create_summary_df(taxon_name = 'blackW', 
                                    timeBin = 'Quat', 
                                    extent = 'Africa',
                                    fold = 'train',
                                    f_class = c('LQP', 'Q'),
                                    beta_values = c(0.025, 0.05, 0.1, 0.25, 0.5, 0.75, 1, 1.5, 2, 2.5, 3, 3.5, 4, 4.5, 5)
)
# find most optimized hyperparameter settings
setwd("~/Library/CloudStorage/Box-Box/ADL_POSTDOC/_Projects/_Wildebeest/WildebeestModels/BlackModels")
blackW_optim <- optimize_maxent_likelihood(summary_df = blackW_summary,   # name of the output summary file
                                           occs = occ_black, # species occurrences
                                           background = back_black, # background data
                                           presenceCol = 19, # `presence` column number
                                           predicCols = c(11:12,17), # column numbers of predictor variables
                                           Boyce = 'cr' # smooth boyce spline method (Liu et al. 2024 Ecography)
)
# this takes ~20 seconds total
blackW_optim |> dplyr::arrange(AICc) |> View()

# the model with LQP feature classes and beta multiplier of 0.1 seems to be the best
# bio10 has nonsensical concave-down behavior
# these models also have greater AICc values
# going back to bio4, bio10, and bio16

# create summary df
blackW_summary <- create_summary_df(taxon_name = 'blackW', 
                                    timeBin = 'Quat', 
                                    extent = 'Africa',
                                    fold = 'train',
                                    f_class = c('LQP', 'Q'),
                                    beta_values = c(0.025, 0.05, 0.1, 0.25, 0.5, 0.75, 1, 1.5, 2, 2.5, 3, 3.5, 4, 4.5, 5)
)
# find most optimized hyperparameter settings
setwd("~/Library/CloudStorage/Box-Box/ADL_POSTDOC/_Projects/_Wildebeest/WildebeestModels/BlackModels")
blackW_optim <- optimize_maxent_likelihood(summary_df = blackW_summary,   # name of the output summary file
                                           occs = occ_black, # species occurrences
                                           background = back_black, # background data
                                           presenceCol = 19, # `presence` column number
                                           predicCols = c(10:11,17), # column numbers of predictor variables
                                           Boyce = 'cr' # smooth boyce spline method (Liu et al. 2024 Ecography)
)
# this takes ~20 seconds total
blackW_optim |> dplyr::arrange(AICc) |> View()

# going with the model with LQP feature classes and beta multiplier of 0.75


##### Running Cross-Validation ---------------------------------------------------------------------------------------------

### LQP_0.75
# create eval object
# default f_class = LQP, so if f_class argument is missing, will assume LQP; only need it there in case NOT LQP
blackW_eval_1 <- create_eval_df(taxon_name = 'blackW', 
                                      timeBin = 'Quat', 
                                      extent = 'Africa', 
                                      beta_values = 0.75, 
                                      f_class = 'LQP')
blackW_eval_1

# make the folders where the maxent objects will go
setwd("~/Library/CloudStorage/Box-Box/ADL_POSTDOC/_Projects/_Wildebeest/WildebeestModels/BlackModels")
make_maxent_folders(blackW_eval_1)

# run eval object
blackW_eval_1 <- maxent_crossval_error(eval_df = blackW_eval_1,
                                             occs = occ_black, # species occurrences
                                             background = back_black, # background data
                                             presenceCol = 19, # `presence` column number
                                             predicCols = c(10:11,17), # column numbers of predictor variables
                                             Boyce = 'cr',  # smooth boyce spline method (Liu et al. 2024 Ecography)
                                             omission_rate = 0.1, # 10% of occurrence data as the threshold
                                             all_background = back_black)

blackW_eval_1$summary
#    taxon timeBin extent   fold features betas   n    thresh  testSens testBoyce
# 1 blackW    Quat Africa test_A      LQP  0.75 133 0.2378496 0.9333333 0.9010454
# 2 blackW    Quat Africa test_B      LQP  0.75 133 0.2696452 0.7666667 0.6556484
# 3 blackW    Quat Africa test_C      LQP  0.75 133 0.2194937 0.9000000 1.0000000
# 4 blackW    Quat Africa test_D      LQP  0.75 133 0.2547604 0.8666667 1.0000000
# 5 blackW    Quat Africa test_E      LQP  0.75 133 0.2562006 0.8666667 1.0000000
#                                                                lambdas
# 1              bio10, bio10^2, bio16^2, bio4^2, bio10*bio4, bio16*bio4
# 2          bio10, bio16^2, bio4^2, bio10*bio16, bio10*bio4, bio16*bio4
# 3              bio10, bio10^2, bio16^2, bio4^2, bio10*bio4, bio16*bio4
# 4 bio10, bio10^2, bio16^2, bio4^2, bio10*bio16, bio10*bio4, bio16*bio4
# 5              bio10, bio10^2, bio16^2, bio4^2, bio10*bio4, bio16*bio4




##### Evaluate Model, Calculate Threshold and Model Variability ------------------------------------------------------------

### LQP 0.1
# evaluate model to calculate weighted suitibility and stdev
blackW_means_1 <- maxent_eval(eval = blackW_eval_1)
black_thresh <- quantile(blackW_means_1$occ$w_mean, 0.1) |> as.numeric()
black_thresh
# [1] 0.2283123
# this is the threshold
# save this value!!

quantile(blackW_means_1$occ$w_mean, seq(0,1,0.05)) |> as.numeric() |> round(3)
#  [1] 0.085 0.190 0.228 0.310 0.365 0.413 0.489 0.534 0.567 0.621 0.645 0.674 0.727
# [14] 0.746 0.806 0.854 0.926 0.975 0.986 0.990 1.000


### plot model stuitability/uncertainty plots
# the only input you need is the object made from the maxent_eval function
# suit_uncert_plot(blackW_means_1)


##### Saving Everything ----------------------------------------------------------------------------------------------------

# null model
write.csv(blackW_null, 'blackW.null.csv', row.names = F)
# model optimized summary object
write.csv(blackW_optim, 'blackW.summary.csv', row.names = F)
# summary evaluation object
write.csv(blackW_eval_LQP0.75$summary, 'blackW.evaluation.csv', row.names=F)


# work on projecting the model to all the extents in another script

##### END --------------------------------------------------------------------------------





##### Projecting Model to all extents --------------------------------------------------------------------------------------

blackW_to_fh <- maxent_everything(eval = blackW_eval_1,
                                  thresh = black_thresh,
                                  everything=j_back,
                                  predic = c(10:11,17))
View(blackW_to_fh)



##### Running Informed MESS Analysis ---------------------------------------------------------------------------------------

# defining the tolerance vector
# 1 = can extrapolate to non-analog conditions
# 0 = cannot extrapolate to non-analog conditions
blackW_tolerance <- c(1,0,   # bio10 - avg temp warmest quarter
                      1,1,   # bio11 - avg temp coldest quarter
                      0,1,   # bio16 - precip wettest quarter
                      1,1)   # bio17 - precip driest quarter


# running mess
blackW_mess <- informed.mess(ref.extent = back_black, 
                             mess.extent = kpgenm.back,
                             coord.cols = 2:3,
                             predic = c(10:11,17),
                             tolerance = blackW_tolerance)





##### Saving Everything ----------------------------------------------------------------------------------------------------

# null model
write.csv(blackW_null, 'Tyrannosaurus.null.csv', row.names = F)
# model optimized summary object
write.csv(blackW_optim, 'Tyrannosaurus.summary.csv', row.names = F)
# model projected to all extents
write.csv(blackW_everything, 'Tyrannosaurus.everything.csv', row.names = F)
# model mess analysis
write.csv(blackW_mess, 'Tyrannosaurus.mess.csv', row.names = F)

# record model threshold

##### END ------------------------------------------------------------------------------------------------------------------

