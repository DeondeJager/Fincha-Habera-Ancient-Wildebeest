

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
back_blue_modern <- fread('back_blue_modern.csv')
back_blue_fossil <- back_fossil[str_detect(back_fossil$ID, pattern='Ct_'),] # Cg indicates Connochaetes taurinus
nrow(back_blue_modern); nrow(back_blue_fossil)
# [1] 36013
# [1] 72145
#   occurrence data
occ_blue_modern <- fread('occ_blue_modern.csv')
occ_blue_fossil <- occ_fossil[str_detect(occ_fossil$ID, pattern='Ct_'),] # Cg indicates Connochaetes taurinus
nrow(occ_blue_modern); nrow(occ_blue_fossil)
# [1] 628
# [1] 66
## fincha habera data
fh_mu <- fread('site_finchahabera_mu.csv')
fh_mu_sd <- fread('site_finchahabera_mu_sd.csv')


### select 10% of the background grid cells to so the ratio of FossilGridCells:ModernGridCells more closely aligns with FossilOccs:ModernOccs
back_blue_fossil_list <- split(back_blue_fossil, f=back_blue_fossil$ID)
length(back_blue_fossil_list)
# [1] 66
for(i in 1:66){
  set.seed(20260202)
  a <- back_blue_fossil_list[[i]]
  a <- a[sample(1:nrow(a), floor(nrow(a)/10), replace = FALSE),]
  back_blue_fossil_list[[i]] <- a
}
back_blue_fossil <- do.call('rbind', back_blue_fossil_list)
nrow(back_blue_fossil)
# [1] 7178

### merging occ data and environmental data together
# black wildebeest 
back_blue <- rbind(back_blue_modern, back_blue_fossil)
occ_blue <- rbind(occ_blue_modern, occ_blue_fossil)


# color for fincha habera
fh_col <- '#E5984E'

PerformanceAnalytics::chart.Correlation(occ_blue[,9:18])



##### Setting up the models ------------------------------------------------------------------------------------------------

# sourcing the R code to run all the functions
setwd("~/Library/CloudStorage/Box-Box/ADL_POSTDOC/_Projects/_Wildebeest/Rcode")
source("NAF_Functions.R")

# working directory for models
setwd("~/Library/CloudStorage/Box-Box/ADL_POSTDOC/_Projects/_Wildebeest/WildebeestModels/BlueModels")

# create null df - empty dataframe that null data will go into later
blue_W_null <- create_null_df(taxon_name = 'blue_W', 
                              timeBin = 'Quat', 
                              extent = 'Africa'
                              )

# create summary df
blue_W_summary <- create_summary_df(taxon_name = 'blue_W', 
                                    timeBin = 'Quat', 
                                    extent = 'Africa',
                                    fold = 'train',
                                    f_class = c('LQP', 'Q'),
                                    beta_values = c(0.025, 0.05, 0.1, 0.25, 0.5, 0.75, 1, 1.5, 2, 2.5, 3, 3.5, 4, 4.5, 5)
                                    )

make_maxent_folders(blue_W_summary)

# run null model
blue_W_null  <- null_aic(null_df = blue_W_null,
                         occs = occ_blue,
                         background = back_blue,
                         first_occ_col = 9)
blue_W_null
#    taxon timeBin extent  fold features betas   n k      ln_L      AIC     AICc Boyce   lambdas
# 1 blue_W    Quat Africa train       NA    NA 694 1 -7407.331 14816.66 14816.67    NA Intercept


##### Optimizing Model Parameters ------------------------------------------------------------------------------------------

colnames(occ_blue)
#  [1] "taxon"    "x"        "y"        "ID"       "prob"     "timeBin"  "extent"   "cellID"   "bio1"     "bio4"    
# [11] "bio10"    "bio11"    "bio12"    "bio13"    "bio14"    "bio15"    "bio16"    "bio17"    "presence" "fold" 

# find most optimized hyperparameter settings
setwd("~/Library/CloudStorage/Box-Box/ADL_POSTDOC/_Projects/_Wildebeest/WildebeestModels/BlueModels")
blue_W_optim1 <- optimize_maxent_likelihood(summary_df = blue_W_summary,   # name of the output summary file
                                           occs = occ_blue, # species occurrences
                                           background = back_blue, # background data
                                           presenceCol = 19, # `presence` column number
                                           predicCols = c(10:11,17:18), # column numbers of predictor variables
                                           Boyce = 'cr' # smooth boyce spline method (Liu et al. 2024 Ecography)
                                           )
Sys.time()
# this takes ~1 minute total
blue_W_optim |> dplyr::arrange(AICc) |> View()

# the model with LQP feature classes and beta multiplier of 0.025 seems to be the best
# the response curves look reasonable!



##### Running Cross-Validation ---------------------------------------------------------------------------------------------

### LQP_0.025
# create eval object
# default f_class = LQP, so if f_class argument is missing, will assume LQP; only need it there in case NOT LQP
blue_W_eval_1 <- create_eval_df(taxon_name = 'blue_W', 
                                       timeBin = 'Quat', 
                                       extent = 'Africa', 
                                       beta_values = 0.025, 
                                       f_class = 'LQP')
blue_W_eval_1

# make the folders where the maxent objects will go
setwd("~/Library/CloudStorage/Box-Box/ADL_POSTDOC/_Projects/_Wildebeest/WildebeestModels/BlueModels")
make_maxent_folders(blue_W_eval_1)

# run eval object
blue_W_eval_1 <- maxent_crossval_error(eval_df = blue_W_eval_1,
                                              occs = occ_blue, # species occurrences
                                              background = back_blue, # background data
                                              presenceCol = 19, # `presence` column number
                                              predicCols = c(10:11,17:18), # column numbers of predictor variables
                                              Boyce = 'cr',  # smooth boyce spline method (Liu et al. 2024 Ecography)
                                              omission_rate = 0.1, # 10% of occurrence data as the threshold
                                              all_background = back_blue)

blue_W_eval_1$summary
#    taxon timeBin extent   fold features betas   n    thresh  testSens testBoyce
# 1 blue_W    Quat Africa test_A      LQP 0.025 557 0.2723273 0.9197080 1.0000000
# 2 blue_W    Quat Africa test_B      LQP 0.025 557 0.2901412 0.8540146 0.9936452
# 3 blue_W    Quat Africa test_C      LQP 0.025 557 0.2715157 0.9124088 0.9915628
# 4 blue_W    Quat Africa test_D      LQP 0.025 557 0.2667596 0.9197080 1.0000000
# 5 blue_W    Quat Africa test_E      LQP 0.025 557 0.3199434 0.8832117 0.9902176
# lambdas
# 1 bio10, bio16, bio17, bio4, bio10^2, bio16^2, bio17^2, bio4^2, bio10*bio16, bio10*bio17, bio10*bio4, bio16*bio17, bio16*bio4, bio17*bio4
# 2 bio10, bio16, bio17, bio4, bio10^2, bio16^2, bio17^2, bio4^2, bio10*bio16, bio10*bio17, bio10*bio4, bio16*bio17, bio16*bio4, bio17*bio4
# 3 bio10, bio16, bio17, bio4, bio10^2, bio16^2, bio17^2, bio4^2, bio10*bio16, bio10*bio17, bio10*bio4, bio16*bio17, bio16*bio4, bio17*bio4
# 4 bio10, bio16, bio17, bio4, bio10^2, bio16^2, bio17^2, bio4^2, bio10*bio16, bio10*bio17, bio10*bio4, bio16*bio17, bio16*bio4, bio17*bio4
# 5 bio10, bio16, bio17, bio4, bio10^2, bio16^2, bio17^2, bio4^2, bio10*bio16, bio10*bio17, bio10*bio4, bio16*bio17, bio16*bio4, bio17*bio4



##### Evaluate Model, Calculate Threshold and Model Variability ------------------------------------------------------------

### LQP 0.1
# evaluate model to calculate weighted suitibility and stdev
blue_W_means_1 <- maxent_eval(eval = blue_W_eval_1)
blue_thresh <- quantile(blue_W_means_1$occ$w_mean, 0.1) |> as.numeric()
blue_thresh
# [1] 0.3022016
# this is the threshold
# save this value!!
quantile(blue_W_means_1$occ$w_mean, seq(0,1,0.05)) |> as.numeric() |> round(3)
#  [1] 0.003 0.094 0.302 0.401 0.467 0.511 0.560 0.596 0.629 0.667 0.711 0.749 0.787
# [14] 0.826 0.853 0.878 0.899 0.914 0.924 0.936 0.955

### plot model stuitability/uncertainty plots
# the only input you need is the object made from the maxent_eval function
# suit_uncert_plot(blue_W_means_1)


##### Saving Everything ----------------------------------------------------------------------------------------------------

# null model
write.csv(blue_W_null, 'blue_W.null.csv', row.names = F)
# model optimized summary object
write.csv(blue_W_optim, 'blue_W.summary.csv', row.names = F)
# summary evaluation object
write.csv(blue_W_eval_LQP0.75$summary, 'blue_W.evaluation.csv', row.names=F)


# work on projecting the model to all the extents in another script

##### END --------------------------------------------------------------------------------





##### Projecting Model to all extents --------------------------------------------------------------------------------------

blue_W_everything_LQP0.1 <- maxent.everything(eval = blue_W_eval_LQP0.1,
                                              thresh = blue_W_thresh_LQP0.1,
                                              everything=kpgenm.back,
                                              predic = c(10:11,17:18))
View(blue_W_everything_LQP0.1)



##### Running Informed MESS Analysis ---------------------------------------------------------------------------------------

# defining the tolerance vector
# 1 = can extrapolate to non-analog conditions
# 0 = cannot extrapolate to non-analog conditions
blue_W_tolerance <- c(1,0,   # bio10 - avg temp warmest quarter
                      1,1,   # bio11 - avg temp coldest quarter
                      0,1,   # bio16 - precip wettest quarter
                      1,1)   # bio17 - precip driest quarter


# running mess
blue_W_mess <- informed.mess(ref.extent = back_blue, 
                             mess.extent = kpgenm.back,
                             coord.cols = 2:3,
                             predic = c(10:11,17:18),
                             tolerance = blue_W_tolerance)





##### Saving Everything ----------------------------------------------------------------------------------------------------

# null model
write.csv(blue_W_null, 'Tyrannosaurus.null.csv', row.names = F)
# model optimized summary object
write.csv(blue_W_optim, 'Tyrannosaurus.summary.csv', row.names = F)
# model projected to all extents
write.csv(blue_W_everything, 'Tyrannosaurus.everything.csv', row.names = F)
# model mess analysis
write.csv(blue_W_mess, 'Tyrannosaurus.mess.csv', row.names = F)

# record model threshold

##### END ------------------------------------------------------------------------------------------------------------------

