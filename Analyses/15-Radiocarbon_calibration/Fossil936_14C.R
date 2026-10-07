##### Calibrating 14C Dates at Fincha Habera -------------------------------

# Nick Freymueller - nicholas.freymueller@adelaide.edu.au
# 6 October 2026

# libraries
library(stringr)
library(purrr)
library(dplyr)
# radiocarbon calibration
library(oxcAAR)
quickSetupOxcal() # ensures your R points to the OxCal software (and downloads it if the first time running)

##### oxcalibrate function -------------------------------

### function that calibrates one or more radiocarbon dates, returning summary stats, probability densities, and diagnostics
## function uses OxCal v4.4.4 (c) Bronk Ramsey (2021)
## function uses the IntCal20 Curve from Reimer et al (2020)
#  function CAN NOT run any other calibration curves (e.g., Marine20) - this is a restriction of the `oxcAAR` package
#  if incorporating reservoir effects, the user MUST apply them first
## the function returns dates in Calibrated years Before Present (Cal yrs BP)
#  this function was designed because `oxcAAR` annoyingly only reports Calibrated years Before Common Era (BCE)
## currently the only functionality in `oxcalibrate` is to run

### ARGUMENTS
## df <- a 3-column data frame (column order MUST be age, error, ID)
## cut_sigma <- threshold (in terms of number of standard deviations)
#               threshold MUST between [0.6744897, 4.417], corresponding to returned probability ranges of (0.5, 0.99999)
#               beyond this specified threshold, the radiocarbon probability is set to zero, then the sequence re-scaled to sum to 1
#               default value is 2.576, representing 99.00% of the probability density
#               if not specified (NULL), function sets cut_sigma == 4.417
#               not recommended to set < 1.96 (< 95%) as this removes too much probability for commonly used analyses
#               not recommended to set > 3.5 (> 99.95%) as this results in many time bins with negligible probability (e.g., p < 1e-4)
#               to figure out a custom probability range, use: `pnorm(cut_sigma) - pnorm(-cut_sigma)`
## agg_years <- optional value to aggregate (i.e., sum) probabilities across (e.g., 30 years)
#               the minimum value is 5 years (this is a restriction of the `oxcAAR` package)
#               will return modified `probabilities` data frames

### RETURNS
## dataframe = a data frame with summary statistics of calibrated dates (median, mean, sd, and ± 1/2/3 sigma probabilities)
#              these dates are in Calibrated years Before Present (Cal yrs BP)
## probabilities = a list of length == nrow(df), where each object is a data.frame
#                  each data.frame is the individual & cumulative probabilities for each date at each year
#                  each data.frame is aggregated into `agg_years` time bins (e.g., 30 years; if provided)
#                  the `CalYearsBP` in each data.frame represents all years concluding at that year
#                  e.g., for agg_years == 30, `44460` represents the probabilities of sum(44460:44489) Cal years BP
## noncalibrated = indicates whether any dates could NOT be calibrated
#                  if all dates could be calibrated, the function returns "All dates could be calibrated."

oxcalibrate <- function(df, cut_sigma = 2.576, agg_years = 30){
     require(oxcAAR)
     require(dplyr)
     require(stringr)
     require(purrr)
     # checks to ensure 
     df <- as.data.frame(df)
     # checking to see if the agg_years 
     if(!is_null(agg_years)){
          if(agg_years%%5 != 0){
               stop('`agg_years` must be divisible by 5.')
          }
     } else {
          agg_years <- 5
     }
     # cut_sigma
     if(!is.numeric(cut_sigma) & !is.null(cut_sigma) ){
          stop('`cut_sigma` must be a positive number or `NULL`.')
     }
     if(cut_sigma < 0.6744897 | cut_sigma > 4.417){
          stop('`cut_sigma` must be ≥ 0.6744897 and ≤ 4.417.')
     }
     if(is.null(cut_sigma) ){
          cut_sigma <- 4.417
     }
     # probability ranges cor cut_sigma
     p_low <- pnorm(-cut_sigma)
     p_high <- pnorm(cut_sigma)
     
     ### this section is based off of a worked example from July Pilowsky
     # similar scripts exist elsewhere on their fossil-database GitHub repo
     # https://github.com/japilo/fossil-database/blob/be9457b4e9cf75a84d2c609ae982aa191bf72a05/Database%20fixes/merge%20fix.Rmd
     # translates the radiocarbon data to Javascript format, then  runs it through OxCal and extracts the mean + sd in Cal yrs BP
     # outputs a 3*nrow(dates) character string with the mean, sd, and ID of the date in Cal yrs BP
     calibrated_ages <- oxcAAR::R_Date(r_dates = df[,1],  # mean 14C years BP 
                                       stds = df[,2],     # uncertainty
                                       names = df[,3]) |> # unique ID for merging back w/original data later
          # useful in case any date can't be converted
          executeOxcalScript() |> readLines() |> 
          # keep the 6 variables measured
          purrr::keep(~any(str_detect(., "likelihood.mean"), str_detect(., "likelihood.sigma"),
                           str_detect(., "R_Date\\("), str_detect(., "likelihood.median"),
                           str_detect(., "likelihood.start"), str_detect(., "likelihood.prob=") ) )
     
     # identifies the rows in df that were able to be calibrated
     # this may fail if dates are very young (<200 years), very old (>50,000), or have non-numeric values for the mean or sd values
     # the magrittr pipe (%>%) is required for piping with a '.' 
     # https://tidyverse.org/blog/2023/04/base-vs-magrittr-pipe/
     complete <- calibrated_ages |> 
          # extract numeric values from the calibrated vector object
          map_chr(str_extract, "\\[[0-9]+\\]") |> map_chr(str_extract, "[0-9]+") |> 
          # make a tibble of the number of vars (mean, sigma, median, R_Date, start, & probability) extracted from each specimen
          as.numeric() |> as_tibble() |> group_by(value) |> summarize(count = n()) |> 
          # only filter to specimens that recorded all 6 variables
          filter(count==6) %>% .$value
     
     # makes a data.frame of everything
     calibrated <- data.frame(ID = calibrated_ages |> purrr::keep(~str_detect(., "comment")) |>
                                   str_split("\\\"") |> map(2) |> str_split(" ") |>
                                   map_chr(1) %>% .[complete],   # unique ID associated with the specimen
                              OxCal_median = calibrated_ages |> purrr::keep(~str_detect(., "likelihood.median")) |>
                                   map(str_split, "=") |> map(1) |> map(2) |> map(str_remove, ";") |> map_dbl(as.numeric) |> 
                                   map_dbl(function(x) {abs(x-1950)}) |> round(),  # OxCal median (Cal Years BP) 
                              OxCal_mean = calibrated_ages |> purrr::keep(~str_detect(., "likelihood.mean")) |>
                                   map(str_split, "=") |> map(1) |> map(2) |> map(str_remove, ";") |> map_dbl(as.numeric) |> 
                                   map_dbl(function(x) {abs(x-1950)}) |> round(),  # OxCal mean (Cal Years BP) 
                              OxCal_error = calibrated_ages |> purrr::keep(~str_detect(., "likelihood.sigma")) |>
                                   map(str_split, "=") |> map(1) |> map(2) |> map(str_remove, ";") |>
                                   map_dbl(as.numeric) |> round(),   # OxCal uncertainty (Cal Years BP)
                              OxCal_start = calibrated_ages |> purrr::keep(~str_detect(., 'likelihood.start')) |>
                                   map(str_split, "=") |> map(1) |> map(2) |> map(str_remove, ";") |> map_dbl(as.numeric) |> 
                                   map_dbl(function(x) {abs(x-1950)}) %>% .[complete] |> ceiling() |> round(),   # start date for OxCal distribution
                              OxCal_prob = calibrated_ages |> purrr::keep(~str_detect(., 'likelihood.prob=')) |>
                                   map(str_split, "=") |> map(1) |> map(2) |> map(str_remove, ";") |>
                                   map_chr(1) %>% .[complete],   # string of OxCal probabilities
                              P_3sigmaMinus = NA, # probability corresponding to -3 sigma (if probability were normal)
                              P_2sigmaMinus = NA, # probability corresponding to -2 sigma (if probability were normal)
                              P_1sigmaMinus = NA, # probability corresponding to -1 sigma (if probability were normal)
                              P_1sigmaPlus = NA,  # probability corresponding to +1 sigma (if probability were normal)
                              P_2sigmaPlus = NA,  # probability corresponding to +2 sigma (if probability were normal)
                              P_3sigmaPlus = NA ) # probability corresponding to +3 sigma (if probability were normal)
     # assigning if 
     # extracts the probabilities from calibrated$OxCal_prob and assigns each of them to a list
     prob_list <- list()
     for(i in 1:length(complete)){
          # extract and reformat probabbilities
          probs <- calibrated$OxCal_prob[i]
          probs <- as.numeric(unlist(strsplit(gsub("\\[|\\]", "", probs),",")))
          # number of probabilities to assign
          n_probs <- length(probs)
          # temporaryy data frame
          temp <- data.frame(CalYearsBP = seq(from=round(calibrated$OxCal_start[i],0), by=-5, length.out=n_probs),
                             Probabilities = round(probs/sum(probs),8),
                             CumProb = NA )
          # cumulative probabilities
          temp$CumProb <- cumsum(temp$Probabilities)
          
          # remove any years that go into negative CalYrsBP
          temp$CalYearsBP[temp$CalYearsBP<0] <- 0
          temp <- unique(temp)
          # assigning confidence interval bounds to calibrated dates
          calibrated$P_3sigmaMinus[i] <- temp[which.min( abs(pnorm(-3)-temp$CumProb) ),1]
          calibrated$P_2sigmaMinus[i] <- temp[which.min( abs(pnorm(-2)-temp$CumProb) ),1]
          calibrated$P_1sigmaMinus[i] <- temp[which.min( abs(pnorm(-1)-temp$CumProb) ),1]
          calibrated$P_1sigmaPlus[i] <- temp[which.min( abs(pnorm(1)-temp$CumProb) ),1]
          calibrated$P_2sigmaPlus[i] <- temp[which.min( abs(pnorm(2)-temp$CumProb) ),1]
          calibrated$P_3sigmaPlus[i] <- temp[which.min( abs(pnorm(3)-temp$CumProb) ),1]
          
          ## truncating radiocarbon probability
          idx_low <- which.min( abs(p_low-temp$CumProb) )
          idx_high <- which.min( abs(p_high-temp$CumProb) )
          temp <- temp[idx_low:idx_high,]
          temp$Probabilities <- temp$Probabilities / sum(temp$Probabilities)
          temp$CumProb <- cumsum(temp$Probabilities)
          
          # re-aggregating data to larger time bins (e.g., 30 years) 
          if(agg_years != 5){
               temp$CalYearsBP <- temp$CalYearsBP - (temp$CalYearsBP %% agg_years)
               temp <- aggregate(temp$Probabilities, by=list(CalYearsBP=temp$CalYearsBP), FUN='sum') |> 
                    dplyr::arrange(-CalYearsBP)
               colnames(temp)[2] <- 'Probabilities'
               temp$CumProb <- cumsum(temp$Probabilities)
          }
          # assigning ID
          temp$ID <- df[i,3]
          
          # assigning temp to prob_list
          prob_list[[i]] <- temp
          
     } # end i loop
     
     # reordering/renaming/rounding columns
     calibrated <- calibrated[,-5:-6]
     colnames(calibrated) <- c("ID", "CalYrBP_median", "CalYrBP_mean", "CalYrBP_sd",
                               "P_3sigmaMinus", "P_2sigmaMinus", "P_1sigmaMinus", "P_1sigmaPlus", "P_2sigmaPlus", "P_3sigmaPlus")
     calibrated[,2:10] <- calibrated[,2:10] |> round(0)
     # renaming calibration curves
     names(prob_list) <- calibrated$ID
     
     # making the output list
     out_list <- list()
     # calibrated dates with metadata
     out_list$dataframe <- calibrated
     # list with PDFs fo each successfully calibrated specimen
     out_list$probabilities <- prob_list
     # vector of rows in the original data frame that could not be calibrated
     noncalibrated <- (1:nrow(df))[-complete]
     if( length(noncalibrated) == 0 ){
          out_list$noncalibrated <- 'All dates could be calibrated.'
     } else {
          out_list$noncalibrated <- noncalibrated
     }
     
     
     return(out_list)
     
}

##### calibrating 14C dates -------------------------------

# fincha habera fossil details
# 37,330 ± 620 BP (UCIAMS-251302)
# specimen is named fossil 936
f936 <- data.frame(age=37330, sd=620, ID='UCIAMS-251302')


### calibrating the fincha habera fossil
# remove the lowest & highest 0.5% probability density
# aggregate 14C probability into 30-year bins to match the climate data

f936_cal <- oxcalibrate(f936)

# checking that all dates could be calibrated - YES!
f936_cal$noncalibrated
# [1] "All dates could be calibrated."

# summary of calibration
f936_cal$dataframe
#              ID CalYrBP_median CalYrBP_mean CalYrBP_sd  P0_5  P2_5   P16   P84 P97_5 P99_5
# 1 UCIAMS-251302          41928      41891          343 42705 42490 42225 41530 41190 40940

f936_cal$probabilities[[1]] |> as.data.table()
#     CalYearsBP Probabilities      CumProb            ID
#          <num>         <num>        <num>        <char>
#  1:      42690  0.0008428347 0.0008428347 UCIAMS-251302
#  2:      42660  0.0015073913 0.0023502260 UCIAMS-251302
#  3:      42630  0.0018467692 0.0041969952 UCIAMS-251302
#  4:      42600  0.0022763920 0.0064733873 UCIAMS-251302
#  5:      42570  0.0028198309 0.0092932182 UCIAMS-251302
#  6:      42540  0.0034971022 0.0127903203 UCIAMS-251302
# ---                        
# 55:      41070  0.0025092353 0.9934890544 UCIAMS-251302
# 56:      41040  0.0020485680 0.9955376225 UCIAMS-251302
# 57:      41010  0.0016670570 0.9972046794 UCIAMS-251302
# 58:      40980  0.0013637023 0.9985683817 UCIAMS-251302
# 59:      40950  0.0011083988 0.9996767805 UCIAMS-251302
# 60:      40920  0.0003232195 1.0000000000 UCIAMS-251302
#     CalYearsBP Probabilities      CumProb            ID
#          <num>         <num>        <num>        <char>
