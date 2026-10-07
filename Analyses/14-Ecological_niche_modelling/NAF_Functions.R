
##### Freymueller Maxent Pipeline -------------------------------------------------------------------------------------

# Custom functions made by Nick Freymueller for his M.S. work
# this pipeline has been generalized to work efficiently for paleo data - especially if . . .
# . . . running multiple species at a time

##### Packages -------------------------------------------------------------------------------------------------------------

# loading in the requisite packages

packs <- c('ade4', 'adehabitatMA', 'adehabitatHR', 'alphahull', 'dismo', 'dplyr', 'ecospat', 'ggplot2', 'jsonlite',
           'kuenm', 'matrixStats', 'predicts', 'raster', 'rgdal', 'rgeos', 'rJava',
           'sf', 'sp', 'splitstackshape', 'tidyr', 'utils', 'wesanderson', 'PerformanceAnalytics', 'SDMTools')

packs <- c('data.table', 'DescTools', 'dplyr', 'enmSdmX', 'ggplot2', 'matrixStats', 'mgcv', 'predicts', 'rJava', 'sf', 'splitstackshape')

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



# packs <- c('alphahull', 'data.table', 'DescTools', 'dplyr', 'ggplot2',  'matrixStats', 'predicts', 
#             'rJava', 'sf', 'splitstackshape', 'tidyr', 'wesanderson', 'PerformanceAnalytics', 'SDMTools')


# optionally writing the package versions as a csv file
# write.csv(packs.df 'package.versions.csv')

# finding where the dismo package is in your computer
# you also need to put a copy of the maxent.jar file there

# test if the maxent.jar file is in the dismo folder
list.files( system.file('java', package = 'dismo') )
# [1] "dismo.jar"  "maxent.jar"


# test if the maxent.jar file is in the predicts package folder
list.files( system.file('java', package = 'predicts'), pattern='jar' )
# [1] "dismo.jar"  "maxent.jar"

##### Spatial Thinning --------------------------------------------------------------

### function that thins a species occurrences to only one occurrence per grid cell per time bin
## can be used in scenarios where a package like spThin may take a long time, or have high variance over multiple iterations

### ARGUMENTS
## taxa_df <- a data frame of species occurrence points
##            this requires at least an x- and y-column, but can also include a "species ID" and a "time bin" column
## rast <- a terra SpatRast with the resolution you want to thin to
## x_col <- column name of the x-coordinate
## y_col <- column name of the y-coordinate
## time <- optional column name of time bins. Only needed if thinning over multiple time bins
## taxa_names <- optional column name of the taxa of interest. Only needed if thinning >1 taxa at a time

thin_one_occ_per_cell <- function(taxa_df, rast, x_col='x', y_col='y', times=NULL, taxa_names=NULL){
  require(data.table)
  require(DescTools)
  require(terra)
  taxa_df <- as.data.table(taxa_df)
  start_occs <- nrow(taxa_df)
  # making vectors of x and y coordinates to snap the occs to
  res_x <- res(rast)[1]
  res_y <- res(rast)[2]
  x_centers <- seq(from=as.numeric(ext(rast)[1]) + (res_x/2),
                   to=as.numeric(ext(rast)[2]) - (res_x/2),
                   by=res_x)
  y_centers <- seq(from=as.numeric(ext(rast)[3]) + (res_y/2),
                   to=as.numeric(ext(rast)[4]) - (res_y/2),
                   by=res_y)
  # making a modified data.frame of occs to be filtered
  if(is.null(times)){
    times <- rep(1, start_occs)
  } else {
    times <- taxa_df[,get(times)]
  }
  if(is.null(taxa_names)){
    taxa_names <- rep(1, start_occs)
  } else {
    taxa_names <- taxa_df[,get(taxa_names)]
  }
  occs1 <- data.frame(x_col=taxa_df[,get(x_col)],
                      y_col=taxa_df[,get(y_col)],
                      times,
                      taxa_names)
  occs1$x_col <- sapply(occs1$x_col, function(x) DescTools::Closest(x_centers, x) )
  occs1$y_col <- sapply(occs1$y_col, function(x) DescTools::Closest(y_centers, x) )
  occs1 <- unique(occs1)
  taxa_df <- taxa_df[as.numeric(row.names(occs1)),]
  final_occs <- nrow(taxa_df)
  row.names(taxa_df) <- 1:final_occs
  writeLines(c(paste0('Unthinned occs: ', start_occs), paste0('Thinned occs: ', final_occs)) )
  return(taxa_df)
}




##### make k-folds --------------------------------------------------------------------------------------

### internal function that divides the data up into k-folds for cross validation
##  the user won't actually interact with this function personally, it's only needed for make_occ_df

### ARGUMENTS
## x <- a data.frame of occurrences from a single taxon in a single area and time bin
## k <- number of folds for cross-validation
## seed <- a seed to be set
make_kfolds <- function(x, k=5, ltrs=letters, seed=NULL){
  if(class(x)[1] != "data.table"){
    class(x) <- "data.table"
  }
  n_occs <- nrow(x) # number of occurrences
  rep_times <- n_occs %/% k   # number of full reps for the rep function
  if(rep_times == 0){
    idx <- 1:n_occs
    folds <- rep('Train', n_occs)
  } else {
    k_bins <- rep_times*k   # number of occs minus any remainder when dividing by k
    remainder <- n_occs - k_bins   # find the remainder
    folds <- rep(ltrs, rep_times)   # make k bins of equal size
    if(!is.null(seed)){
      set.seed(seed)
    } else {
      seed <- Sys.Date() |> str_remove_all('-') |> as.numeric()
      set.seed(seed)
    }
    idx <- sample(1:n_occs, size=k_bins)
    folds <- sample(folds, size=k_bins)
  } # end else 
  
  output <- list(idx=idx, folds=folds)
  # remove seed 
  set.seed(NULL)
  return(output)
}

##### Make Master Files ---------------------------------------------------------

### functions that make occurrence and background data.frames from a provided raster stack

### make_occ_df extracts raster values from grid cells that have occurrences
##  this data frame will be in a format that will make it able to be merged with background data frames . . .
##  .  .  . generated with make_back_df
##  the taxa must already be filtered to a single time bin and Region, but multiple taxa can be ran at the same time

### ARGUMENTS
## rast_stack <- a SpatRaster stack with each layer representing a variable to be used in the modeling
## taxa_df <- a data frame of species occurrences - must have (in this order): 1) species name; 2) x-coordinate; 3) y-coordinate; and 4) ID
##            function will use the first 4 columns for making the data frame 
## x_col <- column name of the x-coordinate in the species data
## y_col <- column name of the y-coordinate in the species data
## ID <- ID of the species. can also be used to average a species' 
## timeBin  <- name of the time period the data is from - to be used for filtering later
## Region <- name of the region the data is from - helpful for filtering later
## probability <- the radiocarbon probability for an occurrence spanning multiple date time bins. if not specified, probability is set to 1
## prob_seed <- the seed to be set with the radiocarbon probability. if not specified, value is set to the current date: Sys.Date() |> str_remove_all('-') |> as.numeric()
## k_folds <- number of folds for cross-validation - to be passed to "make_kfolds" internally
## k_seed <- an optional seed to be sent to "make_kfolds". if not specified, value is set to the current date: Sys.Date() |> str_remove_all('-') |> as.numeric()
## snap_dist <- maximum distance (in km) to snap points in NA-grid cells to the center of the nearest grid cell
## epsg_in <- native crs of the input data (default is longitude/latitude in WGS1984)
## epsg_out <- crs used to project occ. points that are in NA-grid cells 
make_occ_df <- function(rast_stack, taxa_df, x_col = 'x', y_col = 'y', 
                        timeBin = 'time1', extent = 1, 
                        probability = 1, prob_seed = NULL,
                        k_folds = 5, k_seed = NULL, 
                        snap_dist = 50, epsg_in = 4326, epsg_out = NULL ){
  # requiring packages
  require(dplyr)
  require(sf)
  require(terra)
  require(stringr)
  # setting seeds
  seed <- Sys.Date() |> str_remove_all('-') |> as.numeric()
  if(is.null(prob_seed)){
    prob_seed <- seed
  }
  if(is.null(k_seed)){
    k_seed <- seed
  }
  # if probability < 1
  # adding in probability information
  if(probability == 1){
    prob <- probability
  } else if(is.null(probability)){
    prob <- 1
  } else if(probability > 1 | probability < 0){
    stop('`probability` must be between [0,1] ')
  } else {
    prob <- probability
  }
  n_occs <- nrow(taxa_df)
  prob <- data.frame(prob = rep(prob, n_occs))
  # if radiocarbon probabilities are set, forcing k_folds to be 1
  if(prob[1,] != 1){
    k_folds <- 1 
  }
  # ensuring that k_folds is a positive integer between [1,10]
  if( !is.null(k_folds) ){
    if( k_folds < 1){
      k_folds <- 1
      print('Making occurrence data.frame without any k_folds.')
    } else if( k_folds > 10 ){
      k_folds <- 10
      warning('k_folds has been set to a maximum of 10')
    } else if( k_folds%%1 != 0 ){
      warning('k_folds has been rounded to the nearest whole number.')
      k_folds <- round(k_folds,0)
    }
  } else {
    k_folds <- 1
    #print('Making occurrence data.frame without any k_folds.')
  } # end k_folds
  # filtering taxa_df to four columns
  taxa_df <- taxa_df[,1:4]
  # extracting the env-variables to a data.frame
  occs_df <- terra::extract(x=rast_stack, y=taxa_df[,2:3])
  occs_df <- cbind(taxa_df, occs_df[,-1])
  colnames(occs_df)[2:3] <- c(x_col, y_col)
  # checking to see if not all points were matched
  if( nrow(occs_df) != nrow(occs_df[complete.cases(occs_df), ]) ){
    if(is.null(epsg_out)) stop('Points exist in NA grid cells, but epsg_out not provided ')
    # identifying points in need of matching
    n_snap <- nrow(occs_df) - nrow(occs_df[complete.cases(occs_df), ])
    warning(paste0('Snapping ', n_snap, ' points to nearest grid cell based on epsg_out and snap_dist'))
    # extracting and projecting raster cords
    r_coords <- crds(rast_stack, df=TRUE)
    r_coords_sf <- st_as_sf(r_coords, coords=1:2, crs=epsg_in)
    r_coords_sf <- st_transform(r_coords_sf, crs=epsg_out)
    # subsetting non-matched points and re-projecting them
    to_snap <- occs_df[!complete.cases(occs_df), 1:3]
    to_snap_sf <- st_as_sf(to_snap, coords=2:3, crs=epsg_in)
    to_snap_sf <- st_transform(to_snap_sf, crs=epsg_out)
    # finding the nearest grid cell
    near_cell <- st_nearest_feature(to_snap_sf, r_coords_sf)
    r_coords_sf <- r_coords_sf[near_cell,]
    # distances
    dists <- round(as.numeric(st_distance(to_snap_sf, r_coords_sf, by_element = TRUE)/1000), 1)
    to_snap$dist <- dists
    # records to not snap and extract
    to_not_snap <- to_snap[which(to_snap$dist > snap_dist),]
    if(nrow(to_not_snap) > 0) {
      warning(nrow(to_not_snap), ' points were beyond snap_dist and were not matched\nUse row.names(taxa_df) to get indices of matched points')
    }
    # records to snap and extract
    to_snap <- to_snap[which(to_snap$dist <= snap_dist),]
    snap_names <- as.numeric(row.names(to_snap))
    snap_cells <- near_cell[which(to_snap$dist <= snap_dist)]
    # records already snapped and extracted
    id <- as.numeric(row.names(occs_df[complete.cases(occs_df), ]))
    occs_df <- occs_df[complete.cases(occs_df),]
    occs_df$ID <- id
    # extracting records to extract
    to_match <- to_snap
    to_match[,2:3] <- r_coords[snap_cells,]
    to_match <- terra::extract(x=rast_stack, y=to_match[,2:3])[,-1]
    # assigning snapped records to to_snap
    snapped <- cbind(to_snap, to_match)
    snapped$dist <- NULL
    snapped$ID <- snap_names
    # merging snapped back with occs_df
    occs_df <- suppressMessages(full_join(occs_df, snapped))
    occs_df <- occs_df[order(occs_df$ID),]
    row_out <- 1:nrow(occs_df)
    occs_df$ID <- NULL
  } else {row_out <- 1:nrow(occs_df)
  }# closing if-there-were-occs-in-NA-grid-cells
  
  colnames(occs_df)[1] <- 'taxon'
  n_occs <- nrow(occs_df)
  # making the data.frame of the name, TimeBin, extent, and presence columns
  d_cols <- as.data.frame( matrix(1, nrow=n_occs, ncol=3) )
  d_cols[,1] <- timeBin
  d_cols[,2] <- extent
  colnames(d_cols) <- c('timeBin', 'extent', 'presence')
  
  # ID information
  ID <- taxa_df[,4]
  
  # merging everything together
  occs_df <- do.call('cbind', list(occs_df[,1:3], ID, prob, d_cols[,1:2], occs_df[,5:ncol(occs_df)], d_cols[,3] ) )
  colnames(occs_df)[ ncol(occs_df) ] <- 'presence'
  # assigning default k_fold identity - will be altered later if k_folds > 1
  occs_df$fold <- 'train'
  
  # making the cross-validation columns
  if( k_folds > 1 ){
    ltrs <- letters[1:k_folds]
    
    # splitting up the occurrence data.frame by taxon, then 
    occs_list <- base::split(occs_df, occs_df[,1])
    n.taxa <- length(occs_list)
    # assigning k_folds
    # if 
    for(i in 1:n.taxa){
      occ <- occs_list[[i]]
      # relating number of occurrences to k_folds
      mfolds <- make_kfolds(x = occ, k = k_folds, ltrs = LETTERS[1:k_folds])
      occ$fold[mfolds$idx] <- paste0('test_', mfolds$folds)
      # writing back to the list
      occs_list[[i]] <- occ
      
    } # close i loop
    
    # merging all the different taxa back together
    occs_df <- do.call('rbind', occs_list)
    row.names(occs_df) <- row_out
    
  } # closing if( k_folds > 1 )
  
  # output
  return(occs_df)
}



### make_back_df extracts all non-NA grid cells into a data frame
##  this data frame will be in a format that will make it able to be merged with occurrence data frames . . .
##  .  .  . generated with make_occ_df

### ARGUMENTS
## rast_stack <- a SpatRaster stack (using the terra package)
## mask_rast <- a SpatRaster masked to areas you want to extract from - function will set values to 1
## max_cells <- the maximum number of cells you want to extract from - will always select the minimum of c(non-NA cells in mask_rast, max_cells)
## cell_seed <- the seed used to select cells. if not specified, value is set to the current date: Sys.Date() |> str_remove_all('-') |> as.numeric()
## taxa_df <- a data frame of species occurrences - must have (in this order): 1) species name; 2) x-coordinate; 3) y-coordinate; and 4) ID
##            function will use the first 4 columns for making the data frame 
## x_col <- column name of the x-coordinate in the species data
## y_col <- column name of the y-coordinate in the species data
## ID <- ID of the species. can also be used to average a species' 
## timeBin  <- name of the time period the data is from - to be used for filtering later
## extent <- name of the region the data is from - helpful for filtering later

## the function must be ran separately for every combination of extent/time bin you wish to extract for
## if extracting from occurrences that lie entirely within a single time bin (i.e., not-incorporating radiocarbon calibration), . . .
#  . . . set ID equal to a vector of IDs in a data frame
## if extracting from occurrences where the radiocarbon uncertainty spans multiple time bins
make_back_df <- function(rast_stack, mask_rast = NULL, max_cells = NULL, cell_seed = NULL,
                         taxa_df, x_col = 'x', y_col = 'y', ID = 1,
                         timeBin = 1, extent = 1, 
                         probability = 1, prob_seed = NULL){
  # requiring packages
  require(dplyr)
  require(sf)
  require(terra)
  require(stringr)
  # optional masking 
  if(!is.null(mask_rast)){
    mask_rast[!is.na(mask_rast)] <- 1
    rast_stack <- rast_stack*mask_rast
  }
  # number of grid cells - takes the first number of non-NA cells from the first layer of rast_stack
  n_cells <- global(rast_stack, fun='notNA')[1,]
  # if probability < 1
  # adding in probability information
  if(probability == 1){
    prob <- probability
  } else if(is.null(probability)){
    prob <- 1
  } else if(probability > 1 | probability < 0){
    stop('`probability` must be between [0,1] ')
  } else {
    prob <- probability
  }
  prob <- data.frame(prob = rep(prob, n_cells))
  
  # extracting the env-variables to a data.frame
  coords_df <- as.data.frame(rast_stack, xy=TRUE, na.rm=TRUE) # colnumbers c(2:3, 8:10)
  colnames(coords_df)[1:2] <- c(x_col, y_col)
  n_back <- nrow(coords_df)
  # optional if filtering to a certain number of points per occurrence
  if(!is.null(max_cells)){
    if(is.numeric(max_cells) & max_cells > 0){
      max_cells <- round(max_cells, 0)
      max_cells <- min(max_cells, n_back)
    } else {
      stop('`max_cells` must be a positive number.')
    }
    # seed for subsetting background points
    if(!is.null(cell_seed)){
      seed <- cell_seed
    } else {
      seed <- Sys.Date() |> str_remove_all('-') |> as.numeric()
    }
    # filtering background points
    set.seed(seed)
    idx <- sample(1:n_back, max_cells)
    set.seed(NULL)
    coords_df <- coords_df[idx,]
    n_back <- nrow(coords_df)
  }
  # assigning colnames
  colnames(coords_df)[1:2] <- c(x_col, y_col)
  
  # making the data.frame of the name, timeBin, extent, and presence columns
  d_cols <- as.data.frame( matrix(0, nrow=n_back, ncol=4) ) # c(1, 6:7, 11)
  d_cols[,1] <- 'Background'
  d_cols[,2] <- timeBin
  d_cols[,3] <- extent
  colnames(d_cols) <- c('taxon', 'timeBin', 'extent', 'presence')
  d_cols$ID <- rep(ID, n_back) # adding in c(4:5), now c(1, 6:7, 11, 4:5)
  d_cols$fold <- rep(NA, n_back)
  
  # merging everything together
  out <- do.call('cbind', list( d_cols[,1], coords_df[,1:2], d_cols[,5], prob, d_cols[,2:3], coords_df[,3:ncol(coords_df)], d_cols[,c(4,6)] ))
  colnames(out)[1] <- 'taxon'
  colnames(out)[4] <- 'ID'
  colnames(out)[ncol(out)] <- 'fold'
  # output
  return(out)
}


##### aggregate_14C ------------------------------------------------------------------------------------------------------

### maxent and just about every other ENM algorithm can only accept 1/0 data for presence vs background/absences
## this makes fossil occurrences that span multiple time bins (e.g., 30-year climatic normals)  practically challenging to apply in ENMs
## for maxent, there are two ways to properly account for this:
#  1.) bootstrap sample the fossil pseudo-replicates (occurrence + associated background points) relative to their radiocarbon probabilities, run the model many many times, then average
#      given the radiocarbon uncertainty of some fossils (particularly older ones) can span a few thousand years, it would take a LOT of replicates to integrate over this space
#      this is doubly prohibitive because each of these bootstrap replicates needs to be ran for each hyperparameter setting (detault is 30)
#      naturally, this would take a LONG time
#  2.) calculate weighted averages of all pseudo replicates (both for the occurrence and associated background points), then run one model per hyperparameter setting
#      for each 


##### fossil_crossVal ------------------------------------------------------------------------------------------------------

### assigning cross-validation fold information to fossil occurrences
##  to be used for fossil data that spans multiple radiocarbon ages - which can't be . . .
#   . . . assigned cross-validation reps during the intersection with climate data
##  this function re-makes the `fold` column to have cross-validation reps AFTER the data has been processed via `aggregate_14C`
##  function temporally stratifies the data based on `cut_points`
#   e.g., if k_folds == 5, will return 
#   defaults are the boundaries between Early Holocene/Middle Holocene/Late Holocene/Terminal Pleistocene/Last Glacial Maximum
#   can check using cut(x=occs$timeBin, breaks = cut_points)
#   must specify breaks that completely 

### ARGUMENTS
## occs        <- species occurrence data frame generated from `make_occ_df`
## cut_points  <- vector of ages that splits
## k_folds     <- number of folds for cross-validation - to be passed to "make_kfolds" internally
# k_seed       <- occurrence data frame generated from `make_occ_df`

fossil_crossVal <- function(occs, 
                            cut_points = c(0, 4200, 8200, 11700, 21000, 60000), 
                            k_folds = 5, k_seed = NULL){
  cut_points <- sort(cut_points)
  if(min(cut_points) > min(occs$timeBin)){
    cut_points[1] <- min(occs$timeBin)
    warning(paste0('Setting the lowest `cut_points` value to ', min(occs$timeBin), '.' ))
  }
  if(max(cut_points) < max(occs$timeBin) ) {
    cut_points[length(cut_points)] <- max(occs$timeBin)
    warning(paste0('Setting the highest `cut_points` value to ', max(occs$timeBin), '.' ))
  }
  
  # ensuring that k_folds is a positive integer between [1,10]
  if( !is.null(k_folds) ){
    if( k_folds < 1){
      k_folds <- 1
      print('Making occurrence data.frame without any k_folds.')
    } else if( k_folds > 10 ){
      k_folds <- 10
      warning('k_folds has been set to a maximum of 10')
    } else if( k_folds%%1 != 0 ){
      warning('k_folds has been rounded to the nearest whole number.')
      k_folds <- round(k_folds,0)
    }
  } 
  # making a data.frame and lsit to manipulate
  cut_df <- data.frame(idx = 1:nrow(occs), 
                       bin = cut(x=occs$timeBin, breaks = cut_points)
  )
  cut_list <- split(cut_df, cut_df$bin)
  for(i in 1:k_folds){
    occ <- cut_list[[i]]
    occ$fold <- 'train'
    mfolds <- make_kfolds(x = occ, k = k_folds, ltrs = LETTERS[1:k_folds])
    occ$fold[mfolds$idx] <- paste0('test_', mfolds$folds)
    cut_list[[i]] <- occ
  }
  cut_df <- do.call('rbind', cut_list) |> dplyr::arrange(idx)
  cut_df <- cut_df[sort(cut_df$idx),]
  
  occs$fold <- cut_df$fold
  return(occs)
  
  
} # end fossil_crossVal

##### Prep Parameters ------------------------------------------------------------------------------------------------------

## internal function that preps parameters for Maxent models in R with the dismo package
# the user won't interact with this function personally. it is used in optimize_maxent_likelihood and maxent_crossval_error

# this function was made by Xiao Feng
# sourced from https://github.com/shandongfx/workshop_maxent_R/blob/master/code/Appendix2_prepPara.R
#  https://github.com/shandongfx/workshop_maxent_R/blob/master/code

# should be used in concert with Appendix 3 from Feng et al. 2017 (PeerJ Preprints)
# https://doi.org/10.7287/peerj.preprints.3346v1
# manuscript still unpublished as of June 2026
# Appendix 3
# https://github.com/shandongfx/workshop_maxent_R/blob/master/code/Appendix3_maxentParameters_v2.pdf

# some arguments may change, or may/may not be used depending on if you're using raster data vs "samples-with-data" (SWD; column data) 

# A function that implements Maxent parameters using the general R manner
# leave "doclamp" as default - later in the code (internal fxns run), "doclamp" is set to FALSE
prepPara <- function(userfeatures=NULL,             # 41  NULL=autofeature, could be any combination of # c("L", "Q", "H", "P", "T")
                     #     MUST be specified as a single string (e.g., "LQ", "LQP", "LQHPT", etc.)
                     responsecurves=TRUE,           # 1
                     jackknife=TRUE,                # 3
                     outputformat="logistic",       # 4
                     outputfiletype="asc",          # 5
                     projectionlayers=NULL,         # 7
                     randomseed=FALSE,              # 10
                     removeduplicates=TRUE,         # 16
                     betamultiplier=NULL,           # 20, 53-56
                     biasfile=NULL,                 # 22
                     testsamplesfile=NULL,          # 23
                     replicates=1,                  # 24-25
                     replicatetype="crossvalidate", # 24-25
                     writeplotdata=TRUE,            # 37
                     extrapolate=TRUE,              # 39
                     doclamp=TRUE,                  # 42
                     beta_threshold=NULL,           # 20, 53-56
                     beta_categorical=NULL,         # 20, 53-56
                     beta_lqp=NULL,                 # 20, 53-56
                     beta_hinge=NULL,               # 20, 53-56
                     applythresholdrule=NULL        # 60
){
  #20, 29-33, & 41 features, default is autofeature
  if(is.null(userfeatures)){
    args_out <- c("autofeature")
  } else {
    args_out <- c("noautofeature")
    if(grepl("L",userfeatures)) args_out <- c(args_out,"linear") else args_out <- c(args_out,"nolinear")
    if(grepl("Q",userfeatures)) args_out <- c(args_out,"quadratic") else args_out <- c(args_out,"noquadratic")
    if(grepl("H",userfeatures)) args_out <- c(args_out,"hinge") else args_out <- c(args_out,"nohinge")
    if(grepl("P",userfeatures)) args_out <- c(args_out,"product") else args_out <- c(args_out,"noproduct")
    if(grepl("T",userfeatures)) args_out <- c(args_out,"threshold") else args_out <- c(args_out,"nothreshold")
  }
  
  # 1 - generate response curves for each variable
  if(responsecurves) args_out <- c(args_out,"responsecurves") else args_out <- c(args_out,"noresponsecurves")
  
  # 2
  #if(picture) args_out <- c(args_out,"pictures") else args_out <- c(args_out,"nopictures")
  
  # 3 - apply variable jackknife to see how the model changes if that variable is omitted, then if it's the ONLY variable used
  if(jackknife) args_out <- c(args_out,"jackknife") else args_out <- c(args_out,"nojackknife")
  
  # 4 - output format type. choose from c("logistic", "cumulative", "raw")
  args_out <- c(args_out,paste0("outputformat=",outputformat))
  
  # 5 - output file type. choose from c("asc", "mxe", "grd", "bil")
  args_out <- c(args_out,paste0("outputfiletype=",outputfiletype))
  
  # 7 - pathway to projection layers.
  # it seems that the projection layers should be the only files in that folder, just like with the MaxEnt .jar file
  if(!is.null(projectionlayers))    args_out <- c(args_out,paste0("projectionlayers=",projectionlayers))
  
  # 10 - will use different random number generators for selecting training vs testing data and background points (if applicable)
  if(randomseed) args_out <- c(args_out,"randomseed") else args_out <- c(args_out,"norandomseed")
  
  # 16 - remove duplicate coordinates that are in the same grid  - ONLY for raster data, not SWD
  if(removeduplicates) args_out <- c(args_out,"removeduplicates") else args_out <- c(args_out,"noremoveduplicates")
  
  # 20 & 53-56 - various beta (regularization) multipliers to be applied. default = 1
  # 20 applies all parameters by this regularization multiplier.
  # 53-56 can apply uniquely to different feature types
  # check if negative
  betas <- c( betamultiplier,beta_threshold,beta_categorical,beta_lqp,beta_hinge)
  if(! is.null(betas) ){
    for(i in 1:length(betas)){
      if(betas[i] <0) stop("betamultiplier has to be positive")
    }
  }
  if (  !is.null(betamultiplier)  ){
    args_out <- c(args_out,paste0("betamultiplier=",betamultiplier))
  } else {
    if(!is.null(beta_threshold)) args_out <- c(args_out,paste0("beta_threshold=",beta_threshold))
    if(!is.null(beta_categorical)) args_out <- c(args_out,paste0("beta_categorical=",beta_categorical))
    if(!is.null(beta_lqp)) args_out <- c(args_out,paste0("beta_lqp=",beta_lqp))
    if(!is.null(beta_hinge)) args_out <- c(args_out,paste0("beta_hinge=",beta_hinge))
  }
  
  # 22 - pathway to a bias file for selecting background points - ONLY for raster data, not SWD
  if(!is.null(biasfile))    args_out <- c(args_out,paste0("biasfile=",biasfile))
  
  # 23 - pathway to a test data file - can be in csv format
  if(!is.null(testsamplesfile))    args_out <- c(args_out,paste0("testsamplesfile=",testsamplesfile))
  
  # 24 - replicates = number of replicates to run (integer)
  # 25 - replicatetype = what type of replicates to run. choose from c('crossvalidate', 'bootstrap', 'subsample')
  replicates <- as.integer(replicates)
  if(replicates>1 ){
    args_out <- c(args_out,
                  paste0("replicates=",replicates),
                  paste0("replicatetype=",replicatetype) )
  }
  
  # 37 - write output files containing the data used to make response curves
  if(writeplotdata) args_out <- c(args_out,"writeplotdata") else args_out <- c(args_out,"nowriteplotdata")
  
  # 39 - allow extrapolation beyond the limits of the training data
  if(extrapolate) args_out <- c(args_out,"extrapolate") else args_out <- c(args_out,"noextrapolate")
  
  # 42 - apply clamping when projecting
  if(doclamp) args_out <- c(args_out,"doclamp") else args_out <- c(args_out,"nodoclamp")
  
  # 60 - threshold your model to binary 1/0
  # options are: c('Fixed cumulative value 1', 'Fixed cumulative value 5', 'Fixed cumulative value 10', 'Minimum training presence',
  # '10 percentile training presence', 'Equal training sensitivity and specificity', 'Maximum training sensitivity plus specificity').
  if(!is.null(applythresholdrule))    args_out <- c(args_out,paste0("applythresholdrule=",applythresholdrule))
  
  return(args_out)
}

# prepPara()
# [1] "autofeature"           "responsecurves"        "jackknife"             "outputformat=logistic" "outputfiletype=asc"   
# [6] "norandomseed"          "removeduplicates"      "writeplotdata"         "extrapolate"           "doclamp"

##### Create Null Object, Summary Object, Eval Object, and Nested File Structure -------------------------------------------

### functions that create the null data.frame, summary data.frame, and the eval data.frame

### ARGUMENTS 1 - arguments that exist to keep track of everything and do not change how the functions run
## taxon_name <- name of the entity that will get assigned to any null/summary/eval objects
## time_bin <-   name of the time bin that will get assigned to any null/summary/eval objects
## extent <-     name of the extent that will get assigned to any null/summary/eval objects

### ARGUMENTS 2 - arguments that set the model hyper-parameters and will change how the functions run
## cv_runs <-    name(s) of the cross-validation folds. folders will be generated with these names in make_maxent_folders.
#                recommended framework is to treat each cross-validation fold as a letter.
#                In the case of 5-fold cross validation, for running all the data (using the null_aic and optimize_maxent_likelihood . . .
#                . . . functions), you should specify 'abcde'. For running cross validation models with maxent_crossval_error, you . . .
#                . . . should specify c('abcd', 'abce', 'abde', 'acde', 'bcde').
## f_class <-    name(s) of the feature classes used in analysis.
## beta_values <- regularization multipliers used in analysis

# create_summary_df and create_eval_df will create every possible combination of cv_runs, f_class, and beta_values

# function to create the null data.frame

create_null_df <- function(taxon_name, timeBin, extent){
  null_object <- as.data.frame(matrix(NA, nrow=1, ncol=13))
  colnames(null_object) <- c('taxon', 'timeBin', 'extent', 'fold', 'features', 'betas', 'n', 'k', 'ln_L', 'AIC', 'AICc', 'Boyce', 'lambdas')
  null_object$fold <- 'train'
  null_object$taxon <- taxon_name
  null_object$timeBin <- timeBin
  null_object$extent <- extent
  null_object$lambdas <- 'Intercept'
  
  return(null_object)
}



create_summary_df <- function(taxon_name,
                              timeBin,
                              extent,
                              fold = 'train',
                              f_class = c('LQP', 'Q'),
                              beta_values = c(0.025, 0.05, 0.1, 0.25, 0.5, 0.75, 1, 1.5, 2, 2.5, 3, 3.5, 4, 4.5, 5) ){
  
  summary_object <- expand.grid(fold = fold, features = f_class, betas = beta_values, stringsAsFactors = T)
  
  # ifelse functions defining what to do if taxon_name, TimeBin, and extent are not specified
  if( !is.null(taxon_name) ){
    summary_object$taxon <- taxon_name
  } else {
    summary_object$taxon <- 'Taxon1'
  }
  if( !is.null(timeBin) ){
    summary_object$timeBin <- timeBin
  } else {
    summary_object$timeBin <- 'TimeBin1'
  }
  if( !is.null(extent) ){
    summary_object$extent <- extent
  } else {
    summary_object$extent <- 'Extent1'
  }
  
  # re-ordering columns
  summary_object <- summary_object[,c(4:6, 1:3)]
  # adding in all the other parameters
  summary_object$n <- NA   # sample size
  summary_object$k <- NA   # number of non-zero lambdas
  summary_object$ln_L <- NA   # log-likelihood
  summary_object$AIC <- NA   # AIC
  summary_object$AICc <- NA   # AICc corrected for small sample size
  summary_object$Boyce <- NA
  summary_object$lambdas  <- NA   # list of all the non-zero lambdas
  return(summary_object)
}


# function to create the eval data.frame - to be used as an input for maxent_crossval_error
create_eval_df <- function(taxon_name,
                           timeBin,
                           extent,
                           fold = LETTERS[1:5],
                           f_class = 'LQP',
                           beta_values = 1 ){
  # making testfolds
  # tf <- paste0('test_', fold)
  fold <- as.character(fold)
  summary_object <- expand.grid(fold=paste0('test_',fold), features = f_class, betas = beta_values, stringsAsFactors = T)
  
  # ifelse functions defining what to do if taxon_name, TimeBin, and extent are not specified
  if( !is.null(taxon_name) ){
    summary_object$taxon <- taxon_name
  } else {
    summary_object$taxon <- 'Taxon1'
  }
  if( !is.null(timeBin) ){
    summary_object$timeBin <- timeBin
  } else {
    summary_object$timeBin <- 'TimeBin1'
  }
  if( !is.null(extent) ){
    summary_object$extent <- extent
  } else {
    summary_object$extent <- 'Extent1'
  }
  
  # re-ordering columns
  summary_object <- summary_object[,c(4:6, 1:3)]
  # adding in all the other parameters
  summary_object$n <- NA
  summary_object$thresh <- NA
  summary_object$testSens <- NA
  summary_object$testBoyce <- NA
  summary_object$lambdas  <- NA   # list of all the non-zero lambdas
  
  return(summary_object)
}


# function to create the nested file structure that it will use to store the output (including html files) of maxent models 
# summary_eval_df <- the summary_df or the eval_df data.frame.
# the function will use information from the cv_runs, f_class, and beta_values columns to create this nested file structure
# wd <- working directory that the nested file structure is going to be generated in. if running multiple species, it is recommended . . .
# . . . that you make a folder for each species, then run this function in each species folder

make_maxent_folders <- function(summary_eval_df, wd = getwd() ){
  
  # setting the working directory
  setwd( getwd() )
  
  # for loop that goes through the summary/evaluation data.frame and creates the nested file structure
  for( i in 1:nrow(summary_eval_df) ){
    
    # assigning the short names for folders
    fold <- as.character(summary_eval_df$fold[i])
    # if(fold != 'train'){
    #   fold <- paste0('test_', fold)
    # }
    features <- summary_eval_df$features[i]
    betas <- summary_eval_df$betas[i]
    
    ## writing the nested file structure
    # cross-validation folder - either train or test_A:test_E 
    if( !dir.exists( paste(getwd(), fold, sep='/' ) ) ){ 
      writeLines( c('Creating folder:', paste(getwd(), fold, sep='/'), sep='') )
      dir.create( paste(getwd(), fold, sep='/') )
      setwd( paste(getwd(), fold, sep='/') )
    } else {
      setwd( paste(getwd(), fold, sep='/') )
    }
    # feature classes folder
    if( !dir.exists( paste(getwd(), features, sep='/' ) ) ){ 
      writeLines( c('Creating folder:', paste(getwd(), features, sep='/'), sep='') )
      dir.create( paste(getwd(), features, sep='/') )
      setwd( paste(getwd(), features, sep='/') )
    } else {
      setwd( paste(getwd(), features, sep='/') )
    }
    # regularization penalty folder
    if( !dir.exists( paste(getwd(), betas, sep='/' ) ) ){ 
      writeLines( c('Creating folder:', paste(getwd(), paste('beta',  betas, sep='_'),  sep='/'), sep='') )
      dir.create( paste(getwd(), paste('beta',  betas, sep='_'), sep='/') )
      setwd( paste(getwd(), paste('beta',  betas, sep='_'), sep='/') )
    } else {
      setwd( paste(getwd(), paste('beta',  betas, sep='_'), sep='/') )
    }
    
    setwd('../../../')   # go three directories back
    
  } # end i loop
  
}

# setwd("~/Dropbox/SVP_Models/ModelOutput/Tyrano")
# 
# 
# examp <- create_summary_df('Tyrano', 'K', 'Laurimidia')
# make_maxent_folders(examp)
# 
# examp2 <- create_eval_df('Tyrano', 'K', 'Laurimidia')
# make_maxent_folders(examp2)




##### Null AICc ------------------------------------------------------------------------------------------------------------

## function that calculates AIC and AICc values for a null (i.e., intercept-only) model
#  to be compared with the "best" model generated via optimize_maxent_likelihood

### ARGUMENTS
# the arguments are the same as the optimize_maxent_likelihood model, except that null_df should be a data.frame created from the . . .
# . . . create_null_df object
# function will return 
null_aic <- function(null_df, occs, background, first_occ_col){
  require(dplyr)
  
  # number of parameters
  k <- 1
  null_df$k <- k
  
  # # assigning the column number to be sent through maxent
  # col_number <- first_occ_col
  # # filtering the species dataset by col_number and assigning to null_df
  # sp <- occs |> 
  #   dplyr::filter(presence == 1)
  n <- nrow(occs)
  null_df$n <- nrow(occs)
  
  # giving noting the variables/features/hyperparameters with non-zero lambdas
  null_df$lambdas <- 'Intercept'
  
  out_scale <- rep(1, times=n ) / nrow(background)
  
  log_like <- sum(log(out_scale))
  
  null_df$ln_L <- log_like
  
  ### calculating AIC and AICc
  # AIC
  AIC <- 2 - 2*log_like
  null_df$AIC <- AIC
  # AICc for one term null model
  null_df$AICc <- AIC + (2*k^2 + 2*k)/(n - k - 1)
  
  # keep Boyce set at NA. the empirical value for an intercept only model has no standard deviation, and therefore is undefined
  
  return(null_df)
}


##### Optimize Maxent Likelihood -------------------------------------------------------------------------------------------

## function that runs a series of maxent models with varying parameter settings such that the user can . . .
#  .  .  . find one with  optimum settinggsoptimize their
# function will return a filled out summary model
optimize_maxent_likelihood <- function(summary_df,         # summary object that will keep all the model output
                                       # the nested file structure created from `make_maxent_folders` MUST exist
                                       occs,               # species occurrence data frame generated from `make_occ_df`
                                       background,         # background data generated from `make_back_df` (COLUMNS MUST BE IDENTICAL TO occs)
                                       presenceCol,        # column number of the presence column
                                       predicCols,         # column numbers of the predictor variables
                                       Boyce = 'cr',       # smooth boyce spline method (Liu et al. 2024 Ecography)
                                       # choose from c('tp', 'cr', 'bs', 'ps', 'ad', 'm')
                                       # corresponds to 'thin-plate', 
                                       # in practice, selecting 'm' takes longer than model building (~150 occs + 20,000 background)
                                       # all methods broadly agree with each other (SBI within ~0.01)
                                       # 'cr' = cubic-spline
                                       home=getwd(),       # directory where all the models will be ran. 
                                       all_models = TRUE   # do you want to keep all versions of all the models?
                                       # helpful for quickly checking some models, but may consume loads (e.g., >1GB) . . .
                                       # . . . of hard disk space. Recommended to set to FALSE for exploratory analyses.
                                       # if FALSE, function will create a folder called "RunOver" and will . . .
                                       # . . . continuously write-over it for all models, and the only model you see . . .
                                       # . . . at the end will be the last model that was ran
){
  # required packages
  require(mgcv)
  require(dplyr)
  require(splitstackshape)
  # checking boyce
  if(!Boyce %in% c('tp', 'cr', 'bs', 'ps', 'ad', 'm') | length(Boyce) != 1){
    stop('Set `Boyce` to one of: c(`tp`, `cr`, `bs`, `ps`, `ad`, `m`)')
  }
  
  # resetting data to data frames
  occs <- as.data.frame(occs)
  background <- as.data.frame(background)
  
  # Calculating the total number of models
  nmodels <- nrow(summary_df)
  
  # prompting the user if they want to store models in the RAM
  print(paste0('The time is ', round(Sys.time()), '. You are running ', nmodels, ' Maxent models.'))
  
  # setting up the progress bar
  prog <- txtProgressBar(min=0, max=nrow(summary_df), style=3,  char='+')
  
  for(i in 1:nrow(summary_df)){
    
    # setting the working directory for each folder
    if(all_models == TRUE){
      setwd( paste(home, summary_df$fold[i], summary_df$features[i], paste('beta', summary_df$betas[i], sep='_'), sep='/') )
    } else {
      # create the RunOver folder
      if( !dir.exists( paste(home, 'RunOver', sep='/') ) ){
        dir.create( paste(home, 'RunOver', sep='/') )
        setwd( paste(home, 'RunOver', sep='/') )
      } else {
        setwd( paste(home, 'RunOver', sep='/') )
      }
    }
    
    ###########################################################################
    
    ### MaxEnt things happen here
    
    ### preparing the data for the maxent model
    # # filtering the occ object by it's respective cross-validation identity
    # cv_number <- summary_df$cv_num[i]
    # # assigning the column number to be sent through maxent
    # col_number <- first_occ_col + cv_number - 1
    
    # filtering the species dataset by col_number and assigning to summary_df
    sp <- occs |> 
      dplyr::filter(occs[,presenceCol] == 1)
    
    
    n <- nrow(sp)
    summary_df$n[i] <- n
    # bending the occ and background data together
    mx_data <- rbind(sp, background)
    
    # arguments to be passed into maxent
    mx_args <- prepPara(userfeatures = summary_df$features[i], 
                        betamultiplier = summary_df$betas[i], 
                        doclamp = FALSE)
    
    ### running the actual maxent model
    mx_model <- predicts::MaxEnt(p = mx_data[,presenceCol],
                                 x = mx_data[,predicCols],
                                 path =  paste0(getwd()),
                                 args = mx_args
    )
    
    # ### running the actual maxent model
    # mx_model <- dismo::maxent(p = mx_data[,col_number],
    #                           x = mx_data[,predic],
    #                           path =  paste0(getwd()),
    #                           args = prepPara(userfeatures = summary_df[i,6], betamultiplier = summary_df[i,7], doclamp = FALSE)
    # )
    
    ### calculating k and the names of the lambdas
    lambda_file <- as.data.frame(mx_model@lambdas) |> `colnames<-`('lambdas')
    # lambdas data frame
    lambdas_df <- splitstackshape::cSplit(lambda_file, sep=',', splitCols='lambdas')
    colnames(lambdas_df) <- c('feature', 'lambda', 'min', 'max')
    # finding the non-zero lambdas
    non_zero_lambdas <- lambdas_df |> 
      dplyr::filter(!is.na(max)) |> 
      dplyr::filter(lambda != 0)
    class(non_zero_lambdas) <- 'data.frame'
    # assigning the number of parameters
    k <- nrow(non_zero_lambdas)
    # if beta is too high, the model gets over-regularized to the point that all the lambda coefficients get set to 0
    #  this effectively becomes an intercept only model, which is effectively the global mean
    # in this sense, k should get set to 0
    if(k == 0){
      k <- 1
    }
    summary_df$k[i] <- k
    # giving noting the variables/features/hyperparameters with non-zero lambdas
    summary_df$lambdas[i] <- toString(non_zero_lambdas[,1])
    
    ### calculating the log-likelihood, then AIC and AICc
    # if statement calculating if there is an appropriate AIC value
    # e.g., can't fit 4 observations (occurrence points) with 5 variables
    if(n - k < 2){
      summary_df$ln_L[i] <- NA
      summary_df$AIC[i] <- NA
      summary_df$AICc[i] <- NA
    } else {   # if  it is possible to calculate AIC and AICc
      # logistic model output
      # mx.back <- dismo::predict(mx_model, background[,predicCols])
      mx_back <- predicts::predict(mx_model, background[,predicCols])
      # sum of all background point values - mx.back / back_sum should = 1
      back_sum <- sum(mx_back)
      # logistic values of the (k-1)/k occurrences
      # mx_occs <- dismo::predict(mx_model, sp[,predicCols])
      mx_occs <- predicts::predict(mx_model, sp[,predicCols])
      # scaling to make compatible for calculating AIC
      occs_raw <- mx_occs / back_sum
      # log(likelihood)
      log_like  <- sum(log(occs_raw))
      summary_df$ln_L[i]<- log_like
      
      ### calculating AIC and AICc
      # AIC
      AIC <- 2*k - 2*log_like
      summary_df$AIC[i] <- AIC
      # AICc
      summary_df$AICc[i] <- AIC + (2*k^2 + 2*k)/(n - k - 1)
      # Smooth Boyce Index (Liu et al. 2024 Ecography)
      p <- c(mx_occs, mx_back)
      n1 <- length(mx_occs)
      n0 <- length(mx_back)
      prd <- seq(min(p), max(p), length = n0)
      oc <- c(rep(1, n1), rep(0, n0))
      ktry <- 10
      if(Boyce == 'm'){
        md_tp = mgcv::gam(oc ~ s(p, bs = "tp", k = min(ktry, length(unique(p) ) ) ), family = binomial)
        prd_tp = predict(md_tp, newdata = data.frame(p = prd), type = 'response')
        md_cr = mgcv::gam(oc ~ s(p, bs = "cr", k = min(ktry, length(unique(p) ) ) ), family = binomial)
        prd_cr = predict(md_cr, newdata = data.frame(p = prd), type = 'response')
        md_bs = mgcv::gam(oc ~ s(p, bs = "bs", k = min(ktry, length(unique(p) ) ) ), family = binomial)
        prd_bs = predict(md_bs, newdata = data.frame(p = prd), type = 'response')
        md_ps = mgcv::gam(oc ~ s(p, bs = "ps", k = min(ktry, length(unique(p) ) ) ), family = binomial)
        prd_ps = predict(md_ps, newdata = data.frame(p = prd), type = 'response')
        md_ad = mgcv::gam(oc ~ s(p, bs = "ad", k = min(ktry, length(unique(p) ) ) ), family = binomial)
        prd_ad = predict(md_ad, newdata = data.frame(p = prd), type = 'response')
        prd_m = (prd_tp + prd_cr + prd_bs + prd_ps + prd_ad)/5
        boyce <- cor(prd, prd_m, method = "spearman")
      } else {
        md <- mgcv::gam(oc ~ s(p, bs = Boyce, k = min(ktry, length(unique(p) ) ) ), family = binomial)
        prd_new <- predict(md, newdata = data.frame(p = prd), type = 'response')
        boyce <- cor(prd, prd_new, method = "spearman")
      }
      summary_df$Boyce[i] <- boyce
      
    }
    # updating the prograss bar for each run to get an idea of how long things will take
    setTxtProgressBar(prog, i)
    
    
    ###########################################################################
    # returning to the home directory
    setwd(home)
    
  } # close i loop
  
  return(summary_df)
}


# Time difference of 4.658302 mins (when calculating all spline versions of Smooth Boyce Index, then ensembling them together)
# Time difference of 49.71973 secs (when calculating the fastest Smooth Boyce Index - cubic spline)
# saves ~82% of time needed to run the model

##### MaxEnt Cross-Validation Error -----------------------------------------------------------------------------------------------

# function that runs five-fold cross-validation once you've found the optimum hyper-parameters for your maxent model(s)
# functions returns a list containing 3 objects:
# 1.) a filled out eval object
# 2.) an object containing all the maxent models; and . . .
# 3.) a data.frame with dimensions [ 1:nrow(background), 1:nrow(eval_df) ] containing the projections of all the models in
maxent_crossval_error <- function(eval_df,            # object generated from create_eval_df that will keep all the model output
                                  occs,               # species occurrence object
                                  background,         # sampled Background object (columns MUST be identical to columns in occs)
                                  # predic,             # column numbers of the predictor variables
                                  presenceCol,        # column number of the presence column
                                  predicCols,         # column numbers of the predictor variables
                                  Boyce = 'cr',       # smooth boyce spline method (Liu et al. 2024 Ecography)
                                  # choose from c('tp', 'cr', 'bs', 'ps', 'ad', 'm')
                                  # corresponds to 'thin-plate', 
                                  # in practice, selecting 'm' takes longer than model building (~150 occs + 20,000 background)
                                  # all methods broadly agree with each other (SBI within ~0.01)
                                  home=getwd(),       # directory where all the models will be ran
                                  all_models = TRUE,  # do you want to keep all versions of all the models?
                                  # setting to TRUE is helpful for quickly checking some models, but may consume >1GB . . .
                                  # . . . of hard disk space. Recommended to set to FALSE for exploratory analyses.
                                  # if FALSE, function will create a folder called "RunOver" and will . . .
                                  # . . . continuously write-over it for all models, and the only model you see . . .
                                  # . . . at the end will be the last model that was ran
                                  omission_rate=0.1,    # user specified omission rate to be calculated for the threshold
                                  # express as a proportion from 0-1
                                  all_background        # background data from the extents the species exists in
                                  # if using all the potential background points, all_background should == background
                                  ){
  # required packages
  require(mgcv)
  require(dplyr)
  require(splitstackshape)
  
  # checking boyce
  if(!Boyce %in% c('tp', 'cr', 'bs', 'ps', 'ad', 'm') | length(Boyce) != 1){
    stop('Set `Boyce` to one of: c(`tp`, `cr`, `bs`, `ps`, `ad`, `m`)')
  }
  
  # resetting data to data frames
  occs <- as.data.frame(occs)
  background <- as.data.frame(background)
  all_background <- as.data.frame(all_background)
  
  # Calculating the total number of models
  nmodels <- nrow(eval_df)
  
  # stop if an incompatible omission rate is specified
  if(omission_rate < 0 || omission_rate > 1){
    stop('Specify an omission rate between 0-1.')
  }
  
  # prompting the user if they want to store models in the RAM
  print(paste0('The time is ', round(Sys.time()), '. You are running ', nmodels, ' total Maxent models.'))
  
  # setting up the progress bar
  prog <- txtProgressBar(min=0, max=nrow(eval_df), style=3,  char='+')
  
  # setting up various list objects to send output  to
  # list for the maxent models
  maxent_list <- list()
  # list for the testing background data
  test_back_list <- list()
  # list for all the occs
  all_occs_list <- list()
  
  
  
  for(i in 1:nrow(eval_df)){
    
    # arguments
    fold_idx <- eval_df$fold[i]
    features <- eval_df$features[i]
    betas <- eval_df$betas[i]
    
    # setting the working directory for each folder
    if(all_models == TRUE){
      
      # change working directories
      setwd( paste(home, fold_idx, features, paste('beta', betas, sep='_'), sep='/') )
    } else {
      # create the RunOver folder
      if( !dir.exists( paste(home, 'RunOver', sep='/') ) ){
        dir.create( paste(home, 'RunOver', sep='/') )
        setwd( paste(home, 'RunOver', sep='/') )
      } else {
        setwd( paste(home, 'RunOver', sep='/') )
      }
    }
    
    ###########################################################################
    
    ### MaxEnt things happen here
    
    ### preparing the data for the maxent model
    
    # # filtering the occ object by it's respective cross-validation identity
    # cv_number <- eval_df$cv_num[i]
    # # assigning the column number to be sent through maxent
    # col_number <- first_occ_col + cv_number - 1
    
    # filtering the species dataset by col_number and assigning to eval_df
    # sp <- occs |> dplyr::filter(!fold %in% fold_idx)
    occ_train <- occs |> dplyr::filter(!fold %in% fold_idx)
    n <- nrow(occ_train)
    eval_df$n[i] <- n
    
    # generating the testing data
    
    # test.col_number <- first_test_col + cv_number - 1
    # sp_test <- occs |> dplyr::filter(fold %in% fold_idx)
    occ_test <- occs |> dplyr::filter(fold %in% fold_idx)
    
    
    # bending the occ and background data together
    mx_data <- rbind(occ_train, background)
    
    # arguments to be passed into maxent
    mx_args <- prepPara(userfeatures = eval_df$features[i], 
                        betamultiplier = eval_df$betas[i], 
                        doclamp = FALSE)
    
    ### running the actual maxent model
    mx_model <- predicts::MaxEnt(p = mx_data[,presenceCol],
                                 x = mx_data[,predicCols],
                                 path =  paste0(getwd()),
                                 args = mx_args
    )
    
    ### calculating k and the names of the lambdas
    lambda_file <- as.data.frame(mx_model@lambdas) |> `colnames<-`('lambdas')
    # lambdas data frame
    lambdas_df <- suppressWarnings(splitstackshape::cSplit(lambda_file, sep=',', splitCols='lambdas'))
    colnames(lambdas_df) <- c('feature', 'lambda', 'min', 'max')
    # finding the non-zero lambdas
    non_zero_lambdas <- lambdas_df |> 
      dplyr::filter(!is.na(max)) |> 
      dplyr::filter(lambda != 0)
    class(non_zero_lambdas) <- 'data.frame'
    # giving noting the variables/features/hyperparameters with non-zero lambdas
    eval_df$lambdas[i] <- toString( non_zero_lambdas[,1] )
    
    
    ## evaluating the models
    # predicting the training data
    mx_train_occ <-  predicts::predict(mx_model, occ_train[,predicCols])
    # predicting the training background data
    mx_train_back <- predicts::predict(mx_model, background[,predicCols])
    # predicting the testing data
    mx_test_occ <- predicts::predict(mx_model, occ_test[,predicCols])
    # projecting to all extents
    mx_test_back <- predicts::predict(mx_model, all_background[,predicCols] )
    # projecting to all occ points
    mx_all_occs <- predicts::predict(mx_model, occs[,predicCols])
    
    # calculate the threshold
    thresh <- quantile(mx_train_occ, omission_rate) |> as.numeric()
    
    # assigning the threshold value to the output file
    eval_df$thresh[i] <- thresh
    # calculating the test sensitivity
    sens <- length(which(mx_test_occ >= thresh)) / length(mx_test_occ)
    eval_df$testSens[i] <-  sens
    
    
    # Smooth Boyce Index (Liu et al. 2024 Ecography)
    p <- c(mx_test_occ, mx_train_back)
    n1 <- length(mx_test_occ)
    n0 <- length(mx_train_back)
    prd <- seq(min(p), max(p), length = n0)
    oc <- c(rep(1, n1), rep(0, n0))
    ktry <- 10
    if(Boyce == 'm'){
      md_tp = mgcv::gam(oc ~ s(p, bs = "tp", k = min(ktry, length(unique(p) ) ) ), family = binomial)
      prd_tp = predict(md_tp, newdata = data.frame(p = prd), type = 'response')
      md_cr = mgcv::gam(oc ~ s(p, bs = "cr", k = min(ktry, length(unique(p) ) ) ), family = binomial)
      prd_cr = predict(md_cr, newdata = data.frame(p = prd), type = 'response')
      md_bs = mgcv::gam(oc ~ s(p, bs = "bs", k = min(ktry, length(unique(p) ) ) ), family = binomial)
      prd_bs = predict(md_bs, newdata = data.frame(p = prd), type = 'response')
      md_ps = mgcv::gam(oc ~ s(p, bs = "ps", k = min(ktry, length(unique(p) ) ) ), family = binomial)
      prd_ps = predict(md_ps, newdata = data.frame(p = prd), type = 'response')
      md_ad = mgcv::gam(oc ~ s(p, bs = "ad", k = min(ktry, length(unique(p) ) ) ), family = binomial)
      prd_ad = predict(md_ad, newdata = data.frame(p = prd), type = 'response')
      prd_m = (prd_tp + prd_cr + prd_bs + prd_ps + prd_ad)/5
      boyce <- cor(prd, prd_m, method = "spearman")
    } else {
      md <- mgcv::gam(oc ~ s(p, bs = Boyce, k = min(ktry, length(unique(p) ) ) ), family = binomial)
      prd_new <- predict(md, newdata = data.frame(p = prd), type = 'response')
      boyce <- cor(prd, prd_new, method = "spearman")
    }
    eval_df$testBoyce[i] <- boyce
    
    
    
    mx_test_back_df <- as.data.frame( as.matrix(mx_test_back, ncol=1) )
    
    # assigning the objects to the various lists
    maxent_list[[i]] <- mx_model
    test_back_list[[i]] <- mx_test_back
    all_occs_list[[i]] <- mx_all_occs
    
    # updating the prograss bar for each run to get an idea of how long things will take
    setTxtProgressBar(prog, i)
    
    
    ###########################################################################
    # returning to the home directory
    setwd(home)
    
  }   # closes the for loop
  
  # bind the testing background data into a single data
  back_projections <- as.data.frame( do.call('cbind', test_back_list ) )
  colnames(back_projections) <- eval_df$CrossVal
  
  # bind all occs together
  occ_projections <- as.data.frame(do.call('cbind', all_occs_list))
  colnames(occ_projections) <- eval_df$CrossVal
  
  # make a list of the output
  out_list <- list()
  out_list$maxent_models <- maxent_list
  out_list$back_projection <- back_projections
  out_list$occ_projection <- occ_projections
  out_list$summary <- eval_df
  
  return(out_list)
}



##### MaxEnt Evaluation ----------------------------------------------------------------------------------------------------

# function that takes the eval object generated from maxent_crossval_error and calculates the weighted mean and standard deviation
# the weighted mean is calculating my testing sensitivity (1 - omission rate) multiplied by the partial ROC/AUC value.
# this ensures that if a model does not discriminate between presences/non-presences well, it will receive comparatively lower weight
# later package versions will include Boyce index as an additional calibration technique and will offer the user the ability to . . .
# . . . choose which metrics to use for assessing model reliability
# the weighted mean and standard deviation can be plotted to infer model variability/uncertainty

### ARGUMENTS ###
## eval <- eval object generated from `maxent_crossval_error` that will project the model to every grid cell used in the training region
## pROC_error <- the user-specified omission rate.
# (choose from 0.1, 1, 5 - however they almost always end up super correlated with each other)
maxent_eval <- function(eval){
  require(matrixStats)
  
  # model testing sensitivity
  sens <- eval$summary$testSens
  
  # model test boyce index
  boyce <- eval$summary$testBoyce
  
  # weights = sensitivity*AUC_ratio
  weights <- sens * boyce
  weights[is.nan(weights)] <- 0
  
  # making a matrix of the background points
  models_mat <- as.matrix(eval$back_projection)
  # weighted means and standard deviations
  w_mean <- matrixStats::rowWeightedMeans(x=models_mat, w=weights)
  w_sd <- matrixStats::rowWeightedSds(x=models_mat, w=weights)
  # merging the weighted means/sd's of all background points
  back_df <- data.frame(w_mean, w_sd)
  
  # making a matrix of the background points
  occs_mat <- as.matrix(eval$occ_projection)
  # weighted means and standard deviations
  w_mean <- matrixStats::rowWeightedMeans(x=occs_mat, w=weights)
  w_sd <- matrixStats::rowWeightedSds(x=occs_mat, w=weights)
  # merging the weighted means/sd's of all occ points
  occ_df <- data.frame(w_mean, w_sd)
  
  
  # making the output list
  out_list <- list()
  out_list$back <- back_df
  out_list$occ <- occ_df
  out_list$weights <- weights
  
  return(out_list)
  
}


##### maxent_eval_summarize -------------------------------------------------------------------------------------------------

# function that summarizes information across multiple cross-validation reps and find the "best" one using the 
# one-standard error rule
maxent_eval_summarize <- function(eval, omission_rate = 0.1){
  require(dplyr)
  
  # model testing sensitivity
  sens <- eval$summary$testSens
  
  # model test boyce index
  boyce <- eval$summary$testBoyce
  
  # weights = sensitivity*AUC_ratio
  weights <- sens * boyce
  weights[is.nan(weights)] <- 0
  
  eval_sum$weights <- weights
  eval_sum$AUC_ratio <- AUC_ratio
  # grouping by Features and Betas and calculating the performance weight
  eval_group <- eval_sum |> 
    dplyr::group_by(Features, Betas) |> 
    dplyr::summarize(across(c(test_sens, AUC_ratio, weights), mean))
  # calculating the one-standard error rule
  nmodels <- nrow(eval_group)
  SE <- sd(eval_group$weights)/(nmodels^0.5)
  oneMinusSE <- max(eval_group$weights) - SE
  eval_group$standard1mSE <- sapply(eval_group$weights, function(x) abs((x-oneMinusSE)/SE) )
  # returning eval_group
  return(eval_group |> 
           dplyr::arrange(standard1mSE))
}


##### Calculating the MaxEnt Threshold -------------------------------------------------------------------------------------

# function that calculates the threshold value for all the training data

### ARGUMENTS ###
## eval <- eval object generated from maxent_crossval_error that will project the model to every grid cell used in the training region
## occs <- the occurrence data.frame
## predicCols <- column numbers of the predictor variables 
## pROC_error <- the user-specified omission rate.
# (choose from 0.1, 1, 5 - however they almost always end up super correlated with each other)

# function will return the lowest training preference threshold.
# future package versions will give the opportunity to select different thresholds
maxent_thresh1 <- function(eval, occs, predicCols, omission_rate = 0.1){
  require(matrixStats)
  
  n_rows <- nrow(occs)
  n_cols <- nrow(eval$summary)
  
  # making a data.frame of the occurrences
  occ_values <- as.data.frame( matrix(NA, nrow=n_rows, ncol=n_cols ) )
  colnames(occ_values) <- eval$summary$CrossVal
  
  # for loop predicting all the cross validation maxent models
  for(i in 1:nrow(eval$summary) ){
    mx_occs <- predicts::predict(object=eval$maxent_models[[i]], x=occs[,predicCols])
    occ_values[,i] <- mx_occs
  }
  
  # model testing sensitivity
  sens <- eval$summary$testSens
  
  # model test boyce index
  boyce <- eval$summary$testBoyce
  
  # weights = sensitivity*AUC_ratio
  weights <- sens * boyce
  weights[is.nan(weights)] <- 0
  
  # making a matrix of the background points
  models_mat <- as.matrix(occ_values)
  # weighted means
  w_mean <- matrixStats::rowWeightedMeans(x=models_mat, w=weights)
  
  # setting the threshold based on the specified omission rate
  thresh <- quantile(w_mean, omission_rate) |> as.numeric()
  
  return(thresh)
}


##### Projecting the Best Models to All the Extents ------------------------------------------------------------------------

# function to project models to any and all extents you want to
# effectively making new master occurrence and master background files that can easily be saved and re-loaded again
maxent_everything <- function(eval, means, everything, thresh=NULL, predicCols, name='pc'){
  require(matrixStats)
  
  # making a data.frame of the occurrences
  every_value <- as.data.frame( matrix(NA, nrow=nrow(everything), ncol=nrow(eval$summary)) )
  colnames(every_value) <- eval$summary$CrossVal
  
  # for loop predicting all the cross validation maxent models
  for(i in 1:nrow(eval$summary) ){
    mx_everything <- predicts::predict(object=eval$maxent_models[[i]], x=everything[,predicCols])
    every_value[,i] <- mx_everything
  }
  
  # model testing sensitivity
  sens <- eval$summary$testSens
  
  # model test boyce index
  boyce <- eval$summary$testBoyce
  
  # weights = sensitivity*AUC_ratio
  weights <- sens * boyce
  weights[is.nan(weights)] <- 0
  
  # making a matrix of the background points
  models_mat <- as.matrix(every_value)
  # weighted means
  w_mean <- matrixStats::rowWeightedMeans(x=models_mat, w=weights)
  
  # making the output data frame
  out_df <- as.data.frame( matrix(NA, nrow=length(w_mean), ncol=2) )
  colnames(out_df) <- c(paste0(eval$summary[1,1], '_contin'), paste0(eval$summary[1,1], '_thresh') )
  # assigning the weighted means to the output data.frame
  out_df[,1] <- w_mean
  # replicating the weighted means so they can be thresolded
  out_df[,2] <- w_mean
  # calculating the threshold - if no threshold is provided, use the lowest training presence
  if(is.null(thresh)){
    thresh <- min(means$occ$w_mean)
    warning('No threshold provided.\nUsing lowest weighted training presence as the threshold.\nThis value is likely arbitrarily low.')
  }
  # thresolding the values
  out_df[,2][out_df[,2] >= thresh] <- 1
  out_df[,2][out_df[,2] < thresh] <- 0
  
  
  
  return(out_df)
  
}


##### Informed MESS Analysis -----------------------------------------------------------------------------------------------

# function that calculates Multivariate Environmental Suitability Surface (MESS), but allows for some extrapolation . . .
# . . . based on how response curves look.
informed_mess <- function(ref_extent,        # reference data.frame - must contain same column names as 'data'
                          mess_extent,        # data.frame containing the extent to be projected to - ideally this should be all extents
                          predicCols=3:5,    # column numbers of ONLY the final predictor variables
                          tolerance=NULL,    # tolerance vector to extrapolate beyond the limits of the data. see example below
                          tol_quant = 0.05, # minimum tolerance quantile
                          # if tolerance is not specified (default), function will calculate basic MESS
                          keep.layers=FALSE  # decision to retain the MESS values for all variables instead of just the final MESS data
                          # changing to TRUE can help easily identify which variables are causing the MESS value . . .
                          # . . . to be so low (i.e., it may be only 1 variable that is responsible)
                          # can easily re-check this later
){
  # isolating the coordinates and predictor variables
  mess_extent <- as.data.frame(mess_extent)
  ref_extent <- as.data.frame(ref_extent)
  new_vars <- mess_extent[, predicCols]
  ref_vars <- ref_extent[, predicCols]
  # 
  if( !is.null(tolerance) ){
    if(length(predicCols) != length(tolerance)/2 ){
      stop('Ensure that tolerance is 2x as long as `predicCols`')
    }
    # changing tolerance to a data.frame
    tol <- as.data.frame( matrix(tolerance, ncol=2, byrow=TRUE) )
    n_vars <- length(predicCols)
    # pre-changing values outside the range of the input data
    # extra rows
    extra_rows <- matrix(NA, nrow=2, ncol=n_vars)
    colnames(extra_rows) <- colnames(ref_vars)
    ref_vars <- rbind(ref_vars, extra_rows )
    for(i in 1:n_vars){
      # if we can extrapolate beyond the lower end of variable i, allow mess to not 
      if(tol[i,1] == 1){
        min_var <- quantile(ref_vars[,i], tol_quant, na.rm=TRUE) - 0.0001
        # min_var <- min(ref_vars[,i], na.rm=TRUE) - 0.0001
        new_vars[ new_vars[,i] < min_var, i ] <- min_var
        ref_vars[nrow(ref_vars)-1, i] <- min_var
      } else {
        ref_vars[nrow(ref_vars)-1, i] <- min(ref_vars[,i], na.rm=TRUE)
      }
      
      if(tol[i,2] == 1){
        max_var <- quantile(ref_vars[,i], (1 - tol_quant), na.rm=TRUE) + 0.0001
        # max_var <- max(ref_vars[,i], na.rm=TRUE) + 0.0001
        new_vars[ new_vars[,i] > max_var, i ] <- max_var
        ref_vars[nrow(ref_vars), i] <- max_var
      } else {
        ref_vars[nrow(ref_vars), i] <- max(ref_vars[,i], na.rm=TRUE)
      }
      
    } # end for(i in 1:n_vars) loop
    
  }
  
  # running the mess analysis
  mess_vars <- as.data.frame( sapply(1:ncol(new_vars), function(i) .messi3(new_vars[, i], ref_vars[, i])) )
  # re-asigning of the extrapolating points to have a mess value of 0. 
  nref <- nrow(ref_vars)
  fix_vars <- as.data.frame( sapply(1:ncol(new_vars), function(i) .messi3(ref_vars[ (nref-1):nref , i], ref_vars[, i])) )
  for(i in 1:n_vars){
    mess_vars[mess_vars[,i] == fix_vars[1,i],i] <- 0
    mess_vars[mess_vars[,i] == fix_vars[2,i],i] <- 0
    
  }
  
  # making a simple thresholded version of the mess analysis for easy plotting
  final_mess <- as.data.frame( apply(mess_vars, 1, min) )
  colnames(final_mess) <- 'mess_raw'
  mess_thresh <- final_mess$mess_raw
  mess_thresh[mess_thresh > 0] <- 1
  mess_thresh[mess_thresh < 0] <- -1
  final_mess <- cbind(final_mess, mess_thresh)
  
  # deciding whether to keep individual mess layers
  if(keep.layers==TRUE){
    colnames(ref_vars) <- paste0('mess_', colnames(ref_vars))
    final_mess <- cbind(final_mess, ref_vars)
    return(final_mess)
  } else {
    return(final_mess)
  }
  
}

# internal function originally from dismo that runs the actual mess analysis in data.frame format
.messi3 <- function(p,v) { # p=new_vars   v=ref_vars
  # seems 2-3 times faster than messi2
  v <- stats::na.omit(v)
  f <- 100*findInterval(p, sort(v)) / length(v)
  minv <- min(v)
  maxv <- max(v)
  res <- 2*f 
  f[is.na(f)] <- -99
  i <- f>50 & f<100
  res[i] <- 200-res[i]
  
  i <- f==0 
  res[i] <- 100*(p[i]-minv)/(maxv-minv)
  i <- f==100
  res[i] <- 100*(maxv-p[i])/(maxv-minv)
  res
}




##### Informed MESS conceptually -------------------------------------------------------------------------------------------

# principles of informed mess analysis


#            |                                                      #     # we haven't captured the entire response curve for this . . .
#      S   1-|          |                      |                 __ #     # . . . variable. this is a problem when extrapolating . . .
#      u     |                                           ?__-----   #     # . . . to non-analog environments
#      i     |          |                      |    ?___/           #
#      t     |                                  ___/                #     # we're very confident that projecting to environments with . . .
#      i     |          |                      |                    #     # . . . very low temperatures would likely mean low suitability
#      b     |                             ___  ____?____?____?____ #
#      i     |          |              ___/    |                    #     # however, we aren't confident about this at high temperatures.
#      l     |                      __/           \__?              #     # the suitability could tail off at any point, but the model . . .
#      i     |          |        __/           |      \__?          #     # . . . will likely assume ↑temp = ↑suitability forever
#      t     |                 _/                         \_?       #
#      y   0-|  ?___?___|_____/                |              \__?_ #     # to note this in the informed_mess function, you need to . . .
#            |_____________________________________________________ #     # . . . visualize the response curves in the maxent output file
#                       |       Variable       |                    #
#                                  1                                #     # note this in the 'tolerance' argument
#                       |  (e.g., Temperature) |                    #



### the tolerance argument is a vector telling the function what to do with the non-analog environments it encounters
## tolerance is a vector containing 1's and 0's based on what the response curves look like

### specify 1 or 0 for both ends of the response curves for all the final variables in the model
## check the lambdas column in the output summary file for the final variables
# for example, if the lambdas for a madel are: coldmon^2, warmmon^2, & coldmon*warmmon, . . . 
# . . . there are 2 final variables in the model (coldmon and warmmon), so tolerance should be 4 numbers long
# for another example, if the lambdas for a madel are: coldmon^2, warmmon^2, coldmon*warmmon, & seasonality, . . . 
# . . . there are 3 final variables in the model (coldmon, warmmon, and seasonality), so tolerance should be 6 numbers long
# for a final example, if the lambdas for a madel are: coldmon & coldmon^2,, . . . 
# . . . there is 1 final variable in the model (coldmon), so tolerance should be 2 numbers long


### the following example contains three final variables (Temperature, Precip, and Seasonality)
# thus, tolerance should be a vector with length 6
# tolerance <- c(_,_,
#                _,_,
#                _,_)

#            |                                                      #     # because we ARE confident in extrapolating to . . .
#      S   1-|          |                      |                    #     # . . .  low temps, put a 1 for the 1st spot in tolerance
#      u     |                                                      #
#      i     |          |                      |                    #     # tolerance <- c(1,_,
#      t     |    1                                       0         #     #                _,_,
#      i     |          |                      |                    #     #                _,_)
#      b     |                             ___                      #
#      i     |          |              ___/    |                    #     # however, because we're NOT confident in extrapolating . . .
#      l     |                      __/                             #     # . . . to high temps, put a 0 in the 2nd spot in tolerance
#      i     |          |        __/           |                    #
#      t     |                 _/                                   #     # tolerance <- c(1,0,
#      y   0-|          |_____/                |                    #     #                _,_,
#            |_____________________________________________________ #     #                _,_)
#                       |       Variable       |                    #
#                                  1                                #
#                       |  (e.g., Temperature) |                    #



#            |                                                      #     # because we ARE confident in extrapolating to . . .
#      S   1-|          |                      |                    #     # . . . low precip, put a 1 for the 3rd spot in tolerance
#      u     |                                                      #
#      i     |          |                      |                    #     # tolerance <- c(1,0,
#      t     |    1                                       1         #     #                1,_,
#      i     |          |                      |                    #     #                _,_)
#      b     |                                                      #
#      i     |          |          ___         |                    #     # because we ARE ALSO confident about extrapolating to . . .
#      l     |                  __/   \__                           #     # . . . high precip, put a 1 in the 4th spot in tolerance
#      i     |          |    __/         \__   |                    #
#      t     |             _/               \__                     #     # tolerance <- c(1,0,
#      y   0-|          |_/                    |                    #     #                1,1,
#            |_____________________________________________________ #     #                _,_)
#                       |       Variable       |                    #
#                                  2                                #
#                       |    (e.g., Precip)    |                    #



#            |                                                      #     # given that this response curve does NOT turn downwards on . . .
#      S   1-|          |                      |                    #     # . . . the low end, you COULD put a 0 for the 5th spot in . . .
#      u     |                                                      #     # . . . seasonality. HOWEVER, low seasonality is unlikely to . . .
#      i     |          |                      |                    #     # . . . be non-suitable for any species, we can assume that . . .
#      t     |    1                                      1          #     # . . . this extrapolation will NOT be problematic.
#      i     |          |                      |                    #
#      b     |           ___                                        #     # tolerance <- c(1,0,
#      i     |          |   \___               |                    #     #                1,1,
#      l     |                  \___                                #     #                1,_)
#      i     |          |           \__        |                    #
#      t     |                         \_                           #     # because we ARE confident about extrapolating to at high . . .
#      y   0-|          |                \_____|                    #     # . . . seasonality, put a 1 in the 6th spot in tolerance
#            |_____________________________________________________ #
#                       |       Variable       |                    #     # tolerance <- c(1,0,
#                                  3                                #     #                1,1,
#                       | (e.g., Seasonality)  |                    #     #                1,1)


### basically, a 0 tells the code to believe the MESS value that it calculates, whereas a 1 tells the code to . . .
#   . . . change the MESS value because we have confidence that the projection will still be reasonable

### thus, for the previous three response curves, specify the following for mess.tol
# mess.tol <- c(1,0,
#               1,1,
#               1,1)





##### calculate_mop --------------------------------------------------------------------------------------------------------


## calculate_mop calculates Mobility Oriented Parity (MOP; sensu Owens et al. 2013) from a set of reference points

# function identifies env-space that represents novel multivariate combinations relative to a set of reference points
# use this in your analyses if you want to be conservative about projecting ENMs to new areas/time periods


### ARGUMENTS
## ref_extent <- a data.frame or matrix containing the variable values for the reference/training extent, such as . . .
#                . . . any object generated with the create.background.df function
## predic <- the column numbers of the predictor variables in ref_extent
## new_extent <- a data.frame or matrix containing the variable values for the new extent to be projected to, such as . . .
#                . . . any object generated with the create.background.df function
## if new_extent is not provided, function will calculate MOP based on distance from each point in ref_extent to . . .
#  . . . the nearest set of points in ref_extent
## if new_extent is provided, function will calculate MOP based on distance from each point in new_extent to . . .
#  . . . the nearest set of points in ref_extent


## the function returns a data.frame of multivariate distances between a given set of points and other points in the . . .
#  . . . training extent with the same number of rows as either new_extent (if provided) or ref_extent
## distances are multivariate euclidean distances based on a set of reference points that have been standardized (μ=0, σ=1)
## it is recommended that the user run MOP on the ref_extent first to set a threshold for how much . . .
#  . . . environmental difference from ref_extent they are willing to tolerate
## suggested thresholds are near the 
#  note that this threshold is a trade-off. if ref_extent is small (i.e., few grid cells), the user may want to select . . .
#  . . . a higher MOP threshold as there may be greater chance for having holes/gaps in env-space.
#  thus, selecting an artificially low MOP threshold is likely to filter out interpret-able parts of env-space
#  this will also be an issue with the greater number of dimensions ref_extent has.
#  more dimensions = more chance for multivariate holes in env-space. these are often super hard/impossible to detect
## in the output data.frame, the column names are:
# pt1 = distance from a given point to the nearest point in env-space
# pt5 = distance from a given point to the nearest 5 points in env-space
# pt10 = distance from a given point to the nearest 10 points in env-space
# pt1perc = distance from a given point to the nearest 1% of points in env-space
# pt5perc = distance from a given point to the nearest 5% of points in env-space


## worth noting that this will easily blow out R's memory allotment if the number of points is high (>50000) or . . .
#  . . . if the dimensionality is high
# this is because the number of calculations in the distance matrix that is calculated is equal to . . .
# . . . (nrow(ref_extent) + nrow(new_extent))^2
# for example, having 50000 total points results in 2.5B calculations, which takes 3 minutes per dimension and . . .
# . . . will take a stupidly long time to complete
# an example with ~34000 background points took ~13 minutes
# you may see an error like this:
# Error: vector memory exhausted (limit reached?)
# e.g., this failed when there were ~150000 background points over 4 dimensions, even when setting . . .
# . . . R_MAX_VSIZE=100Gb as an environment
# example from https://stackoverflow.com/questions/51295402/r-on-macos-error-vector-memory-exhausted-limit-reached 

calculate_mop <- function(ref_extent, predic, new_extent=NULL){
  # initial condition checking to see if the number of background points is >50000
  if(!is.null(new_extent)){
    total_pts <- nrow(ref_extent) + nrow(new_extent)
  } else {
    total_pts <- nrow(ref_extent)
  }
  if(total_pts >= 50000){
    stop('Calculating ', total_pts, ' points will blow out the memory. Please subset total points to < 50000.')
  }
  # setting start time
  start_time <- Sys.time() |> round()
  message( paste0( 'The time is ', start_time ) )
  # number of variables ran
  n_vars <- length(predic)
  # number of background points
  n_points <- nrow(ref_extent)
  # reducing the data to only the variables of interest for the ref extent
  ref_extent <- ref_extent[,predic]
  # defining mop_extent if new_extent is provided
  # if new_extent is not provided, define mop_extent as ref_extent
  if( !is.null(new_extent) ){
    # number of points in new extent
    n_new_points <- nrow(new_extent)
    # reducing the data to only the variables of interest for the new extent
    new_extent <- new_extent[,predic]
    # ensuring that the columns of ref_extent and new_extent match
    # if so, bind them together (new_extent 1st!). if not, return an error
    if( identical( colnames(ref_extent), colnames(new_extent) ) ){
      mop_extent <- rbind(new_extent, ref_extent)
    } else {
      stop(paste0( ref_extent,' and ', new_extent,' have non-identical columns'  ) )
    }
  } else {
    mop_extent <- ref_extent
    n_new_points <- n_points
  }
  # calculating the mean and sd of each variable
  rescale_df <- as.data.frame( matrix(NA, nrow=2, ncol=n_vars) )
  for( i in 1:n_vars ){
    rescale_df[1,i] <- mean(ref_extent[,i])   # mean of each variable
    rescale_df[2,i] <- sd(ref_extent[,i])   # sd of each variable
  }
  # for loop standardizing each of the variables in both extents
  for (i in 1:ncol(ref_extent) ){
    # center each variable on 0
    mop_extent[,i] <- mop_extent[,i] - rescale_df[1,i]
    # divide each variable by the standard deviation
    mop_extent[,i] <- mop_extent[,i]/rescale_df[2,i]
  }
  # defining the output data.frame that is n.new. points long and measures distances from any point to nearby points
  output_df <- as.data.frame(matrix(NA, nrow=n_new_points, ncol=5))
  colnames(output_df) <- c('pt1',   # distance from each point to the nearest 1 point
                           'pt5',   # average distance from each point to the nearest 5 points
                           'pt10',   # average distance from each point to the nearest 10 points
                           'pt1perc',   # average distance from each point to the nearest 1% of points
                           'pt5perc')   # average distance from each point to the nearest 5% of points
  # defining a distance matrix between every background point and every other background point (for mop_extent)
  start_d <- Sys.time()
  message('Calculating distance matrix')
  d_matrix <- fields::rdist(mop_extent)
  end_d <- Sys.time()
  # if statement contingent on new_extent being provided
  if( is.null(new_extent) ){
    # defining the diagonal and lower half of the matrix as NA values
    # doing this because the distance between points a and b in the matrix is noted in cells [a,b] and [b,a] and . . .
    # . . . we don't want to use that distance twice 
    diag(d_matrix) <- NA
    d_matrix[upper.tri(d_matrix)] <- NA
    # defining the start time
    mop_start <- Sys.time()
    # setting up the progress bar
    message( paste0('Calculating MOP for the reference extent') )
    prog <- txtProgressBar(min=0, max=nrow(output_df), style=3,  char='+')
    # massive for loop calculating the distance from a given point to the nearest X points
    for( i in 1:n_new_points ){
      # the distance from point i to every other point is noted in d_matrix[i,] and d_matrix[,i]
      # filters out distance from point i to any other point and sorts if
      d_point <- sort( c( d_matrix[i,], d_matrix[,i] ) )
      # takes the ranked distance(s) from point i to every other point and averages the closest points
      output_df[i,1] <- d_point[1]   # distance from each point to the nearest 1 point
      output_df[i,2] <- mean(d_point[1:5])   # average distance from each point to the nearest 5 points
      output_df[i,3] <- mean(d_point[1:10])   # average distance from each point to the nearest 10 points
      output_df[i,4] <- mean(d_point[1:floor(n_points/100)])   # average distance from each point to the nearest 1% of points
      output_df[i,5] <- mean(d_point[1:floor(n_points/20)])   # average distance from each point to the nearest 5% of points
      # updating the progress bar for each run to get an idea of how long things will take
      setTxtProgressBar(prog, i)
    }
  } else {
    # filter d_matrix to only distances between new_extent and points in the old extent
    d_matrix <- d_matrix[ -(1:n_new_points) , 1:n_new_points ]
    # defining the start time
    mop_start <- Sys.time()
    # setting up the progress bar
    message( paste0('Calculating MOP for the new extent') )
    prog <- txtProgressBar(min=0, max=nrow(output_df), style=3,  char='+')
    # massive for loop calculating the distance from a given point to the nearest X points
    for( i in 1:n_new_points ){
      # the distance from point i to every other point is noted in d_matrix[i,] and d_matrix[,i]
      # filters out distance from point i to any other point and sorts if
      d_point <- sort(d_matrix[,i])
      # takes the ranked distance(s) from point i to every other point and averages the closest points
      output_df[i,1] <- d_point[1]   # distance from each point to the nearest 1 point
      output_df[i,2] <- mean(d_point[1:5])   # average distance from each point to the nearest 5 points
      output_df[i,3] <- mean(d_point[1:10])   # average distance from each point to the nearest 10 points
      output_df[i,4] <- mean(d_point[1:floor(n_points/100)])   # average distance from each point to the nearest 1% of points
      output_df[i,5] <- mean(d_point[1:floor(n_points/20)])   # average distance from each point to the nearest 5% of points
      # updating the progress bar for each run to get an idea of how long things will take
      setTxtProgressBar(prog, i)
    }
  } # closes ifelse is.null(new_extent) statement
  # returning how long it took to run aspects of the function
  # end time
  end_time <- Sys.time()
  # calculating distance matrix
  message(paste0('Time to calculate distance matrix: ', round(end_d - start_d, 5), ' ', units(end_d - start_d)) )
  # running MOP
  message( paste0('Time to run MOP: ', round(end_time - mop_start, 5), ' ', units(end_time - mop_start)) )
  # total function time
  message( paste0('Time to run function: ', round(end_time - start_time, 5), ' ', units(end_time - start_time)) )
  # return output
  return(output_df)
}


##### Variable Ranges ---------------------------------------------------

# function that returns the minimum and maximum ranges of each variable for a given set of cross-validation models
get_var_ranges <- function(eval){ # eval object generated from `maxent_crossval_error`
  r <- apply(X=eval$maxent_models[[1]]@absence, 2, 'range')
  colnames(r) <- colnames(eval$maxent_models[[1]]@absence)
  rownames(r) <- c('min', 'max')
  return(r)
}


##### Variable Importance ---------------------------------------------------

# function that returns the weighted variable contribution and importance for a given set of cross-validation models
# per Phillips 2006, 
# varContribution depends on the process of how the model is fit. intrepret values of highly-correlated variables with caution
# varImportance is calculated based on randomizing a single variable at a time AFTER the model is made and calculating how much the AUC drops
get_var_contrib_importance <- function(eval, pROC_error = 5){ # eval object generated from `maxent_crossval_error`
  require(matrixStats)
  ### getting model weights
  # extracting model sensitivity
  sens <- eval$summary$test_sens
  # which pROC error amount to use?
  if(pROC_error == 0.1){
    AUC_ratio <- eval$summary$pROC_0.1
  } else if(pROC_error == 1){
    AUC_ratio <- eval$summary$pROC_1
  } else if(pROC_error == 5){
    AUC_ratio <- eval$summary$pROC_5
  } else {
    stop('Select an appropriate partialROC error amount.')
  }
  # weights = sensitivity*AUC_ratio
  weights <- sens * AUC_ratio
  weights[is.nan(weights)] <- 0
  weights[is.na(weights)] <- 0
  # number of cross-validation reps
  n_reps <- length(eval$maxent_models)
  # number of variables
  n_vars <- ncol(eval$maxent_models[[1]]@presence)
  # making data.frames for variables and importance
  contrib_df <- matrix(NA, nrow=n_vars, ncol=n_reps)
  colnames(contrib_df) <- paste0('cv_', 1:n_reps)
  rownames(contrib_df) <- sort(colnames(eval$maxent_models[[1]]@presence)) # maxent sorts variables as strigs, so we have to here also
  impor_df <- contrib_df
  # for loop extracting the variable contribution and importance
  for(i in 1:n_reps){
    contrib_df[,i] <- eval$maxent_models[[i]]@results[grep("contribution", rownames(eval$maxent_models[[i]]@results)), ]
    impor_df[,i] <- eval$maxent_models[[i]]@results[grep("importance", rownames(eval$maxent_models[[i]]@results)), ]
  }
  
  # weighted means and sd
  w_contrib <- data.frame(w_mean = matrixStats::rowWeightedMeans(x=contrib_df, w=weights),
                          w_sd = matrixStats::rowWeightedSds(x=contrib_df, w=weights),
                          row.names = rownames(contrib_df)
  )
  w_impor <- data.frame(w_mean = matrixStats::rowWeightedMeans(x=impor_df, w=weights),
                        w_sd = matrixStats::rowWeightedSds(x=impor_df, w=weights),
                        row.names = rownames(contrib_df)
  )
  
  out_list <- list(w_contrib, w_impor)
  names(out_list) <- c('VarContribution', 'VarImportance')
  return(out_list)
  
}

get_var_contrib_importance1 <- function(eval){ # eval object generated from `maxent_crossval_error`
  require(matrixStats)
  ### getting model weights
  # model testing sensitivity
  sens <- eval$summary$testSens
  # model test boyce index
  boyce <- eval$summary$testBoyce
  # weights = sensitivity*AUC_ratio
  weights <- sens * boyce
  weights[is.nan(weights)] <- 0
  # number of cross-validation reps
  n_reps <- length(eval$maxent_models)
  # number of variables
  n_vars <- ncol(eval$maxent_models[[1]]@presence)
  # making data.frames for variables and importance
  contrib_df <- matrix(NA, nrow=n_vars, ncol=n_reps)
  colnames(contrib_df) <- paste0('cv_', 1:n_reps)
  rownames(contrib_df) <- sort(colnames(eval$maxent_models[[1]]@presence)) # maxent sorts variables as strigs, so we have to here also
  impor_df <- contrib_df
  # for loop extracting the variable contribution and importance
  for(i in 1:n_reps){
    contrib_df[,i] <- eval$maxent_models[[i]]@results[grep("contribution", rownames(eval$maxent_models[[i]]@results)), ]
    impor_df[,i] <- eval$maxent_models[[i]]@results[grep("importance", rownames(eval$maxent_models[[i]]@results)), ]
  }
  
  # weighted means and sd
  w_contrib <- data.frame(w_mean = matrixStats::rowWeightedMeans(x=contrib_df, w=weights),
                          w_sd = matrixStats::rowWeightedSds(x=contrib_df, w=weights),
                          row.names = rownames(contrib_df)
  )
  w_impor <- data.frame(w_mean = matrixStats::rowWeightedMeans(x=impor_df, w=weights),
                        w_sd = matrixStats::rowWeightedSds(x=impor_df, w=weights),
                        row.names = rownames(contrib_df)
  )
  
  out_list <- list(w_contrib, w_impor)
  names(out_list) <- c('varContribution', 'varImportance')
  return(out_list)
  
}



##### Extract Response Curves ---------------------------------------------------

# function will give the response curves for each variable while keeping all other variables at their median values for the training region
get_maxent_response_curves1 <- function(eval, # eval object generated from `maxent_crossval_error`
                                        expand=5, # extrapolation percentage - how far beyond the range of input values to project response curves to
                                        n_gradient=100 # number of points used to calculate the response curve gradient
){
  require(matrixStats)
  ### getting model weights
  # model testing sensitivity
  sens <- eval$summary$testSens
  # model test boyce index
  boyce <- eval$summary$testBoyce
  # weights = sensitivity*AUC_ratio
  weights <- sens * boyce
  weights[is.nan(weights)] <- 0
  ### preparing arguments
  # number of cross-validation reps
  n_reps <- length(eval$maxent_models)
  # make a data.frame of the output
  response_df <- as.data.frame(matrix(NA, nrow=n_gradient, ncol=n_reps))
  colnames(response_df) <- paste0('cv', seq(1:ncol(response_df)))
  # make background points for the gradient
  background <- eval$maxent_models[[1]]@absence
  cn <- colnames(background)
  # number of vars
  n_vars <- ncol(background)
  # structure of the data.frame that will be used to get values for each response curve
  m <- matrix(nrow=1, ncol=ncol(background))
  # defining medians
  m <- as.numeric(apply(background, 2, median)) # assigns median point of response curve
  medians <- as.data.frame(matrix(m, nrow=n_gradient, ncol=length(m), byrow=TRUE))
  colnames(medians) <- cn
  input_range <- as.data.frame(matrix(NA, nrow=n_gradient, ncol=length(m), byrow=TRUE)) # input ranges of response curves
  colnames(input_range) <- cn
  # ranges of input variables
  ranges <- apply(X=background, 2, 'range')
  # creating the output_list
  out_list <- list()
  # filling out input_range
  for(i in 1:length(m)){
    # min and max of the variable
    minmax <- ranges[,i]
    r <- minmax[2]-minmax[1]
    expand1 <- r*round(expand)/100
    # assigning ranges to each variable number
    input_range[,i] <- minmax[1] - expand1 + 0:(n_gradient-1) * (r + 2*expand1)/(n_gradient-1)
  } # close i loop
  
  # nested for loop extracting each variable from each cross-valication rep
  for(var in 1:n_vars){ # going variable by variable, then rep by rep
    # making a list of cross-val reps
    list.reps <- list()
    # for loop analyzing the cross-val reps
    for(rep in 1:n_reps){
      # making a data.frame of medians, but changing the var in question to the values of input_range
      mx <- medians
      mx[,var] <- input_range[,var]
      # projecting maxent model to mx and assigning it to list.reps
      list.reps[[rep]] <- predicts::predict(eval$maxent_models[[rep]], mx)
    }   # closing for(rep in n_reps)
    
    # binding everything together
    reps <- do.call('cbind', list.reps)
    colnames(reps) <- paste0(cn[var], '_cv', 1:n_reps)
    w_mean <- matrixStats::rowWeightedMeans(x=reps, w=weights)
    w_sd <- matrixStats::rowWeightedSds(x=reps, w=weights)
    reps <- as.data.frame(reps)
    reps$w_mean <- w_mean
    reps$w_sd <- w_sd
    colnames(reps) <-  c(paste0(cn[var], '_cv', 1:n_reps), paste0(cn[var], '_mean'), paste0(cn[var], '_sd'))
    reps$input <- input_range[,var]
    
    out_list[[var]] <- reps
  } # closing var loop
  
  names(out_list) <- cn
  
  return(out_list)
}


##### 2-Dimensional Response Curves ---------------------------------------------------


# predicCols <- column numbers of the predictor variables 
# vars2 <- column numbers of the predictor variables 

maxent_2d_respose_curves <- function(eval, # eval object generated from `maxent_crossval_error`
                                     expand=5, # extrapolation percentage - how far beyond the range of input values to project response curves to
                                     n_gradient=50, # number of points used to calculate the response curve gradient
                                     predicCols,
                                     vars2, # which two variables to project to
                                     varlimits = c(0,1,0,1), # the min and max values of the first and second variables in vars2 - expressed in native units
                                     ahull=0.1 # sensitivity of the alpha value in the alpha hull (expressed as a proportion)
){
  ## QC checks
  # to check if the specified `varlimits` is valid
  if(!unique(vars2 %in% predicCols)){
    stop('`vars2` must be a subset of `predicCols`.')
  }
  
  ### getting model weights
  # model testing sensitivity
  sens <- eval$summary$testSens
  # model test boyce index
  boyce <- eval$summary$testBoyce
  # weights = sensitivity*AUC_ratio
  weights <- sens * boyce
  weights[is.nan(weights)] <- 0
  
  
  
}

##### Process Maxent Response Curves ---------------------------------------------------

# function that generates a list of data.frames that represent mean ±1 sd bounding boxes around uncertainties
# input is an object generated by get_maxent_response_curves
process_maxent_response_curves1 <- function(rc){
  # number of variables
  n_vars <- length(rc)
  # making output_list
  output_list <- list()
  # for loop making data.frames of the uncertainties of the variable contributions
  for(i in 1:n_vars){
    
    lower <- data.frame(x=rc[[i]][,grep('input', colnames( rc[[i]] )) ],
                        y=rc[[i]][,grep('mean', colnames( rc[[i]] )) ] - rc[[i]][,grep('sd', colnames( rc[[i]] )) ]
    )
    upper <- data.frame(x=rc[[i]][,grep('input', colnames( rc[[i]] )) ],
                        y=rc[[i]][,grep('mean', colnames( rc[[i]] )) ] + rc[[i]][,grep('sd', colnames( rc[[i]] )) ]
    )
    upper <- upper[rev(rownames(upper)), ]
    rownames(upper) <- NULL
    both <- rbind(lower , upper)
    both[nrow(both) + 1,] <- both[1,]
    rownames(both) <- NULL
    # fixing any values below or above (0,1) to (0,1)
    both$y <- replace(both$y, both$y>1, 1)
    both$y <- replace(both$y, both$y<0, 0)
    
    output_list[[i]] <-  both
  }
  names(output_list) <- names(rc)
  return(output_list)
}


##### suit_uncert_plot -----------------------------------------------------------------------------------------------------

### makes a bivariate plot of weighted suitability and weighted uncertainty (standard deviation)

suit_uncert_plot <- function(means){
  require(ggplot2)
  
  back <- means$back
  occs <- means$occ
  
  ltp <- min(occs$w_mean)
  
  output <- ggplot(data=back, aes(x=w_mean, y=w_sd)) + geom_point(colour='black', size=0.75) +
    geom_point(data=occs, aes(x=w_mean, y=w_sd), size=2, colour='red', shape=18 ) +
    ylim(0, 0.55) + xlim(0,1) + theme_classic() + coord_fixed(1/0.55) +
    annotate(geom='text', x=0.05, y=0.45, label=round(ltp, 4), hjust=0 )
  
  return(output)
  
}



##### END ------------------------------------------------------------------------------------------------------------------




