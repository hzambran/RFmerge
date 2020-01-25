# File RFmerge.R
# Part of the RFmerge R package, https://github.com/hzambran/RFmerge ; 
#                                https://cran.r-project.org/package=RFmerge
# Copyright 2019-2020 Mauricio Zambrano-Bigiarini, Oscar M. Baez-Villanueva
# Distributed under GPL 3 or later
# 
# 
################################################################################
#                                 RFmerge                                      #
################################################################################
# This function creates an improved satellite product by combining one or more #
# satellite datsets with ground-based observations                             #
################################################################################
# Authors : Oscar M. Baez-Villanueva                                           #
#           Mauricio Zambrano-Bigiarini                                        #
#           Juan Diego Giraldo-Osorio                                          #
################################################################################
# Started: Jan 2019 (the individual function)                                  #
# Started: 07-Nov-2019 (the package)                                           #
# Updates: 16-Nov-2019 ; 10-Dec-2019 ; 11-Dec-2019 ; 12-Dec-2019 ; 13-Dec-2019 #
#          14-Dec-2019 ; 17-Dec-2019 ; 20-Dec-2019 ; 23-Dec-2019               #
################################################################################

# 'x'        : zoo object with ground-based values that will be used as the dependent variable to train the RF model.
#              Every column must represent one ground-based station and the codes of the station must be
#              provided as colnames. class(data) must be zoo.
#
# 'metadata' : data.frame with the metadata of the ground-based stations. The first column must have the codesas specified in the 
#              colnames of the data object. The second and third row must have the longitude and latitude.
#              The predefined names of the coordinates are set as "lon" and "lat", respectively.
#
# 'cov'      : covariates used as independent variables. class(cov) must be list. The selected covariates are introduced
#              as a stack file (in the case of the dynamic covariates) or as a raster file (in the case of static covariates
#             , e.g., a digital elevation model). Every covariate must be stored in a different object in the list.
#
# 'mask'     : spatial vector object of the study area. class(mask) must be SpatialPolygonsDataFrame.
#
# 'drty.out' : path to the directory where the final product and the training and evaluation sets will be exported
#
# 'training' : Numerical value: This value indicates the fraction of stations that will be used in the training set. The 
#              valid range is from zero to one. If training = 1, all the stations will be used for training purposes.
#
# 'id'      : name of the column in the metadata where the ID of the stations is stored.
#
# 'lon'      : name of the column in the metadata where the longitude of the stations is stored.
#
# 'lat'      : name of the column in the metadata where the latitude of the stations is stored.
#
# 'seed'     : a single value, interpreted as an integer, or null.
#
# 'ED'       : logical, should the Euclidean distances be computed an used as covariates in the random fores model?. The default value is \code{TRUE}.
#
# 'parallel' :character, indicates how to parallelise \sQuote{hydroPSO} (to be precise, only the evaluation of the objective function \code{fn} is parallelised). Valid values are: \cr
#             -)\kbd{none}: no parallelisation is made (this is the default value)\cr
#             -)\kbd{parallel}: parallel computations for network clusters or machines with multiple cores or CPUs. A \sQuote{FORK} cluster is created with the \code{\link[parallel]{makeForkCluster}} function.  When \code{fn.name="hydromod"} the evaluation of the objective function \code{fn} is done with the \code{\link[parallel]{clusterApply}} function of the \pkg{parallel} package. When \code{fn.name!="hydromod"} the evaluation of the objective function \code{fn} is done with the \code{\link[parallel]{parRapply}} function of the \pkg{parallel} package.\cr
#             -)\kbd{parallelWin}: parallel computations for network clusters or machines with multiple cores or CPUs (this is the only parallel implementation that works on Windows machines). A \sQuote{PSOCK} cluster is created with the \code{\link[parallel]{makeCluster}} function. When \code{fn.name="hydromod"} the evaluation of the objective function \code{fn} is done with the \code{\link[parallel]{clusterApply}} function of the \pkg{parallel} package. When \code{fn.name!="hydromod"} the evaluation of the objective function \code{fn} is done with the \code{\link[parallel]{parRapply}} function of the \pkg{parallel} package. 
# 'par.nnodes': OPTIONAL. Used only when \code{parallel!='none'}. numeric, indicates the number of cores/CPUs to be used in the local multi-core machine, or the number of nodes to be used in the network cluster. By default \code{par.nnodes} is set to the amount of cores detected by the function \code{detectCores()} (\pkg{parallel} package)
#
# 'par.pkgs'  : OPTIONAL. Used only when \code{parallel='parallelWin'}. list of package names (as characters) that need to be loaded on each node for allowing the objective function \code{fn} to be evaluated. By default \code{c("raster", "randomForest", "zoo")}.
#
# 'ntree'     : number of decision trees generated in the Random Forest. The default value is set to 2000. If this value is too low, the prediction may be biased.              
#
# 'write2disk': logical, indicates if the output merged raster layers and the training and avaluation datasets (two files each, one with time series and other with metadata) will be written to the disk. By default \code{write2disk=TRUE}
#
# 'verbose'  : logical, indicates if progress messages are to be printed. By default \code{verbose=TRUE}
#
# 'use.pb'   : logical, indicates if a progress bar should be used to show the progress of the random forest computations (it might reduce a bit the performance of the computations, but it is useful to track if everything is working well). By default \code{use.pb=TRUE}

RFmerge <- function(x, ...) UseMethod("RFmerge")


RFmerge.default <- function(x, metadata, cov, mask, training, 
                            id="id", lat = "lat", lon = "lon", ED = TRUE, 
                            seed = NULL, ntree = 2000, na.action = stats::na.omit,
                            parallel=c("none", "parallel", "parallelWin"),
	                    par.nnodes=parallel::detectCores()-1, 
                            par.pkgs= c("raster", "randomForest", "zoo"), 
                            write2disk=FALSE, drty.out, use.pb=TRUE, verbose=TRUE,
                            ...) {

     # Checking that 'x' is a zoo object
     if ( !is.zoo(x) ) stop("Invalid argument: 'class(x)' must be 'zoo' !!")

     RFmerge.zoo(x=x, metadata=metadata, cov=cov, mask=mask, 
                 training=training, id=id, lon=lon, lat=lat, ED=ED, 
                 seed=seed, ntree=ntree, na.action=na.action, 
                 parallel=parallel, par.nnodes=par.nnodes, par.pkgs=par.pkgs, 
                 write2disk=write2disk, drty.out=drty.out, use.pb=use.pb, 
                 verbose=verbose, ...)

} # 'RFmerge.default' end



RFmerge.zoo <- function(x, metadata, cov, mask, training, 
                        id="id", lat = "lat", lon = "lon", ED = TRUE, 
                        seed = NULL, ntree = 2000, na.action = stats::na.omit,
                        parallel=c("none", "parallel", "parallelWin"),
	                par.nnodes=parallel::detectCores()-1, 
                        par.pkgs= c("raster", "randomForest", "zoo"), 
                        write2disk=FALSE, drty.out, use.pb=TRUE, verbose=TRUE,
                        ...) {

  parallel <- match.arg(parallel)    
  
  # Cheking if the 'x' object provided is a zoo object
  if (!zoo::is.zoo(x) )
    stop("Invalid argument: 'x' must be a zoo object")
  
  # Cheking if the 'metadata' object provided is a data frame
  if (!is.data.frame(metadata) ) {
    stop("Invalid argument: 'metadata' must be a data frame !")
  } else if ( ( !(id %in% colnames(metadata) ) ) | ( !(lat %in% colnames(metadata) ) ) | ( !(lon %in% colnames(metadata) ) ) )
      stop("Invalid argument: 'metadata' must have 'id' (identifier), 'lon' (longitude) and 'lat' (latitude) fields !")

  # Cheking the mask
  mask.crs <- NULL
  if (!missing(mask)) {
    if ( !sf::st_is(mask, c("POLYGON", "MULIPOLYGON")) ) {
      stop("Invalid argument: 'mask' must be a 'sf' (multi)polygon object !!")
    } else mask.crs <- sf::st_crs(mask)[["epsg"]]
  } # IF end
  
  # Cheking if the user specified the seed
  if ( !is.null(seed) ) set.seed(seed)
  
  # Cheking if the value of the 'training' object is within the limits
  if ( (training < 0) | (training > 1) )
    stop("Invalid argument: 'training' must be in [0, 1]")

  # Cheking 'cov'
  cov.crs <- NULL
  if( !is.list(cov) ) { # is it a list object?
    stop("Invalid argument: 'cov' must be a list with all the covariates !")
  } else # does it have names?
      if ( is.null(names(cov)) ) {
        stop("Invalid argument: Please provide names for all the elements in 'cov' !!")
      } else
          if (!raster::compareRaster(cov)) {
            stop("Invalid argument: All the elements in 'cov' must have the same spatial extent, CRS, rotation and geometry !!")
          } else cov.crs <- sf::st_crs(cov[[1]])[["epsg"]]


  # Checking that the CRS of 'mask' and 'cov' are the same, if 'mask' is provided
  if (!missing(mask)) {
    if (mask.crs != cov.crs)
      stop("Invalid argument: 'cov' and 'mask' have different CRS !!")
  } # IF end

  # Cheking that all the time-varying covariates have the same length or they have 1 layer only
  cov.layers <- as.numeric( sapply(cov, raster::nlayers) )
  if ( length( unique(cov.layers) ) > 2 )
    stop("Invalid argument: Check that all the time-varying covariates in 'cov' have the same number of layers !!")
  
  # Create covariates of the same length  
  temp <- which(cov.layers == 1)
  
  set.covariates <- function(cov, temp, cov.layers){
    raster::stack( replicate( max(cov.layers), cov[[temp]] ) )
  }
  
  cov[temp] <- sapply(temp, set.covariates, cov=cov, cov.layers=cov.layers)
  
  # Dividing the Ground-based observations in 'training' and 'validation' samples
  if (verbose) message("[ Creating the training (", round(training*100,2), "%) and evaluation (", round((1-training)*100,2), "%) datasets ... ]" )
  train.sample   <- sort(sample(1:nrow(metadata), round(training * nrow(metadata), 0)))
  
  train.metadata <- metadata[train.sample,]
  train.ts       <- x[, which(names(x) %in% train.metadata[,id])]
  
  eval.metadata  <- metadata[-train.sample,]
  eval.ts        <- x[, which(names(x) %in% eval.metadata[,id])]
  

  # Creating output directories, if necessary
  if (write2disk) {

    if ( !file.exists(drty.out) ) dir.create(drty.out)
  
    # Creating subfolders in the output directory
    merged.drty     <- file.path(drty.out, "RF-MEP")
    grounddata.drty <- file.path(drty.out, "Ground_based_data")
    training.drty   <- file.path(drty.out, "Ground_based_data", "Training")
    evaluation.drty <- file.path(drty.out, "Ground_based_data", "Evaluation")

    if ( !file.exists(grounddata.drty) ) dir.create(grounddata.drty)  
    if ( !file.exists(training.drty) )   dir.create(training.drty)  
    if ( !file.exists(evaluation.drty) ) dir.create(evaluation.drty)  
    if ( !file.exists(merged.drty) )     dir.create(merged.drty)

  } else merged.drty <- ""
  
  # Converting the training metadata into a 'SpatialPointsDataFrame'
  points <- train.metadata
  sp::coordinates(points) <- c(lon, lat)
  
  # If required, computing the euclidean distances
  if (verbose) message("[ Computing the Euclidean distances to each observation of the training set ...]")
  
  lsample <- cov[[1]][[1]]
  if (ED) {
    buff.dist <- as(lsample, "SpatialPixelsDataFrame")
    buff.dist <- .buffer_dist(points, buff.dist, as.factor(1:nrow(data.frame(points))))
    buff.dist <- as(buff.dist, "RasterStack")
  } # IF end

  ########################################################################
  ##                        parallel: start (ini)                        #
  ########################################################################
  if (parallel != "none") {
    
    if ( (parallel=="parallel")  & 
       ( (R.version$os=="mingw32") | (R.version$os=="mingw64") ) )
       stop("[ Fork clusters are not supported on Windows =>  'parallel' can not be set to '", parallel, "' ]")
    
    ifelse(parallel=="parallelWin", parallel.pkg <- "parallel",  parallel.pkg <- parallel)                
    if  ( length(find.package("parallel", quiet=TRUE)) == 0 )  {
            warning("[ Package '", parallel.pkg, "' is not installed =>  parallel='none' ]")
            parallel <- "none"
    }  else { 
      
         if (verbose) message("                               ")
         if (verbose) message("[ Parallel initialization ... ]")
         
         nnodes.pc <- parallel::detectCores()
         if (verbose) message("[ Number of cores/nodes detected: ", nnodes.pc, " ]")
           
         if (write2disk) {
           if ( (parallel=="parallel") | (parallel=="parallelWin") ) {             
              logfile.fname <- paste(file.path(drty.out), "/", "parallel_logfile.txt", sep="") 
              if (file.exists(logfile.fname)) file.remove(logfile.fname)
           } # IF end
         } else 
             if ( (R.version$os=="mingw32") | (R.version$os=="mingw64") ) {
               logfile.fname <- "nul"
             } else logfile.fname <- "/dev/null"
             
         if (is.na(par.nnodes)) {
           par.nnodes <- nnodes.pc
         } else if (par.nnodes > nnodes.pc) {
               warning("[ 'nnodes' > number of detected cores (", par.nnodes, ">", nnodes.pc, ") =>  par.nnodes=", nnodes.pc, " ] !",)
               par.nnodes <- nnodes.pc
           } # ELSE end
           
         if (verbose) message("[ Number of cores/nodes used    : ", par.nnodes, " ]")                 
               
         if (parallel=="parallel") {       
             cl <- parallel::makeForkCluster(nnodes = par.nnodes, outfile=logfile.fname)
         } else if (parallel=="parallelWin") {    
             cl <- parallel::makeCluster(par.nnodes, outfile=logfile.fname)
             pckgFn <- function(packages) {
               for(i in packages) library(i, character.only = TRUE)
             } # 'packFn' END
             parallel::clusterCall(cl, pckgFn, par.pkgs)
             parallel::clusterExport(cl, utils::ls.str(mode="function",envir=.GlobalEnv) )
           } # ELSE end                   
           
       } # ELSE end  
  
  }  # IF end    

  ########################################################################
  ##                        parallel: start (end)                        #
  ########################################################################

  # number of time steps used in the RF procedure
  ntsteps <- raster::nlayers(cov[[1]])
  ldates  <- zoo::index(x)
  
  if (use.pb) pbapply::pboptions(char = "=")

  # Evaluating an R Function 
  if (parallel=="none") {
    if (use.pb) {
      out <- pbapply::pbsapply(X=1:ntsteps, FUN=.lrf, ldates, metadata, id, points, cov, mask, merged.drty, na.action, ntree, ED, train.ts, train.metadata, buff.dist, lsample, write2disk)
    } else out <- sapply(X=1:ntsteps, FUN=.lrf, ldates, metadata, id, points, cov, mask, merged.drty, na.action, ntree, ED, train.ts, train.metadata, buff.dist, lsample, write2disk)
    out <- stack(out)
  } else {
      if (use.pb) {
        out <- pbapply::pbsapply(X=1:ntsteps, FUN=.lrf, ldates, metadata, id, points, cov, mask, merged.drty, na.action, ntree, ED, train.ts, train.metadata, buff.dist, lsample, write2disk, cl=cl, ...)
      } else  out <- parallel::clusterApply(cl= cl, 1:ntsteps, fun=.lrf, ldates, metadata, id, points, cov, mask, merged.drty, na.action, ntree, ED, train.ts, train.metadata, buff.dist, lsample, write2disk, ...)
      out <- stack(out)
      parallel::stopCluster(cl)
      if (verbose) message("[ Parallelisation finished ! ]")
   } # ELSE end    
    

  # If the user wants to write the outputs to the disk
  if ( write2disk) {
    # Exporting the training dataset to the output directory
    fname <- file.path(training.drty, "Training_ts.txt")
    zoo::write.zoo(train.ts, file=fname)

    fname <- file.path(training.drty, "Training_metadata.txt")
    utils::write.table(train.metadata, file=fname, row.names=FALSE, sep=",")
  
    # Exporting the evaluation dataset to the output directory
    if (training < 1) {
      fname <- file.path(evaluation.drty, "Evaluation_ts.txt")
      zoo::write.zoo(eval.ts, file=fname)

      fname <- file.path(evaluation.drty, "Evaluation_metadata.txt")
      utils::write.table(eval.metadata, file=fname, row.names=FALSE, sep=",", quote=FALSE)
    } # IF end
  } # IF end

  # output object
  names(out) <- ldates[1:ntsteps]
  return(out)

} # 'RFmerge.zoo' END



.lrf <- function(day, ldates, metadata, id, points, cov, mask, merged.drty, na.action, ntree, ED, train.ts, train.metadata, buff.dist, lsample, write2disk, ...) {
      
    # Obtaining the ground-based measurements for each day
    nr         <- nrow(train.metadata)
    codes      <- vector("numeric", length= nr)
    obs.values <- codes 
    
    for ( i in 1:nr ) {      
      metadata      <- train.metadata[i,]
      codes[i]      <- as.character(metadata[,id])
      obs.values[i] <- train.ts[day, grep(metadata[,id], names(train.ts))]
    } # FOR end
    
    # Getting the covariates for each day
    ncovs   <- length(cov)
    cov.day <- vector("list", ncovs)
    
    for(i in 1:ncovs)
      cov.day[[i]] <- cov[[i]][[day]]    
    
    cov.day <- raster::stack(cov.day)

    if (ED) cov.day <- raster::stack(cov.day, buff.dist)
    #cov.day <- velox::velox(cov.day)
    
    # Extraction of covariates at the ground-station locations
    extraction     <- data.frame(raster::extract(cov.day, points))
    #extraction     <- data.frame(cov.day$extract(points))
    names(cov.day) <- names(extraction)
     
    # Constructing data frame to train the RF model
    table     <- data.frame(codes, obs.values, extraction)
    table.cov <- table[,2:ncol(table)]
    
    # Initializing 'result' object
    result <- lsample * 0
    
    # Condition: if all the observed values for the training set are zero, the predicted layer will be zero 
    if ( sum(table$obs.values, na.rm = TRUE) > 0) {
      
      # Training of the RF model and generation of the spatial prediction
      rf.model <- randomForest::randomForest(obs.values ~ ., data = table.cov, na.action = na.action, ntree = ntree, ... )
      
      # rf.model <- randomForest::randomForest(obs.values ~ ., data = table.cov, ntree = 2000, na.action = na.omit)
      result   <- raster::predict(cov.day, rf.model)
    }
    
    # Masking the raster layer using the mask of the study area
    if ( !missing(mask) ) 
      result <- raster::mask(result, mask)

    # Exporting the RF-MEP product for each day
    if ( write2disk ) {
      ldate <- paste0("RF-MEP_", ldates[day], ".tif") 
      fname <- file.path(merged.drty, ldate)
      raster::writeRaster(result, fname, format = "GTiff", overwrite = TRUE)
    } # IF end

    return(result)

  } # '.lrf' END



