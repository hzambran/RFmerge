# Merging of satellite datasets with ground observations using Random Forests (RF)

Merging of satellite datasets with ground observations using Random
Forests (RF)

## Usage

``` r
RFmerge(x, ...)

# Default S3 method
RFmerge(x, metadata, cov, mask, training, 
        id="id", lat = "lat", lon = "lon", 
        ED = TRUE, rasterizedED=FALSE, 
        seed = NULL, ntree = 2000, na.action = stats::na.omit,
        parallel=c("none", "parallel", "parallelWin"),
        par.nnodes=parallel::detectCores()-1, 
        par.pkgs= c("terra", "randomForest", "zoo"), write2disk=FALSE,
        drty.out, use.pb=TRUE, verbose=TRUE,...)

# S3 method for class 'zoo'
RFmerge(x, metadata, cov, mask, training, 
        id="id", lat = "lat", lon = "lon", 
        ED = TRUE, rasterizedED=FALSE, 
        seed = NULL, ntree = 2000, na.action = stats::na.omit,
        parallel=c("none", "parallel", "parallelWin"),
        par.nnodes=parallel::detectCores()-1, 
        par.pkgs= c("terra", "randomForest", "zoo"), write2disk=FALSE,
        drty.out, use.pb=TRUE, verbose=TRUE, ...)
```

## Arguments

- x:

  data.frame with the ground-based values that will be used as the
  dependent variable to train the RF model.  
  Every column must represent one ground-based station and the codes of
  the stations must be provided as colnames. class(data) must be zoo.

- metadata:

  data.frame with the metadata of the ground-based stations. At least,
  it MUST have the following 3 columns:  

  -) id: This column stores the unique identifier (ID) of each
  ground-based observation. Default value is "id".  
  -) lat: This column stores the latitude of each ground observation.
  Default value is "lat".  
  -) lon: This column stores the longitude of each ground observation.
  Default value is "lon".

- cov:

  List with all the covariates used as independent variables in the
  Random Forest model. The individual covariates must be SpatRaster
  objects, either when they vary in time (e.g., individual gridded
  precipitation datasets) or does not vary in time (e.g., a digital
  elevation model).  
  All time-varying covariates in `cov` MUST have the same number of
  layers (bands), which is internally checked using the ,
  [`nlyr`](https://rspatial.github.io/terra/reference/dimensions.html)
  function. Covariates that do not change in time (e.g., a DEM) are
  internally transformed into SpatRaster objects with the same number of
  layers as the other time-varying elements in `cov`.

- mask:

  OPTIONAL. If provided, the final merged product mask out all values in
  `cov` outside `mask`.  
  Spatial object (vectorial) with the spatial borders of the study area
  (e.g., catchment, administrative borders). `class(mask)` must be a
  SpatVect object with "POLYGON" or "MULTIPOLYGON" geometry.

- training:

  Numeric indicating the percentage of stations that will be used in the
  training set.  
  The valid range is from zero to one. If `training` = 1, all the
  stations will be used for training purposes.

- id:

  Character, with the name of the column in `metadata` where the
  identification code (ID) of each station is stored.

- lat:

  Character, with the name of the column in `metadata` where the
  latitude of the stations is stored.

- lon:

  Character, with the name of the column in `metadata` where the
  longitude of the stations is stored.

- ED:

  logical, should the Euclidean distances be computed an used as
  covariates in the random forest model?. The default value is `TRUE`.

- rasterizedED:

  logical, should the Euclidean distances between stations and grid
  cells be computed to the actual point coordinate
  (`rasterizedED=FALSE`) or to the rasterized station cell center
  (`rasterizedED=TRUE`). When `rasterizedED=FALSE`, the the
  self-distance between the station and its own grid cell is usually not
  zero. By default, `rasterizedED=FALSE`.

- seed:

  Numeric, indicating a single value, interpreted as an integer, or
  null.

- parallel:

  character, indicates how to parallelise ‘RFmerge’ (to be precise, only
  the evaluation of the objective function `fn` is parallelised). Valid
  values are:  

  -)none: no parallelisation is made (this is the default value)  

  -)parallel: parallel computations for network clusters or machines
  with multiple cores or CPUs. A ‘FORK’ cluster is created with the
  [`makeForkCluster`](https://rdrr.io/r/parallel/makeCluster.html)
  function.  

  -)parallelWin: parallel computations for network clusters or machines
  with multiple cores or CPUs (this is the only parallel implementation
  that works on Windows machines). A ‘PSOCK’ cluster is created with the
  [`makeCluster`](https://rdrr.io/r/parallel/makeCluster.html) function.

- par.nnodes:

  OPTIONAL. Used only when `parallel!='none'`  

  numeric, indicates the number of cores/CPUs to be used in the local
  multi-core machine, or the number of nodes to be used in the network
  cluster.  

  By default `par.nnodes` is set to the amount of cores detected by the
  function `detectCores()` (parallel package)

- par.pkgs:

  OPTIONAL. Used only when `parallel='parallelWin'`  

  list of package names (as characters) that need to be loaded on each
  node for allowing the objective function `fn` to be evaluated. By
  default `c("terra", "randomForest", "zoo")`.

- ntree:

  Numeric indicating the maximum number trees to grow in the Random
  Forest algorithm. The default value is set to 2000. This should not be
  set to too small a number, to ensure that every input row gets
  predicted at least a few times. If this value is too low, the
  prediction may be biased.

- na.action:

  A function to specify the action to be taken if NAs are found. (NOTE:
  If given, this argument must be named.)

- write2disk:

  logical, indicates if the output merged SpatRaster layers and the
  training and evaluation datasets (two files each, one with time series
  and other with metadata) will be written to the disk. By default
  `write2disk=FALSE`

- drty.out:

  Character with the full path to the directory where the final merged
  product will be exported as well as the training and evaluation
  datasets. Only used when `write2disk=TRUE`

- use.pb:

  logical, indicates if a progress bar should be used to show the
  progress of the random forest computations (it might reduce a bit the
  performance of the computations, but it is useful to track if
  everything is working well). By default `use.pb=TRUE`

- verbose:

  logical, indicates if progress messages are to be printed. By default
  `verbose=TRUE`

- ...:

  further arguments to be passed to the low level function
  randomForest.default.

## Value

It returns a `SpatRaster` object with as many layers as time steps are
present in `x`. Each one of the layers in the output object has the same
spatial resolution and spatial extent as the `cov` argument.

## References

Baez-Villanueva, O. M.; Zambrano-Bigiarini, M.; Beck, H.; McNamara, I.;
Ribbe, L.; Nauditt, A.; Birkel, C.; Verbist, K.; Giraldo-Osorio, J.D.;
Thinh, N.X. (2020). RF-MEP: a novel Random Forest method for merging
gridded precipitation products and ground-based measurements. *Remote
Sensing of Environment*, 239, 111610.
[doi:10.1016/j.rse.2019.111606](https://doi.org/10.1016/j.rse.2019.111606)
.  

Hengl, T.; Nussbaum, M.; Wright, M. N.; Heuvelink, G. B.; Gräler, B.
(2018). Random forest as a generic framework for predictive modeling of
spatial and spatio-temporal variables. *PeerJ*, 6, e5518.
[doi:10.7717/peerj.5518](https://doi.org/10.7717/peerj.5518) .

## Author

Oscar M. Baez-Villanueva, <oscar.baezvillanueva@ugent.be>  
Mauricio Zambrano-Bigiarini, <mzb.devel@gmail>  
Juan D. Giraldo-Osorio, <j.giraldoo@javeriana.edu.co>

## See also

[`terra`](https://rspatial.github.io/terra/reference/terra-package.html),
[`rast`](https://rspatial.github.io/terra/reference/rast.html),
[`nlyr`](https://rspatial.github.io/terra/reference/dimensions.html),
[`resample`](https://rspatial.github.io/terra/reference/resample.html),
[`rotate`](https://rspatial.github.io/terra/reference/rotate.html),
[`crop`](https://rspatial.github.io/terra/reference/crop.html).

## Examples

``` r
library(terra)

data(ValparaisoPPts)    
data(ValparaisoPPgis) 
   
ValparaisoSHP.fname <- system.file("extdata/ValparaisoSHP.shp",package="RFmerge")
ValparaisoSHP       <- terra::vect(ValparaisoSHP.fname)

chirps.fname   <- system.file("extdata/CHIRPS5km.tif",package="RFmerge")
prsnncdr.fname <- system.file("extdata/PERSIANNcdr5km.tif",package="RFmerge")
dem.fname      <- system.file("extdata/ValparaisoDEM5km.tif",package="RFmerge")

CHIRPS5km        <- rast(chirps.fname)
PERSIANNcdr5km   <- rast(prsnncdr.fname)
ValparaisoDEM5km <- rast(dem.fname)

covariates <- list(chirps=CHIRPS5km, persianncdr=PERSIANNcdr5km, 
                   dem=ValparaisoDEM5km)

# \donttest{

# The following code assumes that the region is small enough for neglecting
# the impact of computing Euclidean distances in geographical coordinates.
# If this is not the case, please read the vignette 'Tutorial for merging 
# satellite-based precipitation datasets with ground observations using RFmerge' 


# without using parallelisation
rfmep <- RFmerge(x=ValparaisoPPts, metadata=ValparaisoPPgis, cov=covariates, 
                id="Code", lat="lat", lon="lon", mask=ValparaisoSHP, training=1)
#> [ Creating the training (100%) and evaluation (0%) datasets ... ]
#> Warning: [vect] guessed crs
#> [ Computing the Euclidean distances to each observation of the training set ...]


# Detecting if your OS is Windows or GNU/Linux, 
# and setting the 'parallel' argument accordingly:
onWin <- ( (R.version$os=="mingw32") | (R.version$os=="mingw64") )
ifelse(onWin, parallel <- "parallelWin", parallel <- "parallel")
#> [1] "parallel"

#Using parallelisation, with a maximum number of nodes/cores to be used equal to 2:
par.nnodes <- min(parallel::detectCores()-1, 2)
rfmep <- RFmerge(x=ValparaisoPPts, metadata=ValparaisoPPgis, cov=covariates,
                 id="Code", lat="lat", lon="lon", mask=ValparaisoSHP, 
                 training=0.8, parallel=parallel, par.nnodes=par.nnodes)
#> [ Creating the training (80%) and evaluation (20%) datasets ... ]
#> Warning: [vect] guessed crs
#> [ Computing the Euclidean distances to each observation of the training set ...]
#> [ Parallelisation finished ! ]
# }
```
