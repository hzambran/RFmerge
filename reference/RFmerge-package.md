# Merging of Satellite Datasets with Ground Observations using Random Forests

S3 implementation of the Random Forest MErging Procedure (RF-MEP), which
combines two or more satellite-based datasets (e.g., precipitation
products, topography) with ground observations to produce a new dataset
with improved spatio-temporal distribution of the target field. In
particular, this package was developed to merge different
Satellite-based Rainfall Estimates (SREs) with measurements from rain
gauges, in order to obtain a new precipitation dataset where the time
series in the rain gauges are used to correct different types of errors
present in the SREs. However, this package might be used to merge other
hydrological/environmental satellite fields with point observations. For
details, see Baez-Villanueva et al. (2020)
\<doi:10.1016/j.rse.2019.111606\>. Bugs / comments / questions /
collaboration of any kind are very welcomed.

## Details

|  |  |
|----|----|
| Package: | RFmerge |
| Type: | Package |
| Version: | 0.3-3 |
| Date: | 2026-05-07 |
| License: | GPL \>= 3 |
| LazyLoad: | yes |
| Packaged: | Thu May 7 09:41:19 -04 2026 ; MZB |
| BuiltUnder: | R version 4.6.0 (2026-04-24) – "Because it was There" ; aarch64-apple-darwin23 |

## Author

Mauricio Zambrano-Bigiarini, Oscar M. Baez-Villanueva  

Maintainer: Mauricio Zambrano-Bigiarini \<mzb.devel@gmail\>

## References

Baez-Villanueva, O. M.; Zambrano-Bigiarini, M.; Beck, H.; McNamara, I.;
Ribbe, L.; Nauditt, A.; Birkel, C.; Verbist, K.; Giraldo-Osorio, J.D.;
Thinh, N.X. (2020). RF-MEP: a novel Random Forest method for merging
gridded precipitation products and ground-based measurements, Remote
Sensing of Environment, 239, 111610.
[doi:10.1016/j.rse.2019.111606](https://doi.org/10.1016/j.rse.2019.111606)
. \<https://doi.org/10.1016/j.rse.2019.111606\>.  

Hengl, T., Nussbaum, M., Wright, M. N., Heuvelink, G. B., & Gräler, B.
(2018). Random forest as a generic framework for predictive modeling of
spatial and spatio-temporal variables. PeerJ, 6, e5518.
[doi:10.7717/peerj.5518](https://doi.org/10.7717/peerj.5518) .

## See also

<https://cran.r-project.org/package=terra>.  
<https://cran.r-project.org/package=hydroGOF>.

## Examples

``` r
library(terra)
#> terra 1.9.27

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
