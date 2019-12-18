# RFmerge
[![CRAN](http://www.r-pkg.org/badges/version/RFmerge)](https://cran.r-project.org/package=RFmerge) [![status](https://tinyverse.netlify.com/badge/RFmerge)](https://CRAN.R-project.org/package=RFmerge)  [![License](https://img.shields.io/badge/license-GPL%20%28%3E=%203%29-lightgrey.svg?style=flat)](http://www.gnu.org/licenses/gpl-3.0.html) [![monthly](http://cranlogs.r-pkg.org/badges/RFmerge)](https://www.rpackages.io/package/RFmerge) [![total](http://cranlogs.r-pkg.org/badges/grand-total/RFmerge)](https://www.rpackages.io/package/RFmerge)

RFmerge provides an S3 implementation of the Random Forest MErging Procedure (RF-MEP), which combines two or more satellite-based datasets (e.g., precipiation products, topography) with ground observations to produce a new dataset with improved spatio-temporal distribution of the target field. 

In particular, this package was developed to merge different Satellite-based Rainfall Estimates (SREs) with measurements from rain gauges, in order to obtain a new precipitation dataset where the time series in the rain gauges are used to correct different types of errors present in the SREs. However, this package might be used to merge other hydrological/environmental satellite fields with point observations. 

Bugs / comments / questions / collaboration of any kind are very welcomed.


## Installation
Installing the latest stable version from [CRAN](https://CRAN.R-project.org/package=RFmerge):
```{r}
install.packages("RFmerge")
```

Alternatively, you can also try the under-development version from [Github](https://github.com/hzambran/RFmerge):
```{r}
if (!require(devtools)) install.packages("devtools")
library(devtools)
install_github("hzambran/RFmerge")
```


## A simple first application:

### Loading required packages:

```{r Loading_other_pks, eval = TRUE, message=FALSE}
library(zoo)
library(hydroTSM)
library(hydroGOF)
library(sf)
library(raster)
library(RFmerge)
```

### Loading times series and metadata of ground observations:

   
```{r Loading_GroundObservarions, eval = TRUE}
data(ValparaisoPPts)
data(ValparaisoPPgis) 
data(ValparaisoSHP)
```

```{r SpatialMetadata}
stations <- ValparaisoPPgis
( stations <- st_as_sf(stations, coords = c('lon', 'lat'), crs = 4326) )
```

### Loading satellite-based datasets:
   
```{r LoadingSatelliteData, eval = TRUE}
chirps.fname   <- system.file("extdata/CHIRPS5km.tif",package="RFmerge")
prsnncdr.fname <- system.file("extdata/PERSIANNcdr5km.tif",package="RFmerge")
dem.fname      <- system.file("extdata/ValparaisoDEM5km.tif",package="RFmerge")

CHIRPS5km        <- brick(chirps.fname)
PERSIANNcdr5km   <- brick(prsnncdr.fname)
ValparaisoDEM5km <- raster(dem.fname)
```

### Reprojecting input datsets

Reprojecting the input datsets from geographic coordinates into WGS 84 / UTM zone 19S (EPSG:32719):

```{r ReprojectingRasters}
utmz19s.p4s <- CRS("+init=epsg:32719") # WGS 84 / UTM zone 19S

CHIRPS5km.utm        <- projectRaster(from=CHIRPS5km, crs=utmz19s.p4s)
PERSIANNcdr5km.utm   <- projectRaster(from=PERSIANNcdr5km, crs=utmz19s.p4s)
ValparaisoDEM5km.utm <- projectRaster(from=ValparaisoDEM5km, crs=utmz19s.p4s)

stations.utm <- sf::st_transform(stations, crs=32719) # for 'sf' objects

ValparaisoSHP.utm <- sf::st_transform(ValparaisoSHP, crs=32719)
```

Creating a new object with the spatial location of the ground observations in  WGS 84 / UTM zone 19S:

```{r FinalMEtadata}
st.coords <- st_coordinates(stations.utm)
lon       <- st.coords[, "X"]
lat       <- st.coords[, "Y"]

ValparaisoPPgis.utm <- data.frame(ID=stations.utm[["Code"]], lon=lon, lat=lat)
```

### Covariates

Raster covariates to be used in `RFmerge`:

```{r CovariatesCreation}
covariates.utm <- list(chirps=CHIRPS5km.utm, persianncdr=PERSIANNcdr5km.utm, 
                   dem=ValparaisoDEM5km.utm)
```

### Running `RFmerge` 

Runing `RFmerge` with parallelisation in GNU/Linux machines:

```{r RFmergeWithLinuxParallelisation, eval = TRUE}
drty.out <- "~/Test.par"
rfmep <- RFmerge(x=ValparaisoPPts, metadata=ValparaisoPPgis.utm, cov=covariates.utm,
                 id="ID", lat="lat", lon="lon",  
                 mask=ValparaisoSHP.utm, drty.out = drty.out, training=0.8,
                 parallel="parallel")
```



## Reporting bugs, requesting new features

If you find an error in some function, or want to report a typo in the documentation, or to request a new feature (and wish it be implemented :) you can do it [here](https://github.com/hzambran/RFmerge/issues)


## Citation 
```{r}
citation("RFmerge")
```

To cite RFmerge in publications use:

> Zambrano-Bigiarini, M.; Baez-Villanueva, O.M., Giraldo-Osorio, J. RFmerge: Merging of Satellite Datasets with Ground Observations using Random Forests. R package version 0.1-0. URL https://hzambran.github.io/RFmerge/. DOI:10.5281/zenodo.1287350

> Baez-Villanueva, O. M.; Zambrano-Bigiarini, M.; Beck, H.; McNamara, I.; Ribbe, L.; Nauditt, A.; Birkel, C.; Verbist, K.; Giraldo-Osorio, J.D.; Thinh, N.X. (2019). RF-MEP: a novel Random Forest method for merging gridded precipitation products and ground-based measurements, Remote Sensing of Environment, (accepted).


BibTeX entries for LaTeX users are:

> @Manual{Zambrano-Bigiarini+al2019-RFmerge_pkg,
>     title = {RFmerge: Merging of Satellite Datasets with Ground Observations using Random Forests},
>     author = {{Mauricio Zambrano-Bigiarini} and {Oscar M. Baez-Villanueva} and {Juan Giraldo-Osorio}},
>     note = {R package version 0.1 . doi:10.5281/zenodo.1287350},
>     url = {https://CRAN.R-project.org/package=RFmerge},
>   }

> @Article{BaezVillanueva+al2019-RFmerge_article,
>     title = {RF-MEP: a novel Random Forest method for merging gridded precipitation products and ground-based measurements},
>     journal = {Remote Sensing of Environment},
>     author = {{Baez-Villanueva} and O. M. and {Zambrano-Bigiarini} and {M.} and {Beck} and {H.} and {McNamara} and {I.} and {Ribbe} and {L.} and {Nauditt} and {A.} and {Birkel} and {C.} and {Verbist} and {K.} and {Giraldo-Osorio} and {J.D.} and {Thinh} and {N.X.}},
>     year = {2019},
>   }

## Vignette 
[Here](https://github.com/hzambran/RFmerge/blob/master/vignettes/RFmerge-RainfallExample.pdf) you can find an introductory vignette showing the use of `RFmege` to create an improved precipitation dataset by combining the satellite-based CHIRPSv2 and PERSIANN-CDR precipitation products, elevations from a DEM and rainfall observations recorded in rain gauges.



## Related Material 

* *A novel methodology for merging different gridded precipitation products and ground-based measurements* (**EGU-2019**)  [abstract](https://meetingorganizer.copernicus.org/EGU2019/EGU2019-10659.pdf). EGU General Assembly 2019. Wien, Austria. (oral presentation). HS7.2 Precipitation Modelling: uncertainty, variability, assimilation, ensemble simulation and downscaling. Abstract EGU2019-10659


## See Also 

* [hydroTSM: Time Series Management, Analysis and Interpolation for Hydrological Modelling](https://github.com/hzambran/hydroTSM).

* [hydroGOF: Goodness-of-fit functions for comparison of simulated and observed hydrological time series](https://github.com/hzambran/hydroGOF).

