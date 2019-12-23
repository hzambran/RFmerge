%% File RFmerge-package.Rd
%% Part of the RFmerge R package, https://github.com/hzambran/RFmerge ; 
%%                                http://cran.r-project.org/web/packages/RFmerge/
%% Copyright 2019-2019 Mauricio Zambrano-Bigiarini
%% Distributed under GPL 3 or later

\name{RFmerge-package}
\Rdversion{1.1}
\alias{RFmerge-package}
\docType{package}
\title{
Management, analysis, interpolation and plot of hydrological time series, with focus on hydrological modelling
}
\description{
S3 functions for management, analysis, interpolation and plotting of time series used in hydrology and related environmental sciences. In particular, this package is highly oriented to hydrological modelling tasks. 
The focus of this package has been put in providing a collection of tools useful for the daily work of hydrologists (although an effort was made to optimise each function as much as possible, functionality has had priority over speed). Bugs / comments / questions / collaboration of any kind are very welcomed, and in particular, datasets that can be included in this package for academic purposes.
}
\details{
\tabular{ll}{
Package: \tab RFmerge\cr
Type: \tab Package\cr
Version: \tab 0.1-3\cr
Date: \tab 2019-12-23\cr
License: \tab GPL >= 3 \cr
LazyLoad: \tab yes\cr
Packaged: \tab Mon Dec 23 12:20:24 -03 2019; MZB \cr
BuiltUnder: \tab R version 3.6.2 (2019-12-12) -- "Dark and Stormy Night" ; x86_64-pc-linux-gnu (64-bit) \cr
}
%%~~ An overview of how to use the package, including the most important functions ~~
}
\author{
Mauricio Zambrano-Bigiarini, Oscar M. Baez-Villanueva \cr

Maintainer: Mauricio Zambrano-Bigiarini <mzb.devel@gmail>
}
%%\references{
%%~~ Literature or other references for background information ~~
%%}

\keyword{ package }
\seealso{
\url{https://github.com/hzambran/hydroGOF}. \cr
\url{https://github.com/hzambran/hydroPSO}. \cr
}
\examples{
library(rgdal)
library(raster)

data(ValparaisoPPts)    
data(ValparaisoPPgis) 
data(ValparaisoSHP)    

chirps.fname   <- system.file("extdata/CHIRPS5km.tif",package="RFmerge")
prsnncdr.fname <- system.file("extdata/PERSIANNcdr5km.tif",package="RFmerge")
dem.fname      <- system.file("extdata/ValparaisoDEM5km.tif",package="RFmerge")

CHIRPS5km        <- brick(chirps.fname)
PERSIANNcdr5km   <- brick(prsnncdr.fname)
ValparaisoDEM5km <- raster(dem.fname)

covariates <- list(chirps=CHIRPS5km, persianncdr=PERSIANNcdr5km, 
                   dem=ValparaisoDEM5km)

\donttest{

# The following code assumes that the region is small enough for neglecting
# the impact of computing Euclidean distances in geographical coordinates.
# If this is not the case, please read the vignette 'Tutorial for merging 
# satellite-based precipitation datasets with ground observations using RFmerge' 

# without using parallelisation
RFmerge(x=ValparaisoPPts, metadata=ValparaisoPPgis, cov=covariates, 
        id="Code", lat="lat", lon="lon", mask=ValparaisoSHP, training=1)

# Limiting the maximum number of cores to be used in this toy example to 2
par.nnodes <- max(parallel::detectCores()-1, 2)

# using parallelisation in GNU/Linux machines
RFmerge(x=ValparaisoPPts, metadata=ValparaisoPPgis, cov=covariates,  
        id="Code", lat="lat", lon="lon", mask=ValparaisoSHP, training=1,
        parallel="parallel", par.nnodes=par.nodes)


# using parallelisation in Window$ machines
RFmerge(x=ValparaisoPPts, metadata=ValparaisoPPgis, cov=covariates,  
        id="Code", lat="lat", lon="lon", mask=ValparaisoSHP, training=1,
        parallel="parallelWin", par.nnodes=par.nodes)
}

}
