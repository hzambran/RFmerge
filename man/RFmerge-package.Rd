%% File RFmerge-package.Rd
%% Part of the RFmerge R package, https://github.com/hzambran/RFmerge ; 
%%                                https://cran.r-project.org/package=RFmerge
%% Copyright 2019-2020 Mauricio Zambrano-Bigiarini, Oscar M. Baez-Villanueva
%% Distributed under GPL 3 or later

\name{RFmerge-package}
\Rdversion{1.1}
\alias{RFmerge-package}
\docType{package}
\title{
Merging of Satellite Datasets with Ground Observations using Random Forests
}
\description{
S3 implementation of the Random Forest MErging Procedure (RF-MEP), which combines two or more satellite-based datasets (e.g., precipitation products, topography) with ground observations to produce a new dataset with improved spatio-temporal distribution of the target field. In particular, this package was developed to merge different Satellite-based Rainfall Estimates (SREs) with measurements from rain gauges, in order to obtain a new precipitation dataset where the time series in the rain gauges are used to correct different types of errors present in the SREs. However, this package might be used to merge other hydrological/environmental satellite fields with point observations. For details, see Baez-Villanueva et al. (2020) <doi:10.1016/j.rse.2019.111606>. Bugs / comments / questions / collaboration of any kind are very welcomed.
}
\details{
\tabular{ll}{
Package: \tab RFmerge\cr
Type: \tab Package\cr
Version: \tab 0.1-5\cr
Date: \tab 2020-01-07\cr
License: \tab GPL >= 3 \cr
LazyLoad: \tab yes\cr
Packaged: \tab Tue Jan  7 21:17:04 -03 2020; MZB \cr
BuiltUnder: \tab R version 3.6.2 (2019-12-12) -- "Dark and Stormy Night" ; x86_64-pc-linux-gnu (64-bit) \cr
}
%%~~ An overview of how to use the package, including the most important functions ~~
}
\author{
Mauricio Zambrano-Bigiarini, Oscar M. Baez-Villanueva \cr

Maintainer: Mauricio Zambrano-Bigiarini <mzb.devel@gmail>
}
\references{
\cite{Baez-Villanueva, O. M.; Zambrano-Bigiarini, M.; Beck, H.; McNamara, I.; Ribbe, L.; Nauditt, A.; Birkel, C.; Verbist, K.; Giraldo-Osorio, J.D.; Thinh, N.X. (2020). RF-MEP: a novel Random Forest method for merging gridded precipitation products and ground-based measurements, Remote Sensing of Environment, 239, 111610. \doi{10.1016/j.rse.2019.111606}. <https://authors.elsevier.com/c/1aKrd7qzSnJWL>}. \cr

\cite{Hengl, T., Nussbaum, M., Wright, M. N., Heuvelink, G. B., & Gr\"{a}ler, B. (2018). Random forest as a generic framework for predictive modeling of spatial and spatio-temporal variables. PeerJ, 6, e5518}.
}

\keyword{ package }
\seealso{
\url{https://cran.r-project.org/package=raster}. \cr
\url{https://cran.r-project.org/package=hydroGOF}.
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
rfmep <- RFmerge(x=ValparaisoPPts, metadata=ValparaisoPPgis, cov=covariates, 
                id="Code", lat="lat", lon="lon", mask=ValparaisoSHP, training=1)

# Limiting the maximum number of cores to be used in this toy example to 2
par.nnodes <- min(parallel::detectCores()-1, 2)

# using parallelisation in GNU/Linux machines
rfmep <- RFmerge(x=ValparaisoPPts, metadata=ValparaisoPPgis, cov=covariates,  
                id="Code", lat="lat", lon="lon", mask=ValparaisoSHP, training=1,
                parallel="parallel", par.nnodes=par.nnodes)


# using parallelisation in Window$ machines
rfmep <- RFmerge(x=ValparaisoPPts, metadata=ValparaisoPPgis, cov=covariates,  
                id="Code", lat="lat", lon="lon", mask=ValparaisoSHP, training=1,
                parallel="parallelWin", par.nnodes=par.nnodes)
}

}
