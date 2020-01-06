%%File RFmerge.Rd
%% Part of the RFmerge R package, https://github.com/hzambran/RFmerge ; 
%%                                https://cran.r-project.org/package=RFmerge
%% Copyright 2019-2020 Mauricio Zambrano-Bigiarini, Oscar M. Baez-Villanueva
%% Distributed under GPL 3 or later

\name{RFmerge}
\Rdversion{1.1}
\alias{RFmerge}
\alias{RFmerge.default}
\alias{RFmerge.zoo}


%- Also NEED an '\alias' for EACH other topic documented here.
\title{
Merging of satellite datasets with ground observations using Random Forests (RF)
}
\description{
Merging of satellite datasets with ground observations using Random Forests (RF)
}
\usage{
RFmerge(x, ...)

\method{RFmerge}{default}(x, metadata, cov, mask, training, 
        id="id", lat = "lat", lon = "lon", ED = TRUE, 
        seed = NULL, ntree = 2000, na.action = stats::na.omit,
        parallel=c("none", "parallel", "parallelWin"),
        par.nnodes=parallel::detectCores()-1, 
        par.pkgs= c("raster", "randomForest", "zoo"), write2disk=FALSE,
        drty.out, use.pb=TRUE, verbose=TRUE,...)

\method{RFmerge}{zoo}(x, metadata, cov, mask, training, 
        id="id", lat = "lat", lon = "lon", ED = TRUE, 
        seed = NULL, ntree = 2000, na.action = stats::na.omit,
        parallel=c("none", "parallel", "parallelWin"),
        par.nnodes=parallel::detectCores()-1, 
        par.pkgs= c("raster", "randomForest", "zoo"), write2disk=FALSE,
        drty.out, use.pb=TRUE, verbose=TRUE, ...)


}
%- maybe also 'usage' for other objects documented here.
\arguments{
   \item{x}{
data.frame witht the ground-based values that will be used as the dependent variable to train the RF model. \cr
Every column must represent one ground-based station and the codes of the stations must be provided as colnames. class(data) must be zoo.
}  
    \item{metadata}{
data.frame with the metadata of the ground-based stations.  At least, it MUST have the following 3 columns: \cr
-) \kbd{id}: This column stores the unique identifier (ID) of each ground-based observation. Default value is \kbd{"id"}. \cr
-) \kbd{lat}: This column stores the latitude of ech ground observation. Default value is \kbd{"lat"}. \cr
-) \kbd{lon}: This column stores the longitude of ech ground observation. Default value is \kbd{"lon"}.
}
    \item{cov}{
List with all the covariates used as independent variables in the Random Forest model. The individual covariates can be a \kbd{RasterStack} or \kbd{RasterBrick} object when they vary in time, or they can be a single \kbd{RasterLayer} object when they do not change in time (e.g., a digital elevation model). \cr
All time-varying covariates in \code{cov} MUST have the same number of layers (bands). Covariates that do not change in time (e.g., a DEM) are internally trasnformed into \kbd{RasterStack} or \kbd{RasterBrick} objects with the same number of layers as the other time-varying elements in \code{cov}
}
    \item{mask}{
OPTIONAL. If provided, the final merged product maks out all values in \code{cov} outside \code{mask}. \cr
Spatial object (vectorial) with the spatial borders of the study area (e.g., catchment, administrative borders). \code{class(mask)} must be \kbd{SpatialPolygonsDataFrame} or \kbd{sf} with "POLYGON" or "MULTIPOLYGON" geometry.
}
   \item{training}{
Numeric indicating the percentage of stations that will be used in the training set. \cr
The valid range is from zero to one. If \code{training} = 1, all the stations will be used for training purposes.
}
   \item{id}{
Character, with the name of the column in \code{metadata} where the identification code (ID) of each station is stored.
}
   \item{lat}{
Character, with the name of the column in \code{metadata} where the latitude of the stations is stored.
}
   \item{lon}{
Character, with the name of the column in \code{metadata} where the longitude of the stations is stored.
}
   \item{ED}{
logical, should the Euclidean distances be computed an used as covariates in the random fores model?. The default value is \code{TRUE}.
}
   \item{seed}{
Numeric, indicating  a single value, interpreted as an integer, or null.
}
  \item{parallel}{
character, indicates how to parallelise \sQuote{RFmerge} (to be precise, only the evaluation of the objective function \code{fn} is parallelised). Valid values are: \cr
-)\kbd{none}: no parallelisation is made (this is the default value)\cr
-)\kbd{parallel}: parallel computations for network clusters or machines with multiple cores or CPUs. A \sQuote{FORK} cluster is created with the \code{\link[parallel]{makeForkCluster}} function.  When \code{fn.name="hydromod"} the evaluation of the objective function \code{fn} is done with the \code{\link[parallel]{clusterApply}} function of the \pkg{parallel} package. When \code{fn.name!="hydromod"} the evaluation of the objective function \code{fn} is done with the \code{\link[parallel]{parRapply}} function of the \pkg{parallel} package.\cr
-)\kbd{parallelWin}: parallel computations for network clusters or machines with multiple cores or CPUs (this is the only parallel implementation that works on Windows machines). A \sQuote{PSOCK} cluster is created with the \code{\link[parallel]{makeCluster}} function. When \code{fn.name="hydromod"} the evaluation of the objective function \code{fn} is done with the \code{\link[parallel]{clusterApply}} function of the \pkg{parallel} package. When \code{fn.name!="hydromod"} the evaluation of the objective function \code{fn} is done with the \code{\link[parallel]{parRapply}} function of the \pkg{parallel} package. 
} 
  \item{par.nnodes}{
OPTIONAL. Used only when \code{parallel!='none'} \cr
numeric, indicates the number of cores/CPUs to be used in the local multi-core machine, or the number of nodes to be used in the network cluster. \cr
By default \code{par.nnodes} is set to the amount of cores detected by the function \code{detectCores()} (\pkg{parallel} package)
} 
  \item{par.pkgs}{
OPTIONAL. Used only when \code{parallel='parallelWin'} \cr
list of package names (as characters) that need to be loaded on each node for allowing the objective function \code{fn} to be evaluated. By default \code{c("raster", "randomForest", "zoo")}.
} 
   \item{ntree}{
Numeric indicating the maximum number trees to grow in the Random Forest algorithm. The default value is set to 2000. This should not be set to too small a number, to ensure that every input row gets predicted at least a few times. If this value is too low, the prediction may be biased.
}   
   \item{na.action}{
A function to specify the action to be taken if NAs are found.  (NOTE: If given, this argument must be named.)
}   
  \item{write2disk}{
logical, indicates if the output merged raster layers and the training and avaluation datasets (two files each, one with time series and other with metadata) will be written to the disk. By default \code{write2disk=TRUE}
}
   \item{drty.out}{
Character with the full path to the directory where the final merged product will be exported as well as the training and evaluation datasets. Only used when \code{write2disk=TRUE}
}

  \item{use.pb}{
logical, indicates if a progress bar should be used to show the progress of the random forest computations (it might reduce a bit the performance of the computations, but it is useful to track if everything is working well). By default \code{use.pb=TRUE}
} 
  \item{verbose}{
logical, indicates if progress messages are to be printed. By default \code{verbose=TRUE}
} 
  \item{\dots}{
further arguments to be passed to the low level function \kbd{randomForest.default}.
}
}
%%\details{
%%  ~~ If necessary, more details than the description above ~~
%%}
\value{
It returns a \code{RasterStack} object with as many layers as time steps are present in \code{x}. Each one of the layers in the output object has the same spatial resolution and spatial extent as the \code{cov} argument.
}
\references{
\cite{Baez-Villanueva, O. M.; Zambrano-Bigiarini, M.; Beck, H.; McNamara, I.; Ribbe, L.; Nauditt, A.; Birkel, C.; Verbist, K.; Giraldo-Osorio, J.D.; Thinh, N.X. (2020). RF-MEP: a novel Random Forest method for merging gridded precipitation products and ground-based measurements, Remote Sensing of Environment, 239, 111610. \doi{10.1016/j.rse.2019.111606}. <https://authors.elsevier.com/c/1aKrd7qzSnJWL>}. \cr

\cite{Hengl, T., Nussbaum, M., Wright, M. N., Heuvelink, G. B., & Gr\"{a}ler, B. (2018). Random forest as a generic framework for predictive modeling of spatial and spatio-temporal variables. PeerJ, 6, e5518}.
}
\author{
Oscar M. Baez-Villanueva, \email{obaezvil@th-koeln.de} \cr
Mauricio Zambrano-Bigiarini, \email{mzb.devel@gmail}   \cr
Juan D. Giraldo-Osorio, \email{j.giraldoo@javeriana.edu.co}
}
%%\note{
%%  ~~further notes~~
%%}

%% ~Make other sections like Warning with \section{Warning }{....} ~

\seealso{
\code{\link[raster]{raster}}, \code{\link[raster]{stack}}, \code{\link[raster]{brick}}, \code{\link[raster]{resample}}, \code{\link[raster]{rotate}}, \code{\link[raster]{crop}}.
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
% Add one or more standard keywords, see file 'KEYWORDS' in the
% R documentation directory.
\keyword{manip}
%%\keyword{ ~kwd2 }
