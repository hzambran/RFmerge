%%File ValparaisoPPts.Rd
%% Part of the RFmerge R package, https://github.com/hzambran/RFmerge ; 
%%                                https://cran.r-project.org/package=RFmerge
%% Copyright 2019-2020 Mauricio Zambrano-Bigiarini, Oscar M. Baez-Villanueva
%% Distributed under GPL 3 or later

\name{ValparaisoSHP}
\docType{data}
\alias{ValparaisoSHP}
\title{
Administrative borders of the Valparaiso Region (Chile)
}
\description{
Simple feature polygon collection with 1 feature and 10 fields
}
\usage{
data(ValparaisoSHP)
}

\details{
Projection: EPSG:4326
}

\format{
\link[sp]{SpatialPolygonsDataFrame-class}. \cr

The fields stored in the \code{@data} slot of this object are: \cr 

-) \var{NOM_REG}: Name of the administrative region \cr
-) \var{NOM_PROV}: Name of the administrative province \cr
-) \var{COD_COM}: ID of the administrative comune \cr 
-) \var{NOM_COM}: Name of the administrative comune \cr
-) \var{COD_REGI}: Numeric code of the administrative region \cr 
-) \var{SUPERFICIE}: Spatial area within the administrative borders of the region, [km2] \cr 
-) \var{POBLAC02}: Probably, it corresponds to the number of inhabitants of the region according to the 2002 census. \cr 
-) \var{POBL2010}: Probably, it corresponds to the number of inhabitants of the region according to the 2010 census. \cr 
-) \var{SHAPE_Leng}: Total length of the administrative border of the region, [m] \cr 
-) \var{SHAPE_Area}: Spatial area within the administrative borders of the region, [m2]

}
\source{
Originally downloaded from the Biblioteca del Congreso Nacional de Chile and then reprojected into geographic coordiantes (EPSG:4326). Last accessed [March 2016]). \cr
These data are intended to be used for research purposes only, being distributed in the hope that it will be useful, but WITHOUT ANY WARRANTY.
}

\note{
The orginal file is not longer available at he Biblioteca del Congreso Nacional de Chile
}

%%\references{
%% ~put references to the literature/web site here ~
%%}
\keyword{datasets}
