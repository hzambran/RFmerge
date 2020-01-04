%%File ValparaisoPPgis.Rd
%% Part of the RFmerge R package, https://github.com/hzambran/RFmerge ; 
%%                                https://cran.r-project.org/package=RFmerge
%% Copyright 2019-2020 Mauricio Zambrano-Bigiarini, Oscar M. Baez-Villanueva
%% Distributed under GPL 3 or later

\name{ValparaisoPPgis}
\docType{data}
\alias{ValparaisoPPgis}
\title{
Spatial location of rain gauges in the Valparaiso region (Chile)
}
\description{
Spatial location of the 34 rain gauges with daily precipitation for the Valparaiso region (dataset \var{'ValparaisoPPts'}), Chile, with more than 70\% of days with information (without missing values) 
}
\usage{
data(ValparaisoPPgis)
}

\details{
Projection: EPSG:4326
}
\format{
A data.frame with seven fields: \cr
*) 'ID             : identifier of each gauging station. \cr
*) 'STATION_NAME'  : name of the gauging station. \cr
*) 'lon'           : easting coordinate of the gauging station, EPSG:4326. \cr
*) 'lat'           : northing coordinate of the gauging station, EPSG:4326. \cr 
*) 'ELEVATION'     : elevation of the gauging station, [m a.s.l.]. \cr
*) 'BASIN_ID'      : identifier of the subbasin in which the gauging station s located. \cr
*) 'BASIN_NAME'    : name of the subbasin in which the gauging station s located.
}
\source{
Downloaded ('Red de Control Meteorologico') from the web site of the Confederacion Hidrografica del Ebro (CHE) \url{http://www.chebro.es/} (original link http://oph.chebro.es/ContenidoCartoClimatologia.htm, last accessed [March 2008]), and then the name of 7 selected fields were translated into English language.  \cr

These data are intended to be used for research purposes only, being distributed in the hope that it will be useful, but WITHOUT ANY WARRANTY.
}

%%\references{
%% ~put references to the literature/web site here ~
%%}
\keyword{datasets}
