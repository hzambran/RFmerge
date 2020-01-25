%%File ValparaisoPPts.Rd
%% Part of the RFmerge R package, https://github.com/hzambran/RFmerge ; 
%%                                https://cran.r-project.org/package=RFmerge
%% Copyright 2019-2020 Mauricio Zambrano-Bigiarini, Oscar M. Baez-Villanueva
%% Distributed under GPL 3 or later

\name{ValparaisoPPts}
\docType{data}
\alias{ValparaisoPPts}
\title{
Daily Precipitation Time Series for Valparaiso Region (Chile)
}
\description{
Daily time series for the year 1983 on 34 rain gauges of the Valparaiso region (Chile), with more than 90\% of days with information (without missing values) 

}
\usage{
data(ValparaisoPPts)
}

\details{
Daily time series of ground-based daily precipitation for 1900-2018 were downloaded from a dataset of 816 rain gauges from the Center of Climate and Resilience Research (CR2; \url{http://www.cr2.cl/datos-de-precipitacion/}). \cr
The 34 stations in this dataset wer selected because they had less than 10\% of missing values in year 1983.
}
\format{
A zoo object with 34 columns (one for each rain gauge) and 365 rows (one for each day in 1983). \code{colnames(ValparaisoPPts)} must coincide with the \var{ID} column in \var{ValparaisoPPgis}.
}
\source{
The \href{http://www.cr2.cl/datos-de-precipitacion/}{CR2 dataset} unifies individual datasets provided by Direcci\'{o}n General de Aguas (DGA) and Direcci\'{o}n Meteorol\'{o}gica de Chile (DMC), the Chilean water and meteorological agencies, respectively. \cr
These data are intended to be used for research purposes only, being distributed in the hope that it will be useful, but WITHOUT ANY WARRANTY.
}

%%\references{
%% ~put references to the literature/web site here ~
%%}
\keyword{datasets}
