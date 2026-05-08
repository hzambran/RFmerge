# Spatial location of rain gauges in the Valparaiso region (Chile)

Spatial location of the 34 rain gauges with daily precipitation for the
Valparaiso region (dataset `'ValparaisoPPts'`), Chile, with more than
70% of days with information (without missing values)

## Usage

``` r
data(ValparaisoPPgis)
```

## Details

Projection: EPSG:4326

## Format

A data.frame with seven fields:  
\*) 'ID : identifier of each gauging station.  
\*) 'STATION_NAME' : name of the gauging station.  
\*) 'lon' : easting coordinate of the gauging station, EPSG:4326.  
\*) 'lat' : northing coordinate of the gauging station, EPSG:4326.  
\*) 'ELEVATION' : elevation of the gauging station, \[m a.s.l.\].  
\*) 'BASIN_ID' : identifier of the subbasin in which the gauging station
s located.  
\*) 'BASIN_NAME' : name of the subbasin in which the gauging station s
located.

## Source

Downloaded ('Red de Control Meteorologico') from the web site of the
Confederacion Hidrografica del Ebro (CHE) <http://www.chebro.es/>
(original link http://oph.chebro.es/ContenidoCartoClimatologia.htm, last
accessed \[March 2008\]), and then the name of 7 selected fields were
translated into English language.  

These data are intended to be used for research purposes only, being
distributed in the hope that it will be useful, but WITHOUT ANY
WARRANTY.
