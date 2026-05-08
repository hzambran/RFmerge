# Daily Precipitation Time Series for Valparaiso Region (Chile)

Daily time series for the year 1983 on 34 rain gauges of the Valparaiso
region (Chile), with more than 90% of days with information (without
missing values)

## Usage

``` r
data(ValparaisoPPts)
```

## Details

Daily time series of ground-based daily precipitation for 1900-2018 were
downloaded from a dataset of 816 rain gauges from the Center of Climate
and Resilience Research (CR2;
<https://www.cr2.cl/datos-de-precipitacion/>).  
The 34 stations in this dataset wer selected because they had less than
10

## Format

A zoo object with 34 columns (one for each rain gauge) and 365 rows (one
for each day in 1983). `colnames(ValparaisoPPts)` must coincide with the
`ID` column in `ValparaisoPPgis`.

## Source

The [CR2 dataset](https://www.cr2.cl/datos-de-precipitacion/) unifies
individual datasets provided by Dirección General de Aguas (DGA) and
Dirección Meteorológica de Chile (DMC), the Chilean water and
meteorological agencies, respectively.  
These data are intended to be used for research purposes only, being
distributed in the hope that it will be useful, but WITHOUT ANY
WARRANTY.
