# RFmerge
[![Research software impact](http://depsy.org/api/package/cran/hydroTSM/badge.svg)](http://depsy.org/package/r/hydroTSM) [![Build Status](https://travis-ci.org/hzambran/hydroTSM.svg?branch=master)](https://travis-ci.org/hzambran/hydroTSM)

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

## Reporting bugs, requesting new features

If you find an error in some function, or want to report a typo in the documentation, or to request a new feature (and wish it be implemented :) you can do it [here](https://github.com/hzambran/hydroTSM/issues)


## Citation 
```{r}
citation("RFmerge")
```

To cite RFmerge in publications use:

>  Mauricio Zambrano-Bigiarini. hydroTSM: Time Series Management, Analysis and Interpolation for Hydrological Modelling. R package version 0.5-1. URL https://hzambran.github.io/hydroTSM/. DOI:10.5281/zenodo.839864

> Baez-Villanueva, O. M.; Zambrano-Bigiarini, M.; Beck, H.; McNamara, I.; Ribbe, L.; Nauditt, A.; Birkel, C.; Verbist, K.; Giraldo-Osorio, J.D.; Thinh, N.X. (2019). RF-MEP: a novel Random Forest method for merging gridded precipitation products and ground-based measurements, Remote Sensing of Environment, (accepted).

A BibTeX entry for LaTeX users is

>  @Manual{hydroTSM,  
>    title = {hydroTSM: Time Series Management, Analysis and Interpolation for Hydrological Modelling},  
>    author = {{Mauricio Zambrano-Bigiarini}},  
>    note = {R package version 0.5-1},  
>    url = {https://hzambran.github.io/hydroTSM/},  
>    doi = {10.5281/zenodo.839864},  
>  }


## Vignette 
[Here](https://cran.r-project.org/web/packages/hydroTSM/vignettes/hydroTSM_Vignette-knitr.pdf) you can find an introductory vignette showing the use of several hydroTSM functions.



## Related Material 

* *R: a statistical environment for hydrological analysis* (**EGU-2010**)  [abstract](http://meetingorganizer.copernicus.org/EGU2010/EGU2010-13008.pdf), [poster](http://www.slideshare.net/hzambran/egu2010-ra-statisticalenvironmentfordoinghydrologicalanalysis-9095709).

* *Using R for analysing spatio-temporal datasets: a satellite-based precipitation case study* (**EGU-2017**) [abstract](http://meetingorganizer.copernicus.org/EGU2017/EGU2017-18343.pdf), [poster](https://doi.org/10.5281/zenodo.570145).



## See Also 

* [hydroGOF: Goodness-of-fit functions for comparison of simulated and observed hydrological time series](https://github.com/hzambran/hydroGOF).

* [hydroPSO: Model-independent Particle Swarm Optimisation (PSO) for environmental/hydrological models](https://github.com/hzambran/hydroPSO).
