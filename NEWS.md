NEWS/ChangeLog for RFmerge
--------------------------

# Changes in version 0.3-0   01-May-2026

        -) The parallel option is working now, thank to the asssistance of Codex (ChatGPT).
        -) Several samll enhancements and quality control changes.


# Changes in version 0.2-0   26-Jul-2023 (not public)
        o The package was completely checked and modified when deemed appropiate in order **to make it work with the 'terra' pacakge instead of the superseeded 'raster' package**. This was motivated by a tehnical course OScr Baez-Villanueva an I gave to SISSA in Buenos Aires (Argentina). However, the parallel option did not work due to problems with the use of SpatRast and SpatVec objects in the terra package. For this reason, this version remained only on Github, without being submitted to CRAN.  


# Changes in version 0.1-10  21-May-2020
        o Updated version after CRAN reported an unnecessary 'RFmerge-internal.Rd' file, which was removed for this release.


# Changes in version 0.1-9   12-May-2020
        o Updated version after CRAN reported an error in the vignette ("Invalid argument: 'cov' and 'mask' have different CRS !!"). this error was previously undetected because the package was tested with R-devel in combination with GEOS 3.5.1, GDAL 2.2.2, PROJ 4.9.2, while CRAN uses GEOS 3.8.1, GDAL 3.0.4, PROJ 7.0.0.


# Changes in version 0.1-8   28-Apr-2020 (not public)
        o Updated version after CRAN reported a problem with the vignette when tested against R-devel.



# Changes in version 0.1-7   27-Apr-2020 (not public)
        o 'rfmerge': sf::st_crs(x])[["epsg"]] was changed to 'sf::st_crs(x)$proj4string' for identifying the projection of a vectorial or raster object. This was due to a report from Edzer Pebema mentioning that the old call used a legacy 'crc' structure, which was changed in the new verison of the 'sf' package.
       

# Changes in version 0.1-6   07-Jan-2020
        o 'rmarkdown' was added to 'Suggests' after failing a new CRAN test.
        o the minimal vignette was modified to automatically set the 'parallel' argumento to 'paralell' or 'parallelWin' depending whether the OS is based on GNU/Linux or Windows, respectively.


# Changes in version 0.1-5   04-Jan-2020
        o 'RFmerge-package.Rd': now the example correctly uses 2 cores as maximum when package is tested with the "--run-donttest" option


# Changes in version 0.1-4   04-Jan-2020
        o DOI was added to the CITATION file and all the references within the package (the article was available online in January 2nd).
        o A maximum of two cores is used in all examples (fixed, with respect to the previous version).
        o Improved vignettes (minor details).


# Changes in version 0.1-3   20-Dec-2019
        o 'RFmerge.R': new argument 'write2disk' to allow the user to choose whether the output merged layers should be written to disk or not. Default value is 'FALSE'. 
        o 'RFmerge.R': 'installed.packages()' was replaced by 'find.package', for increasing efficiency to find out if a named package is installed or not.
        o Examples are now run within a 'donttest{}' piece of code instead of a '\dontrun{}' one.
        o A maximum of two cores is used in all examples, even if the number of cores available is much larger than that (required by CRAN)
        o RFmerge.Rd: '\value' tag was added, explaining the returning object of the function.
        o Hengl et al. (2018) was added as a reference' did not match any file(s) known to git.


# Changes in version 0.1-2   18-Dec-2019
        o The package was built using the '--compact-vignettes="gs+qpdf"'  to avoid the following WARNIG from CRAN: "'gs+qpdf' made some significant size reductions".


# Changes in version 0.1-1   18-Dec-2019
        o Full PDF vignette was replaced by a  minimal RMarkdown version to avoid dependencies from rgdal, hydroTSM, and hydroGOF.


# Changes in version 0.1-0   15-Dec-2019
        o First release
        o It depends on R >= 3.5.0 because serialized objects (.rds data files) in serialize/load version 3 cannot be read in older versions
        o RMarkdown vignette was replaced by a sttic PDF  one to avoid dependencies from rgdal, hydroTSM, and hydroGOF.

