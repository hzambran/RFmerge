# File RFmerge.R
# Part of the RFmerge R package, https://github.com/hzambran/RFmerge ; 
#                                https://cran.r-project.org/package=RFmerge
# Copyright 2019-2020 Mauricio Zambrano-Bigiarini, Oscar M. Baez-Villanueva
# Distributed under GPL 3 or later

.onAttach <- function(libname, pkgname) {

  packageStartupMessage("(C) 2019 M. Zambrano-Bigiarini and Oscar M. Baez-Villanueva (GPL >=3 license)\n",
                         "Type 'citation('RFmerge')' to see how to cite this package")
  invisible()
    
} # '.onAttach' END

