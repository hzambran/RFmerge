# All the code of this function was adapted from the function 'buffer.dist.R' in package GSIF (https://cran.r-project.org/package=GSIF); developed by Tomislav Hengl, Bas Kempen and Gerard Heuvelink with license GPL>=2; and then slighly modified.
# This was the only function needed from that package, but the orignal funcion only worked with the old raster package. Threfore, the original function was slighly modified
# and used in the 'RFmerge' package just to avoid all the dependencies required to install GSIF package
# 

# .buffer_dist <- function(observations, predictionDomain, classes, width, ...){

#   if(missing(width)) 
#     width <- sqrt(sp::areaSpatialGrid(predictionDomain)) 

#   if(!length(classes)==length(observations))
#     stop("Length of 'observations' and 'classes' does not match.")

#   ## remove classes without any points:
#   xg <- summary(classes, maxsum=length(levels(classes)))
#   selg.levs = attr(xg, "names")[xg > 0]

#   if(length(selg.levs)<length(levels(classes))){
#     fclasses <- as.factor(classes)
#     fclasses[which(!fclasses %in% selg.levs)] <- NA
#     classes <- droplevels(fclasses)
#   }

#   ## derive buffer distances
#   s <- list(NULL)
#   for(i in 1:length(levels(classes))){
#     #s[[i]] <- terra::distance(terra::rasterize(observations[which(classes==levels(classes)[i]),1]@coords, y=terra::rast(predictionDomain)), width=width, ...)
#     s[[i]] <- raster::distance(raster::rasterize(observations[which(classes==levels(classes)[i]),1]@coords, y=raster::raster(predictionDomain)), width=width)
#   }

#   s <- s[sapply(s, function(x){!is.null(x)})]
#   s <- terra::rast(s)
#   s <- as(s, "SpatialPixelsDataFrame")
#   s <- s[predictionDomain@grid.index,]

#   return(s)

# } # '.buffer_dist' END


# obs.t <- vect(observations)
# cov1 <- lsample
# crs(obs.t) <- crs(cov1)

.buffer_dist <- function(observations, predictionDomain) {

 s2 <- list(NULL)
  for(i in 1:length(levels(classes))){
    s2[[i]] <- terra::distance(cov1, obs.t[i], rasterize=TRUE)
  }

} # '.buffer_dist' END

#s[[1]] - raster::raster(s2[[1]])

