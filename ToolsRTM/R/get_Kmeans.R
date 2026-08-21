

#' K-means Classification on a RasterStack
#' This function performs a k-means unsupervised classification of a stack of rasters. 
#' Number of clusters can be specified, as well as number of iterations and starting sets. An optional geographic weighting system can be turned on that constrains clusters to a geographic area, by including coordinates in the clustering. 
#' All variables are normalized before clustering is performed. 
#'https://rdrr.io/github/ozjimbob/ecbtools/man/raster.kmeans.html
#'
#' @param stack A RasterStack object, or a string pointing to the directory where all raster layers are stored. All rasters should have the same extent, resolution and coordinate system.  
#' @param k Number of clusters to classify to. 
#' @param iter.max Maximum number of iterations allowed
#' @param nstart  Number of random sets to be chosen. 
#' @param geo True/False - should geographic weighting be used? 
#' @param geo.weight A weighting multiplier indicating the strength of geographic weighting relative to the other variables. A value of 1 gives equal weight.

#'
#' @return Kmenans
#' @export
#'
#' 
getKmeans <- function(stack,k=12,iter.max=100,nstart=10,geo=T,geo.weight=1){
  
  # NB (raster->terra migration): "RasterStack" is a raster-package S4 class
  # that no longer applies to a terra SpatRaster -- accept a SpatRaster
  # directly instead (inherits("SpatRaster")), keeping the folder-of-tifs
  # path unchanged.
  if(inherits(stack, "SpatRaster")){
    stk=stack
  }else{
    stl=list.files(stack,full.names = TRUE,include.dirs = FALSE,pattern = "tif")
    stk=terra::rast(stl)
  }


  # NB: terra::as.data.frame() defaults to na.rm = TRUE (drops NA cells),
  # unlike raster::as.data.frame() which kept every cell by default -- the
  # xyFromCell() call right below assumes oDF has exactly ncell(stk) rows,
  # one per cell index 1:ncell in order, so na.rm must stay FALSE here to
  # preserve that row<->cell alignment (same reasoning applies below wherever
  # oDF/oDFk/of are indexed positionally against cell numbers).
  oDF=terra::as.data.frame(stk, na.rm = FALSE)
  oo=terra::xyFromCell(stk,1:(terra::ncol(stk) * terra::nrow(stk)))
  if(geo == T){
    oDF$x=oo[,1]
    oDF$y=oo[,2]
  }
  
  for(idx in 1:length(oDF)){
    oDF[,idx]=normalize(oDF[,idx])
  }
  
  if(geo == T){
    oDF$x = oDF$x * geo.weight
    oDF$y = oDF$y * geo.weight
  }
  
  oDFk=oDF
  oDFk$idx=1:length(oDFk[,1])
  oDFi=subset(oDFk,complete.cases(oDFk))
  oDF=subset(oDF,complete.cases(oDF))
  
  E <- kmeans(oDF, k, iter.max = iter.max, nstart = nstart)
  
  oDFi$cluster=E$cluster
  
  of=plyr::join(oDFk,oDFi,by="idx",type="left")
  
  # NB (raster->terra migration): a terra SpatRaster has no @ncols/@nrows/
  # @crs/@extent S4 slots (those were raster-package Raster* class slots) --
  # use the equivalent accessor functions instead.
  classify.m <- matrix(of$cluster, nrow=terra::nrow(stk),ncol=terra::ncol(stk), byrow=TRUE)
  raster.classify <- terra::rast(classify.m, crs=terra::crs(stk),
                                 extent=terra::ext(terra::xmin(stk), terra::xmax(stk), terra::ymin(stk), terra::ymax(stk)))
  return(raster.classify)
}

#' Normalize
#'
#' @param x 
#' @param low 
#' @param high 
#' Rescale a numeric vector to a given range
#'
#' @param x numeric vector to rescale.
#' @param low numeric. Lower bound of the target range. Default 0.
#' @param high numeric. Upper bound of the target range. Default 1.
#'
#' @return normalize values
#' @export
#'
normalize <- function(x,low=0,high=1){
  low+(x-min(x,na.rm=T))*(high-low)/(max(x,na.rm=T)-min(x,na.rm=T))
}
