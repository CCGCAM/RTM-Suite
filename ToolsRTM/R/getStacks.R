
#' @title Get TIFF from a folder and generate the stack. 
#'
#' @param rasterFiles  path with the Netcdf
#' @param frequency Daily
#' @param bands  a vector with the names of the inputs of the NetCDF
#' @param output  path of the outputs
#' 
#' @description
#' This function is for Sentinel-2 data using the names B1,B2 ...and SCL
#' 
#' @return a stack
#' @export
#' @examples
#' \dontrun{
#' rasterFiles <- "path/to/netcdf/files"
#' frequency <- "Daily"
#' bands <- c("B1", "B2", "B3", "SCL")
#' output <- "path/to/output/folder"
#' stack <- getStacks(rasterFiles, frequency, bands, output)
#' }
#'

getStacks<-function(rasterFiles=NULL, frequency='Daily', bands=NULL,output=NULL){
  options(warn=-1) ###avoid warnings
  #getwd()
  files = list.files(rasterFiles,pattern="day_*", full.names=F)
  dates<-as.Date(substr(files,5,14))
  dates.unique<-unique(dates)



  #### create the Stack of images in all files for "
  progress_bar = txtProgressBar(min=0, max=length(dates.unique), style = 3, char="=")
  for (k in c(1:length(dates.unique))){
    setTxtProgressBar(progress_bar, k)
    path_out=paste(rasterFiles,'/day-',dates.unique[k],'',sep='')
    files = list.files(path_out,pattern="*.tif$", full.names=TRUE)
    #print(files)
    #rs <- stack(files)
    ## by order B1-B2-B3-B4-B5-B6-B7-B8-B8A-B9-B11-B12
    # Extract the pattern from each file name
    pattern <- stringr::str_extract(string = files, pattern = paste(bands, collapse = "|"))

    # Sort the files based on the pattern vector
    sorted_files <- files[order(match(pattern, bands))]

    # Print the sorted file names
    #print(sorted_files)
    # Load the first file to initialize the raster stack
    rs <- raster(sorted_files[1])

    # Loop over the remaining files and add them to the raster stack
    for (i in 2:length(sorted_files)) {
      rs <- stack(rs, raster(sorted_files[i]))
    }
    names(rs) <- bands
    #plot(rs)
    path_out<-paste(output,'/Stacks',sep="")
    ifelse(!dir.exists(path_out), dir.create(path_out), FALSE)
    raster.file<-paste(path_out,'/SE2A-',dates[k],'',sep="")
    writeRaster(rs, raster.file, "GTiff", overwrite=TRUE,bylayer=F)


}
  close(progress_bar)
  m_final<-message('Stacks files were generated sucessfully')
  return(m_final)


}
