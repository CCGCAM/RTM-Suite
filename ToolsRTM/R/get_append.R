#' Append  tables
#'
#' @param input  a list of tables
#'
#' @return append table
#' @export 
#'
#' 
get_append<-function(input = NULL){
  
  file_list = input
  
  for (file in file_list){
    
    # if the merged dataset doesn't exist, create it
    if (!exists("db")){
      db <- read.table(file, header=TRUE, sep=",")
    }
    
    # if the merged dataset does exist, append to it
    if (exists("db")){
      temp_ <-read.table(file, header=TRUE, sep=",")
      db<-rbind(db, temp_)
      rm(temp_)
    }
  }
  return(db)
}

