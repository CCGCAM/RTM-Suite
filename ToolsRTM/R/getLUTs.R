
#' Get LUTs for radiative transfer models
#'
#' @param inputs  a LUT table with variables with specific ranges.
#' @param nLUT  the number of rows for the LUT
#' @param dependencies nLUT  the number of rows for the LUT
#' @param setseed setseed  a seed number to contral random process. By default is FALSE or null. Please add the number.
#' @return a dataframe with all parameters
#' @export
#'
#'
getLUTs<-function(inputs=NULL, nLUT=100, dependencies='Car',setseed = 123){

  if (!(any(c("data.frame", "matrix") %in%  class(inputs)))) {
    stop("Data needs to be in correct format. data.frame or matrix are permitted.")
  }


  if ( is.null(setseed) | setseed ==FALSE ){
    set.seed(Sys.time())
  } else {
    set.seed(setseed)
  }


  if (is.null(dependencies)){

  } else {

    trait_dep = dependencies
  }


  if (is.null(nLUT)){
    message('number of varitions for each input is fixed to 100')
    nLUT = 100
  }

  var.list = list()
  inputs_names<- names(inputs)

  for (i in  c(1:length(inputs_names))) {

    trait <- inputs_names[i]
    table.sb<-inputs[trait]
    var.list[[trait]] <- stats::runif(nLUT,min = min(table.sb),max=max(table.sb))
  }
  if (trait_dep %in% names(inputs)){
    cat(' get correlation with chlorohyll ....')
    var.list[[trait_dep]] <- ToolsRTM::correlatedValue(x=var.list[['Cab']]/4, r=.8)
  }  else {
    cat(' LUT without correaltion with Cab....')
  }
  LUT <- as.data.frame(do.call(cbind,var.list ))

  return(LUT)

}
