#' getPredicts
#'
#' @param model a keras model
#' @param type.model keras model type; the options are: 'CNN'; 'Hidden-layers'
#' @param data the dataset in matrix format for predicting depVar
#' @param data.trans PCA and normalized, options are avalaible
#' @param data.Y  the Y variable: Avalaible options are: 'NULL' in real cases or a vector or matrix or data.frame
#' @param depVar variable name 
#' @param scaler.depVar Scaler for Y 
#' @return predictions
#' @export 
#'
#' 
getPredicts<-function(model=NULL, type.model='CNN',data=NULL,data.trans=NULL ,
                      data.Y=NULL, depVar='Cab', scaler.depVar = NULL) {
  

  if (is.null(data.trans)){
    stop('PLease insert a tranformation for applying: PCA or preProcess  ...')
  
  }
  
  # 
  # if (is.null(scaler.depVar)){
  #   stop('PLease insert the scaler for the model  ...')
  #   
  # }
  
  if (is.null(scaler.depVar)){
    depVar.trans=FALSE
  } else {
    depVar.trans=TRUE
  }
  
  
  if (is.null(depVar)){
    stop('PLease insert the name variable to predict ...')
  }
  
  if (type.model == "CNN" | type.model == "Hidden-layers"){
    #cat('get predictions .....')
  } else{
    stop('PLease insert a type model: CNN; Hidden-layers  or All options are avalaible')
  }
  if (is.null(model)){
    stop('PLease insert a model for predicting  ...')
  }
  if (is.null(model)){
    stop('PLease insert a model for predicting  ...')
  }
  if (class(data)[1] != "matrix" | is.null(data)){
    stop('PLease insert a matrix with inputs ...')
  } 
  
  if (!inherits(data.Y, "numeric")){
   
    data.Y.to<-data.Y[,depVar]
  } else{
    data.Y.to <- data.Y
  }
  
  
  if (type.model == 'Hidden-layers'){
    
    
    pred.model.3hlayers <- as.numeric(predict(model, data))
    pred.model.3hlayers<-as.matrix(pred.model.3hlayers)
    colnames(pred.model.3hlayers)<-depVar
    if (depVar.trans == TRUE){
      
      pred.model.3hlayers <- getReverse.trans(preProc=scaler.depVar,data=pred.model.3hlayers)

    }
    if (is.null(data.Y.to)){
        df.plot<-data.frame(data,predicted.3hlayer=pred.model.3hlayers[,1],predicted.cnn=NA)
    } else {
          df.plot<-data.frame(depVar=data.Y.to,predicted.3hlayer=pred.model.3hlayers[,1],predicted.cnn=NA)
    }
    
    
  } else if (type.model == 'CNN'){
    
    data.cnn <- array_reshape(data, c(nrow(data), ncol(data), 1))
    pred.model.cnn <- predict(model,data.cnn)
    pred.model.cnn<-as.matrix(pred.model.cnn)
    colnames(pred.model.cnn)<-depVar
    if (depVar.trans == TRUE){
      pred.model.cnn <- getReverse.trans(preProc=scaler.depVar,data=pred.model.cnn)
    }
    if (is.null(data.Y.to)){
      df.plot<-data.frame(data,predicted.3hlayer=NA,predicted.cnn=pred.model.cnn[,1])
    } else {
      df.plot<-data.frame(depVar=data.Y.to,predicted.3hlayer=NA,predicted.cnn=pred.model.cnn[,1])
    }
  
    
  } else if (type.model == 'All'){
    
    pred.model.3hlayers <- predict(model, data)
    pred.model.3hlayers<-as.matrix(pred.model.3hlayers)
    colnames(pred.model.3hlayers)<-depVar
    data.cnn <- array_reshape(data, c(nrow(data), ncol(data), 1))
    pred.model.cnn <- predict(model,data.cnn)
    pred.model.cnn<-as.matrix(pred.model.cnn)
    colnames(pred.model.cnn)<-depVar
    
    if (depVar.trans == TRUE){
      pred.model.3hlayers <- getReverse.trans(preProc=scaler.depVar,data=pred.model.3hlayers)
      pred.model.cnn <- getReverse.trans(preProc=scaler.depVar,data=pred.model.cnn)
    }
    if (is.null(data.Y.to)){
      df.plot<-data.frame(data,predicted.3hlayer=pred.model.3hlayers,predicted.cnn=pred.model.cnn)
      inputs_depVar= names(df.plot)[grep("pred.", names(df.plot))]
      df.plot$depVar_pred_avg=rowMeans(df.plot[,c(inputs_depVar)], na.rm=TRUE)
      df.plot$depVar_pred_sd<-apply(df.plot[,c(inputs_depVar)], 1, sd)
    } else {
      df.plot<-data.frame(depVar=data.Y.to,predicted.3hlayer=pred.model.3hlayers,predicted.cnn=pred.model.cnn)
      inputs_depVar= names(df.plot)[grep("pred.", names(df.plot))]
      df.plot$depVar_pred_avg=rowMeans(df.plot[,c(inputs_depVar)], na.rm=TRUE)
      df.plot$depVar_pred_sd<-apply(df.plot[,c(inputs_depVar)], 1, sd)
    }
   
  
  }
  
  
  return(df.plot)
  
}

#' Reverse a caret preProcess transformation and round the result
#'
#' @param preProc a "preProcess" object created by \code{caret::preProcess()}, describing the transformation to reverse.
#' @param data a matrix. Transformed data to convert back to its original scale.
#' @param digits integer. Number of decimal places to round the reversed values to. Default 10.

#' @return reverse trans
#' @export
#'

getReverse.trans <- function(preProc=NULL, data=NULL, digits = 10) {
  
  stopifnot(class(preProc) == "preProcess")
  stopifnot(class(data)[1] == "matrix")
  
  if (is.null(preProc)){
    stop('PLease insert the scaler for the model  ...')
    
  }
  if (is.null(data)){
    stop('PLease insert a matrix to reverse the applied transfromation  ...')
    
  }
  if (is.null(data)){
    stop('PLease insert a matrix to reverse the applied transfromation  ...')
    
  }
  df.to<-as.data.frame(data)
  
  nc <- ncol(df.to); nr <- nrow(df.to)
  df.mean <- t(replicate(nr, preProc$mean))
  df.std <- t(replicate(nr, preProc$std))
  boundaries  <- preProc$rangeBounds
  df.max <- t(replicate(nr, preProc$ranges[2,]))
  df.min <- t(replicate(nr, preProc$ranges[1,]))
  
  if(sum(!is.na(match(c("center", "scale"), 
                      names(preProc$method)))) == 2) {
    df.transformed <- df.to * df.std + df.mean
    
  } else if(sum(!is.na(match("center", 
                             names(preProc$method)))) == 1) {
    df.transformed <- df.to + df.mean
    
  } else if(sum(!is.na(match("scale", 
                             names(preProc$method)))) == 1) {
    df.transformed <- df.to * df.std
  } else {
    df.transformed <- (df.to-boundaries[1])/(boundaries[2]-boundaries[1])*(df.max - df.min) + df.min
  }
  
  return(round(df.transformed, digits))
}




