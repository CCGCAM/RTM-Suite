#' getSplitData for ML models with no messages
#'
#' @param data A dataframe with inputs and variable to predic
#' @param inputs  the a vector with variables names for training and testing the model
#' @param depVar name of the variable to predict
#' @param data.trans a data.transformation method, options are: 'PCA','preProcess',
#' @param method.preProcess the data.transformation method for preProcess, data.transformation are: 'Normalize', 'YeoJohnson','BoxCox', Standarize', 'Center','Scale', and 'PCA'
#' @param depVar.trans a boolean variable for applying data transformation in Y variable, options are: TRUE or FALSE.
#' @param prop.split a vector with proportion for spliting the dataset. prop.split =c(0.8,02) will be used as default.
#'
#' @return the training and testing dataset with plot and scalers
#' @export
#'
#' 
#' 

getSplitData_noMessages<-function(data=NULL, depVar='Cab',inputs=NULL, 
                       data.trans=NULL,method.preProcess=NULL,depVar.trans=NULL,
                       prop.split=NULL) {

  
  
  if (is.null(data.trans)){
    stop('PLease insert a transformation: PCA, preProcess are avalaible.')
  }
  if (data.trans == 'PCA'){
    #message('PCA method:')
  } else if (data.trans == 'preProcess'){
    
    #message('preProcess method:')
  }
  
  if (is.null(depVar.trans)){
    depVar.Trans=FALSE
  }
  
  if (is.null(depVar)){
    stop('PLease insert a vector with the variable to predict ...')
  }

  if (is.null(inputs)){
    stop('PLease insert a vector with the inputs  ...')
  }
  
  if (is.null(method.preProcess)){
    method.preProcess = 'Normalize'
  }
  
  if (is.null(method.preProcess) & (data.trans == 'preProcess')){
    
    method.preProcess = 'Normalize'
  }
  
  if (data.trans == 'PCA'){
  
    method.preProcess = 'NA'
  }
  
  methods<-c('Normalize', 'YeoJohnson','BoxCox', 'Standarize', 'Center','Scale', 'PCA')
  
  if ((method.preProcess %in% methods) == TRUE){
   
    #message(paste(method.preProcess,' will be used for splitting the data',sep=''))
    
  } else {
   
    #message(paste('PCA method will be used for splitting the data',sep=''))
  }
  
  if (!inherits(data, 'data.frame') | is.null(data)){
    
    stop('PLease insert a dataframe with the depend variable and inputs ...')
    
  } else {
    
   inputs_to<-which(names(data)  %in% inputs)
   inputs_to<-names(data)[inputs_to]
   data<-data[,c(depVar,inputs_to)]
   ## Nplit the data
   
   if (is.null(prop.split)){
     split.train <- 0.8
     split.val <- 0.2
     #message('split dataset: 80% for training and 20% for testing ...')
   } else{
     split.train <- prop.split[1]
     split.val <- prop.split[2]
     #message(paste('split dataset: ',round(split.train*100),'% for training and ',round(split.val*100),'% for testing ...',sep=''))
   }
   ind <- sample(2, nrow(data), replace=TRUE, prob = c(split.train,split.val))
   inputs.to.include = names(data)[-1] ## for rfl bands
   #cat('inputs are: ',inputs.to.include)
   #cat('/')
   # change the data to matrix
   LUT.keras<-data[,c(depVar,inputs.to.include)]
   #skimr::skim(LUT.keras)
   
   ### Clean data
   #lapply(LUT.keras, function(x) sum(is.na(x))) |> str()
   LUT.keras <- na.omit(LUT.keras)
   names_to<-c(depVar,inputs.to.include)
  }
  
 

  
 
  
  if (data.trans == 'preProcess'){
    # calculate the pre-process parameters from the dataset
    LUT.to <-cbind(LUT.keras[,depVar],LUT.keras[,inputs.to.include])
    colnames(LUT.to)<-c(depVar,inputs.to.include)
    
    data.Xtrain_ <- LUT.to[ind==1, 2:dim(LUT.to)[2]]
    data.Xval_ <- LUT.to[ind==2, 2:dim(LUT.to)[2]]
    dim(data.Xval_)
    
    data.Ytrain <- round(LUT.to[ind==1, 1],4)
    length(data.Ytrain)
    data.Yval <- round(LUT.to[ind==2, 1],4)
    length(data.Yval)
    
    if (method.preProcess == 'Normalize'){
      #print('Normalize')
      preprocess.scalar <- caret::preProcess(data.Xtrain_, rangeBounds = c(0,1),method=c('range'))
     
    } else if (method.preProcess == 'YeoJohnson'){
      #print('YeoJohnson')
      preprocess.scalar <- caret::preProcess(data.Xtrain_, method=c('YeoJohnson'))
    } else if (method.preProcess == 'BoxCox'){
      #print('BoxCox')
      preprocess.scalar <- caret::preProcess(data.Xtrain_, method=c('BoxCox'))
    } else if (method.preProcess == 'Standarize'){
      #print('Standarize')
      preprocess.scalar <- caret::preProcess(data.Xtrain_, method=c('center', 'scale'))
    } else if (method.preProcess == 'Center'){
      #print('Center')
      preprocess.scalar <- caret::preProcess(data.Xtrain_, method=c('center'))
    } else if (method.preProcess == 'Scale'){
      #print('Scale')
      preprocess.scalar <- caret::preProcess(data.Xtrain_, method=c('scale'))
    } else if (method.preProcess == 'PCA'){
      #print('PCA')
      preprocess.scalar <- caret::preProcess(data.Xtrain_, method=c('center', 'scale', 'pca'))
     
    }
    data.Xtrain<-as.matrix(predict(preprocess.scalar, data.Xtrain_))
    
    data.Xval<-as.matrix(predict(preprocess.scalar, data.Xval_))
    
    M <- cor(data.Xtrain)
    # tl.cex=0.5 is only legible up to a few dozen variables, and past ~80
    # variables (e.g. ~230 PRISMA bands) the tiles themselves become too
    # small to show any visible colour/structure even with labels dropped --
    # at that point the plot is just noise, so skip it entirely rather than
    # render something unreadable.
    n.vars <- ncol(M)
    if (n.vars <= 80) {
      plot.cor<-corrplot::corrplot(M, method = 'square', order = 'FPC',
                                   type = 'lower',
                                   tl.cex=min(0.5, 15 / n.vars), tl.pos = 'lt',
                                   tl.col = 'black', cl.ratio=0.4,
                                   diag = FALSE)
    } else {
      plot.cor <- NULL
    }


    
    split.database<- list('Xtrain' = data.Xtrain,'Ytrain' =data.Ytrain,
                          'Xval' =data.Xval,'Yval' = data.Yval,
                          'Scalar.train' = preprocess.scalar,
                          'Scalar.Ytrain' = NA,
                          'plot.train' = plot.cor)
    
    if (depVar.trans == TRUE){
      
      data.Ytrain.to<-as.matrix(data.Ytrain)
      colnames(data.Ytrain.to)<-depVar
      preprocess.scalarY <- caret::preProcess(data.Ytrain.to,rangeBounds = c(0,1),method=c('range'))
      
      data.Ytrain <- predict(preprocess.scalarY, data.Ytrain.to)
      LUT.to.Train <-cbind(data.Ytrain.to, data.Ytrain,data.Xtrain)
      colnames(LUT.to.Train)<-c(depVar,paste(depVar,'_Tra',sep=''),colnames(data.Xtrain))
      LUT.to.Train<-as.matrix(LUT.to.Train)
      
      data.Ytrain <- round(LUT.to.Train[,2],4)
      
      data.Yval.to<-as.matrix(data.Yval)
      colnames(data.Yval.to)<-depVar
      
      data.Yval <- round(predict(preprocess.scalarY, data.Yval.to),4)
      LUT.to.Val <-cbind(data.Yval.to, data.Yval,data.Xval)
      colnames(LUT.to.Val)<-c(depVar,paste(depVar,'_Tra',sep=''),colnames(data.Xval))
      
      data.Yval <- round(LUT.to.Val[,2],4)
      
      
      split.database<- list('Xtrain' = data.Xtrain,'Ytrain' =data.Ytrain,
                            'Xval' =data.Xval,'Yval' = data.Yval,
                            'Scalar.train' = preprocess.scalar,
                            'Scalar.Ytrain' = preprocess.scalarY,
                            'plot.train' = plot.cor)
      
    } 
    
    
  } else if (data.trans == 'PCA' ){
    
  
    LUT.to <-cbind(LUT.keras[,depVar],LUT.keras[,inputs.to.include])
    colnames(LUT.to)<-c(depVar,inputs.to.include)
    
    data.Xtrain_ <- LUT.to[ind==1, 2:dim(LUT.to)[2]]
    data.Xval_ <- LUT.to[ind==2, 2:dim(LUT.to)[2]]
    
    data.Ytrain <- round(LUT.to[ind==1, 2],4)
    data.Yval <- round(LUT.to[ind==2, 2],4)
    
    
    LUT.to.xtrain <-prcomp(data.Xtrain_, scale = TRUE)
    LUT.to.xVal <-prcomp(data.Xval_, scale = TRUE)
    
    plot.pca.1<-factoextra::fviz_eig(LUT.to.xtrain,)
    print(plot.pca.1)
    
    plot.pca.2<-factoextra::fviz_pca_var(LUT.to.xtrain,
                                         col.var = 'contrib', # Color by contributions to the PC
                                         gradient.cols = c('#00AFBB', '#E7B800', '#FC4E07'),
                                         repel = TRUE)     # Avoid text overlapping)

    
    Cumm_pca<-summary(LUT.to.xtrain)$importance[3,]
    max.pca<-which(Cumm_pca  > 0.99)
    
    LUT.to.xtrain <-LUT.to.xtrain$x
    LUT.to.xtrain <- LUT.to.xtrain[,1:max.pca[1]]
    LUT.to.xVal <- LUT.to.xVal$x[,1:max.pca[1]]
    
    data.Xtrain <- LUT.to.xtrain
    data.Xval <- LUT.to.xVal
    
    split.database<- list('Xtrain' = data.Xtrain,'Ytrain' =data.Ytrain,
                          'Xval' =data.Xval,'Yval' = data.Yval,
                          'Scalar.train' = NA,
                          'Scalar.Ytrain' = NA,
                          'plot.train' = plot.pca.2)
    
    if (depVar.trans == TRUE){
      
      data.Ytrain.to<-as.matrix(data.Ytrain)
      colnames(data.Ytrain.to)<-depVar
      preprocess.scalarY <- caret::preProcess(data.Ytrain.to,method=c('range'))
      
      data.Ytrain <- predict(preprocess.scalarY, data.Ytrain.to)
      LUT.to.Train <-cbind(data.Ytrain.to, data.Ytrain,LUT.to.xtrain)
      colnames(LUT.to.Train)<-c(depVar,paste(depVar,'_Tra',sep=''),colnames(LUT.to.xtrain))
      LUT.to.Train<-as.matrix(LUT.to.Train)
      
      data.Ytrain <- round(LUT.to.Train[,2],4)
      
      data.Yval.to<-as.matrix(data.Yval)
      colnames(data.Yval.to)<-depVar
      
      data.Yval <- round(predict(preprocess.scalarY, data.Yval.to),4)
      LUT.to.Val <-cbind(data.Yval.to, data.Yval,LUT.to.xVal)
      colnames(LUT.to.Val)<-c(depVar,paste(depVar,'_Tra',sep=''),colnames(LUT.to.xVal))
      
      data.Yval <- round(LUT.to.Val[,2],4)
      
      split.database<- list('Xtrain' = data.Xtrain,'Ytrain' =data.Ytrain,
                            'Xval' =data.Xval,'Yval' = data.Yval,
                            'Scalar.train' = NA,
                            'Scalar.Ytrain' = preprocess.scalarY,
                            'plot.train' = plot.pca.2)
      
    } 
    
    
  }
  
return(split.database)

}
  
  

  
 
  