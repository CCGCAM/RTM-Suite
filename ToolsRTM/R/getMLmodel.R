
#' getMLmodel with a prefixed configuration
#'
#' @param dataset a dataframe
#' @param depVar name of the variable to predict
#' @param model a ML model. options are: 'CNN','Hidden-layers',
#' @param optimizer the optimizer for the model. options are: 'adam','adadelta','adagrad', 'adamax', 'nadam', 'msprop', 'sgd'
#' @param batch.size batch size used for each epoch. By default is 125
#' @param n.epochs  number of epoch. By default is 100
#' @param save.model a boolean variable for saving ML model, options are: TRUE or FALSE. if TRUE, please use path.model to give a folder for the model
#' @param path.model a path for saving the models. By default path.model ='Models'
#' @param prop.split a vector with proportion for spliting the dataset. prop.split =c(0.8,02) will be used as default.
#' @param data.trans a data.transformation method, options are: 'PCA','preProcess',
#' @param method.preProcess the data.transformation method for preProcess, data.transformation are: 'Normalize', 'YeoJohnson','BoxCox', Standarize', 'Center','Scale', and 'PCA'
#' @param depVar.trans a boolean variable for applying data transformation in Y variable, options are: TRUE or FALSE.
#'
#' @return a list with models and plots
#' @export
#'
#'


getMLmodel<-function(dataset=NULL, depVar='Cab',model='CNN',optimizer='adam',
                     batch.size=125,n.epochs=10, save.model=T, path.model=NULL,
                     prop.split=c(0.8,0.2),
                     data.trans='preProcess',method.preProcess='Normalize',
                     depVar.trans=FALSE) {

  if (!requireNamespace("keras", quietly = TRUE)) {
    stop("getMLmodel() requires the 'keras' package. Install it with install.packages(\"keras\") (and run keras::install_keras() once) and try again.", call. = FALSE)
  }

  stopifnot(class(dataset) == 'data.frame')

  stopifnot(model != 'CNN' | model != 'Hidden-layers')

  if (is.null(depVar)){
    stop('DepVar to predict should be indicated.')
  }


  if (!(depVar %in% colnames(dataset))){
    stop('inputs to predict should be in dataset, please check it ...')
  }
  if (is.null(depVar.trans)){
    depVar.trans=FALSE
  }


  inputs_<-colnames(dataset)
  if(any(inputs_ %in% depVar)){
    #print('DepVar is present')
    ##remove from dat
    inputs_<-inputs_[! inputs_ %in% depVar]
    ## reorder and put DepVar first in the vector.
    inputs_<-c(depVar,inputs_)
  } else{
    inputs_<-c(depVar,inputs_)
  }
  #print(inputs_)
  if (save.model == TRUE){
    if (is.null(path.model)){
      # output folder
      path.model='Models/'
      ifelse(!dir.exists(path.model), dir.create(path.model), FALSE)
      message(paste('model will save in ',path.model,' ',sep=''))
    } else{
      # output folder
      path.model=path.model
      ifelse(!dir.exists(path.model), dir.create(path.model), FALSE)
      message(paste('model will save in ',path.model,' ',sep=''))
    }


  }

  split.data<-ToolsRTM::getSplitData(data=dataset[,inputs_], depVar=depVar,inputs=inputs_[-1],
                                     data.trans=data.trans,prop.split=prop.split,method.preProcess=method.preProcess,
                                     depVar.trans=depVar.trans)
  ##########################################################################################
  ##### Parameters for the models
  ##########################################################################################

  callbacks = list(keras::callback_early_stopping(monitor = "val_loss",mode='min',patience = 5, restore_best_weights = TRUE))

  if (is.null(n.epochs)){
    n.epochs = 100
  } else {
    n.epochs = n.epochs
  }

  if (is.null(batch.size)){
    batch.size = 125
  } else {
    batch.size = batch.size
  }

  if (is.null(optimizer)) {
    optimizer = 'adam'
  }

  if (optimizer == 'adam') {
    opt = keras::optimizer_adam(learning_rate = 0.0001,beta_1 = 0.9, beta_2 = 0.999)

  } else if (optimizer == 'adadelta') {
    opt = keras::optimizer_adadelta(learning_rate = 1, rho = 0.95, epsilon = NULL, decay = 0)
    
  } else if (optimizer =='adagrad') {
    opt = keras::optimizer_adagrad(learning_rate = 0.01, epsilon = NULL, decay=0)
    
  } else if (optimizer =='adamax'){
    opt = keras::optimizer_adamax( learning_rate = 0.002, beta_1 = 0.9,beta_2 = 0.999)
     
  } else if (optimizer =='nadam'){
    opt = keras::optimizer_nadam(learning_rate = 0.002, beta_1 = 0.9,  beta_2 = 0.999, epsilon = NULL, schedule_decay = 0.004)
    
  } else if (optimizer =='rmsprop'){
    opt = keras::optimizer_rmsprop(learning_rate = 0.001, rho = 0.9, epsilon = NULL, decay = 0)
    
  } else if (optimizer =='sgd'){
    opt = keras::optimizer_sgd(learning_rate = 0.01, momentum = 0, decay = 0,  nesterov = FALSE)
  } else {
    stop("Unknown optimizer: '", optimizer, "'. Supported: 'adam', 'adadelta', 'adagrad', 'adamax', 'nadam', 'rmsprop', 'sgd'.")
  }
  stats<-list()
  ###########################################################################
  ############ step to improve the model
  ###########################################################################
  ## 1.Reduce your learning rate to a very small number like 0.001 or even 0.0001.
  ## 2.Provide more data.
  ## 3.Set Dropout rates to a number like 0.2. Keep them uniform across the network.
  ## 4.Try decreasing the batch size.
  ###########################################################################


  ##########################################################################################
  ##### Run the models
  ##########################################################################################

  if (model == 'Hidden-layers') {

    ##############################################################################################################################
    # sequential ML model ---
    ##############################################################################################################################


    data.Xtrain.reshape <- keras::array_reshape(split.data[['Xtrain']], c(nrow(split.data[['Xtrain']]), ncol(split.data[['Xtrain']])))
    dim(data.Xtrain.reshape)
    # Create the configuration for the model
    model.dML <- keras::keras_model_sequential() |>
      keras::layer_dense(units=64,activation="relu",
                  input_shape=c(dim(data.Xtrain.reshape)[2]))  |>
      #keras::layer_dropout(rate=0.1) |>
      keras::layer_dense(units = 32, activation = 'relu') |>
      keras::layer_dropout(rate=0.1) |>
      keras::layer_dense(units = 16, activation = 'relu') |>
      #keras::layer_dropout(rate=0.1) |>
      keras::layer_dense(units = 1, activation = 'relu')

    # Compile the configuration for the model
    model.dML |> keras::compile(loss = "mse",
                               optimizer = opt,#'adam',#,get_optimizer(),#"adam", #'sgd' can also be used
                               metrics = list("mean_absolute_error"))
    model.dML |> summary()


    # fit the configuration for the model
    history.model.dML<- model.dML |> keras::fit(data.Xtrain.reshape, split.data[['Ytrain']],
                                               epochs = n.epochs, batch_size = batch.size, verbose=1,shuffle=F,callbacks =callbacks,
                                               validation_split = 0.2)

    # evaluate the configuration for the model
    model.dML |> keras::evaluate(split.data[['Xval']], split.data[['Yval']])

    stats[['Hidden-layers']] <- model.dML |> keras::evaluate(split.data[['Xval']], split.data[['Yval']])

    # Save the model
    if (save.model == TRUE){
    model.dML |> keras::save_model_hdf5(paste(path.model,'/model_3hlayers_for_',depVar,'.hdf5',sep=''))
    model.dML |> keras::save_model_hdf5(paste(path.model,'Model-3hlayers-for-',depVar,'-',method.preProcess,'.hdf5',sep=''))
    saveRDS(split.data[['Scalar.train']], file = paste(path.model,'1-ScalerX-Model-3hlayers-for-',depVar,'-',method.preProcess,'.rds',sep=''))
      
    }

  } else if (model == 'CNN'){

    ##############################################################################################################################
    # CNN ML model ---
    ##############################################################################################################################

    #Reshaping the data for CNN
    ### These dimensions don't look correct; switch ncol() with nrow()

    data.Xtrain.CNN <- keras::array_reshape(split.data[['Xtrain']], c(nrow(split.data[['Xtrain']]), ncol(split.data[['Xtrain']]), 1))
    data.Xval.CNN <- keras::array_reshape(split.data[['Xval']], c(nrow(split.data[['Xval']]), ncol(split.data[['Xval']]), 1))

    # Create the configuration for the model
    dataset.dim=data.Xtrain.CNN
    model.dML <- keras::keras_model_sequential() |>
      keras::layer_conv_1d(filters=64, kernel_size=4, activation="relu",
                    ### The input shape doesn't look correct; instead of
                    ### `c(ncol(dataTrain_x), nrow(dataTrain_x))` (354, 13)
                    ### I believe you want `dim(dataTest_x)` (13, 1)
                    input_shape=c(ncol(dataset.dim), 1)) |>
      keras::layer_max_pooling_1d(pool_size=2) |>
      keras::layer_conv_1d(filters=32, kernel_size=2, activation="relu") |>
      keras::layer_max_pooling_1d(pool_size=2) |>
      #keras::layer_dropout(rate=0.4) |>
      keras::layer_flatten() |>
      keras::layer_dense(units=16, activation="relu") |>
      keras::layer_dropout(rate=0.1) |>
      keras::layer_dense(units=1, activation="relu")

    # Compile the configuration for the model

    model.dML |> keras::compile(loss = "mse",
                          optimizer =  opt,#"adam", #'sgd' can also be used
                          metrics = list("mean_absolute_error")
    )

    model.dML |> summary()
    # fit the configuration for the model
    history.model.dML <- model.dML |> keras::fit(data.Xtrain.CNN, split.data[['Ytrain']],
                                     epochs = n.epochs,batch_size = batch.size, verbose=1,shuffle=F,
                                     callbacks = callbacks,
                                     validation_split = 0.2)
    # evaluate the configuration for the model
    stats[['CNN-model']] <-model.dML |> keras::evaluate(data.Xval.CNN, split.data[['Yval']])
    # Save the model
    if (save.model == TRUE){
      model.dML |> keras::save_model_hdf5(paste(path.model,'Model-CNN-for-',depVar,'-',method.preProcess,'.hdf5',sep=''))
      saveRDS(split.data[['Scalar.train']], file = paste(path.model,'1-ScalerX-Model-CNN-for-',depVar,'-',method.preProcess,'.rds',sep=''))
    }
  }

  if (exists('model.dML')){

    if (depVar.trans == FALSE) {
      df.plot<-ToolsRTM::getPredicts(model=model.dML, type.model=model,
                                     data=split.data[['Xval']], data.trans=data.trans,
                                     data.Y=split.data[['Yval']],
                                     depVar=depVar)
      df.plot<- df.plot[, colSums(is.na(df.plot)) != nrow(df.plot)]
      colnames(df.plot) <- c(depVar,paste(depVar,'.predicted',sep=''))
    } else {
      scaler.depVar = split.data[['Scalar.Ytrain']]
      df.plot<-ToolsRTM::getPredicts(model=model.dML, type.model=model,
                                   data=split.data[['Xval']], data.trans=data.trans,
                                   data.Y=split.data[['Yval']],
                                   depVar=depVar, scaler.depVar= scaler.depVar)
      df.plot <- df.plot[, colSums(is.na(df.plot)) != nrow(df.plot)]
      colnames(df.plot) <- c(depVar,paste(depVar,'.predicted',sep=''))
      df.plot[,1] <-ToolsRTM::getReverse.trans(preProc=scaler.depVar,data=as.matrix(df.plot[,depVar]))

    }

  }
    #############################################################################################################################
    #	Plot in Validation dataset  -----
    ##############################################################################################################################

    axis_x<-expr(paste('measured ', !!depVar,sep=''))# axis x
    axis_y<-expr(paste('predicted ', !!depVar,sep=''))# axis y
    max_ <- max(df.plot[,depVar])
    depVar.predicted =paste(depVar,'.predicted',sep='')

    #### scatterplot for ANN model

    r2.stats<-round(cor(df.plot[,depVar],df.plot[,depVar.predicted],use='pairwise.complete.obs')^2,2)
    rmse.stats<-round(ToolsRTM::RMSE(df.plot[,depVar],df.plot[,depVar.predicted]),2)
    mylabel.r = bquote(bold(z: r)^2 == .(format(r2.stats, digits = 3)))

    statsLabel = paste0("Deep ML model; r2 = ", round(r2.stats,2), ", RMSE = ",round(rmse.stats,2))

    scatter.model<-  ggplot(df.plot, aes_string(x=depVar, y=depVar.predicted)) +
      geom_point(alpha=0.4,shape = 16,aes(), size=1.5) +  #geom_smooth(method=lm, formula = 'y ~ x', se=F,lty=2, lwd=1) +
      theme_bw() + xlim(0,max_) + ylim(0,max_)  + ggtitle('ANN model')  +

      scale_color_gradient(low = "#0091ff", high = "#f0650e") +
      theme(legend.position="bottom",
            plot.title = element_text(hjust = 0.5, size=10,face="bold"),
            axis.title = element_text(face="bold", size=10),
            axis.title.x = element_text(face="bold", size=12),
            axis.title.y = element_text(face="bold", size=12),
            axis.text.y=element_text(hjust = 0.5, size=10,face="bold"),
            axis.text.x=element_text(hjust = 0.5, size=10,face="bold"),
            legend.title=element_blank()) +
      labs(title = statsLabel, x=axis_x,y=axis_y, size=10,face="bold") +
      stat_smooth(method = "lm",formula = y ~ x,geom = "smooth",col='red',lty=2,se=T)
    
    # Save the scatterplot
    
    if (save.model == TRUE){
      ggsave(paste(path.model,'1-ScatterPLot-',model,'-for_',depVar,'-',method.preProcess,'.png',sep=''),
             width = 10, height = 10,  dpi = 300,units = "cm")
    }
    


  if (depVar.trans == FALSE) {
    model.outputs<- list('model' = model.dML,'history' =history.model.dML,
                        'stats' = stats, 'plot.val' =scatter.model,
                       'Scalar.train' = split.data[['Scalar.train']],
                       'Scalar.Ytrain' = NA,
                       'plot.train' = split.data[['plot.train']])
  } else {
    model.outputs<- list('model' = model.dML,'history' =history.model.dML,
                         'stats' = stats, 'plot.val' =scatter.model,
                         'Scalar.train' = split.data[['Scalar.train']],
                         'Scalar.Ytrain' = split.data[['Scalar.Ytrain']],
                         'plot.train' = split.data[['plot.train']])
  }

  return(model.outputs)
}


