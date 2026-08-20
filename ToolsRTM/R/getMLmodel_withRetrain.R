
#' getMLmodel is a function for retrain a deep model with a prefixed configuration
#'
#' @param dataset a dataframe
#' @param depVar name of the variable to predict
#' @param model a ML model. options are: 'CNN','Hidden-layers',
#' @param optimizer the optimizer for the model. options are: 'adam','adadelta','adagrad', 'adamax', 'nadam', 'msprop', 'sgd'
#' @param n.times  number of times to repaet the model. By default is 1
#' @param n.neurons  number of neurons to use by the the model. By default is 128 (only for hidden layer model is implemented)
#' @param n.layers  number of hidden-layer to use by the the model. By default is 4 (only for hidden layer model is implemented))

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

getMLmodel.withRetrain<-function(dataset=NULL, depVar='Cab',model='CNN',optimizer='adam',
                     n.times=NULL,
                     n.neurons=128,n.layers=4,
                     batch.size=125,n.epochs=100,
                     save.model=T, path.model=NULL,
                     prop.split=c(0.8,0.2),
                     data.trans='preProcess',method.preProcess='Normalize',
                     depVar.trans=FALSE) {



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

  if (is.null(n.times)) {
    N.times = 1
  } else {
    N.times = n.times
  }

  if (is.null(n.neurons)) {
    neurons = 128
  } else {
    neurons = n.neurons
  }

  if (is.null(n.layers)) {
    n.layers = 4
  } else {
    n.layers = n.layers
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


  ##########################################################################################
  ##### Parameters for the models
  ##########################################################################################

  callbacks_ = callback_early_stopping(monitor = 'val_loss', mode='min',patience = 5,restore_best_weights = TRUE)

  if (is.null(n.epochs)){
    n.epochs = 100
  } else {
    n.epochs = n.epochs
  }

  if (is.null(batch.size)){
    batch.size = 32
  } else {
    batch.size = batch.size
  }

  if (is.null(optimizer)) {
    optimizer = 'adam'
  }

  # BUG (fixed): opt/opt.retrain used to be built ONCE here, before the
  # `for (i.times in 1:N.times)` retrain loop below, and then reused
  # (via `optimizer = opt` in the compile() calls) to compile a fresh
  # keras model on every loop iteration. TF/Keras optimizers track
  # variables by reference from the model they were first compiled
  # against, so reusing the same optimizer object across DIFFERENT model
  # instances crashes on the 2nd+ iteration with:
  #   KeyError: 'The optimizer cannot recognize variable conv1d_3/kernel:0 ...'
  # (reproduced directly via Scripts/3-Inversion.R's CNN retrain loop,
  # n.times = 3 -- iteration 1 succeeds, iteration 2 crashes). Fixed by
  # building a brand-new optimizer instance inside the loop, once per
  # retrain iteration, instead of once for the whole function call.
  make_optimizer <- function(optimizer, retrain = FALSE) {
    if (optimizer == 'adam') {
      if (retrain) optimizer_adam(learning_rate = 0.00001, beta_1 = 0.9, beta_2 = 0.999)
      else optimizer_adam(learning_rate = 0.0001, beta_1 = 0.9, beta_2 = 0.999)

    } else if (optimizer == 'adadelta') {
      if (retrain) optimizer_adadelta(learning_rate = 0.1, rho = 0.95, epsilon = NULL, decay = 0)
      else optimizer_adadelta(learning_rate = 1, rho = 0.95, epsilon = NULL, decay = 0)

    } else if (optimizer =='adagrad') {
      if (retrain) optimizer_adagrad(learning_rate = 0.001, epsilon = NULL, decay=0)
      else optimizer_adagrad(learning_rate = 0.01, epsilon = NULL, decay=0)

    } else if (optimizer =='adamax'){
      if (retrain) optimizer_adamax( learning_rate = 0.0002, beta_1 = 0.9,beta_2 = 0.999)
      else optimizer_adamax( learning_rate = 0.002, beta_1 = 0.9,beta_2 = 0.999)

    } else if (optimizer =='nadam'){
      if (retrain) optimizer_nadam(learning_rate = 0.0002, beta_1 = 0.9,  beta_2 = 0.999, epsilon = NULL, schedule_decay = 0.004)
      else optimizer_nadam(learning_rate = 0.002, beta_1 = 0.9,  beta_2 = 0.999, epsilon = NULL, schedule_decay = 0.004)

    } else if (optimizer =='rmsprop'){
      if (retrain) optimizer_rmsprop(learning_rate = 0.0001, rho = 0.9, epsilon = NULL, decay = 0)
      else optimizer_rmsprop(learning_rate = 0.001, rho = 0.9, epsilon = NULL, decay = 0)

    } else if (optimizer =='sgd'){
      if (retrain) optimizer_sgd(learning_rate = 0.001, momentum = 0, decay = 0,  nesterov = FALSE)
      else optimizer_sgd(learning_rate = 0.01, momentum = 0, decay = 0,  nesterov = FALSE)
    } else {
      stop("Unknown optimizer: '", optimizer, "'. Supported: 'adam', 'adadelta', 'adagrad', 'adamax', 'nadam', 'rmsprop', 'sgd'.")
    }
  }
  # validate the optimizer name once up front (was previously validated as a
  # side effect of the if/else chain above) without instantiating an unused
  # optimizer object
  if (!(optimizer %in% c('adam', 'adadelta', 'adagrad', 'adamax', 'nadam', 'rmsprop', 'sgd'))) {
    stop("Unknown optimizer: '", optimizer, "'. Supported: 'adam', 'adadelta', 'adagrad', 'adamax', 'nadam', 'rmsprop', 'sgd'.")
  }

  ###########################################################################
  ############ step to improve the model
  ###########################################################################
  ## 1.Reduce your learning rate to a very small number like 0.001 or even 0.0001.
  ## 2.Provide more data.
  ## 3.Set Dropout rates to a number like 0.2. Keep them uniform across the network.
  ## 4.Try decreasing the batch size.
  ###########################################################################

  stats<-list()
  preds.model <- list()
  preds.model.retrain <- list()
  scatters <- list()
  models.keras <-list()
  history.keras <-list()
  Scalar.to  <-list()
  plot.cor.to <-list()

  #barProgress <- txtProgressBar(min = 1, max = N.times, style = 3)

  for (i.times in c(1:N.times)){

    message(paste(model,' with N Time :',i.times,sep=''))
    #print(i.times)
    #setTxtProgressBar(barProgress, i.times)
    # fresh optimizer instance for THIS iteration's model -- see note above
    # make_optimizer() definition for why these must not be reused across
    # iterations.
    opt <- make_optimizer(optimizer, retrain = FALSE)
    opt.retrain <- make_optimizer(optimizer, retrain = TRUE)
    ##########################################################################################
    ##### Split the Dataset
    ##########################################################################################
    set.seed(1234+i.times)
    split.data<-ToolsRTM::getSplitData_noMessages(data=dataset[,inputs_], depVar=depVar,inputs=inputs_[-1],
                                       data.trans=data.trans,prop.split=prop.split,method.preProcess=method.preProcess,
                                       depVar.trans=depVar.trans)
    scaler.train <-  split.data[['Scalar.train']]
    ##########################################################################################
    ##### Run the models
    ##########################################################################################

    if (model == 'Hidden-layers') {

      ##############################################################################################################################
      # sequential ML model ---
      ##############################################################################################################################


      data.Xtrain.reshape <- array_reshape(split.data[['Xtrain']], c(nrow(split.data[['Xtrain']]), ncol(split.data[['Xtrain']])))
      dim(data.Xtrain.reshape)
      # Create the configuration for the model
      n.units= neurons
      n.layers = n.layers

      create_model <- function(neurons, n.layers) {
        model <- keras_model_sequential()

        # Input layer
        model <- model |>
          layer_dense(units = neurons, activation = "relu", input_shape = c(dim(data.Xtrain.reshape)[2]))

        # Hidden layers
        for (i in 1:n.layers) {
          units <- max(1, round(neurons / (2^i), 0)) # Calculate units dynamically
          model <- model |>
            #layer_dropout(rate=0.1) |>
            layer_dense(units = units, activation = 'relu')
        }

        # Output layer
        model <- model |>
          layer_dense(units = 1, activation = 'relu')

        return(model)
      }
      model.dML <- create_model(neurons = n.units, n.layers = n.layers)

      # Compile the configuration for the model
      model.dML |> compile(loss = "mse",
                            optimizer = opt,#'adam',#,get_optimizer(),#"adam", #'sgd' can also be used
                            metrics = list("mean_absolute_error"))
      model.dML |> summary()

      # fit the configuration for the model
      history.model.dML<- model.dML |> fit(data.Xtrain.reshape, split.data[['Ytrain']],
                                            epochs = n.epochs, batch_size = batch.size, verbose=1,shuffle=F,callbacks =callbacks_,
                                            validation_split = 0.2)

      # evaluate the configuration for the model
      model.dML |> evaluate(split.data[['Xval']], split.data[['Yval']])

      ### skill scores
      #stats[['Hidden-layers']] <- model.dML |> evaluate(split.data[['Xval']], split.data[['Yval']])

      if (depVar.trans == FALSE) {

        df.val<-ToolsRTM::getPredicts(model=model.dML, type.model='Hidden-layers',
                                       data=split.data[['Xval']], data.trans=data.trans,
                                       data.Y=split.data[['Yval']],
                                       depVar=depVar)
        df.val<- df.val[, colSums(is.na(df.val)) != nrow(df.val)]
        colnames(df.val) <- c(depVar,paste(depVar,'.predicted',sep=''))
      } else {
        scaler.depVar = split.data[['Scalar.Ytrain']]
        df.val<-ToolsRTM::getPredicts(model=model.dML, type.model='Hidden-layers',
                                       data=split.data[['Xval']], data.trans=data.trans,
                                       data.Y=split.data[['Yval']],
                                       depVar=depVar, scaler.depVar= scaler.depVar)
        df.val <- df.val[, colSums(is.na(df.val)) != nrow(df.val)]
        colnames(df.val) <- c(depVar,paste(depVar,'.predicted',sep=''))
        df.val[,1] <-ToolsRTM::getReverse.trans(preProc=scaler.depVar,data=as.matrix(df.val[,depVar]))

      }

      ### skill scores
      table.stats <- list()
      table.stats['VarDep'] <- depVar
      table.stats['Model'] <- 'Hidden-layers'
      table.stats['Trans'] <- method.preProcess
      table.stats['db'] <-'Testing'
      table.stats['i.Time'] <- i.times
      table.stats['MAE'] <-  round(MLmetrics::MAE(df.val[,depVar], df.val[,paste(depVar,'.predicted',sep='')]),3)
      table.stats['RMSE'] <-  round(MLmetrics::RMSE(df.val[,depVar], df.val[,paste(depVar,'.predicted',sep='')]),3)
      table.stats['R2'] <-  round(MLmetrics::R2_Score(df.val[,depVar], df.val[,paste(depVar,'.predicted',sep='')]),3)
      table.stats<-data.frame(do.call(cbind,table.stats))
      #print(table.stats)

      #############################################################################################################################
      #	Update model wit new predictions  -----
      ##############################################################################################################################

      data.to.retrain <- dataset[,inputs_]
      #data.to.retrain <-na.omit(data.to.retrain)


      dim(data.to.retrain)
      split.retrain<-ToolsRTM::getSplitData_noMessages(data=data.to.retrain, depVar=depVar,inputs=inputs_[-1],
                                            data.trans=data.trans,prop.split=c(0.95,0.05),method.preProcess=method.preProcess,
                                            depVar.trans=depVar.trans)
      scaler.train.retrain <-  split.retrain[['Scalar.train']]

      data.Xtrain.reshape <- array_reshape(split.retrain[['Xtrain']], c(nrow(split.retrain[['Xtrain']]), ncol(split.retrain[['Xtrain']])))
      dim(data.Xtrain.reshape)

      # fit the configuration for the model
      history.model.dML<- model.dML |> fit(data.Xtrain.reshape, split.retrain[['Ytrain']],
                                            epochs = floor(n.epochs/2), batch_size = floor(batch.size/2), verbose=1,shuffle=F,callbacks =callbacks_,
                                            validation_split = 0.2)

      # evaluate the configuration for the model
      model.dML |> evaluate(split.retrain[['Xval']], split.retrain[['Yval']])

      ### skill scores
      #stats[['Hidden-layers']] <- model.dML |> evaluate(split.retrain[['Xval']], split.retrain[['Yval']])

      if (depVar.trans == FALSE) {

        df.val.retrain<-ToolsRTM::getPredicts(model=model.dML, type.model='Hidden-layers',
                                      data=split.data[['Xval']], data.trans=data.trans,
                                      data.Y=split.data[['Yval']],
                                      depVar=depVar)
        df.val.retrain<- df.val.retrain[, colSums(is.na(df.val.retrain)) != nrow(df.val.retrain)]
        colnames(df.val.retrain) <- c(depVar,paste(depVar,'.predicted',sep=''))
      } else {
        scaler.depVar = split.retrain[['Scalar.Ytrain']]
        df.val.retrain<-ToolsRTM::getPredicts(model=model.dML, type.model='Hidden-layers',
                                      data=split.data[['Xval']], data.trans=data.trans,
                                      data.Y=split.data[['Yval']],
                                      depVar=depVar, scaler.depVar= scaler.depVar)
        df.val.retrain <- df.val.retrain[, colSums(is.na(df.val.retrain)) != nrow(df.val.retrain)]
        colnames(df.val.retrain) <- c(depVar,paste(depVar,'.predicted',sep=''))
        df.val.retrain[,1] <-ToolsRTM::getReverse.trans(preProc=scaler.depVar,data=as.matrix(df.val.retrain[,depVar]))

      }


      table.stats.r<-list()
      table.stats.r['VarDep'] <- depVar
      table.stats.r['Model'] <- 'Hidden-layers'
      table.stats.r['Trans'] <- method.preProcess
      table.stats.r['db'] <-'Val.retrain'
      table.stats.r['i.Time'] <- i.times

      table.stats.r['MAE'] <-  round(MLmetrics::MAE(df.val.retrain[,depVar], df.val.retrain[,paste(depVar,'.predicted',sep='')]),3)
      table.stats.r['RMSE'] <-  round(MLmetrics::RMSE(df.val.retrain[,depVar], df.val.retrain[,paste(depVar,'.predicted',sep='')]),3)
      table.stats.r['R2'] <-  round(MLmetrics::R2_Score(df.val.retrain[,depVar], df.val.retrain[,paste(depVar,'.predicted',sep='')]),3)

      table.stats.r<-data.frame(do.call(cbind,table.stats.r))

      table.stats.to.export<-rbind(table.stats,table.stats.r)

      # Save the model
      if (save.model == TRUE){
        model.dML |> save_model_hdf5(paste(path.model,'Model-3hlayers-for-',depVar,'-',method.preProcess,'-',i.times,'.h5',sep=''))
        model.dML |> save_model_weights_hdf5(paste(path.model,'Model-3hlayers-for-',depVar,'-',method.preProcess,'-',i.times,'-weights.h5',sep=''))

        saveRDS(split.retrain[['Scalar.train']], file = paste(path.model,'1-ScalerX-Model-3hlayers-for-',depVar,'-',method.preProcess,'-',i.times,'.rds',sep=''))
        write.table(table.stats.to.export, file = paste(path.model,'1-Statistcal_scores_for_Model-3hlayers-for-',depVar,'-',method.preProcess,'-',i.times,'.csv',sep=''),sep=',',row.names = F)

      }

    } else if (model == 'CNN'){

      ##############################################################################################################################
      # CNN ML model ---
      ##############################################################################################################################

      #Reshaping the data for CNN
      ### These dimensions don't look correct; switch ncol() with nrow()

      data.Xtrain.CNN <- array_reshape(split.data[['Xtrain']], c(nrow(split.data[['Xtrain']]), ncol(split.data[['Xtrain']]), 1))
      data.Xval.CNN <- array_reshape(split.data[['Xval']], c(nrow(split.data[['Xval']]), ncol(split.data[['Xval']]), 1))

      # Create the configuration for the model
      dataset.dim=data.Xtrain.CNN
      n.units = neurons

      #:::::::::::::::::::::::::::::::::::::::::::::::::::::::::::::::::
      # CNN model by function
      #:::::::::::::::::::::::::::::::::::::::::::::::::::::::::::::::::

      create_cnn_model <- function(neurons, n.layers) {
        model <- keras_model_sequential()

        # Input layer
        model <- model |>
          layer_conv_1d(filters = neurons, kernel_size = 4, activation = "relu", input_shape = c(ncol(dataset.dim), 1)) |>
          layer_max_pooling_1d(pool_size = 2) |>
          layer_conv_1d(filters = n.units, kernel_size = 2, activation = "relu") |>
          layer_flatten()
        # Subsequent hidden-layers
        for (i in 1:(n.layers - 1)) {
          units <- max(1, round(neurons / (2^i), 0)) # Calculate filters dynamically
          model <- model |>
            #layer_dropout(rate=0.1) |>
            layer_dense(units = units, activation = 'relu')
        }

        # Dense layers
        model <- model |>
          layer_dense(units = 1, activation = "relu")

        return(model)
      }
      model.dML <- create_cnn_model(neurons = n.units, n.layers = n.layers)
      # Summary the CNN model
       summary(model.dML)
      # Compile the configuration for the model

      model.dML |> compile(loss = "mse",
                            optimizer =  opt,#"adam", #'sgd' can also be used
                            metrics = list("mean_absolute_error"))

      #model.dML |> summary()
      # fit the configuration for the model
      history.model.dML <- model.dML |> fit(data.Xtrain.CNN, split.data[['Ytrain']],
                                             epochs = n.epochs,batch_size = batch.size, verbose=1,shuffle=F,
                                             callbacks = callbacks_,
                                             validation_split = 0.2)
      # evaluate the configuration for the model
      #stats[['CNN-model']] <-model.dML |> evaluate(data.Xval.CNN, split.data[['Yval']])


      if (depVar.trans == FALSE) {

        df.val<-ToolsRTM::getPredicts(model=model.dML, type.model='CNN',
                                      data=split.data[['Xval']], data.trans=data.trans,
                                      data.Y=split.data[['Yval']],
                                      depVar=depVar)
        df.val<- df.val[, colSums(is.na(df.val)) != nrow(df.val)]
        colnames(df.val) <- c(depVar,paste(depVar,'.predicted',sep=''))
      } else {
        scaler.depVar = split.data[['Scalar.Ytrain']]
        df.val<-ToolsRTM::getPredicts(model=model.dML, type.model='CNN',
                                      data=split.data[['Xval']], data.trans=data.trans,
                                      data.Y=split.data[['Yval']],
                                      depVar=depVar, scaler.depVar= scaler.depVar)
        df.val <- df.val[, colSums(is.na(df.val)) != nrow(df.val)]
        colnames(df.val) <- c(depVar,paste(depVar,'.predicted',sep=''))
        df.val[,1] <-ToolsRTM::getReverse.trans(preProc=scaler.depVar,data=as.matrix(df.val[,depVar]))

      }

      ### skill scores
      table.stats <- list()
      table.stats['VarDep'] <- depVar
      table.stats['Model'] <- 'CNN'
      table.stats['Trans'] <- method.preProcess
      table.stats['db'] <-'Testing'
      table.stats['i.Time'] <- i.times
      table.stats['MAE'] <-  round(MLmetrics::MAE(df.val[,depVar], df.val[,paste(depVar,'.predicted',sep='')]),3)
      table.stats['RMSE'] <-  round(MLmetrics::RMSE(df.val[,depVar], df.val[,paste(depVar,'.predicted',sep='')]),3)
      table.stats['R2'] <-  round(MLmetrics::R2_Score(df.val[,depVar], df.val[,paste(depVar,'.predicted',sep='')]),3)
      table.stats<-data.frame(do.call(cbind,table.stats))
      #print(table.stats)
      #############################################################################################################################
      #	Update model wit new predictions  -----
      ##############################################################################################################################

      data.to.retrain <- dataset[,inputs_]
      #data.to.retrain <-na.omit(data.to.retrain)
      dim(data.to.retrain)
      split.retrain<-ToolsRTM::getSplitData_noMessages(data=data.to.retrain, depVar=depVar,inputs=inputs_[-1],
                                            data.trans=data.trans,prop.split=c(0.90,0.1),method.preProcess=method.preProcess,
                                            depVar.trans=depVar.trans)

      scaler.train.retrain <-  split.retrain[['Scalar.train']]

      data.Xtrain.reshape <- array_reshape(split.retrain[['Xtrain']], c(nrow(split.retrain[['Xtrain']]), ncol(split.retrain[['Xtrain']]), 1))
      data.Xval.reshape <- array_reshape(split.retrain[['Xval']], c(nrow(split.retrain[['Xval']]), ncol(split.retrain[['Xval']]), 1))
      dim(data.Xtrain.reshape)

      # fit the configuration for the model
      history.model.dML<- model.dML |> fit(data.Xtrain.reshape, split.retrain[['Ytrain']],
                                            epochs = floor(n.epochs/2), batch_size = floor(batch.size/2), verbose=1,shuffle=F,callbacks =callbacks_,
                                            validation_split = 0.1)

      # evaluate the configuration for the model
      model.dML |> evaluate(data.Xval.reshape, split.retrain[['Yval']])

      ### skill scores
      #stats[['Hidden-layers']] <- model.dML |> evaluate(split.retrain[['Xval']], split.retrain[['Yval']])

      if (depVar.trans == FALSE) {

        df.val.retrain<-ToolsRTM::getPredicts(model=model.dML, type.model='CNN',
                                              data=split.data[['Xval']], data.trans=data.trans,
                                              data.Y=split.data[['Yval']],
                                              depVar=depVar)
        df.val.retrain<- df.val.retrain[, colSums(is.na(df.val.retrain)) != nrow(df.val.retrain)]
        colnames(df.val.retrain) <- c(depVar,paste(depVar,'.predicted',sep=''))
      } else {
        scaler.depVar = split.retrain[['Scalar.Ytrain']]
        df.val.retrain<-ToolsRTM::getPredicts(model=model.dML, type.model='CNN',
                                              data=split.data[['Xval']], data.trans=data.trans,
                                              data.Y=split.data[['Yval']],
                                              depVar=depVar, scaler.depVar= scaler.depVar)
        df.val.retrain <- df.val.retrain[, colSums(is.na(df.val.retrain)) != nrow(df.val.retrain)]
        colnames(df.val.retrain) <- c(depVar,paste(depVar,'.predicted',sep=''))
        df.val.retrain[,1] <-ToolsRTM::getReverse.trans(preProc=scaler.depVar,data=as.matrix(df.val.retrain[,depVar]))

      }

      # skill scores
      table.stats.r<-list()
      table.stats.r['VarDep'] <- depVar
      table.stats.r['Model'] <- 'CNN'
      table.stats.r['Trans'] <- method.preProcess
      table.stats.r['db'] <-'Val.retrain'
      table.stats.r['i.Time'] <- i.times

      table.stats.r['MAE'] <-  round(MLmetrics::MAE(df.val.retrain[,depVar], df.val.retrain[,paste(depVar,'.predicted',sep='')]),3)
      table.stats.r['RMSE'] <-  round(MLmetrics::RMSE(df.val.retrain[,depVar], df.val.retrain[,paste(depVar,'.predicted',sep='')]),3)
      table.stats.r['R2'] <-  round(MLmetrics::R2_Score(df.val.retrain[,depVar], df.val.retrain[,paste(depVar,'.predicted',sep='')]),3)

      table.stats.r<-data.frame(do.call(cbind,table.stats.r))

      table.stats.to.export<-rbind(table.stats,table.stats.r)


      # Save the model
      if (save.model == TRUE){
        model.dML |> save_model_hdf5(paste(path.model,'Model-CNN-for-',depVar,'-',method.preProcess,'-',i.times,'.h5',sep=''))
        model.dML |> save_model_weights_hdf5(paste(path.model,'Model-CNN-for-',depVar,'-',method.preProcess,'-',i.times,'-weights.h5',sep=''))
        saveRDS(split.retrain[['Scalar.train']], file = paste(path.model,'1-ScalerX-Model-CNN-for-',depVar,'-',method.preProcess,'-',i.times,'.rds',sep=''))

        write.table(table.stats.to.export, file = paste(path.model,'1-Statistcal_scores_for_Model-CNN-for-',depVar,'-',method.preProcess,'-',i.times,'.csv',sep=''),sep=',',row.names = F)


      }
    }
    ##############################################################################################################################
    # Save table with skill scores for times  ---
    ##############################################################################################################################

    table.stats.to.save <- rbind(table.stats,table.stats.r)
    print(table.stats.to.save)
    stats[[i.times]] <-table.stats.to.save

    y.predicted = paste(depVar,'.predicted',sep='')
    preds.model[[i.times]] <-df.val.retrain[,c(depVar,y.predicted)]
    preds.model.retrain[[i.times]] <-df.val.retrain[,c(depVar,y.predicted)]



    ##############################################################################################################################
    # Plot predicted DepVar by model at Testing dataset  ---
    ##############################################################################################################################

    axis_x<-expr(paste('measured ', !!depVar,sep=''))# axis x
    axis_y<-expr(paste('predicted ', !!depVar,sep=''))# axis y
    max_ <- max(df.val.retrain[,depVar])

    #### scatterplot for ANN model

    r2.SE2<-round(cor(df.val.retrain[,depVar],df.val.retrain[,y.predicted],use='pairwise.complete.obs')^2,2)
    rmse.se<-round(ToolsRTM::RMSE(df.val.retrain[,depVar],df.val.retrain[,y.predicted]),2)
    mylabel.r.se2 = bquote(bold(SE2: r)^2 == .(format(r2.SE2, digits = 3)))

    statsLabel = paste0("r2 = ", round(r2.SE2,2), ", RMSE = ",round(rmse.se,2))

    scatter.model <-  ggplot(df.val.retrain, aes_string(x=depVar, y=y.predicted)) +
      geom_point(alpha=0.1,shape = 16,aes(), size=1.5) +  #geom_smooth(method=lm, formula = 'y ~ x', se=F,lty=2, lwd=1) +
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
      labs(title = statsLabel, x=axis_x,y=axis_y, size=8,face="bold") +
      stat_smooth(method = "lm",formula = y ~ x,geom = "smooth",col='#C1CDC1',lty=2,se=T)
    print(scatter.model)

    # Save the scatterplot
    if (save.model == TRUE){
      ggsave(paste(path.model,'1-ScatterPLot-',model,'-for_',depVar,'-',method.preProcess,'-',i.times,'.png',sep=''),
             width = 10, height = 10,  dpi = 300,units = "cm")
    }


    ####

    scatters[[i.times]] <- scatter.model
    models.keras[[i.times]]<- model.dML
    history.keras[[i.times]] <- history.model.dML

    Scalar.train<-split.retrain[['Scalar.train']]
    Scalar.to[[i.times]] <- Scalar.train

    plot.cor <-split.retrain[['plot.train']]

    plot.cor.to[[i.times]] <- plot.cor



  }
  stats.to.export<-data.frame(do.call(rbind, stats))


  #preds.model.to.export<-data.frame(do.call(cbind, preds.model))
  #colnames(preds.model.to.export) <-c(paste(depVar,'model',method.preProcess,c(1:N.times),sep='.'))

  #preds.model.retrain.to.export<-data.frame(do.call(cbind, preds.model.retrain))
  #colnames(preds.model.retrain.to.export) <-c(paste(depVar,'retrain',method.preProcess,c(1:N.times),sep='.'))

  if (depVar.trans == FALSE) {
    model.outputs<- list('model' = models.keras,'history' =history.keras,
                         'stats' = stats.to.export,
                         'Scalar.train' = Scalar.to,
                         'Scalar.Ytrain' = NA,
                         'plot.val' =scatters,
                         'plot.cor' = plot.cor.to)
  }
  else {
    model.outputs<- list('model' = models.keras, 'history' =history.keras,
                         'stats' = stats.to.export,
                         'Scalar.train' = Scalar.to,
                         'Scalar.Ytrain' = split.retrain[['Scalar.Ytrain']],
                         'plot.val' = scatters,
                         'plot.cor' = plot.cor.to)
  }

  return(model.outputs)
  #close(barProgress)
}



