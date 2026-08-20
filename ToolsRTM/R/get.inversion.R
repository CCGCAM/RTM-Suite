#' get.inversion of the plant traits using different machine learning models
#'
#' @param data The data frame containing the variables.
#' @param depVar The dependent variable.
#' @param inputs The independent variables.
#' @param algorithm The type of machine learning model. Options: "PLSR" (Partial Least Squares Regression), "SVM" (Support Vector Machine),
#' "RF" (Random Forest); "NN" (Neural Network); "GB" (Gradient boosting); "xGB" (eXtreme Gradient Boosting (XGBoost) with linear base learners);
#' 'Bayesian' ( Bayesian Additive Regression Trees); 'AdaBag' ( Bagged AdaBoost); "qLASSO" (Quantile Regression with LASSO penalty); "RVM" (Relevance Vector Machines (RVM) with linear kernel);
#' 'BRNN' (Bayesian Regularized Neural Networks); "Ensemble" (Stacking Ensemble models)
#' Default is "PLSR".
#' @param method.resampling The resampling method for controlling tht ML: Options are: "boot" (Bootstrapping); "boot632" (Bootstrapping-632);
#' "optimism_boot"; "boot_all"; "cv" (cross-Validation); "repeatedcv" (repeats k-fold cross-validation with 3 times);
#' "LOOCV" (Leave-One-Out Cross-Validation with 3 times); "LGOCV"
#' @param n.cores The number of cores
#' @param ... Additional arguments (currently unused, reserved for future extensions).

#' @param seed The seed for reproducibility. Default is 123.
#' @param n.samples A integer with the number of samples used for tunning search (nsample/2) and create the ML model (n.sample)
#' @param save.model Logical indicating whether to save the trained models. Default is FALSE.
#' @param save.path Path to save the trained models. Required if save_models is TRUE.

#' @return A list containing predictions, statistics, and plots for the specified machine learning model.
#' @export
#'
#' @examples
#' \dontrun{
#' get.inversion(data = my_data, depVar = "Cab", inputs = c("NDVI", "TCARI"), ML = "SVM", seed = 123)
#' }

get.inversion <- function(data, depVar, inputs, algorithm='PLSR',method.resampling=NULL,n.cores=NULL,
                          seed=123, n.samples=500, save.model = FALSE, save.path = NULL, ...) {

  #https://topepo.github.io/caret/available-models.html
  #https://topepo.github.io/caret/model-training-and-tuning.html#model-training-and-parameter-tuning
  #https://topepo.github.io/caret/train-models-by-tag.html#gaussian-process
  # Set the seed for reproducibility
  if(is.null(seed)) {
    set.seed(as.numeric(Sys.time()) %% 10000)  # Default seed value
  } else{
    set.seed(seed)
  }


  # Set the n.cores
  # Was: both branches called parallel::makeCluster(...) here regardless of
  # the n.cores argument -- assigning a *cluster object* to n.cores (not a
  # number), which then got leaked (never stopped, since every algorithm
  # branch below immediately overwrote its own `clusters <- makeCluster(n.cores)`
  # using that leaked cluster object as the `spec` argument). n.cores is now
  # kept as the actual number of workers requested, honored below: passing
  # n.cores = 1 skips cluster creation entirely and runs sequentially.
  if(is.null(n.cores)) {
    n.cores <- max(1, parallel::detectCores() - 2)
  }

  # Set the seed the resampling method in the Model
  if(is.null(method.resampling)) {
    method.resampling <- 'boot'  # Default method
  } else{
    method.resampling <- method.resampling
  }


  # Reduce sample size for tuning if n.samples is provided
  if (!is.null(n.samples) && n.samples <= nrow(data)) {
    rows.r <- sample(nrow(data), n.samples)
  } else {
    stop("n.samples must be less than or equal to the number of rows in data.")
  }

  #Subset data for modeling
  data <- data[rows.r, ]

  # Create train-test split indices based on the 'depVar' column
  indices <- caret::createDataPartition(data[,depVar], p = 0.7, list=F)

  if (length(indices) == nrow(data)){
    indices <- sample(1:nrow(data), size = round(0.7 * nrow(data)))
  }

  df.train <- data[indices, ]  # Training set
  df.test <- data[-indices, ]   # Testing set



  # Force to run PLSR if ML is null or empty
  if(is.null(algorithm) || algorithm == "") {
    algorithm <- "PLSR"
  }

  if(algorithm == "PLSR") {

    print('processing hybrid approach using PLSR ...')

    fmla.n <- as.formula(paste(depVar," ~ ", paste(inputs, collapse= "+")))

    # Parallelize tuning process
    clusters <- if (n.cores > 1) parallel::makeCluster(n.cores) else NULL
    if (!is.null(clusters)) doParallel::registerDoParallel(clusters) else foreach::registerDoSEQ()

    myfolds <- caret::createMultiFolds(df.train[,depVar], k = 10, times = 10)
    control <- caret::trainControl(method= method.resampling, index = myfolds, selectionFunction = "oneSE")

    # Train PLS model
    model_ <- caret::train(fmla.n, data = df.train[,],
                   method = "pls",
                   metric = "RMSE",
                   tuneLength = 20,
                   trControl = control,
                   preProc = c("zv","center","scale"))

    k=model_$bestTune

    if (!is.null(clusters)) parallel::stopCluster(clusters)
    # estimate variable importance
    importance <- caret::varImp(model_, scale=FALSE)

  } else if(algorithm == "SVM") {

    print('processing hybrid approach using Support Vector Machine ...')

    fmla.n <- as.formula(paste(depVar," ~ ", paste(inputs, collapse= "+")))

    # Parallelize tuning process
    clusters <- if (n.cores > 1) parallel::makeCluster(n.cores) else NULL
    if (!is.null(clusters)) doParallel::registerDoParallel(clusters) else foreach::registerDoSEQ()

    tobj2 <- e1071::tune.svm(fmla.n, data = df.train[,], sampling = "fix",
                             gamma = 2^c(-10, -8, -6, -4),  # search space for gamma
                             cost = 2^c(-5, -3, -1, 1),     # search space for cost
                             tunecontrol =  e1071::tune.control(cross = 5))  # number of cross-validation folds
    cc <- as.numeric(tobj2$best.parameters[2])
    gg <- as.numeric(tobj2$best.parameters[1])

    # SVM
    model_ <- e1071::svm(fmla.n, kernel = "radial", data = df.train[, ], gamma = gg, cost = cc,
                     type = "eps-regression", probability = FALSE)

    if (!is.null(clusters)) parallel::stopCluster(clusters)

    # estimate variable importance
    # Get the coefficients
    coefficients <- model_$coefs
    # Calculate the magnitude of coefficients
    importance <- apply(coefficients, 1, function(x) sqrt(sum(x^2)))
    importance <- NA

  } else if(algorithm == "RF") {

    print('processing hybrid approach using Random Forest ...')

    fmla.n <- as.formula(paste(depVar," ~ ", paste(inputs, collapse= "+")))

    # Parallelize tuning process
    clusters <- if (n.cores > 1) parallel::makeCluster(n.cores) else NULL
    if (!is.null(clusters)) doParallel::registerDoParallel(clusters) else foreach::registerDoSEQ()

    # Define tuning grid with reduced search space
    mtry2 <- randomForest::tuneRF(df.train[, inputs], y = df.train[, depVar],
                    ntreeTry = ncol(df.train[, inputs])/3, stepFactor = 1.5, improve = 0.01,
                    trace = F, plot = F)
    best.m <- mtry2[mtry2[, 2] == min(mtry2[, 2]), 1]
    metric <- "RMSE"
    tunegrid <- expand.grid(.mtry = best.m)
    n.trees <- round(ncol(df.train[, inputs])/3,0)

    # Define training control
    fit.control <- caret::trainControl(method = method.resampling, number = 3,
                 search = "grid", repeats = 3, allowParallel = TRUE)
    # Random Forest
    model_ <- caret::train(fmla.n, data = df.train[, ], method = "rf", metric = metric,
                           trControl = fit.control, verbose=F, tuneGrid = tunegrid)

    if (!is.null(clusters)) parallel::stopCluster(clusters)

    # estimate variable importance
    importance <- caret::varImp(model_, scale=FALSE)

  } else if (algorithm == 'GB'){
    print('processing hybrid approach using Gradient Boosting ...')

    fmla.n <- as.formula(paste(depVar," ~ ", paste(inputs, collapse= "+")))

    # Parallelize tuning process
    clusters <- if (n.cores > 1) parallel::makeCluster(n.cores) else NULL
    if (!is.null(clusters)) doParallel::registerDoParallel(clusters) else foreach::registerDoSEQ()

    # Define training control
    fit.control <- caret::trainControl(method=method.resampling, allowParallel=T,
                                       returnResamp = "all",
                                       savePredictions = "all",
                                       number=3, repeats=3, search="random")
    # Define tuning grid with reduced search space
    tune.grid <- expand.grid(shrinkage = seq(0.1, 1, by = 0.3),
                             interaction.depth = c(1, 5),
                             n.minobsinnode = c(2, 5),
                             n.trees = c(100, 300, 1000))

    ##Gradient Boosting
    model_<- caret::train(fmla.n, data = df.train[,],  method = "gbm", metric='RMSE',
                         # preProc = c('center', 'scale','BoxCox', 'YeoJohnson', 'expoTrans', 'ica'),
                         trControl = fit.control, tuneGrid =tune.grid,verbose = F)

    if (!is.null(clusters)) parallel::stopCluster(clusters)
    # estimate variable importance
    importance <- NA #caret::varImp(model_, scale=FALSE)

  }  else if (algorithm == 'NN'){

    print('processing hybrid approach using a simple Neural Network')

    fmla.n <- as.formula(paste(depVar," ~ ", paste(inputs, collapse= "+")))

    # Parallelize tuning process
    clusters <- if (n.cores > 1) parallel::makeCluster(n.cores) else NULL
    if (!is.null(clusters)) doParallel::registerDoParallel(clusters) else foreach::registerDoSEQ()

    # Define training control
    # Was: index = caret::createFolds(df.train[,inputs], 5) -- createFolds()
    # expects a single outcome vector to stratify on, not the full ~thousands-
    # of-columns predictor matrix (`inputs`). That produced degenerate/garbage
    # fold assignments, which is why every resample's Rsquared came back NA
    # ("Something is wrong; all the Rsquared metric values are missing").
    fit.control <- caret::trainControl(method=method.resampling, allowParallel=T,
                                       number=3, repeats=3, search="random",
                                       index = caret::createFolds(df.train[,depVar], 5),
                                      # sampling ='smote',
                                     # na.action = na.pass,
                                       returnResamp = "all",
                                       savePredictions = "all")
    # Define tuning grid with reduced search space
    nnet.grid <- expand.grid(.decay = seq(0,0.1,by=0.01),
                             .size = seq(1,10,by=1))#, .bag=F) for vNNet

    ##Neural-Network Model
    model_<- caret::train(fmla.n, data = df.train[,],
                             method = "nnet", repeats = 1, trControl = fit.control,
                             preProc = c("center", "scale",'BoxCox', 'YeoJohnson'),
                             cross=10,trace=F, ##remove message with Trace=False
                             threshold = 0.3,
                             na.action  = na.pass,
                             maxit = 1000, linout = 1,
                             metric='Rsquared', tuneGrid = nnet.grid,
                             verbose = FALSE)

    size<-getElement(model_,"bestTune")$size
    decay<-getElement(model_,"bestTune")$decay

    if (!is.null(clusters)) parallel::stopCluster(clusters)
    # estimate variable importance
    importance <- caret::varImp(model_, scale=FALSE)


  }  else if (algorithm == 'Bayesian') {
    print('processing Bayesian Additive Regression Trees ...')

    fmla.n <- as.formula(paste(depVar," ~ ", paste(inputs, collapse= "+")))

    # Parallelize tuning process
    clusters <- if (n.cores > 1) parallel::makeCluster(n.cores) else NULL
    if (!is.null(clusters)) doParallel::registerDoParallel(clusters) else foreach::registerDoSEQ()

   # Define training control
    fit.control <- caret::trainControl(method=method.resampling, allowParallel=F,
                                       returnResamp = "all",
                                       savePredictions = "all",
                                       number=3, repeats=3, search="random")
    # Define tuning grid
    tune.grid <- expand.grid(num_trees = c(50),
                             k = c(0.5, 1),
                             alpha = c(0.95, 0.99),
                             beta = c(1.0, 2.0),
                             nu = c(3,5,7))

    ## Bayesian Generalized Linear Model
    model_ <- caret::train(fmla.n, data = df.train[,], method = "bartMachine", metric='RMSE',
                           trControl = fit.control, tuneGrid = tune.grid,verbose = FALSE)

    if (!is.null(clusters)) parallel::stopCluster(clusters)
    # estimate variable importance
    importance <- caret::varImp(model_, scale=FALSE)

  } else if (algorithm == 'AdaBag') {
    print('processing Bagged AdaBoost ...')

    # Define formula
    fmla.n <- as.formula(paste(depVar," ~ ", paste(inputs, collapse= "+")))

    # Parallelize tuning process
    clusters <- if (n.cores > 1) parallel::makeCluster(n.cores) else NULL
    if (!is.null(clusters)) doParallel::registerDoParallel(clusters) else foreach::registerDoSEQ()

    # Define training control
    fit.control <- caret::trainControl(method = method.resampling,
                                       number = 3, repeats = 3,
                                       search = "random",
                                       allowParallel = TRUE)

    # Define tuning grid
    tune.grid <- expand.grid(mfinal = c(50, 100, 200),
                             maxdepth = c(1, 3, 5))

    ## Bagged AdaBoost
    model_ <- caret::train(fmla.n, data = df.train[,],
                           method = "AdaBag",
                           trControl = fit.control,
                           tuneGrid = tune.grid,
                           metric = 'RMSE',  # Use RMSE as the evaluation metric
                           verbose = FALSE)

    if (!is.null(clusters)) parallel::stopCluster(clusters)

    # estimate variable importance
    importance <- caret::varImp(model_, scale=FALSE)

  } else if (algorithm == 'BRNN') {
    print('processing Bayesian Regularized Neural Networks ...')

    # Define formula
    fmla.n <- as.formula(paste(depVar, " ~ ", paste(inputs, collapse = "+")))

    # Parallelize tuning process
    clusters <- if (n.cores > 1) parallel::makeCluster(n.cores) else NULL
    if (!is.null(clusters)) doParallel::registerDoParallel(clusters) else foreach::registerDoSEQ()

    # Define training control
    fit.control <- caret::trainControl(method = method.resampling,
                                       number = 3, repeats = 3,
                                       search = "random",
                                       allowParallel = TRUE)

    # Define tuning grid
    tune.grid <- expand.grid(neurons = c(5, 10, 20, 30))

    ## Bayesian Regularized Neural Networks (BRNN)
    model_ <- caret::train(fmla.n, data = df.train[,],
                           method = "brnn",
                           trControl = fit.control,
                           tuneGrid = tune.grid,
                           metric = 'RMSE',  # Use RMSE as the evaluation metric
                           verbose = FALSE)

    if (!is.null(clusters)) parallel::stopCluster(clusters)

    # estimate variable importance
    importance <- caret::varImp(model_, scale=FALSE)

  } else if (algorithm == 'xGB') {
    print('processing eXtreme Gradient Boosting (XGBoost) with linear base learners ...')

    # Define formula
    fmla.n <- as.formula(paste(depVar, " ~ ", paste(inputs, collapse = "+")))

    # Parallelize tuning process
    clusters <- if (n.cores > 1) parallel::makeCluster(n.cores) else NULL
    if (!is.null(clusters)) doParallel::registerDoParallel(clusters) else foreach::registerDoSEQ()

    # Define training control
    fit.control <- caret::trainControl(method = method.resampling,
                                       number = 3, repeats = 3,
                                       search = "random",
                                       allowParallel = TRUE)

    # Define tuning grid
    tune.grid <- expand.grid(nrounds = c(50, 100, 200),
                             lambda = c(0, 0.01, 0.1),
                             alpha = c(0, 0.01, 0.1),
                             eta = c(0.01, 0.05, 0.1))

    ## eXtreme Gradient Boosting (XGBoost) with linear base learners
    model_ <- caret::train(fmla.n, data = df.train[,],
                           method = "xgbLinear",
                           trControl = fit.control,
                           tuneGrid = tune.grid,
                           metric = 'RMSE',  # Use RMSE as the evaluation metric
                           verbose = FALSE)

    if (!is.null(clusters)) parallel::stopCluster(clusters)

    # estimate variable importance
    importance <- NA


  } else if (algorithm == 'RVM') {
    print('processing Relevance Vector Machines (RVM) with linear kernel ...')

    # Define formula
    fmla.n <- as.formula(paste(depVar, " ~ ", paste(inputs, collapse = "+")))

    # Parallelize training process if applicable
    clusters <- if (n.cores > 1) parallel::makeCluster(n.cores) else NULL
    if (!is.null(clusters)) doParallel::registerDoParallel(clusters) else foreach::registerDoSEQ()

    # Define training control
    fit.control <- caret::trainControl(method = method.resampling,
                                       number = 3, repeats = 3,
                                       search = "random",
                                       allowParallel = TRUE)

    ## Relevance Vector Machines (RVM) with linear kernel
    model_ <- caret::train(fmla.n, data = df.train[,],
                           method = "rvmLinear",
                           trControl = fit.control,
                           metric = 'RMSE',  # Use RMSE as the evaluation metric
                           verbose = FALSE)



    if (!is.null(clusters)) parallel::stopCluster(clusters)
    # estimate variable importance
    importance <- caret::varImp(model_, scale=FALSE)


  } else if (algorithm == 'qLASSO') {
    print('processing Quantile Regression with LASSO penalty ...')

    # Define formula
    fmla.n <- as.formula(paste(depVar, " ~ ", paste(inputs, collapse = "+")))

    # Parallelize training process if applicable
    clusters <- if (n.cores > 1) parallel::makeCluster(n.cores) else NULL
    if (!is.null(clusters)) doParallel::registerDoParallel(clusters) else foreach::registerDoSEQ()

    # Define training control
    fit.control <- caret::trainControl(method = method.resampling,
                                       number = 3, repeats = 3,
                                       search = "random",
                                       allowParallel = TRUE)

    # Define tuning grid
    tune.grid <- expand.grid(lambda = c(0.01, 0.1, 1, 10))

    ## Quantile Regression with LASSO penalty
    model_ <- caret::train(fmla.n, data = df.train[,],
                           method = "rqlasso",
                           trControl = fit.control,
                           tuneGrid = tune.grid,
                           metric = 'RMSE',  # Use RMSE as the evaluation metric
                           verbose = FALSE)

    if (!is.null(clusters)) parallel::stopCluster(clusters)

    # estimate variable importance
    importance <- caret::varImp(model_, scale=FALSE)

  } else if (algorithm == 'Ensemble') {
    print('processing Ensemble approach by stacking 3 models ....')
    print('list of models: SVM,Gradient Boosting and Neural Network ')


    # Parallelize tuning process
    clusters <- if (n.cores > 1) parallel::makeCluster(n.cores) else NULL
    if (!is.null(clusters)) doParallel::registerDoParallel(clusters) else foreach::registerDoSEQ()

    algorithmList <- c('gbm', #Gradient-boosted machines
                       'svmRadial', #SVM with RBF Kernel
                       'nnet') #neural network

    ### tuning parameters for each model

    # gmb
    tune.grid.gbm <- expand.grid(shrinkage = seq(0.1, 1, by = 0.3),
                                 interaction.depth = c(1, 5),
                                 n.minobsinnode = c(2, 5),
                                 n.trees = c(100, 300, 1000))
    # svm
    tuneGrid.svm = expand.grid(C = c(2^(1:4)),sigma=c(2^(-4:1)))
    # nnet
    tune.grid.nne <- expand.grid(.decay = seq(0.001,0.2,by=0.01), .size = seq(1,10,by=1))

    ### Ensemble mddels

    models.ensemble=list(gbm=caretEnsemble::caretModelSpec(method="gbm",  #metric='MAE',
                                            preProc = c("center", "scale",'BoxCox', 'YeoJohnson'),
                                            tuneGrid=tune.grid.gbm),
                         svmRadial=caretEnsemble::caretModelSpec(method="svmRadial", # metric='MAE',
                                                  preProc = c("center", "scale",'BoxCox', 'YeoJohnson'),
                                                  tuneGrid=tuneGrid.svm,#tuneLength=10,
                                                  threshold = 0.3),
                         #tuneGrid=tuneGrid.svm),
                         nnet=caretEnsemble::caretModelSpec(method="nnet", #metric='MAE',
                                             preProc = c("center", "scale",'BoxCox', 'YeoJohnson'),
                                             threshold = 0.3,
                                             tuneGrid=tune.grid.nne,
                                             maxit = 1000, linout = 1, trace=F))

    fit.control<- caret::trainControl(method=method.resampling,
                               number=3,
                               savePredictions=TRUE,
                               #index = caret::createFolds(data.train[,input], 5),
                               repeats=3,
                               search = "random")

    models <- caretEnsemble::caretList(fmla.n, data = df.train[,],
                                       trControl=fit.control,
                                       verbose=FALSE,
                                       tuneList = models.ensemble,
                                       methodList=algorithmList)
    # Combine Predictions from multiple models
    stackControl <- caret::trainControl(method=method.resampling,
                                 number=3,
                                 repeats=3,
                                 index = caret::createFolds(df.train[,], 5),
                                 savePredictions = "all",
                                 search = "random")

    # Ensemble the predictions of `models` to form a new combined prediction based on glm
    model_ <- caretEnsemble::caretStack(models, method="glm", trControl=stackControl)

    if (!is.null(clusters)) parallel::stopCluster(clusters)

    # estimate variable importance
    importance <- NA

  } else {
    stop("Invalid model type. Choose from 'PLSR', 'SVM', 'RF', 'GB',  'xGB', 'NN', 'qLASSO','Bayesian','AdaBag','BRNN', 'RVM', or 'Ensemble.")
  }



  pred.train <- stats::predict(object = model_, df.train[, c(depVar, inputs)])
  pred.test <- stats::predict(object = model_, df.test[, c(depVar, inputs)])

  # Prediction
  n.rows_.train <- as.numeric(names(pred.train))
  n.rows_.test <- as.numeric(names(pred.test))

  # Statistics
  stats <- get.stats(model = model_, train = df.train[, c(depVar, inputs)],
                     test = df.test[, c(depVar, inputs)], var = depVar)

  # Plot
  plot.result <- get.plot.ML(model = model_, df.train[, c(depVar, inputs)],
                             df.test[, c(depVar, inputs)], depVar)

  # Save the model, statistics, and plot if the save.model flag is TRUE
  if (save.model) {
    # Save the model as an RDS file
    saveRDS(model_, file = file.path(save.path, paste0('model-', algorithm, '.rds')))

    # Write the statistics as a CSV file
    write.csv(stats, file = file.path(save.path, paste0('stats-', algorithm, '.csv')), row.names = FALSE)

    # Save the plot as a PNG file
    ggplot2::ggsave(filename = file.path(save.path, paste0('Plot-', algorithm, '.png')), plot = plot.result)
  }


  results <- list(
    model.label = algorithm,
    model= model_,
    predictions = list(train = pred.train, test = pred.test),
    statistics = stats,
    plot = plot.result,
    importance = importance
  )

  return(results)
}
