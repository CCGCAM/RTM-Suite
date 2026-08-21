
#' Inversion of plant traits using ML models (Ensemble variant)
#'
#' @param LUT  Dataset with inputs and Bands
#' @param input variable to estimate
#' @param split  ratio between 0 and 1 for splitting the dataset in training and testing
#' @param setseed  set random number
#' @param pattern  Please indicate the number of bands with same pattern'B'.
#' @param collinearity  collinearity-and-stepwise-vif-selection or CARS method implemmented, options='VIF' and 'CARS'
#' @param Field.data  dataframe with observations
#' @param acron acronynm for the observation measure: e.g., Cab_obsrv, where acron='_observ' and Cab has same name as input
#' @export
#'
#' @return A list: \code{model} (the fitted ensemble model), \code{Stats} (accuracy
#'   statistics on the held-out test split), and \code{Plot} (a ggplot scatter of
#'   predicted vs. measured/simulated values for \code{input}).
hybrid_inversionE<-function(LUT=NULL,input=NULL,split=0.8,setseed=NULL,
                           collinearity=NULL, pattern=NULL,
                           Field.data=NULL, acron=NULL){
  if   (!is.null(Field.data)){
    if (is.null(acron)){
      message('please add the acronym for the ground data')
      stop()
    }
  }
   # models used for generating the ensemble predictios

  options(warn=-1) ###avoid warnings
  dataset=LUT
  xi_y1=split
  ### partition datasets
  index <- as.vector(caret::createDataPartition(dataset[,input], p = xi_y1, list = F))
  
  data.train<-dataset[index,]
  data.test<-dataset[-index,]
  if (is.null(pattern)){
    keep.variables<-names(dataset)
    input_order = grep(colnames(dataset),pattern=input,fixed = TRUE)
    keep.variables[-input_order]
    fmla <- as.formula(paste(input," ~ ", paste(keep.variables, collapse= "+")))
    
  }
  if (is.null(collinearity)) {
    
    keep.variables<-names(dataset[,grep(colnames(dataset),pattern=pattern,fixed = TRUE)])
    fmla <- as.formula(paste(input," ~ ", paste(keep.variables, collapse= "+")))
    
  } else if(collinearity == 'VIF') {
      ## require fmsb package
      inputs_names<-names(dataset[,grep(colnames(dataset),pattern=pattern,fixed = TRUE)])
      keep.variables<-getVIF(dataset[,inputs_names],thresh=10,trace=T)  
      fmla <- as.formula(paste(input," ~ ", paste(keep.variables, collapse= "+")))
 
  } else if(collinearity == 'CARS') {
  ## require pracma and pls packages
    inputs_names<-names(dataset[,grep(colnames(dataset),pattern=pattern,fixed = TRUE)])
    carspls_selection<-carspls(data.train[,inputs_names],y=data.train[,input],nLV=5,fold=10,scale.pretreat=1,iteration=100,PartitionType="interleaved")
    keep.variables<-colnames(data.train[,inputs_names])[c(carspls_selection$SelectedVariables)]
    fmla <- as.formula(paste(input," ~ ", paste(keep.variables, collapse= "+")))
  
  }
  
set.seed(setseed)
list.models<-list
message('Step-1: processing hybrid inversion using SVM model ....')

clusters <- max(1, parallel::detectCores()-1)
if (nzchar(Sys.getenv("_R_CHECK_LIMIT_CORES_"))) clusters <- min(clusters, 2L)
cl <- parallel::makePSOCKcluster(clusters)
doParallel::registerDoParallel(cl)
set.seed(setseed)
## Support Vector Machines (SVM)
tobj <- e1071::tune.svm(fmla, data = data.train, gamma = 2^(-4:1), cost = 2^(1:4),
                 tunecontrol=e1071::tune.control(cross=10), # by default scale=T
                 parallel.cores =cl )

cc <- as.numeric(tobj$best.parameters[2]) ##cost
gg <- as.numeric(tobj$best.parameters[1]) ### gamma

list.models[['svm']]<- e1071::svm(fmla, kernel="radial", data = data.train, gamma= gg, cost= cc,
            cross=10,  # k-fold de 10 reduce overffiting
            probability=F,  ## probabilidad
            fitted=T) ### predichos
parallel::stopCluster(cl)

message('Step-2: processing hybrid approach using Random Forest model')
clusters <- max(1, parallel::detectCores()-1)
if (nzchar(Sys.getenv("_R_CHECK_LIMIT_CORES_"))) clusters <- min(clusters, 2L)
cl <- parallel::makePSOCKcluster(clusters)
set.seed(setseed)
##random forest model RF
#mtry: Number of variables randomly sampled as candidates at each split.
#ntree: Number of trees to grow.
# Random Search
fit.control <- caret::trainControl(method="repeatedcv", allowParallel=T,
                       returnResamp = "all",
                       savePredictions = "all",
                      number=10, repeats=3, search="random")

tune.grid <- expand.grid(.mtry = c(1: dim(data.train[,keep.variables])[2]-1))
#tune.grid <- expand.grid(.mtry=sqrt(ncol(data.train))*2)

list.models[['rf']] <- caret::train(fmla, data=data.train, method="rf", metric='RMSE', 
               # preProc = c('center', 'scale','BoxCox', 'YeoJohnson', 'expoTrans'),
                trControl=fit.control, tuneGrid=tune.grid, ntree=500)
## for ntree search in caret package
## not implemented for the computing time
# #results_rf <- list()
# for (ntree in c(100,200,500,600,800, 1000, 2000)) {
#   set.seed(setseed)
#   m_trees <- train(fmla,data = data_train,method = "rf",metric='RMSE',
#                       preProcess = c("center", "scale"),
#                       trControl=fit.control, tuneGrid=tune.grid, 
#                        nodesize = 14,
#                        maxnodes = 24,
#                        ntree = ntree)
#   i_model <- toString(ntree)
#   results_rf[[i_model]] <- m_trees
# }
# results_tree <- resamples(store_maxtrees)
# summary(results_tree)

parallel::stopCluster(cl)

message('Step-3: processing hybrid approach using Gradient Boosting')
clusters <- max(1, parallel::detectCores()-1)
if (nzchar(Sys.getenv("_R_CHECK_LIMIT_CORES_"))) clusters <- min(clusters, 2L)
cl <- parallel::makePSOCKcluster(clusters)
##Gradient Boosting
set.seed(setseed)
fit.control <- caret::trainControl(method="repeatedcv", allowParallel=T,
                     returnResamp = "all",
                     savePredictions = "all",
                    number=10, repeats=3, search="random")
tune.grid <- expand.grid(shrinkage = seq(0.1, 1, by = 0.2), 
                  interaction.depth = c(1, 3, 7, 10),
                  n.minobsinnode = c(2, 5, 10),
                  n.trees = c(100, 300, 500, 1000))
list.models[['gbm']]<- caret::train(fmla, data = data.train,  method = "gbm", metric='RMSE',
             # preProc = c('center', 'scale','BoxCox', 'YeoJohnson', 'expoTrans', 'ica'),
              trControl = fit.control, tuneGrid =tune.grid, verbose = FALSE)


parallel::stopCluster(cl)

message('Step-4: processing hybrid approach using a nnet')
clusters <- max(1, parallel::detectCores()-1)
if (nzchar(Sys.getenv("_R_CHECK_LIMIT_CORES_"))) clusters <- min(clusters, 2L)
cl <- parallel::makePSOCKcluster(clusters)
##nnet Model 
set.seed(setseed)
fit.control <- caret::trainControl(method="repeatedcv", allowParallel=T,
                                   number=10, repeats=3, search="random",
                                   index = caret::createFolds(data.train[,input], 5),
                                   returnResamp = "all",
                                   savePredictions = "all")


nnet.grid <- expand.grid(.decay = seq(0,0.1,by=0.01), .size = seq(1,10,by=1))
nnet.fit <- caret::train(fmla, data = data.train,method = "nnet",  trControl = fit.control,
                   preProc = c("center", "scale"),
                    cross=10,
                    threshold = 0.3,
                    metric='Rsquared',maxit = 100, tuneGrid = nnet.grid,linout=TRUE,
                         verbose = FALSE) 
size<-getElement(nnet.fit,"bestTune")$size
decay<-getElement(nnet.fit,"bestTune")$decay

#nne model with the best tune parameters for 500 iterations
model <- nnet::nnet(fmla,data=data.train,size=size,decay=decay,trace=F,linout=TRUE,skip=T,
                    maxit=200)
best.value <- model$value
value <- NULL
progress_bar = txtProgressBar(min=0, max=1000, style = 3, char="=")
for(i in 1:1000)
{
  aux.nnet <- nnet::nnet(fmla,data=data.train,size=size,decay=decay,trace=F,linout=TRUE,skip=F,
                         maxit=200)
  value[i] <- aux.nnet$value
  if(aux.nnet$value < best.value)
  {
    model <- aux.nnet
    best.value <- model$value
  }
  setTxtProgressBar(progress_bar, value = i)
}
close(progress_bar)
list.models[['nnet']] = model
parallel::stopCluster(cl)

message('Step 5: processing hybrid approach using Ensemble by stacking approach ....')
message('list of models: SVM,Gradient Boosting and Neural Network ')

clusters <- max(1, parallel::detectCores()-1)
if (nzchar(Sys.getenv("_R_CHECK_LIMIT_CORES_"))) clusters <- min(clusters, 2L)
cl <- parallel::makePSOCKcluster(clusters)
set.seed(setseed)

algorithmList <- c('gbm', #Gradient-boosted machines
                   'svmRadial', #SVM with RBF Kernel
                   'nnet') #neural network
### tuning parameters for each model

# gmb
tune.grid.gbm <- expand.grid(shrinkage = seq(0.1, 1, by = 0.2), 
                         interaction.depth = c(1, 3, 7, 10),
                         n.minobsinnode = c(2, 5, 10),
                         n.trees = c(100, 300, 500, 1000))
# svm
tuneGrid.svm = expand.grid(C = c(2^(1:4)),sigma=c(2^(-4:1)))
# nnet
tune.grid.nne <- expand.grid(.decay = seq(0.001,0.2,by=0.01), .size = seq(1,10,by=1))

### tunning in a list
#in caret PreProc can be:
#BoxCox, YeoJohnson, expoTrans, invHyperbolicSine, center, scale, range, 
#nnImpute, bagImpute, medianImpute, pca, ica, spatialSign, ignore, keep, 
#remove, zv, nzv, conditionalX, corr
models.ensemble=list(gbm=caretEnsemble::caretModelSpec(method="gbm",  #metric='MAE', 
                            #   preProc = c("center", "scale"),
                                 tuneGrid=tune.grid.gbm),
              svmRadial=caretEnsemble::caretModelSpec(method="svmRadial", # metric='MAE',
                         # preProc = c("center","scale"),
                          tuneGrid=tuneGrid.svm,#tuneLength=10,
                           threshold = 0.3),
                           #tuneGrid=tuneGrid.svm),
              nnet=caretEnsemble::caretModelSpec(method="nnet", #metric='MAE',
                               #   preProc = c("center", "scale"),
                                  threshold = 0.3,
                                  tuneGrid=tune.grid.nne,
                                  #tuneGrid = tune.grid.nne,
                                  linout=TRUE,
                                  trace=FALSE))
fit.control<- caret::trainControl(method="repeatedcv", 
                          number=10, 
                          savePredictions=TRUE,
                          #index = caret::createFolds(data.train[,input], 5),
                          repeats=3,
                          search = "random")

models <- caretEnsemble::caretList(fmla, data = data.train, 
                                   trControl=fit.control,
                                   verbose=FALSE,
                                   tuneList = models.ensemble,
                                   methodList=algorithmList) 
# Combine Predictions from multiple models
set.seed(setseed)
stackControl <- caret::trainControl(method="repeatedcv", 
                             number=10, 
                             repeats=3,
                             index = caret::createFolds(data.train[,input], 5),
                             savePredictions = "all",
                             search = "random")

# Ensemble the predictions of `models` to form a new combined prediction based on glm
list.models[['ensemble']] <- caretEnsemble::caretStack(models, method="glm", trControl=stackControl)



  

## Predictions on train and test
pred.train<-c(predict(object = model,data.train))
pred.test<-c(predict(object = model,data.test))

#### Skill scores for training data
r2.train<-round(cor(pred.train,data.train[input],use='pairwise.complete.obs')^2,2)
rmse.train<-round(ToolsRTM::RMSE(pred.train,data.train[,input]),2)
mae.train<-round(ToolsRTM::MAE(pred.train,data.train[,input]),2)
#### Skill scores for testing data
r2.test<-round(cor(pred.test,data.test[,input],use='pairwise.complete.obs')^2,2)
rmse.test<-round(ToolsRTM::RMSE(pred.test,data.test[,input]),2)
mae.test<-round(ToolsRTM::MAE(pred.test,data.test[,input]),2)
## add skill scores in table
stats<-data.frame(r2=c(r2.train,r2.test),
                  rmse=c(rmse.train,rmse.test),
                  mae=c(mae.train,mae.test))
### Some input for scatter-plots

mylabel.r.test = bquote(bold(r)^2 == .(format(stats[2,1], digits = 3)))
mylabel.rmse.test = bquote(bold(rmse) == .(format(stats[2,2], digits = 3)))

statsLabel = paste0("r2 = ", round(stats[2,1],2), ", RMSE = ", round(stats[2,2],4))

## save results in data frame for plotting  
data.plot<-list()
data.plot$input<-c(data.test[,input])
data.plot$pred<-pred.test

data.plot<-data.frame(do.call(cbind,data.plot))

axis_x<-bquote(bold(.(input)['measured'])) # axis x
axis_y<-bquote(bold(.(input)['predicted']))# axis y

scatter_plot<-ggplot(data.plot, aes(y=pred, x=input)) +
  geom_point(alpha=0.6) + geom_smooth(method=lm, formula = 'y ~ x') + theme_bw()+
  geom_abline(intercept = 0, slope = 1,linetype="dashed", size=0.5,color='gray')+
  coord_fixed(ratio = 1,xlim = c(0, max(data.plot$input)), ylim = c(0, max(data.plot$pred))) +
  xlab(axis_x) + ylab(axis_y) + ggtitle(statsLabel) 

if (IncludeModel == T) {
}
  
  

if (is.null(Field.data)) {
  message('no field data is added ....')
  #################

  Hybrid = list('model'=model,'Stats'=stats,'Plot'=scatter_plot)
  
  return(Hybrid)
  
} else {
 
  predict.obs<-c(predict(object = model,Field.data[,keep.variables]))
  Field.data$pred<-predict.obs
  names_f<-names(Field.data)
  colnames(Field.data)<-c(names(Field.data)[1:(length(names_f)-1)],paste(input,'_pred',sep=''))
  #### Skill scores for training data
  r2.obs<-round(cor(predict.obs,Field.data[,paste(input, '_obsv',sep='')],use='pairwise.complete.obs')^2,2)
  rmse.obs<-round(ToolsRTM::RMSE(predict.obs,Field.data[,paste(input,acron,sep='')]),2)
  mae.obs<-round(ToolsRTM::MAE(predict.obs,Field.data[,paste(input,acron,sep='')]),2)
  
  stats<-data.frame(r2=c(r2.train,r2.test,r2.obs),
                    rmse=c(rmse.train,rmse.test,rmse.obs),
                    mae=c(mae.train,mae.test,mae.obs))
  
  mylabel.r.test = bquote(bold(r)^2 == .(format(stats[3,1], digits = 3)))
  mylabel.rmse.test = bquote(bold(rmse) == .(format(stats[3,2], digits = 3)))
  statsLabel = paste0("r2 = ", round(stats[3,1],2), ", RMSE = ", round(stats[3,2],4))
  ## save results in data frame for plotting  
  data.plot<-list()
  data.plot$input<-c(Field.data[,paste(input,acron,sep='')])
  data.plot$pred<-predict.obs
  data.plot<-data.frame(do.call(cbind,data.plot))
  
  axis_x<-bquote(bold(.(input)['measured at Field level'])) # axis x
  axis_y<-bquote(bold(.(input)['predicted']))# axis y
  
  scatter_obs<-ggplot(data.plot, aes(y=pred, x=input)) +
    geom_point(alpha=0.6) + geom_smooth(method=lm, formula = 'y ~ x') + theme_bw()+
    geom_abline(intercept = 0, slope = 1,linetype="dashed", size=0.5,color='gray')+
    coord_fixed(ratio = 1,xlim = c(0, max(data.plot$input)), ylim = c(0, max(data.plot$pred))) +
    xlab(axis_x) + ylab(axis_y) + ggtitle(statsLabel) 
  
  
  
  
  Hybrid = list('model'=model,'Stats'=stats,'Plot'=scatter_plot,'Plot_field'=scatter_obs, 'Field.pred'=Field.data)
  return(Hybrid)
  
}
  
}

# MAE() and RMSE() are defined once, in hybrid_inversion.R — removed the
# duplicate copies that used to live here (identical code, was exported
# twice under the same names, which is invalid for a single package).
