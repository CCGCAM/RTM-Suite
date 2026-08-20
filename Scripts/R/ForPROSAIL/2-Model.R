
rm(list=ls())

##############################################################################################################################
#	0. load main Libraries   -----
##############################################################################################################################

if (!require("RColorBrewer")) { install.packages("RColorBrewer"); require("RColorBrewer") }  ### colors
if (!require("caret")) { install.packages("caret"); require("caret") }  ### findCorrelation()

if (!require("dplyr")) { install.packages("dplyr"); require("dplyr") }  ###
if (!require("tidyverse")) { install.packages("tidyverse"); require("tidyverse") }  ###

if (!require("keras")) { install.packages("keras"); require("keras") }  ###
if (!require("tensorflow")) { install.packages("tensorflow"); require("tensorflow") }  ###
# NOTE: hsdar / fastICA / GGally / corrgram / skimr were required here but never
# actually used anywhere in this script -- dropped to cut install weight. hsdar
# in particular was archived from CRAN on 2023-07-06 and can no longer be
# installed via install.packages() on current R versions anyway.

if (!requireNamespace("ToolsRTM", quietly = TRUE)) {
  remotes::install_gitlab("caminoccg/toolsrtm")
}
library(ToolsRTM)

out_root <- "outs/ForPROSAIL"  # outputs always go under <repo root>/outs/, never inside Scripts/

##############################################################################################################################
#	1. Open  LUT with correlations  -----
##############################################################################################################################

lut_path <- file.path(out_root, "LUTs/LUT_PROSAIL_with_n20k_filtered.csv")
if (!file.exists(lut_path)) {
  stop("LUT not found at '", lut_path, "'. Run 1-GetSimulationsLUTs.R then 1-ReduceLUTs.R first ",
       "(both need Carlos' private field data, see Scripts/ForPROSAIL/data/README.md).")
}
LUT = read.table(lut_path, sep=',', header = T,fileEncoding="UTF-8-BOM")
nSamples =dim(LUT)[1]

print(nSamples)

LUT <-ToolsRTM::getIndices(LUT[,1:188], pattern.rfl = 'R.', spectral.domain = 'VNIR')
head(LUT)
LUT[,'CBC'] <- LUT[,'CBC'] * 1000

inputs.bands = names(LUT)[grep("R.", names(LUT))]
inputs.bands =inputs.bands[seq(1,length(inputs.bands),10)]
inputs.bands = inputs.bands[1:17]



field_data_path <- file.path(out_root, "FieldData/2019_2020RawData_Vcmax_withRFL_resampled.csv")
if (!file.exists(field_data_path)) {
  stop("Field data not found at '", field_data_path, "'. This is produced by running ",
       "1-GetSimulationsLUTs.R through step 9, which itself needs Carlos' private 2019/2020 ",
       "Vcmax field campaign data (see Scripts/ForPROSAIL/data/README.md).")
}
data.field = read.table(field_data_path, sep=',', header = T)

data.field <-ToolsRTM::getIndices(data.field[,1:173], pattern.rfl = 'R.', spectral.domain = 'VNIR')
head(data.field)

spectral.indices<-names(data.field)[174:dim(data.field)[2]]

##############################################################################################################################
#	 2. Configuration the inputs for the Model  -----
##############################################################################################################################

##############################################################################################################################
### avalaible options for getsplit function
##############################################################################################################################

data.trans <- c('PCA','preProcess')
depVar <- c('Cab','EWT','Car','Anth','LMA','CBC','LIDFa','LAI')
i=1 # to write model
method.preProcess <- c('Normalize', 'YeoJohnson','BoxCox', 'Standarize', 'Center','Scale', 'PCA',NULL)
depVar.trans=c(TRUE,FALSE)

inputs_ = c(inputs.bands,spectral.indices)


split.data<-ToolsRTM::getSplitData(data=LUT, depVar='Cab',inputs=inputs_,
                                   data.trans='preProcess',prop.split=c(0.8,0.2),method.preProcess='Normalize',
                                   depVar.trans=FALSE)
colnames(split.data[['Xtrain']])
names(split.data)
scaler.train <-  split.data[['Scalar.train']]
scaler.Ytrain <-  split.data[['Scalar.Ytrain']]

plot.train<-split.data[['plot.train']]



##############################################################################################################################
#	 3. Tensorflow and Keras Models  -----
##############################################################################################################################
n.epochs=200
batch.size=125
stats<-list()

# output folder
paths.models <- paste0(out_root, "/Models/")
ifelse(!dir.exists(paths.models), dir.create(paths.models, recursive = TRUE), FALSE)

##############################################################################################################################
# 3.1. CNN model ---
##############################################################################################################################

#Reshaping the data for CNN
### These dimensions don't look correct; switch ncol() with nrow()
data.Xtrain.CNN <- array_reshape(split.data[['Xtrain']], c(nrow(split.data[['Xtrain']]), ncol(split.data[['Xtrain']]), 1))
data.Xval.CNN <- array_reshape(split.data[['Xval']], c(nrow(split.data[['Xval']]), ncol(split.data[['Xval']]), 1))

# Create the configuration for the model

dataset.dim=data.Xtrain.CNN
model.cnn <- keras_model_sequential() %>%
  layer_conv_1d(filters=64, kernel_size=4, activation="relu",
                ### The input shape doesn't look correct; instead of
                ### `c(ncol(dataTrain_x), nrow(dataTrain_x))` (354, 13)
                ### I believe you want `dim(dataTest_x)` (13, 1)
                input_shape=c(ncol(dataset.dim), 1)) %>%
  layer_max_pooling_1d(pool_size=2) %>%
  layer_conv_1d(filters=32, kernel_size=2, activation="relu") %>%
  layer_max_pooling_1d(pool_size=2) %>%
  layer_dropout(rate=0.4) %>%
  layer_flatten() %>%
  layer_dense(units=16, activation="relu") %>%
  layer_dropout(rate=0.1) %>%
  layer_dense(units=1, activation="relu")

# Compile the configuration for the model

model.cnn %>% compile(loss = "mse",
  optimizer =  "adam", #'sgd' can also be used
  metrics = list("mean_absolute_error")
)

model.cnn %>% summary()
callbacks = list(callback_early_stopping(monitor = "loss", patience = 75, restore_best_weights = TRUE))
# fit the configuration for the model
history.cnn <- model.cnn %>% fit(data.Xtrain.CNN, split.data[['Ytrain']],
                         epochs = n.epochs,batch_size = 50, verbose=1,shuffle=F,
                         callbacks = callbacks,
                         validation_split = 0.2)

# evaluate the configuration for the model


stats[['CNN-model']] <-model.cnn %>% evaluate(data.Xval.CNN, split.data[['Yval']])
print(stats)

# Save the model
model.cnn %>% save_model_hdf5(paste(paths.models,'model_cnn_for_',depVar[i],'.hdf5',sep=''))
# Restore the weights
#model.cnn <- load_model_hdf5('Tables/Models/model_cnn.hdf5')


##############################################################################################################################
# 3.2. Deep ML with 3-hidden layers ---
##############################################################################################################################

ntrain = dim( split.data[['Xtrain']])[1]
batch.size = 125
step_epoch <- ntrain %/% batch.size

data.Xtrain.reshape <- array_reshape(split.data[['Xtrain']], c(nrow(split.data[['Xtrain']]), ncol(split.data[['Xtrain']])))
dim(data.Xtrain.reshape)

lr_schedule <- keras::learning_rate_schedule_inverse_time_decay( initial_learning_rate=0.001,
  decay_steps = step_epoch * 1000,decay_rate = 1, staircase = FALSE)

get_optimizer <- function() {
  optimizer_adam(lr_schedule)
}


## define Hidden layer
## training Data Samples/Factor * (Input Neurons + Output Neurons)
## where factor is 1 (A factor can take a value between 1 and 10)
ndim=dim(data.Xtrain.reshape)[2]
# Create the configuration for the model
model.3hlayers <- keras_model_sequential() %>%
  layer_dense(units=64,activation="relu",
                input_shape=c(dim(data.Xtrain.reshape)[2]))  %>%
  #layer_dropout(rate=0.1) %>%
  layer_dense(units = 32, activation = 'relu') %>%
  layer_dropout(rate=0.1) %>%
  layer_dense(units = 16, activation = 'relu') %>%
  #layer_dropout(rate=0.1) %>%
  layer_dense(units = 1, activation = 'relu')

# Compile the configuration for the model
model.3hlayers %>% compile(loss = "mse",
                  optimizer = 'adam',#,get_optimizer(),#"adam", #'sgd' can also be used
                  metrics = list("mean_absolute_error"))
model.3hlayers %>% summary()

callbacks = list(callback_early_stopping(monitor = "loss", patience = 75, restore_best_weights = TRUE))

# fit the configuration for the model
history.3hlayers <- model.3hlayers %>% fit(data.Xtrain.reshape, split.data[['Ytrain']],
  epochs = n.epochs, batch_size = batch.size, verbose=1,shuffle=F,callbacks =callbacks,
  validation_split = 0.2)

# evaluate the configuration for the model
model.3hlayers %>% evaluate(split.data[['Xval']], split.data[['Yval']])

stats[['3hlayer-model']] <- model.3hlayers %>% evaluate(split.data[['Xval']], split.data[['Yval']])
print(stats)
# Save the model
model.3hlayers %>% save_model_hdf5(paste(paths.models,'/model_3hlayers_for_',depVar[i],'.hdf5',sep=''))


# Restore the weights
#model.cnn <- load_model_hdf5('Tables/Models/model_cnn.hdf5')
##############################################################################################################################
# 3.3.  Predictions  X Val ---
##############################################################################################################################

Ytrans=FALSE
preds.3hl<-ToolsRTM::getPredicts(model=model.3hlayers, type.model='Hidden-layers',
                                 data=split.data[['Xval']], data.trans='preProcess',
                                 data.Y=split.data[['Yval']],
                                 depVar=depVar[i])#scaler.depVar=scaler.Ytrain)

preds.cnn<- ToolsRTM::getPredicts(model=model.cnn, type.model='CNN',
                                  data=split.data[['Xval']], data.trans='preProcess',
                                  data.Y=split.data[['Yval']],
                                  depVar=depVar[i])#scaler.depVar=scaler.Ytrain)


df.plot<-cbind(preds.3hl[,c('depVar','predicted.3hlayer')],preds.cnn[,'predicted.cnn'])
colnames(df.plot)<-c('depVar','predicted.3hlayer','predicted.cnn')
input_preds = names(df.plot)[grep("pred.", names(df.plot))]
df.plot[,'depVar.pred.avg']=rowMeans(df.plot[,c(input_preds)], na.rm=TRUE)
df.plot[,'depVar.pred.std']<-apply(df.plot[,c(input_preds)], 1, sd)
head(df.plot)


if (Ytrans == TRUE){
  Cab_to<-getReverse.trans(preProc=scaler.Ytrain,data=as.matrix(df.plot[,'depVar']))
  df.plot<-cbind(Cab_to,df.plot)
  colnames(df.plot)<-c('depVar','depVar.trans','predicted.3hlayer','predicted.cnn','depVar.pred.avg','depVar.pred.std')
  head(df.plot)
}
df.plot<-subset(df.plot, predicted.3hlayer <=70)
df.plot<-subset(df.plot, predicted.cnn <=70)
dim(df.plot)

MLmetrics::MAPE(df.plot[,'depVar'], df.plot[,'predicted.3hlayer'])
MLmetrics::MAPE(df.plot[,'depVar'], df.plot[,'predicted.cnn'])
MLmetrics::MAPE(df.plot[,'depVar'], df.plot[,'depVar.pred.avg'])


#############################################################################################################################
#	4. Plot vs Val dataset  -----
##############################################################################################################################

axis_x<-expr(paste('measured ', !!depVar[i],sep=''))# axis x
axis_y<-expr(paste('predicted ', !!depVar[i],sep=''))# axis y
max_ <- max(df.plot[,'depVar'])
#### scatterplot for ANN model

r2.SE2<-round(cor(df.plot[,'depVar'],df.plot[,'predicted.3hlayer'],use='pairwise.complete.obs')^2,2)
rmse.se<-round(ToolsRTM::RMSE(df.plot[,'depVar'],df.plot[,'predicted.3hlayer']),2)
mylabel.r.se2 = bquote(bold(SE2: r)^2 == .(format(r2.SE2, digits = 3)))

statsLabel = paste0("ANN model; r2 = ", round(r2.SE2,2), ", RMSE = ",round(rmse.se,2))

scatter.ann <-  ggplot(df.plot, aes_string(x='depVar', y='predicted.3hlayer')) +
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

print(scatter.ann)

paths.plots <- paste0(out_root, "/DeepML/")
ifelse(!dir.exists(paths.plots), dir.create(paths.plots, recursive = TRUE), FALSE)


ggsave(paste(paths.plots,'1-Model_predicted.3hlayer_for_',depVar[i],'.png',sep=''),
       width = 10, height = 10,  dpi = 300,units = "cm")


##############################################################################################################################
# 4.  Predictions  in observations  ---
##############################################################################################################################

data.field.to.predict <- as.matrix(predict(scaler.train, data.field[,inputs_]))
data.field.to.predict <- array_reshape(data.field.to.predict, c(nrow(data.field.to.predict), ncol(data.field.to.predict)))
dim(data.field.to.predict)
preds.field.3hl<-ToolsRTM::getPredicts(model=model.3hlayers, type.model='Hidden-layers',
                                       data=data.field.to.predict, data.trans='preProcess',
                                       data.Y=NULL,
                                       depVar=depVar[i])#scaler.depVar=scaler.Ytrain)


preds.field.cnn<-ToolsRTM::getPredicts(model=model.cnn, type.model='CNN',
                                       data=data.field.to.predict, data.trans='preProcess',
                                       data.Y=NULL,
                                       depVar=depVar[i])#,scaler.depVar=scaler.Ytrain)

df.field.plot<-cbind(data.field,preds.field.3hl[,c('predicted.3hlayer')],preds.field.cnn[,'predicted.cnn'])
colnames(df.field.plot)<-c(names(data.field),'predicted.3hlayer','predicted.cnn')
input_preds = names(df.field.plot)[grep("pred.", names(df.field.plot))]
df.field.plot[,'depVar.pred.avg']=rowMeans(df.field.plot[,c(input_preds)], na.rm=TRUE)
df.field.plot[,'depVar.pred.std']<-apply(df.field.plot[,c(input_preds)], 1, sd)
head(df.field.plot)


#############################################################################################################################
#	2.4 Plot vs testing dataset  -----
##############################################################################################################################

axis_x<-bquote(bold('measured Vcmax')) # axis x
axis_y<-bquote(bold('predicted Cab')) # axis x
paths.plots <- paste0(out_root, "/FieldData_plots/")
ifelse(!dir.exists(paths.plots), dir.create(paths.plots, recursive = TRUE), FALSE)


################################################################
#### scatterplot predcited Cab vs Vcmax
################################################################

r2.S<-round(cor(df.field.plot[,'Vcmax'],df.field.plot[,'predicted.cnn'],use='pairwise.complete.obs')^2,2)
mylabel.r.se2 = bquote(bold(r)^2 == .(format(r2.S, digits = 3)))
statsLabel = paste0("CNN-model; r2 = ", round(r2.S,2))

scatter.field <-  ggplot(df.field.plot, aes_string(x='Vcmax', y='predicted.cnn')) +
  geom_point(alpha=0.4,shape = 16,aes(), size=1.5) +  #geom_smooth(method=lm, formula = 'y ~ x', se=F,lty=2, lwd=1) +
  theme_bw()  + xlim(30,60) + ylim(30,60) +
  scale_color_gradient(low = "#0091ff", high = "#f0650e") +
  theme(legend.position="bottom",
        plot.title = element_text(hjust = 0.5, size=10,face="bold"),
        axis.title = element_text(face="bold", size=10),
        axis.text.y=element_text(hjust = 0.5, size=10,face="bold"),
        axis.text.x=element_text(hjust = 0.5, size=10,face="bold"),
        legend.title=element_blank()) +
  labs(title = statsLabel, x=axis_x,y=axis_y, size=10,face="bold") +
  stat_smooth(method = "lm",formula = y ~ x,geom = "smooth",col='red',lty=2,se=T)

print(scatter.field)


ggsave(paste(paths.plots,'1-ScatterPlot_Vcmax_vs_predicted_',depVar[i],'.png',sep=''),plot=scatter.field,
       width = 10, height = 10,  dpi = 300,units = "cm")

################################################################
#### scatter plot predicted Cab vs TCARI
################################################################

axis_x<-bquote(bold('TCARI')) # axis x
axis_y<-bquote(bold('predicted Cab')) # axis x

r2.S<-round(cor(df.field.plot[,'TCARI'],df.field.plot[,'predicted.cnn'],use='pairwise.complete.obs')^2,2)
mylabel.r.se2 = bquote(bold(r)^2 == .(format(r2.S, digits = 3)))
statsLabel = paste0("CNN-model; r2 = ", round(r2.S,2))

scatter.indices <- ggplot(df.field.plot, aes_string(x='TCARI', y='predicted.cnn')) +
  geom_point(alpha=0.4,shape = 16,aes(), size=1.5) +  #geom_smooth(method=lm, formula = 'y ~ x', se=F,lty=2, lwd=1) +
  theme_bw() +
  scale_color_gradient(low = "#0091ff", high = "#f0650e") +
  theme(legend.position="bottom",
        plot.title = element_text(hjust = 0.5, size=10,face="bold"),
        axis.title = element_text(face="bold", size=10),
        axis.text.y=element_text(hjust = 0.5, size=10,face="bold"),
        axis.text.x=element_text(hjust = 0.5, size=10,face="bold"),
        legend.title=element_blank()) +
  labs(title = statsLabel, x=axis_x,y=axis_y, size=10,face="bold") +
  stat_smooth(method = "lm",formula = y ~ x,geom = "smooth",col='red',lty=2,se=T)

print(scatter.indices)

ggsave(paste(paths.plots,'1-ScatterPlot_TCARI_vs_predicted_',depVar[i],'.png',sep=''),plot=scatter.indices,
       width = 10, height = 10,  dpi = 300,units = "cm")



################################################################
#### scatter plot Vcmax vs Vcmax
################################################################

axis_x<-bquote(bold('TCARI')) # axis x
axis_y<-bquote(bold('predicted Cab')) # axis x
#fit logistic regression model
model.vcmax <- glm(Vcmax ~ depVar.pred.avg + CUR + TCARI, data=df.field.plot)
df.field.plot[,'Vcmax.predicted'] <-predict(model.vcmax, df.field.plot)
#summary(model.gpp)

