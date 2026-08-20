
rm(list=ls())

##############################################################################################################################
#	0. load main Libraries   -----
##############################################################################################################################

if (!require("RColorBrewer")) { install.packages("RColorBrewer"); require("RColorBrewer") }  ### colors

# My packages in R
if (!require("ToolsRTM")) { install.packages("ToolsRTM"); require("ToolsRTM") }
if (!require("caret")) { install.packages("caret"); require("caret") }  ### findCorrelation()

if (!require("dplyr")) { install.packages("dplyr"); require("dplyr") }  ###
if (!require("tidyverse")) { install.packages("tidyverse"); require("tidyverse") }  ###

if (!require("keras")) { install.packages("keras"); require("keras") }  ###
if (!require("tensorflow")) { install.packages("tensorflow"); require("tensorflow") }  ###
# NOTE: hsdar / fastICA / GGally / corrgram / skimr were required here but never
# actually used anywhere in this script -- dropped to cut install weight. hsdar
# in particular was archived from CRAN on 2023-07-06 and can no longer be
# installed via install.packages() on current R versions anyway.

use.dev.source <- TRUE  # TRUE = load ToolsRTM from local source (ToolsRTM/R); FALSE = use the installed package
if (use.dev.source) devtools::load_all("../../../ToolsRTM/R")

out_root <- "../../../outs/ForPROSAIL"  # outputs always go under <repo root>/outs/, never inside Scripts/

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

inputs.bands = names(LUT)[grepl("R.", fixed=T,names(LUT))] ## same for grep
print(inputs.bands)
#inputs.bands =inputs.bands[seq(1,length(inputs.bands),10)]

##############################################################################################################################
#	2. Open  Field Data   -----
##############################################################################################################################

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
#	 3. Prepare the data for modelling  -----
##############################################################################################################################

depVar <- c('Cab','Car','Anth','CBC','LIDFa','LAI')

inputs.to = c(inputs.bands,spectral.indices)

LUT.to.reduce <- LUT[,inputs.to]
dim(LUT.to.reduce)
cor.LUT<-cor(LUT.to.reduce)

#the index of the columns to be removed because they have a high correlation
index <- findCorrelation(cor.LUT, cutoff=0.99)

#the name of the columns chosen above
to.be.removed <- colnames(cor.LUT)[index]
print(to.be.removed)
#now go back to df and use to_be_removed to subset the original df
LUT.to.reduce<-LUT.to.reduce[!names(LUT.to.reduce) %in% to.be.removed]

inputs.to = names(LUT.to.reduce)
print(inputs.to)

##############################################################################################################################
#	 4. Run Deep ML models  -----
##############################################################################################################################

## folders for outputs

## for plots (comparison with Field data)
paths.plots <- paste0(out_root, "/FieldData_plots/")
ifelse(!dir.exists(paths.plots), dir.create(paths.plots, recursive = TRUE), FALSE)
## for plots (comparison with simulations)
paths.outs <- paste0(out_root, "/DeepML/")
ifelse(!dir.exists(paths.outs), dir.create(paths.outs, recursive = TRUE), FALSE)

## Model

deepML <- 'Hidden-layers' ### options: 'CNN','Hidden-layers'
transformation ='preProcess' ## options: 'PCA','preProcess'
method.preProcess = 'Normalize' ## options: 'Normalize', 'YeoJohnson','BoxCox', Standarize', 'Center','Scale', and 'PCA'
optimizer.ml ='adamax' ### options: 'adam','adadelta','adagrad', 'adamax', 'nadam', 'msprop', 'sgd'
n.epochs=10
batch.size = 250

save.field<-list()


for (i in c(1:length(depVar))){
  print(depVar[i])
  print(summary(LUT[,depVar[i]]))

  ##############################################################################################################################
  ## 4.1. this function split and run a depp ML in a single step. For that function used a fixed structure (check code: 2.Model.R)
  ##############################################################################################################################

  paths.models <- paste0(out_root, "/Models/")
  ifelse(!dir.exists(paths.models), dir.create(paths.models, recursive = TRUE), FALSE)

  model.prep<-getMLmodel(dataset=LUT[,c(depVar[i],inputs.to)], depVar=depVar[i],model=deepML,
                         optimizer = optimizer.ml,
                         batch.size=batch.size,n.epochs=n.epochs, save.model=T, path.model=paths.models,
                         prop.split=c(0.8,0.2),
                         data.trans=transformation,method.preProcess=method.preProcess,
                         depVar.trans=FALSE)

  plot.val <- model.prep[['plot.val']]
  print(plot.val)
  ggsave(paste(paths.outs,'1-Model_predicted-',depVar[i],'-',deepML,'.png',sep=''),
         plot=plot.val,width = 10, height = 10,  dpi = 300,units = "cm")

  scaler.train = model.prep[['Scalar.train']]

  ##############################################################################################################################
  #	 2.2. Exploring dataset  -----
  ##############################################################################################################################

  data.field.to.predict <- as.matrix(predict(scaler.train, data.field[,inputs.to]))

  data.field.to.predict <- array_reshape(data.field.to.predict, c(nrow(data.field.to.predict), ncol(data.field.to.predict)))
  dim(data.field.to.predict)
  colnames(data.field.to.predict)<- names(data.field[,inputs.to])
  df.plot<-ToolsRTM::getPredicts(model=model.prep[['model']], type.model=deepML,
                                 data=data.field.to.predict, data.trans=transformation,
                                 data.Y=NULL,
                                 depVar=depVar[i])

  df.plot<- df.plot[, colSums(is.na(df.plot)) != nrow(df.plot)]

  colnames(df.plot) <- c(names(data.field[,inputs.to]),paste(depVar[i],'.predicted',sep=''))

  df.plot<-cbind(data.field[,1:6],df.plot)

  save.field[[i]]<- df.plot[,c(paste(depVar[i],'.predicted',sep=''))]

  ################################################################
  #### scatter plot predicted depVar vs Vcmax
  ################################################################



  depVar.predicted =paste(depVar[i],'.predicted',sep='')


  ################################################################
  #### scatterplot predcited Cab vs Vcmax
  ################################################################
  axis_x<-bquote(bold('measured Vcmax')) # axis X
  axis_y<-expr(paste('predicted ', !!depVar[i],sep=''))# axis y

  r2.S<-round(cor(df.plot[,'Vcmax'],df.plot[,depVar.predicted],use='pairwise.complete.obs')^2,2)
  mylabel.r.se2 = bquote(bold(r)^2 == .(format(r2.S, digits = 3)))
  statsLabel = paste0("Deep ML model; r2 = ", round(r2.S,2))

  scatter.vcmax <-  ggplot(df.plot, aes_string(x='Vcmax', y=depVar.predicted)) +
    geom_point(alpha=0.4,shape = 16,aes(), size=1.5) +
    theme_bw()  + xlim(30,90) +
    scale_color_gradient(low = "#0091ff", high = "#f0650e") +
    theme(legend.position="bottom",
          plot.title = element_text(hjust = 0.5, size=10,face="bold"),
          axis.title = element_text(face="bold", size=10),
          axis.text.y=element_text(hjust = 0.5, size=10,face="bold"),
          axis.text.x=element_text(hjust = 0.5, size=10,face="bold"),
          legend.title=element_blank()) +
    labs(title = statsLabel, x=axis_x,y=axis_y, size=10,face="bold") +
    stat_smooth(method = "lm",formula = y ~ x,geom = "smooth",col='red',lty=2,se=T)

  print(scatter.vcmax)
  ggsave(paste(paths.plots,'1-ScatterPlot_Vcmax_vs_predicted-',depVar[i],'-',deepML,'.png',sep=''),
         plot=scatter.vcmax,width = 10, height = 10,  dpi = 300,units = "cm")


  ################################################################
  #### scatter plot predicted VarDep vs spectral indices
  ################################################################


  if (depVar[i] == 'LAI'){
    inputs.to_<-inputs.to[ inputs.to %in% c('MCARI','MCARI2','SIFe1','MSAVI')]

  } else if (depVar[i] == 'Cab') {
    inputs.to_<-inputs.to[ inputs.to %in% c('TCARI','T.O','PSSRb','GM2')]
  } else if (depVar[i] == 'Car') {
    inputs.to_<-inputs.to[ inputs.to %in% c('DCabxc','CAR','SIFe1','GM2')]
  }  else if (depVar[i] == 'Anth') {
    inputs.to_<-inputs.to[ inputs.to %in% c('CRI700m','LIC1','PRIM2','PRI515')]
  } else if (depVar[i] == 'CBC') {
    inputs.to_<-inputs.to[ inputs.to %in% c('PSSRa','PRI515','MSR','CUR')]
  } else if (depVar[i] == 'LIDFa') {
    inputs.to_<-inputs.to[ inputs.to %in% c('TVI','SIFe1','SIFe1','SIFe4')]
  }

  if (!(is.null(inputs.to_[1]))){
    print(inputs.to_[1])
    axis_x<-expr(paste('predicted ', !!inputs.to_[1],sep=''))# axis y
    axis_y<-expr(paste('predicted ', !!depVar[i],sep=''))# axis y
    depVar.predicted =paste(depVar[i],'.predicted',sep='')

    r2.S<-round(cor(df.plot[,c(inputs.to_[1])],df.plot[,depVar.predicted],use='pairwise.complete.obs')^2,2)
    mylabel.r.se2 = bquote(bold(r)^2 == .(format(r2.S, digits = 3)))
    statsLabel = paste0("Deep ML model; r2 = ", round(r2.S,2))

    scatter.indices <-  ggplot(df.plot, aes_string(x=inputs.to_[1], y=depVar.predicted)) +
      geom_point(alpha=0.4,shape = 16,aes(), size=1.5) +  #geom_smooth(method=lm, formula = 'y ~ x', se=F,lty=2, lwd=1) +
      theme_bw()  +
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
    ggsave(paste(paths.plots,'2-ScatterPlot-',inputs.to_[1],'_vs_predicted-',depVar[i],'-',deepML,'.png',sep=''),
           plot=scatter.indices,width = 10, height = 10,  dpi = 300,units = "cm")
  }


}


data.to.predict<-data.frame(do.call(cbind,  save.field))
colnames(data.to.predict)<-paste(depVar,'.predicted',sep='')

df.plot<-cbind(data.field[,1:6],data.to.predict)


#fit logistic regression model
model.vcmax<- glm(Vcmax ~ Cab.predicted  + Car.predicted   + LIDFa.predicted  + LAI.predicted , data=df.plot)

df.plot[,'Vcmax.pred']<-predict(model.vcmax, df.plot)


r2.S<-round(cor(df.plot[,'Vcmax'],df.plot[,'Vcmax.pred'],use='pairwise.complete.obs')^2,2)
mylabel.r.se2 = bquote(bold(r)^2 == .(format(r2.S, digits = 3)))
statsLabel = paste0("Deep ML model; r2 = ", round(r2.S,2))

ggplot(df.plot, aes_string(x='Vcmax', y='Vcmax.pred')) +
  geom_point(alpha=0.4,shape = 16,aes(), size=1.5) +  #geom_smooth(method=lm, formula = 'y ~ x', se=F,lty=2, lwd=1) +
  theme_bw()  +
  scale_color_gradient(low = "#0091ff", high = "#f0650e") +
  theme(legend.position="bottom",
        plot.title = element_text(hjust = 0.5, size=10,face="bold"),
        axis.title = element_text(face="bold", size=10),
        axis.text.y=element_text(hjust = 0.5, size=10,face="bold"),
        axis.text.x=element_text(hjust = 0.5, size=10,face="bold"),
        legend.title=element_blank()) +
  labs(title = statsLabel, size=10,face="bold") +
  stat_smooth(method = "lm",formula = y ~ x,geom = "smooth",col='red',lty=2,se=T)



