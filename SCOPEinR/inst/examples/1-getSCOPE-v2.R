
#::::::::::::::::::::::::::::::::::::::::::::::::::::::::::::::::::::::::::::::::::::::::::::::::::::::::::::
### Running SCOPE (Soil Canopy Observation, Photochemistry and Energy balance model)
### SCOPE model v.2.1
### Main Author: Carlos Camino
#::::::::::::::::::::::::::::::::::::::::::::::::::::::::::::::::::::::::::::::::::::::::::::::::::::::::::::

rm(list= ls())

if (!require('parallel')) { install.packages('parallel'); require('parallel') }  ### Paralell
if (!require('doParallel')) { install.packages('doParallel'); require('doParallel') }  ### Paralell foreach and caret
if (!require('ToolsRTM')) { install.packages('ToolsRTM'); require('ToolsRTM') }  ###
if (!require('SCOPEinR')) { install.packages('SCOPEinR'); require('SCOPEinR') }  ###

if (!require('dplyr')) { install.packages('dplyr'); require('dplyr') }  ### db
if (!require('tidyr')) { install.packages('tidyr'); require('tidyr') }  ### db
if (!require('ggplot2')) { install.packages('ggplot2'); require('ggplot2') }  ### Plots



if (!require('signal')) { install.packages('signal'); require('signal') }  ### SCOPEinR and ToolsRTM

if (!require('pracma')) { install.packages('pracma'); require('pracma') }  ### SCOPEinR and ToolsRTM
if (!require('expint')) { install.packages('expint'); require('expint') }  ### SCOPEinR and ToolsRTM
if (!require('copula')) { install.packages('copula'); require('copula') }  ### SCOPEinR and ToolsRTM
if (!require('progress')) { install.packages('progress'); require('progress') }  ### SCOPEinR and ToolsRTM
if (!require('data.table')) { install.packages('data.table'); require('data.table') }  ### SCOPEinR and ToolsRTM

#::::::::::::::::::::::::::::::::::::::::::::::::::::::::::::::::::::::::::::::::::::::::::::::::::::::::::::
### 0.Get options for SCOPE model
#::::::::::::::::::::::::::::::::::::::::::::::::::::::::::::::::::::::::::::::::::::::::::::::::::::::::::::

#Get options for running SCOPE mode
table.with.opts<-read.table('input/setoptions.csv',header=T, sep=',')
paths.inputs <- 'Tables/FieldData/v2024/Processed/'

df.field = read.table(paste(paths.inputs,'4-GeckoRfl_average_with_2019_2020RawData_Vcmaxwith_FDL2.csv',sep=''), sep=',',dec = '.', header = T)
summary(df.field$LMA) /10000
summary(df.field$Vcmax25)



##::::::::::::::::::::::::::::::::::::::::::::::::::::::::::::::::::::::::::::::::::::::::::::::::::::::::::::
# 1.Get LUT ----
#::::::::::::::::::::::::::::::::::::::::::::::::::::::::::::::::::::::::::::::::::::::::::::::::::::::::::::

n.samples =20000
n_seed<-round(runif(1,1,n.samples),0)

inputLUT=read.table('input/inputs_SCOPE.csv',header=T, sep=',')
LUT <-getLUT.SCOPE(inputLUT=inputLUT,nLUT=n.samples)
# Parameters to be modified:

### biochemical compounds
LUT$Cab = stats::runif(n.samples,min = 5,max=70)
LUT$Car = stats::runif(n.samples,min = 1,max=15)
LUT$Anth = stats::runif(n.samples,min = 0,max=5)
LUT$Cbrown = stats::runif(n.samples,min = 0,max=1)
LUT$Cx = stats::runif(n.samples,min = 0,max=1)
LUT$LMA = stats::runif(n.samples,min = 0.003723 ,max=0.006104 )


### Photosyntetic capacity
LUT$Vcmax25 = stats::runif(n.samples,min = 5 ,max=90 )


LUT$EWT = 0.01
### Structural and viewing angles
LUT$tts = stats::runif(n.samples,min = 15,max=75)

LUT$tto = stats::runif(n.samples,min = 0,max=30)
LUT$hspot = stats::runif(n.samples,min = 0,max=1)
LUT$LAI = stats::runif(n.samples,min = 0.05 ,max=6)

LUT.lidfs<-getCor(n_inputs = 2,setseed = n_seed,distribution = 'Uniform',nLUT = n.samples, rho = 0.20,
                  Varnames = c('LIDFa','LIDFb'),MinRange = c(-0.5,-0.5), MaxRange = c(0.2,0.2))
LUT$LIDFa<-LUT.lidfs$LUT$LIDFa
LUT$LIDFb<-LUT.lidfs$LUT$LIDFb

dim(LUT)

plot(LUT$LIDFa, LUT$LIDFb)

##::::::::::::::::::::::::::::::::::::::::::::::::::::::::::::::::::::::::::::::::::::::::::::::::::::::::::::
# 2. Make simulations  in paralllel ----
#::::::::::::::::::::::::::::::::::::::::::::::::::::::::::::::::::::::::::::::::::::::::::::::::::::::::::::
n.rows <- nrow(LUT)
print(n.rows)
chunk_size <- n.rows /5

for (i in c(1:5)){
  #print(i)
  start.row <- ((i - 1) * chunk_size) +1
  #print(start.row)
  end.row <-  min(i * chunk_size, n.rows)
  #print(end.row)
  LUT_ <- LUT[start.row:end.row,]
  # Get random indices to subset your LUT data
  #n.samples =100
  #n.randoms <- sample(1:nrow(LUT), n.samples)

  db.sims <-SCOPEinR::get.SCOPE.parallel(LUT=LUT_, n.LUT=chunk_size,options.SCOPE=table.with.opts,optipar=SCOPEinR::optipar2017.ProspectD,
                                         leaf.model='fluspect-CX',canopy.model='fourSAIL', parallel = T,
                                         get.outputs='ALL', get.plots = F, get.csv =T)
}



##### Put the simulations in a
db.rfl <-getCSV(path.out='outs',n.folders=5,  files.names='All')

##::::::::::::::::::::::::::::::::::::::::::::::::::::::::::::::::::::::::::::::::::::::::::::::::::::::::::::
# 3. Make saddtional plots  ----
#::::::::::::::::::::::::::::::::::::::::::::::::::::::::::::::::::::::::::::::::::::::::::::::::::::::::::::

# Specify the folder path
folder_path <- "outs/"
# Get a list of all subdirectories
subdirectories <- list.dirs(folder_path, recursive = FALSE)
print(subdirectories)
last_subdirectory <- subdirectories[length(subdirectories)]
### optiops are 'fluorescence'; reflectance, radiance
traits <- c('Vcmax25','Anth')
get.SCOPE.plots(path.files=last_subdirectory, plant.trait=traits, get.plots='reflectance')
get.SCOPE.plots(path.files=last_subdirectory, plant.trait=traits, get.plots='fluorescence')
get.SCOPE.plots(path.files=last_subdirectory, plant.trait=traits, get.plots='radiance')






