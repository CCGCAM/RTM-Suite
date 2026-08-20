

###########################################################################################################
### Running SCOPE (Soil Canopy Observation, Photochemistry and Energy balance model)
### SCOPE model v.2.1
### Main Author: Carlos Camino
###########################################################################################################

rm(list= ls())

if (!require("parallel")) { install.packages("parallel"); require("parallel") }  ### Paralell
if (!require("doParallel")) { install.packages("doParallel"); require("doParallel") }  ### Paralell foreach and caret
if (!require("ToolsRTM")) { install.packages("ToolsRTM"); require("ToolsRTM") }  ### Paralell foreach and caret
if (!require("SCOPEinR")) { install.packages("SCOPEinR"); require("SCOPEinR") }  ### Paralell foreach and caret

if (!require("dplyr")) { install.packages("dplyr"); require("dplyr") }  ### Paralell
if (!require("tidyr")) { install.packages("tidyr"); require("tidyr") }  ### Paralell
if (!require("ggplot2")) { install.packages("ggplot2"); require(ggplot2) }  ### Paralell

##################################################################################
### 1.Get options for SCOPE model
##################################################################################

#Get options for running SCOPE mode
table.with.opts<-read.table('inst/input/setoptions.csv',header=T, sep=',')


##################################################################################
### 2.Get LUT with Inputs needed for SCOPE model
##################################################################################

N.Samples =1
file.LUT = c('LUT','default')
file.LUT = file.LUT[2]

start_time <- Sys.time()
sim.scope<-list()
if (file.LUT == 'default'){
  Table.LUT=read.table('inst/input/LUT_input.csv',header=T,sep=',')

  db.sim <-SCOPEinR::get.SCOPE(LUT=Table.LUT[1,],options.SCOPE=table.with.opts,
                                   optipar=SCOPEinR::optipar2021.Pro.CX,
                                   leaf.model='fluspect-CX',canopy.model='fourSAIL',
                                   get.outputs = 'ALL', get.plots = F)


} else{

  inputLUT=read.table('input/inputs_SCOPE.csv',header=T, sep=',')
  Table.LUT <-getLUT.SCOPE(inputLUT=inputLUT,nLUT=N.Samples)

  db.sim <-get.SCOPE(LUT=Table.LUT, n.LUT = N.Samples,options.SCOPE=table.with.opts,
                     optipar=SCOPEinR::optipar2021.Pro.CX,
                     leaf.model='fluspect-CX',canopy.model='fourSAIL',
                     get.outputs = 'ALL', get.plots = F)




}
execution_time <-end_time <- Sys.time()
print(end_time - start_time)
print(paste("Execution time:", execution_time))

## get main inputs
get.SCOPE.outputs(data.sim = db.sim, N.sims=N.Samples,LUT=Table.LUT, path.out ='outs/',
                  get.more.inputs=c('refl','lidf','LIDFb','Ft_Fo','rdo'),
                  get.plots=T)

##############################################################################################################################
# 2.Get Simulations  in parallell ----
##############################################################################################################################

## choose number of processors/cores
no_cores <- detectCores() - 2
cl <- makeCluster(no_cores)
registerDoParallel(cl)

start_time <- Sys.time()


N.Samples =10

inputLUT=read.table('input/inputs_SCOPE.csv',header=T, sep=',')
Main.LUT <-getLUT.SCOPE(inputLUT=inputLUT,nLUT=N.Samples)

n_seed<-round(runif(1,1,N.Samples),0)

LUT.lidfs<-getCor(n_inputs = 2,setseed = n_seed,distribution = 'Uniform',nLUT = N.Samples, rho = 0.20,
                Varnames = c('LIDFa','LIDFb'),MinRage = c(-0.5,-0.5), MaxRange = c(0.2,0.2))
Main.LUT$LIDFa<-LUT.lidfs$LUT$LIDFa
Main.LUT$LIDFb<-LUT.lidfs$LUT$LIDFb

dim(Main.LUT)
table.with.opts<-read.table('input/setoptions.csv',header=T, sep=',')
sim.scope<-list()
db.sims.paralell<-foreach::foreach(i=1:N.Samples, .packages = c("SCOPEinR"), # .combine = list not working
                                   .export = c('Main.LUT','table.with.opts')) %dopar% {

  db.sim<- SCOPEinR::get.SCOPE.ind(LUT=Main.LUT[i,],options.SCOPE=table.with.opts,
                              optipar=SCOPEinR::optipar2021.Pro.CX,
                              leaf.model='fluspect-CX',canopy.model='fourSAIL',
                              get.outputs = 'Main', get.plots = F)
  sim.scope[[i]]<-db.sim



} ##end paralle

stopCluster(cl)
execution_time <-end_time <- Sys.time()
print(end_time - start_time)
print(paste("Execution time:", execution_time))

length(db.sims.paralell)


get.SCOPE.outputs(data.sim = db.sims.paralell, N.sims=N.Samples,LUT=Main.LUT,
                             path.out ='outs/',get.plots=T)


# Specify the folder path
folder_path <- "outs/"
# Get a list of all subdirectories
subdirectories <- list.dirs(folder_path, recursive = FALSE)
print(subdirectories)

### optiops are 'fluorescence'; reflectance, radiance
traits <- c('Vcmax25','EWT','Anth')
get.SCOPE.plots(path.files=subdirectories[8], plant.trait=traits, get.plots='reflectance')
get.SCOPE.plots(path.files=subdirectories[6], plant.trait=traits, get.plots='fluorescence')

get.SCOPE.plots(path.files=subdirectories[6], plant.trait=traits, get.plots='radiance')






