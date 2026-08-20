
#::::::::::::::::::::::::::::::::::::::::::::::::::::::::::::::::::::::::::::::::::::::::::::::::::::::::::::
### Running SCOPE (Soil Canopy Observation, Photochemistry and Energy balance model)
### SCOPE model v.2.1
### Main Author: Carlos Camino
#::::::::::::::::::::::::::::::::::::::::::::::::::::::::::::::::::::::::::::::::::::::::::::::::::::::::::::


rm(list= ls())

if (!require("parallel")) { install.packages("parallel"); require("parallel") }  ### Paralell
if (!require("doParallel")) { install.packages("doParallel"); require("doParallel") }  ### Paralell foreach and caret
if (!require("ToolsRTM")) { install.packages("ToolsRTM"); require("ToolsRTM") }  ### Paralell foreach and caret
if (!require("SCOPEinR")) { install.packages("SCOPEinR"); require("SCOPEinR") }  ### Paralell foreach and caret

if (!require("dplyr")) { install.packages("dplyr"); require("dplyr") }  ### Paralell
if (!require("tidyr")) { install.packages("tidyr"); require("tidyr") }  ### Paralell
if (!require("ggplot2")) { install.packages("ggplot2"); require(ggplot2) }  ### Paralell

#::::::::::::::::::::::::::::::::::::::::::::::::::::::::::::::::::::::::::::::::::::::::::::::::::::::::::::
### 0.Get options for SCOPE model
#::::::::::::::::::::::::::::::::::::::::::::::::::::::::::::::::::::::::::::::::::::::::::::::::::::::::::::

#Get options for running SCOPE mode
table.with.opts<-read.table('examples/Vcmax/data/setoptions.csv',header=T, sep=',')
paths.inputs <- 'examples/Vcmax/data/'

df.field = read.table(paste(paths.inputs,'4-GeckoRfl_average_with_2019_2020RawData_Vcmaxwith_FDL2.csv',sep=''), sep=',',dec = '.', header = T)
summary(df.field$LMA) /10000
summary(df.field$Vcmax25)


##::::::::::::::::::::::::::::::::::::::::::::::::::::::::::::::::::::::::::::::::::::::::::::::::::::::::::::
# 1.Get LUT ----
#::::::::::::::::::::::::::::::::::::::::::::::::::::::::::::::::::::::::::::::::::::::::::::::::::::::::::::


n.samples =20000
n_seed<-round(runif(1,1,n.samples),0)

inputLUT=read.table('examples/Vcmax/data/inputs_SCOPE.csv',header=T, sep=',')
LUT <-getLUT.SCOPE(inputLUT=inputLUT,nLUT=n.samples)

#Main Parameter to estimate GPP
# LUT$Rdparam == 0


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

LUT$alpha <-5
LUT$EWT = 0.01
### Structural and viewing angles
LUT$tts = stats::runif(n.samples,min = 15,max=75)

LUT$tto = stats::runif(n.samples,min = 0,max=30)
LUT$hspot = stats::runif(n.samples,min = 0,max=1)
LUT$LAI = stats::runif(n.samples,min = 0.05 ,max=6)

LUT.lidfs<-getCor(n_inputs = 2,setseed = n_seed,distribution = 'Uniform',nLUT = n.samples, rho = 0.5,
                  Varnames = c('LIDFa','LIDFb'),MinRange = c(-0.5,-0.5), MaxRange = c(0.5,0.5))
LUT$LIDFa<-LUT.lidfs$LUT$LIDFa
LUT$LIDFb<-LUT.lidfs$LUT$LIDFb

plot(LUT$LIDFb,LUT$LIDFa)
dim(LUT)

##::::::::::::::::::::::::::::::::::::::::::::::::::::::::::::::::::::::::::::::::::::::::::::::::::::::::::::
# 3 .Get Simulations  in parallell ----
#::::::::::::::::::::::::::::::::::::::::::::::::::::::::::::::::::::::::::::::::::::::::::::::::::::::::::::
set.seed(3425)

n.LUT<-10
rows_to_select <- sample(nrow(LUT), n.LUT)
print(rows_to_select)
LUT_ <- LUT[rows_to_select,]




db.sims <-get.SCOPE.parallel(LUT=LUT_,options.SCOPE=table.with.opts,optipar=SCOPEinR::optipar2017.ProspectD,
                             leaf.model='fluspect-CX',canopy.model='fourSAIL', parallel = T,
                             get.outputs='Main', get.plots = F, get.csv =T, n.cores = 8)

matrix.reflapp <- sapply(seq_along(db.sims), function(i) db.sims[[i]]$data.rad$reflapp[1:2001])
# Assuming matrix.refl is your matrix
matplot(wave_,matrix.reflapp, type = "l", lty = 1, col = 1:ncol(matrix.reflapp), xlab = "", ylab = "Reflectance", main = "SCOPE v 2.1")
legend("topright", legend = 1:ncol(matrix.reflapp), col = 1:ncol(matrix.reflapp), lty = 1, title = "")


get.csv(path.out = 'out',n.folder=4, files.names = 'All')

# Specify the folder path
folder_path <- "outs/"
# Get a list of all subdirectories
subdirectories <- list.dirs(folder_path, recursive = FALSE)
# Get the last subdirectory
last_subdirectory <- tail(subdirectories, n = 1)

### optiops are 'fluorescence'; reflectance, radiance
traits <- c('Vcmax25','Car','Anth')
get.SCOPE.plots(path.files=last_subdirectory, plant.trait=traits, get.plots='reflectance')
get.SCOPE.plots(path.files=last_subdirectory, plant.trait=traits, get.plots='fluorescence')

get.SCOPE.plots(path.files=last_subdirectory, plant.trait=traits, get.plots='radiance')

matrix.reflapp <- sapply(seq_along(db.sims), function(i) db.sims[[i]]$data.rad$reflapp[1:2001])
# Assuming matrix.refl is your matrix
matplot(wave_,matrix.reflapp, type = "l", lty = 1, col = 1:ncol(matrix.refl), xlab = "", ylab = "Reflectance", main = "SCOPE v 2.1")
legend("topright", legend = 1:ncol(matrix.refl), col = 1:ncol(matrix.refl), lty = 1, title = "")





##::::::::::::::::::::::::::::::::::::::::::::::::::::::::::::::::::::::::::::::::::::::::::::::::::::::::::::
# 2. Make simulations  ----
#::::::::::::::::::::::::::::::::::::::::::::::::::::::::::::::::::::::::::::::::::::::::::::::::::::::::::::


start_time <- Sys.time()
set.seed(1236)
rows_to_select <- sample(nrow(LUT), n.LUT)
LUT_ <- LUT[rows_to_select,]

pb <- progress::progress_bar$new( format = "[:bar] :percent ETA: :eta", total = 10)

db.sims <- list()
for (i in c(1:10)){

  db.sims[[i]] <-get.SCOPE.ind(LUT=LUT_[i,],options.SCOPE=table.with.opts,
                       optipar=SCOPEinR::optipar2021.Pro.CX,
                       leaf.model='fluspect-CX',canopy.model='fourSAIL',
                       get.outputs = 'ALL', get.plots = F)
  # Update progress bar
  #pb$tick()
}

wave_<-db.sims[[1]]$data.spectral$wlS[1:2001]
matrix.reflapp <- sapply(seq_along(db.sims), function(i) db.sims[[i]]$data.rad$reflapp[1:2001])
# Assuming matrix.refl is your matrix
matplot(wave_,matrix.reflapp, type = "l", lty = 1, col = 1:ncol(matrix.reflapp), xlab = "", ylab = "Reflectance", main = "SCOPE v 2.1")
legend("topright", legend = 1:ncol(matrix.refl), col = 1:ncol(matrix.refl), lty = 1, title = "")



## get main inputs
get.SCOPE.outputs(data.sim = db.sims, N.sims=10,LUT=LUT[1:10,], path.out ='outs/',
                  get.outputs='ALL',
                  #get.more.inputs=c('refl','lidf','LIDFb','Ft_Fo','rdo'),
                  get.plots=T)


##::::::::::::::::::::::::::::::::::::::::::::::::::::::::::::::::::::::::::::::::::::::::::::::::::::::::::::
# 2. Make simulations  ----
#::::::::::::::::::::::::::::::::::::::::::::::::::::::::::::::::::::::::::::::::::::::::::::::::::::::::::::


db.sims<- SCOPEinR::get.SCOPE.ind(LUT=LUT[1:10,],options.SCOPE=table.with.opts,
                                 optipar=SCOPEinR::optipar2021.Pro.CX,
                                 leaf.model='fluspect-CX',canopy.model='fourSAIL',
                                 get.outputs = 'ALL', get.plots = F)


db.sims[[1]]


##::::::::::::::::::::::::::::::::::::::::::::::::::::::::::::::::::::::::::::::::::::::::::::::::::::::::::::
# 2. Make simulations  ----
#::::::::::::::::::::::::::::::::::::::::::::::::::::::::::::::::::::::::::::::::::::::::::::::::::::::::::::


start_time <- Sys.time()
sims <- get.SCOPE(LUT, n.LUT=10,options.SCOPE=table.with.opts,optipar=SCOPEinR::optipar2017.ProspectD,
                 leaf.model='fluspect-CX',canopy.model='fourSAIL',
                 get.outputs='Main', get.plots = F)
end_time <- Sys.time()
print(end_time - start_time)


##::::::::::::::::::::::::::::::::::::::::::::::::::::::::::::::::::::::::::::::::::::::::::::::::::::::::::::
# 4. Make simulations with chunks  ----
#::::::::::::::::::::::::::::::::::::::::::::::::::::::::::::::::::::::::::::::::::::::::::::::::::::::::::::


# Define function to process each chunk of LUT data
process_chunk <- function(chunk) {
  SCOPEinR::get.SCOPE(chunk, n.LUT = nrow(chunk), options.SCOPE = table.with.opts,
            optipar = SCOPEinR::optipar2017.ProspectD,
            leaf.model = 'fluspect-CX', canopy.model = 'fourSAIL',
            get.outputs = 'Main', get.plots = FALSE)
}
# Split LUT into chunks
num_cores <- detectCores()-1 # Get the number of available cores

rows.r <- sample(nrow(LUT), 1000)

chunks <- split(LUT[rows.r,], sort(rep_len(1:num_cores, length.out = nrow(LUT[rows.r,]))))
# Initialize progress bar
pb <- progress_bar$new(total = length(chunks), format = "[:bar] :percent eta: :eta")

# Initialize parallel processing
cl <- makeCluster(num_cores)  # Create a cluster using all available cores
registerDoParallel(cl)

start_time <- Sys.time()
# Parallel execution
sims <- foreach(chunk = chunks, .combine = c) %dopar% {
  pb$tick()  # Increment progress bar
  process_chunk(chunk)
}
end_time <- Sys.time()
print(end_time - start_time)
# Clean up parallel resources
stopCluster(cl)







