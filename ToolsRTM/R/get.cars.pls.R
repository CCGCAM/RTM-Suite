#' Perform CARS-PLS analysis.
#'
#' This function performs Competitive Adaptive Reweighted Sampling (CARS) for partial least squares (PLS) analysis.
#'
#' @param X The sample matrix, where samples are in rows and variables are in columns.
#' @param y The response variable.
#' @param nLV The number of latent variables in PLS. Default is 2.
#' @param fold The number of segments for cross-validation. Default is 10.
#' @param scale.pretreat Whether to scale the variables. 1 for scaling, 0 for no scaling (only centered). Default is 1.
#' @param iteration The number of Monte Carlo samplings in CARS. Default is 50.
#' @param PartitionType The partition type for cross-validation: "random", "consecutive", or "interleaved". Default is "interleaved".
#' @author Ported by Carlos Camino; orignal code oin matlab by Yizeng Liang, and Hongdong Li
#' @return A list containing the results of the CARSPLS analysis.
#' @export
#'
#' @examples
#' \dontrun{
#' # X: sample matrix (rows = samples, columns = variables); y: response vector
#' get.cars.pls(X, y)
#' get.cars.pls(X, y, nLV = 3, fold = 5, scale.pretreat = 0)
#' }
get.cars.pls<-function(X,y,nLV=2,fold=10,scale.pretreat=1,iteration=50,PartitionType="interleaved") {


Num_row<-nrow(X)
Num_col<-ncol(X)

NewOrder<-order(y)
X<-X[NewOrder,]
y<-y[NewOrder]

data.CARS<-data.frame(indepdent=X,response=y)  #+++ Convert into dataframe
RMSECV<-rep(0,iteration)
NumLV<-rep(0,iteration)
VarIndex<-1:Num_col

#Some predifined data for CARS
subsetVariable<-1:Num_col     #elected variables for PLS modelling
Coef<-matrix(rep(0,Num_col*iteration),Num_col) #Coefficient matrix in CARS
Nvar<-rep(0,iteration)
ycal<-y

#Parameter of exponentially decreasing function.
ratio0<-1
ratio1<-2/Num_col
b=log(ratio0/ratio1)/(iteration-1)
a=ratio0*exp(b)
#Main Loop for CARS Algorithm
for (iter in 1:iteration){
     Xcal<-X[,subsetVariable]
     #+++ PLS modeling
     data.CARS.cal<-data.frame(indepdent=Xcal,response=ycal)  #Convert into dataframe
     nLV<-min(c(nLV,dim(Xcal)))
     ncomp=nLV
     ncomp
     if (scale.pretreat==1) {plsr.fit<- mvr(response~.,ncomp,data = data.CARS.cal,method = "simpls",scale=TRUE)}
     else{plsr.fit<- mvr(response~.,ncomp,data = data.CARS.cal,method = "simpls",scale=FALSE)}

     #Model selection by cross validation
     CV<-crossval(plsr.fit, segments =fold,data=data.CARS.cal,segment.type=PartitionType)

     validation<-CV$validation
     PRESS<-validation$PRESS
     RMSECV.temp=sqrt(PRESS/Num_row)
     RMSECV[iter]<-min(RMSECV.temp)
     NumLV[iter]<-which.min(RMSECV.temp)


     #xtract the coefficients and store them into the matrix: Coef.
     coef0<-rep(0,Num_col)
     coef.iter<-plsr.fit$coefficients[,,nLV]
     coef0[subsetVariable]<-coef.iter
     Coef[,iter]<-coef0

     #Weights of each variable
     weight<-abs(coef0)
     weight.order<-order(weight,decreasing=TRUE)

     #Calculate the ratio of variables to be retained by EDF in CARS.
     ratioVariable<-a*exp(-b*(iter+1))
     Nvar[iter]<-length(which(coef0!=0))
     K<-ceiling(Num_col*ratioVariable)

     #Eliminate the variables of small regression coefficients by force.
     weight[weight.order[K+1:Num_col]]<-0

     #Retained variables
     subsetVariable<-which(weight!=0)

     #screen print
     screen.output<-paste("The",iter,"th CARS-PLS iteration finished.")
     print(screen.output)
}

MinError<-min(RMSECV)
OPT.iter<-which(RMSECV==MinError)
OPT.iter<-OPT.iter[length(OPT.iter)]
min.RMSECV=min(RMSECV)
SelectedVariables<-which(Coef[,OPT.iter]!=0)
CARS<-list(Coef=Coef,Nvar=Nvar,RMSECV=RMSECV,
      NumLV=NumLV,Optimal.iteration=OPT.iter,MinError=MinError,
      SelectedVariables=SelectedVariables)
return(CARS)
}






