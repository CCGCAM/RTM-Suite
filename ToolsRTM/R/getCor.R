
#' Generate multiple mutually-correlated variables at a target correlation level
#'
#' Builds \code{n_inputs} variables, each \code{nLUT} values long, all
#' pairwise-correlated at (approximately) \code{rho}, then rescales each to
#' its own \code{[MinRange, MaxRange]} bound. Useful for LUTs where several
#' traits should co-vary (e.g. LAI and canopy height) rather than being
#' drawn independently.
#'
#' @param n_inputs integer. How many variables to generate.
#' @param nLUT integer. How many values per variable (LUT rows). Default 100.
#' @param distribution character. \code{"Uniform"} or \code{"Normal"}. Default \code{"Uniform"}.
#' @param setseed integer. Random seed. Default 123.
#' @param rho numeric (-1 to 1). Target pairwise correlation between every pair of variables.
#' @param Varnames character vector, length \code{n_inputs}. Column names for the output;
#'   if \code{NULL}, columns are named \code{Var_1}, \code{Var_2}, ...
#' @param MinRange numeric vector, length \code{n_inputs}. Minimum of each variable's output range.
#' @param MaxRange numeric vector, length \code{n_inputs}. Maximum of each variable's output range.
#' @return A list: \code{LUT} (data.frame, \code{nLUT} rows x \code{n_inputs} columns) and
#'   \code{Covarianza} (the realized correlation matrix of \code{LUT}, for checking how close
#'   the sample came to the requested \code{rho}).
#' @export
#'
#' @examples
#' out <- getCor(n_inputs = 2, nLUT = 200, distribution = "Uniform", rho = 0.7,
#'                Varnames = c("LAI", "Height"), MinRange = c(0.5, 2), MaxRange = c(7, 30))
#' cor(out$LUT$LAI, out$LUT$Height)  # close to 0.7
getCor<-function(n_inputs=NULL, nLUT=100,  distribution = 'Uniform',setseed = 123, rho=NULL,
                         Varnames = NULL,
                         MinRange = NULL,
                         MaxRange = NULL){
  
  if (is.null(rho)){
    message('Please insert r value, only valid for Uniform distribution')
    stop()
  } 
  
  if (is.null(n_inputs)){
    message('Please insert number of variables that you need correlate')
    stop()
  } 
  
  if (is.null(distribution)){
    message('Please indicate the distribution function for all correlated inputs ...')
  } 
  if (is.null(n_inputs)){
    message('Please insert length of the variable')
    stop()
  }
  if (is.null(setseed)){
    set.seed(1234)
  }
  
  if (length(MinRange) != length(MinRange) ){
    message('Min and Max vectors should have same lengths')
    stop()
  }
  if (is.null(MinRange) | is.null(MinRange)){
    message('please some inputs is missing')
    stop()
  }
  #############################################################
  ### same variables
  ### list of functions
  list_distributions = list(function(L)rnorm(L,0,2), #normal
                            function(L)runif(L,0,1), #uniform
                            function(L) qpois(L,7), # poisson
                            function(L) round(qnorm(L,100,10)), ## normal
                            function(L) qnorm(L,-100,1)) ## normal

  n=n_inputs
  n.samples = nLUT
  #############################################################
  
  if ( distribution == 'Normal'){
    message('Generating a Normal distribution for all correlated inputs ...')
    
    # Was: Sigma <- t(A) %*% A for a fully random A -- ignored the `rho`
    # argument entirely, so the requested correlation had no effect for
    # distribution = 'Normal' (only the Uniform branch below used rho).
    # Build Sigma with off-diagonal = rho (like the Uniform branch does)
    # instead, so both branches honor the same `rho` argument.
    if (is.null(rho)) rho <- 0
    Sigma <- matrix(rho, nrow = n, ncol = n)
    diag(Sigma) <- 1
    mu <- rep(0, n)

    df.matrix <- MASS::mvrnorm(n.samples, mu = mu, Sigma = Sigma )  # from Mass package
   
    for (i in c(1:dim(df.matrix)[2])){
      
      df.matrix[,i]<-scales::rescale(df.matrix[,i], to = c(MinRange[i], MaxRange[i]))   
      df.matrix[,i]<-(jitter(df.matrix[,i], factor=2, amount = NULL))
       
    }
    
    #summary(df.matrix)
    # Calculate kernel density estimate
    #df.matrix.kde <- MASS::kde2d(df.matrix[,1], df.matrix[,n], n = n.samples/0.4)   # from MASS package
    
    # Contour plot overlayed on heat map image of results
    #print(graphics::image(df.matrix.kde))       # from base graphics package
    #print(contour(df.matrix.kde, add = TRUE)  )
    
    df.export<-data.frame(df.matrix)
    if ( is.null(Varnames)){
      colnames(df.export)<-paste('Var_',c(1:n),sep = '')
    } else {
      colnames(df.export)<-Varnames
    }
    M <-cor(df.matrix)
    #print(M)
    

  } 
  
  if ( distribution == 'Uniform'){
    message('Generating a Uniform distribution for all correlated inputs ...')
    #norm.cop <- copula::normalCopula(rho,dim=n); 
    #df.uniform <- copula::rCopula(n.samples, norm.cop)
    
    # Generating uncorrelated uniform random variables
    df.uniform <- matrix(runif(n.samples * n), ncol = n)
    # If rho is provided, make them correlated
    if (!is.null(rho) && rho != 0) {
      Sigma <- matrix(rho, nrow = n, ncol = n)
      diag(Sigma) <- 1
      df.uniform <- t(t(chol(Sigma)) %*% t(df.uniform))
    }
    

    ## convert to uniform
    #df.uniform = pnorm(df.uniform) 
    #df.uniform = sapply(1:n, FUN = function(i) list_distributions[[2]](df.uniform[,i]))
    #pCopula(as.matrix(df.uniform),norm.cop)
    
    for (i in c(1:dim(df.uniform)[2])){
      
      df.uniform[,i]<-scales::rescale(df.uniform[,i], to = c(MinRange[i], MaxRange[i]))   
      df.uniform[,i]<-(jitter(df.uniform[,i], factor=2, amount = NULL))
      
    }
    
    df.export<-as.data.frame(df.uniform)
    if ( is.null(Varnames)){
      colnames(df.export)<-paste('Var_',c(1:n),sep = '')
    } else {
      colnames(df.export)<-Varnames
    }
    
    # Calculate kernel density estimate
    #mvn.kde <- MASS::kde2d(df.uniform[,1], df.uniform[,n], n = n.samples/0.4)  
    #print(graphics::image(mvn.kde))       # from base graphics package
    #print(contour(mvn.kde, add = TRUE))
    
    M <-cor(df.export)
    #print(M)
  }
  
  LUTdata = list('LUT'=df.export,'Covarianza'=M)
  return(LUTdata)
  
}
