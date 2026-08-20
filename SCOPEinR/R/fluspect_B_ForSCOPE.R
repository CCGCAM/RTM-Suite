#' Leaf FLUSPECT-B model for SCOPE
#'
#' \code{getFluspect.B.SCOPE} calculates reflectance and transmittance spectra of a leaf using FLUSPECT-B,
#' calculates reflectance and transmittance spectra of a leaf using FLUSPECT-B,
#'plus four excitation-fluorescence matrices
#' @param inputsLeaf  a LUT for FLUXspect
#' @param inputsOptipar  internal parameters
#' @param version  Leaf model uses for Kcal: 'D' refers as PROSPECT-D; 'PRO', refers as PROSPECT-PRO
#'
#' @return  a list which contains:
#' * refl (reflectance)
#' * tran (transmittance)
#' * Mb (backward scattering fluorescence matrix, I for PSI and II for PSII)
#' * Mf (forward scattering fluorescence matrix,  I for PSI and II for PSII)
#' @export

#' @author Wout Verhoef, Christiaan van der Tol, Joris Timmermans, Nastassia Vilfan (Original version in Matlab)
#' @author Carlos Camino (Ported version into R)
#'
# Authors: Wout Verhoef, Christiaan van der Tol, Joris Timmermans,
# Date: 2007
# Update from PROSPECT to FLUSPECT: January 2011 (CvdT)
#
#     Nov 2012 (CvdT) Output EF-matrices separately for PSI and PSII
#   31 Jan 2013 (WV)   Adapt to SCOPE v_1.40, using structures for I/O
#   30 May 2013 (WV)   Repair bug in s for non-conservative scattering
#   24 Nov 2013 (WV)   Simplified doubling routine
#   25 Nov 2013 (WV)   Restored piece of code that takes final refl and
#                      tran outputs as a basis for the doubling routine
#   03 Dec 2013 (WV)   Major upgrade. Border interfaces are removed before
#                      the fluorescence calculation and later added again
#   23 Dec 2013 (WV)   Correct a problem with N = 1 when calculating k
#                      and s; a test on a = Inf was included
#   01 Apr 2014 (WV)   Add carotenoid concentration (Cca and Kca)
#   19 Jan 2015 (WV)   First beta version for simulation of PRI effect
#   20 Jan 2021 (CvdT) Include PROSPECT-PRO coefficients
#   01 Apr 2014 (WV)   Add carotenoid concentration (Cca and Kca)
#   19 Jan 2015 (WV)   First beta version for simulation of PRI effect
#'
#' @examples
#' \dontrun{
#' inputs <- ToolsRTM::inputsRTMs
#' LUT <- as.data.frame(ToolsRTM::getLUT(inputs = inputs, nLUT = 1, setseed = 1234))
#' sim <- getFluspect.B.SCOPE(inputsLeaf = LUT, inputsOptipar = SCOPEinR::optipar, version = 'D')
#' }

getFluspect.B.SCOPE<-function(inputsLeaf,inputsOptipar, version = 'D' )  {

  if (is.null(version)){
    version = 'D'
  }
  if ("Prot" %in% colnames(inputsLeaf)){
    version = 'Cx'
  }  else {
    version = 'D'
  }

  # fixed parameters for the fluorescence module
  ndub = 15;           # number of doublings applied
  # Fluspect parameters
  #inputsLeaf=LUT.flus[1,]

  #define alll inputs in the models. retreived from LUT tables
  Cab = inputsLeaf[,'Cab'];Car = inputsLeaf[,'Car'];EWT = inputsLeaf[,'EWT']; LMA = inputsLeaf[,'LMA'];
  Cs = inputsLeaf[,'Cs']; N = inputsLeaf[,'N']; fqe_ = inputsLeaf[,'fqe'];
  fqe <- c(fqe_/5, fqe_)
  Cx = inputsLeaf[,'Cx']  # Violaxanthin - Zeaxanthin transition status [0-1]
  if (version == 'Cx') {
    Prot = inputsLeaf[,'Prot']; CBC = inputsLeaf[,'CBC'];
    Anth =inputsLeaf[,'Anth']

    if (LMA > 0 & (Prot > 0 | CBC > 0)) {
      LMA <- 0
    }
  }


  # If inputsOptipar is not present, print a message
  if (missing(inputsOptipar)) {

    ### define inputsOptipar
    optipar_default <- SCOPEinR::optipar2017.ProspectD

    nr = optipar_default[['nr']]; Kdm = optipar_default[['Kdm']];Kab = optipar_default[['Kab']];

    Kca = optipar_default[['Kca']]; Kw = optipar_default[['Kw']]; Ks = optipar_default[['Ks']];
    Kant = optipar_default[['Kant']];
    Kcbc = optipar_default[['Kcbc']]; Kp = optipar_default[['Kp']];
    phiI = optipar_default[['phiI']]; phiII = optipar_default[['phiII']]
    phi = optipar_default[['phi']]
    if (Cx == -999) {
      # Use old Kca spectrum if this is given as input
      Kca <-optipar_default[['Kca']]
    } else {
      # Otherwise make linear combination based on Cx
      # For Cx going from 0 to 1 we go from Viola to Zea
      Kca <- (1 - Cx) *optipar_default[['KcaV']] + Cx * optipar_default[['KcaZ']]
    }

  } else {   # If inputsOptipar is present, continue with the function
    nr = inputsOptipar[['nr']]; Kdm = inputsOptipar[['Kdm']];Kab = inputsOptipar[['Kab']];
    Kca = inputsOptipar[['Kca']];Kw = inputsOptipar[['Kw']]; Ks = inputsOptipar[['Ks']];
    Kant = inputsOptipar[['Kant']];
    Kcbc = inputsOptipar[['Kcbc']]; Kp = inputsOptipar[['Kp']];
    phiI = inputsOptipar[['phiI']]; phiII = inputsOptipar[['phiII']];
    phi = inputsOptipar[['phi']];
    if (Cx == -999) {
      # Use old Kca spectrum if this is given as input
      Kca <-inputsOptipar[['Kca']]
    } else {
      # Otherwise make linear combination based on Cx
      # For Cx going from 0 to 1 we go from Viola to Zea
      Kca <- (1 - Cx) *inputsOptipar[['KcaV']] + Cx * inputsOptipar[['KcaZ']]
    }


  }





  if (version == 'D') {
    #### PROSPECT D calculations
    Kall = (Cab * Kab + Car * Kca + LMA * Kdm + EWT * Kw  + Cs * Ks) / N  # Compact leaf layer


  } else {
    #### PROSPECT  PRO calculations
    Kall = (Cab * Kab + Car * Kca + LMA * Kdm + EWT * Kw  + Cs * Ks + Anth * Kant + Prot * Kp + CBC * Kcbc) / N  # Compact leaf layer


  }


  j <- which(Kall > 0)  # Non-conservative scattering (normal case)
  t1 <- (1 - Kall) * exp(-Kall)
  t2 <- Kall^2 * expint::expint(Kall)
  tau <- rep(1, length(t1))
  tau[j] <- t1[j] + t2[j]
  kChlrel <- rep(0, length(t1))
  kChlrel[j] <- Cab * Kab[j] / (Kall[j] * N)

  talf <- ToolsRTM::calctav(59, nr)
  ralf <- 1 - talf

  t12 <- ToolsRTM::calctav(90, nr)
  r12 <- 1 - t12

  t21 <- t12 / (nr ^ 2)
  r21 <- 1 - t21

  # top surface side
  denom <- 1 - r21 * r21 * tau * tau
  Ta <- talf * tau * t21 / denom
  Ra <- ralf + r21 * tau * Ta

  # bottom surface side
  t <- t12 * tau * t21 / denom
  r <- r12 + r21 * tau * t

  # Stokes equations to compute properties of next N-1 layers (N real)
  # Normal case

  D <- sqrt((1 + r + t) * (1 + r - t) * (1 - r + t) * (1 - r - t))
  rq <- r^2
  tq <- t^2
  a <- (1 + rq - tq + D) / (2 * r)
  b <- (1 - rq + tq + D) / (2 * t)

  bNm1 <- b^(N - 1)
  bN2 <- bNm1^2
  a2 <- a^2
  denom <- a2 * bN2 - 1
  Rsub <- a * (bN2 - 1) / denom
  Tsub <- bNm1 * (a2 - 1) / denom

  # Case of zero absorption
  j <- which(r + t >= 1)
  Tsub[j] <- t[j] / (t[j] + (1 - t[j]) * (N - 1))
  Rsub[j] <- 1 - Tsub[j]

  # Reflectance and transmittance of the leaf: combine top layer with next N-1 layers
  denom <- 1 - Rsub * r
  tran <- Ta * Tsub / denom
  refl <- Ra + Ta * Rsub * t / denom

  # From here a new path is taken: The doubling method used to calculate
  #  fluoresence is now only applied to the part of the leaf where absorption
  # takes place, that is, the part exclusive of the leaf-air interfaces. The
  #  reflectance (rho) and transmittance (tau) of this part of the leaf are
  # now determined by "subtracting" the interfaces

  Rb  = (refl - ralf) / (talf * t21 + (refl - ralf) * r21)  # Remove the top interface
  Z   = tran * (1 - Rb * r21) / (talf * t21)              # Derive Z from the transmittance

  rho = (Rb - r21 * Z^2)/(1-(r21 * Z)^2)   # Reflectance and transmittance
  tau = (1 - Rb * r21) / ( 1 - (r21 * Z)^2) * Z   # of the leaf mesophyll layer
  t   = tau
  r <- pmax(rho, 0)  # Avoid negative r

  # Derive Kubelka-Munk s and k
  I_rt <- (r + t) < 1
  D <- sqrt((1 + r[I_rt] + t[I_rt]) *
              (1 + r[I_rt] - t[I_rt]) *
              (1 - r[I_rt] + t[I_rt]) *
              (1 - r[I_rt] - t[I_rt]))
  a <- (1 + r[I_rt]^2 - t[I_rt]^2 + D) / (2 * r[I_rt])
  b <- (1 - r[I_rt]^2 + t[I_rt]^2 + D) / (2 * t[I_rt])
  a[!I_rt] <- 1
  b[!I_rt] <- 1

  s <- r / t
  I_a <- (a > 1 & a != Inf)

  I_na <- is.na(s) | is.na(a) | is.na(b)
  s[I_a & !I_na] <- 2 * a[I_a & !I_na] / (a[I_a & !I_na]^2 - 1) * log(b[I_a & !I_na])

  #s[I_a] <- 2 * a[I_a] / (a[I_a]^2 - 1) * log(b[I_a])

  k <- log(b)
  I_na <- is.na(k) | is.na(a) | is.na(b)
  k[I_a & !I_na] <- (a[I_a & !I_na] - 1) / (a[I_a & !I_na] + 1) * log(b[I_a & !I_na])

  #k[I_a] <- (a[I_a] - 1) / (a[I_a] + 1) * log(b[I_a])
  kChl <- kChlrel * k

  #  Fluorescence of the leaf mesophyll layer
  if (fqe_ > 0) {

    ## get bands from an object
    spectral <- SCOPEinR::define.bands()

    wle <- spectral[['wlE']]  # excitation wavelengths, transpose to column
    #wle <- matrix(wle, nrow = length(wle), ncol = 1)
    wlf <- spectral[['wlF']]  # fluorescence wavelengths, transpose to column
    #wlf <- matrix(wlf, nrow = length(wlf), ncol = 1)
    wlp <- c(spectral[['wlP']])     # PROSPECT wavelengths, kept as a row vector

    minwle <- min(wle)
    maxwle <- max(wle)
    minwlf <- min(wlf)
    maxwlf <- max(wlf)

    # indices of wle and wlf within wlp over PROSPECT model (2100 wavelengths)
    Iwle <- which(wlp >= minwle & wlp <= maxwle)
    Iwlf <- which(wlp >= minwlf & wlp <= maxwlf)
    eps <- 2^(-ndub)

    # initialisations
    te <- 1 - (k[Iwle] + s[Iwle]) * eps
    tf <- 1 - (k[Iwlf] + s[Iwlf]) * eps
    ### for range between 400-750
    re <- s[Iwle] * eps
    ## for range 640 and 840
    rf <- s[Iwlf] * eps


    # Calculate sigmoid matrix
   sigmoid <- 1/(1+ outer(exp(-wlf/10), exp((wle)/10))) #

   # sigmoid <- matrix(0, nrow = length(wlf), ncol = length(wle))  # initialize sigmoid matrix with zeros
   #
   #  for (i in 1:length(wlf)) {
   #    for (j in 1:length(wle)) {
   #      sigmoid[i, j] <- 1 / (1 + exp(-wlf[i] / 10) * exp(wle[j] / 10))  # compute sigmoid element-wise
   #    }
   #  }
    dim(sigmoid)


    MfI  <- MbI  <- fqe[1] * ((0.5 * phiI[Iwlf]) * eps) %*% t(as.matrix(kChl[Iwle])) * sigmoid
    MfII <- MbII <- fqe[2] * ((0.5 * phiII[Iwlf]) * eps) %*% t(as.matrix(kChl[Iwle])) * sigmoid


    Ih <-  matrix(1, nrow = 1, ncol = length(te))#rep(1, length(te))  # row of ones
    Iv <- matrix(1, nrow = length(tf), ncol = 1) #rep(1, length(tf))  # column of ones

    # Doubling routine
    for (i in 1:ndub) {
      xe <- te / (1 - re * re)
      ten <- te * xe
      ren <- re * (1 + ten)

      xf <- tf / (1 - rf * rf)
      tfn <- tf * xf
      rfn <- rf * (1 + tfn)


      A11  <- xf %*% Ih + Iv %*% xe;
      A12 <- (xf %*% t(xe)) * (rf %*% Ih + Iv %*% re)
      A21  <- 1+(xf %*% t(xe)) * (1 + rf %*% t(re))
      A22 <- (xf * rf) %*% Ih + Iv %*% (xe * re)



      MfnI  <- MfI  * A11 + MbI  * A12
      MbnI  <- MbI  * A21 + MfI  * A22
      MfnII <- MfII * A11 + MbII * A12
      MbnII <- MbII * A21 + MfII * A22

      te <- ten
      re <- ren
      tf <- tfn
      rf <- rfn

      MfI <- MfnI
      MbI <- MbnI
      MfII <- MfnII
      MbII <- MbnII

    } ## end ndum

    # Here we add the leaf-air interfaces again for obtaining the final
    # leaf level fluorescences.

    g1 <- MbI
    g2 <- MbII
    f1 <- MfI
    f2 <- MfII

    Rb <- rho + tau^2 * r21 / (1 - rho * r21)

    Xe <- Iv %*% (talf[Iwle] / (1- r21[Iwle] * Rb[Iwle]))
    Xf <- t21[Iwlf] / (1 - r21[Iwlf] * Rb[Iwlf]) %*% Ih
    Ye <- Iv %*% (tau[Iwle] * r21[Iwle] / (1 - rho[Iwle] * r21[Iwle]))
    Yf <- tau[Iwlf] * r21[Iwlf] / (1 - rho[Iwlf] * r21[Iwlf]) %*% Ih

    A <- Xe * (1 + Ye * Yf) * Xf
    B <- Xe * (Ye + Yf) * Xf

    g1n <- A * g1 + B * f1
    f1n <- A * f1 + B * g1
    g2n <- A * g2 + B * f2
    f2n <- A * f2 + B * g2

    # outputs:
    # refl          reflectance
    # tran          transmittance
    # Mb            backward scattering fluorescence matrix, I for PSI and II for PSII
    # Mf            forward scattering fluorescence matrix,  I for PSI and II for PSII

  LRT<- list(lambda =spectral[['wlP']], refl = refl, tran= tran, kChlrel = kChlrel,
             MbI  = g1n,  MbII = g2n,
             MfI  = f1n, MfII = f2n)
  return(LRT)
  }

}
