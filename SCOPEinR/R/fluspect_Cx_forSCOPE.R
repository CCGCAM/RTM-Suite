
#' Leaf FLUSPECT-B-Cx model for SCOPE
#'
#' \code{getFluspect.Cx.SCOPE} calculates reflectance and transmittance spectra of a leaf using FLUSPECT-B,
#' calculates reflectance and transmittance spectra of a leaf using FLUSPECT-B,
#'plus four excitation-fluorescence matrices
#'
#' @author Wout Verhoef, Christiaan van der Tol, Joris Timmermans, Nastassia Vilfan (Original version in Matlab)
#' @author Carlos Camino (Ported version into R)
#'
#' Matlab version: Date: 2007-2020
#'
#' Ported in R: Jun-2023
#'
#' @param inputsLeaf  a LUT with main plnt traits for Fluspect-B Cx
#' @param inputsOptipar  internal parameters
#' @param version  Leaf model uses for Kcal: 'Cx'
#' @param step a step for wavelengths for getting a matrix of Mf and Mb, by default uses a step of 1 nm, This step return matrix of 211x351
#' the SCOPE model uses a step of 5 for getting matrix of 53x71
#' @return  a list which contains:
#' * refl (reflectance)
#' * tran (transmittance)
#' * Mb (backward scattering fluorescence matrix, I for PSI and II for PSII)
#' * Mf (forward scattering fluorescence matrix,  I for PSI and II for PSII)
#' @export
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
#'
#' @examples
#' \dontrun{
#' inputs <- ToolsRTM::inputsRTMs
#' LUT <- as.data.frame(ToolsRTM::getLUT(inputs = inputs, nLUT = 1, setseed = 1234))
#' sim <- getFluspect.Cx.SCOPE(inputsLeaf = LUT, inputsOptipar = SCOPEinR::optipar, version = 'Cx')
#' }

getFluspect.Cx.SCOPE<-function(inputsLeaf,inputsOptipar, version = 'Cx', step)  {

  if (is.null(version)){
    version = 'SCOPE'
  }

  # fixed parameters for the fluorescence module
  ndub = 15;           # number of doublings applied
  if (missing(step)){
    step= 5;
  } else{
    step = step
  }




  # Fluspect parameters
  #inputsLeaf=LUT[1,]
  #inputsOptipar <- ToolsRTM::optipar


  #define alll inputs in the models. retreived from LUT tables
  if  (version ==  'Cx'){
    Cab = inputsLeaf[,'Cab'];Car = inputsLeaf[,'Car'];EWT = inputsLeaf[,'EWT']; LMA = inputsLeaf[,'LMA'];
    Cs = inputsLeaf[,'Cs']; N = inputsLeaf[,'N']; fqe_ = inputsLeaf[,'fqe'];
    fqe <- c(fqe_/5, fqe_)
    Cx = inputsLeaf[,'Cx']  # Violaxanthin - Zeaxanthin transition status [0-1]

    Prot = inputsLeaf[,'Prot']; CBC = inputsLeaf[,'CBC'];
    Anth =inputsLeaf[,'Anth']
  } else if (version ==  'SCOPE') {

    Cab = inputsLeaf[['Cab']];Car = inputsLeaf[['Car']];EWT = inputsLeaf[['EWT']]; LMA = inputsLeaf[['LMA']];
    Cs = inputsLeaf[['Cs']]; N = inputsLeaf[['N']]; fqe_ = inputsLeaf[['fqe']];
    fqe <- c(fqe_/5, fqe_)
    Cx = inputsLeaf[['Cx']]  # Violaxanthin - Zeaxanthin transition status [0-1]

    Prot = inputsLeaf[['Prot']]; CBC = inputsLeaf[['CBC']];
    Anth =inputsLeaf[['Anth']]
  }



  if (LMA > 0 & (Prot > 0 | CBC > 0)) {
      LMA <- 0
      }


  # If inputsOptipar is not present, print a message
  if (missing(inputsOptipar)) {

    ### define inputsOptipar
    optipar_default <- SCOPEinR::optipar2017.ProspectD

    nr = optipar_default[['nr']]; Kdm = optipar_default[['Kdm']]; Kab = optipar_default[['Kab']];

    Kca = optipar_default[['Kca']]; Kw = optipar_default[['Kw']]; Ks = optipar_default[['Ks']];
    Kant = optipar_default[['Kant']];
    Kcbc =optipar_default[['Kcbc']]; Kp = optipar_default[['Kp']];
    phiI =optipar_default[['phiI']]; phiII = optipar_default[['phiII']];
    phi =optipar_default[['phi']];

    if (Cx == -999 | is.na(inputsLeaf[['Cx']])) {
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


    if (Cx == -999 | is.na(inputsLeaf[['Cx']])) {
      # Use old Kca spectrum if this is given as input
      Kca <- inputsOptipar[['Kca']]
    } else {
      # Otherwise make linear combination based on Cx
      # For Cx going from 0 to 1 we go from Viola to Zea
      Kca <- (1 - Cx) * inputsOptipar[['KcaV']] + Cx * inputsOptipar[['KcaZ']]
    }


  }



  #### PROSPECT  PRO calculations
  Kall = (Cab * Kab + Car * Kca + LMA * Kdm + EWT * Kw  + Cs * Ks + Anth * Kant + Prot * Kp + CBC * Kcbc) / N;   # Compact leaf layer
  #plot(Kall,type='l',col='red')
  j <- which(Kall > 0)  # Non-conservative scattering (normal case)
  t1 <- (1 - Kall) * exp(-Kall)
  t2 <- Kall^2 * expint::expint(Kall)
  tau <- rep(1, length(t1))
  tau[j] <- t1[j] + t2[j]
  kChlrel <- rep(0, length(t1))
  kCarrel <- rep(0, length(t1))
  kChlrel[j] <- Cab * Kab[j] / (Kall[j] * N)
  kCarrel[j] <- Car * Kca[j] / (Kall[j] * N)

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
  #plot(kChl,type='l',col='red', lwd=2)
  #  Fluorescence of the leaf mesophyll layer
  # Fluorescence part is skipped for fqe = 0

  if (fqe_ > 0) {

    ## get bands from an object
    spectral <- SCOPEinR::define.bands()
    wlp <- c(spectral[['wlP']])     # PROSPECT wavelengths, kept as a row vector
    if (step == 1){
      wle <- spectral[['wlE']]  # excitation wavelengths, transpose to column
      wlf <- spectral[['wlF']]  # fluorescence wavelengths, transpose to column
    } else {
      wle <- seq(min(spectral[['wlE']]),max(spectral[['wlE']]),step)

      wlf <- seq(min(spectral[['wlF']]),max(spectral[['wlF']]),step-1) # fluorescence wavelengths, transpose to column
    }

    ### get position of wavelengths in the full spectrum (400:2400)
    Iwle <- match(wle, wlp)
    Iwlf <- match(wlf, wlp)

    #### news
    k_iwle = signal::interp1(wlp,k,wle)
    s_iwle = signal::interp1(wlp,s,wle)
    kChl_iwle = signal::interp1(wlp,kChl,wle);
    r21_iwle = signal::interp1(wlp,r21,wle);
    rho_iwle = signal::interp1(wlp,rho,wle);
    tau_iwle = signal::interp1(wlp,tau,wle);
    talf_iwle = signal::interp1(wlp,talf,wle);
    #plot(talf_iwle,type='l',col='red', lwd=2)
    ####
    eps <- 2^(-ndub)
    ####
    # initialisations
    te <- 1 - (k_iwle+ s_iwle) * eps
    tf <- 1 - (k[Iwlf] + s[Iwlf]) * eps

    ### for range between 400-750
    re <- s_iwle * eps
    ## for range 640 and 840
    rf <- s[Iwlf] * eps

    #plot(rf*1e5,type='l',col='red', lwd=2)


    # Calculate sigmoid matrix
    #  Matlab code --> sigmoid= 1 / ( 1 + exp( -wlf / 10)* exp(wle'/10) )  # matrix computed as an outproduct

    #sigmoid_term.a <- matrix(1 + exp( -wlf / 10),ncol=1)
    #sigmoid_term.b <- matrix(exp((wle)/10),nrow = 1)

    #plot(sigmoid_term.a,type='l',col='red', lwd=2)
    #plot(sigmoid[1,],type='l',col='red', lwd=2)

    sigmoid <- 1/(1+ outer(exp(-wlf/10), exp((wle)/10))) #

    #plot(sigmoid[,1],type='l',col='red', lwd=2)
    #phiI <-SCOPEinR::optipar2017.ProspectD$phi
    phi_plot <-phi[Iwlf]


    #plot(phi_plot,type='l',col='red', lwd=2)

    Mf  <- step * fqe[2] * ((0.5 * phi_plot) * eps) %*% t(kChl_iwle) * sigmoid
    #plot(Mf[2,],type='l',col='red', lwd=2)


    #plot(Mf[,1],type='l')
    Mb <- Mf


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

      A11  <- xf %*% Ih + Iv %*% t(xe);
      A12 <- (xf %*% t(xe)) * (rf %*% Ih + Iv %*% t(re))
      A21  <- 1 + (xf %*% t(xe)) * ( 1 + rf %*% t(re))
      A22 <- (xf * rf) %*% Ih + Iv %*% t((xe * re))

      Mfn <- Mf * A11 + Mb * A12;
      Mbn <- Mb * A21 + Mf * A22;

      te <- ten
      re <- ren
      tf <- tfn
      rf <- rfn

      Mf <- Mfn
      Mb <- Mbn


    } ## end ndum


    # Here we add the leaf-air interfaces again for obtaining the final
    # leaf level fluorescences.
    g = Mb; f = Mf;

    Rb <- rho + tau^2 * r21 / (1 - rho * r21)
    Rb_iwle = signal::interp1(wlp,Rb,wle);

    Xe <- Iv %*% t(talf_iwle / (1-r21_iwle * Rb_iwle))



    Xf <- t21[Iwlf]/(1 - r21[Iwlf] * Rb[Iwlf]) %*% Ih
    Ye <- Iv %*% t(tau_iwle * r21_iwle /(1-rho_iwle * r21_iwle))

    Yf <- tau[Iwlf]*r21[Iwlf]/(1-rho[Iwlf] * r21[Iwlf]) %*% Ih

    A <- Xe * (1 + Ye * Yf) * Xf
    B <- Xe * (Ye + Yf) * Xf


    gn <- A * g + B * f
    fn <- A * f + B * g


    # outputs:
    # refl          reflectance
    # tran          transmittance
    # Mb            backward scattering fluorescence matrix, I for PSI and II for PSII
    # Mf            forward scattering fluorescence matrix,  I for PSI and II for PSII
    # kChlrel Relative portion of chlorophyll contribution to reflecntace / transmittance in the spectral range, 400 to 2400 nm

    LRT<- list(lambda =spectral[['wlP']], refl = refl, tran= tran, kChlrel = kChlrel,  kCarrel = kCarrel,
               Mb  = gn,
               Mf  = fn)
    return(LRT)

  }

}
