

##############################################################################
### Biochemical model  get.biochemical.MD12
##############################################################################

#' MD12 algorithm for the computation of fluorescence yield
#'
#' @param ps numeric. Photochemical yield (fraction of absorbed light used in photochemistry).
#' @param Ja numeric. Actual electron transport rate.
#' @param Jms numeric. Maximum (light-saturated) electron transport rate.
#' @param kps numeric. Rate constant for photochemistry.
#' @param kf numeric. Rate constant for fluorescence.
#' @param kds numeric. Rate constant for (light-independent) thermal dissipation.
#' @param kDs numeric. Rate constant for (light-dependent, sustained) thermal dissipation.
#'
#' @return fs: PSII fluorescence yield.
#' @export
#'
#' @author 	Christiaan van der Tol (Original version in Matlab)
#' @author Carlos Camino (Ported version into R)
#'
#' @examples
#' get.MD12(ps = 0.3, Ja = 40, Jms = 120, kps = 4, kf = 0.05, kds = 0.95, kDs = 0)
get.MD12 <- function(ps, Ja, Jms, kps, kf, kds, kDs) {

  fs1 <- ps * (kf / kps) / (1 - Ja / Jms)   # [E/E]   PSII fluorescence yield under CO2-limited conditions

  par1 <- kps / (kps - kds)                # [E/E]   empirical parameter in the relationship under light-limited conditions
  par2 <- par1 * (kf + kDs + kds) / kf     # [E/E]   empirical parameter in the relationship under light-limited conditions
  fs2 <- (par1 - ps) / par2                # [E/E]   PSII fluorescence yield under light-limited conditions

  fs <- pmin(fs1, fs2)                     # [E/E]   PSII fluorescence yield

  return(fs)

} ### end of function biochemical




#' @title get.BallBerry
#' \code{get.BallBerry}   get BerryBall value from Berry Model
#'
#' @param Cs is CO2 at leaf surface
#' @param RH is relative humidity
#' @param A is Net assimilation in 'same units of CO2 as Cs'/m2/s
#' @param BallBerrySlope parameter for  BallBerry model
#' @param BallBerry0 parameter for  BallBerry model
#' @param minCi minimum Ci as a fraction of Cs (in case RH is very low?)
#' @param Ci_input will use only for Ci == Ci_input
#'
#' @return A list with gs (stomatal conductance) and Ci (intercellular CO2 concentration).
#' @export
#'
#' @author 	Christiaan van der Tol (Original version in Matlab)
#' @author Carlos Camino (Ported version into R)
#'
#' @examples
#' get.BallBerry(Cs = 400, RH = 0.7, A = 15, BallBerrySlope = 9, BallBerry0 = 0.01, minCi = 0.3)
get.BallBerry<-function(Cs, RH, A, BallBerrySlope, BallBerry0, minCi, Ci_input){


  if (missing(Ci_input) == FALSE) {
    # Ci is given: try to compute gs
    Ci <- Ci_input
    gs <- NULL
    if (!is.null(A)) {
      gs <- get.gsFun(Cs, RH, A, BallBerrySlope, BallBerry0)
    }
  } else if (all(BallBerry0 == 0) || is.null(A)) {
    # EXPLANATION:   *at equilibrium* CO2_in = CO2_out => A = gs(Cs - Ci) [1]
    #  so Ci = Cs - A/gs (at equilibrium)                                 [2]
    #  Ball-Berry suggest: gs = m (A RH)/Cs + b   (also at equilib., see Leuning 1990)
    #  if b = 0 we can rearrange B-B for the second term in [2]:  A/gs = Cs/(m RH)
    #  Substituting into [2]
    #  Ci = Cs - Cs/(m RH) = Cs ( 1- 1/(m RH)  [ the 1.6 converts from CO2- to H2O-diffusion ]
    Ci <- pmax(minCi * Cs, Cs * (1 - 1.6 / (BallBerrySlope * RH)))
    gs <- NULL
  } else {

    # EXPLANATION:
    #  if b > 0  Ci = Cs( 1 - 1/(m RH + b Cs/A) )
    # if we use Leuning 1990, Ci = Cs - (Cs - Gamma)/(m RH + b(Cs - Gamma)/A)  [see def of Gamma, above]
    # note: the original B-B units are A: umol/m2/s, ci ppm (umol/mol), RH (unitless)
    #   Cs input was ppm but was multiplied by ppm2bar above, so multiply A by ppm2bar to put them on the same scale.
    #  don't let gs go below its minimum value (i.e. when A goes negative)

    gs <- get.gsFun(Cs, RH, A, BallBerrySlope, BallBerry0)
    Ci <- pmax(minCi * Cs, Cs - 1.6 * A / gs)
  }

  return(list(gs = gs,Ci = Ci))
}



#' @title get.gsFun
#' \code{get.gsFun}   get stomatal conductance
#'
#' @param Cs numeric. CO2 concentration at the leaf surface.
#' @param RH numeric. Relative humidity (0-1).
#' @param A numeric. Net assimilation rate, in the same CO2 units as Cs, per m2 per s.
#' @param BallBerrySlope numeric. Slope parameter of the Ball-Berry stomatal conductance model.
#' @param BallBerry0 numeric. Intercept (minimum conductance) parameter of the Ball-Berry model.
#'
#' @return gs: stomatal conductance.
#' @export
#'
#' @author 	Christiaan van der Tol (Original version in Matlab)
#' @author Carlos Camino (Ported version into R)
#'
#' @examples
#' get.gsFun(Cs = 400, RH = 0.7, A = 15, BallBerrySlope = 9, BallBerry0 = 0.01)
get.gsFun<- function(Cs, RH, A, BallBerrySlope, BallBerry0){
  # add in a bit just to avoid div zero. 1 ppm = 1e-6
  #(note since A < 0 if Cs ==0, it gives a small gs rather than maximal gs
  gs = pmax(BallBerry0,  BallBerrySlope * A * RH / (Cs + 1e-9)  + BallBerry0)
  # clean it up:
  # gs[ Cs == 0 ] <- max_gs  # eliminate infinities
  gs[is.nan(Cs) | is.infinite(Cs)] <- NaN

  return(gs)
}


##############################################################################
### Biochemical model  get.biochemical
##############################################################################


#' @title get.Fluorescence.model
#' \code{get.Fluorescence.model} Fluorescence model
#'
#' @param ps numeric. Photochemical yield (fraction of absorbed light used in photochemistry).
#' @param x numeric. Normalized (0-1) intermediate variable driving non-photochemical quenching, passed in to avoid recomputing it internally.
#' @param Kp numeric. Rate constant for photochemistry.
#' @param Kf numeric. Rate constant for fluorescence.
#' @param Kd numeric. Rate constant for (light-independent) thermal dissipation.
#' @param Knparams numeric vector of length 3: Kno (max NPQ rate constant), alpha, beta - shape parameters of the non-photochemical quenching response.
#'
#' @return A list with eta, qE, qQ, fs, fo, fm, fo0, fm0, and Kn - fluorescence yields and quenching parameters.
#' @export
#'
#' @examples
#' get.Fluorescence.model(ps = 0.3, x = 0.5, Kp = 4, Kf = 0.05, Kd = 0.95,
#'                         Knparams = c(2.48, 2.83, 0.114))
#'
get.Fluorescence.model <- function(ps, x, Kp, Kf, Kd, Knparams) {
  # note: x isn't strictly needed as an input parameter but it avoids code-duplication (of po0) and it's inherent risks.

  Kno <- Knparams[1]
  alpha <- Knparams[2]
  beta <- Knparams[3]

  x_alpha <- exp(log(x) * alpha)  # this is the most expensive operation in this fn; doing it twice almost doubles the time spent here (MATLAB 2013b doesn't optimize the duplicate code)
  Kn <- Kno * (1 + beta) * x_alpha / (beta + x_alpha)

  fo0 <- Kf / (Kf + Kp + Kd)        # dark-adapted fluorescence yield Fo,0
  fo <- Kf / (Kf + Kp + Kd + Kn)     # light-adapted fluorescence yield in the dark Fo
  fm <- Kf / (Kf + Kd + Kn)     # light-adapted fluorescence yield Fm
  fm0 <- Kf / (Kf + Kd)        # dark-adapted fluorescence yield Fm
  fs <- fm * (1 - ps)            # steady-state (light-adapted) yield Ft (aka Fs)
  eta <- fs / fo0
  qQ <- 1 - (fs - fo) / (fm - fo)    # photochemical quenching
  qE <- 1 - (fm - fo) / (fm0 - fo0)  # non-photochemical quenching

  return(list(eta=eta, qE=qE, qQ=qQ, fs=fs, fo=fo, fm=fm, fo0=fo0, fm0=fm0, Kn=Kn))
}



#' @title get.high.temp.inhibtionC3
#' \code{get.high.temp.inhibtionC3} High Temperature Inhibition Function:The following function pertains to C3 photosynthesis
#'
#' @param Tref numeric. Reference temperature, in Kelvin.
#' @param R numeric. Universal gas constant.
#' @param T numeric. Leaf temperature, in Kelvin.
#' @param deltaS numeric. Entropy term of the high-temperature inhibition function.
#' @param deltaHd numeric. Deactivation energy (enthalpy of high-temperature inhibition).
#'
#' @return fHTv: high-temperature inhibition factor applied to C3 photosynthesis rates.
#' @export
#'
#' @examples
#' get.high.temp.inhibtionC3(Tref = 298.15, R = 8.314, T = 305, deltaS = 650, deltaHd = 200000)
get.high.temp.inhibtionC3 <- function(Tref, R, T, deltaS, deltaHd) {
  # High Temperature Inhibition Function
  hightempfunc_num <- (1 + exp((Tref * deltaS - deltaHd) / (Tref * R)))
  hightempfunc_deno <- (1 + exp((deltaS * T - deltaHd) / (R * T)))
  fHTv <- hightempfunc_num / hightempfunc_deno
  return(fHTv)
}


#' @title get.temperature.functionC3
#' \code{get.temperature.functionC3} Temperature Correction Functions:The following function pertains to C3 photosynthesis
#'
#' @param Tref numeric. Reference temperature, in Kelvin.
#' @param R numeric. Universal gas constant.
#' @param Temp numeric. Leaf temperature, in Kelvin.
#' @param deltaHa numeric. Activation energy of the temperature-dependent process.
#'
#' @return fTv: Arrhenius-type temperature correction factor.
#' @export
#'
#' @examples
#' get.temperature.functionC3(Tref = 298.15, R = 8.314, Temp = 305, deltaHa = 65330)
get.temperature.functionC3<-function(Tref,R,Temp,deltaHa){
  # Temperature function
  tempfunc1 <- (1 - Tref / Temp);
  fTv <- exp(deltaHa / (Tref * R) * tempfunc1);
  return(fTv)
}


#' @title get.Ci.next
#' \code{get.Ci.next} Function to calculate the difference between "guessed" Ci (Ci_in) and Ci computed using BB after computing A
#' Test-function for iteration (note that it assigns A in the function's context.)
#' As with the next section, this code can be read as if the function body executed at this point. (if iteration was used).
#' In other words, A is assigned at this point in the file (when iterating).
#'
#' @param Ci_in numeric. The "guessed" intercellular CO2 concentration to test.
#' @param Cs numeric. CO2 concentration at the leaf surface.
#' @param RH numeric. Relative humidity (0-1).
#' @param minCi numeric. Minimum Ci as a fraction of Cs.
#' @param BallBerrySlope numeric. Slope parameter of the Ball-Berry stomatal conductance model.
#' @param BallBerry0 numeric. Intercept (minimum conductance) parameter of the Ball-Berry model.
#' @param A_fun function. Assimilation function, called as \code{A_fun(Cs)} to compute net assimilation.
#' @param ppm2bar numeric. Conversion factor from ppm to bar for CO2 partial pressure.
#'
#' @return A list with err (difference between guessed and Ball-Berry-computed Ci) and Ci_out (the Ball-Berry result).
#' @export
#'
#' @examples
#' \dontrun{
#' get.Ci.next(Ci_in = 280, Cs = 400, RH = 0.7, minCi = 0.3, BallBerrySlope = 9,
#'             BallBerry0 = 0.01, A_fun = function(Cs) list(A = 15), ppm2bar = 1e-6)
#' }
get.Ci.next <- function(Ci_in, Cs, RH, minCi, BallBerrySlope, BallBerry0, A_fun, ppm2bar) {
  # compute A
  Assimilation.variables <- A_fun(Cs)
  if (is.null(Assimilation.variables$A)){
    A_bar =NULL
  } else {
    A_bar <- Assimilation.variables$A * ppm2bar
  }


  # compute Ci using Ball-Berry model
  Ci_out <- get.BallBerry(Cs, RH, A_bar, BallBerrySlope, BallBerry0, minCi)

  # calculate the error
  err <- Ci_out[[2]] - Ci_in

  return(list(err=err, Ci_out= Ci_out))
}





#'  quadratic formula, root of least magnitude
#'
#' @param a numeric. Quadratic coefficient (a in a*x^2 + b*x + c = 0).
#' @param b numeric. Linear coefficient.
#' @param c numeric. Constant term.
#' @param dsign numeric. Sign of the discriminant term to select which root to return (+1 or -1).
#'
#' @return The selected root of least magnitude.
#' @export
#'
#' @examples
#' sel_root(a = 1, b = -3, c = 2, dsign = 1)
sel_root <- function(a, b, c, dsign) {
  # sel_root - select a root based on the fourth arg (dsign <- discriminant sign)
  # for the eqn ax^2 + bx + c,
  # if dsign is:
  #    -1, 0: choose the smaller root
  #    +1: choose the larger root
  # NOTE: technically, we should check a, but in biochemical, a is always > 0

  if (a == 0) {
    x <- -c / b  # note: this works because 'a' is a scalar parameter!
  } else {
    if (any(dsign == 0)) {
      dsign[dsign == 0] <- -1  # technically, dsign==0 iff b <- c <- 0, so this isn't strictly necessary except, possibly for ill-formed cases
    }
    # disc_root <- sqrt(b.^2 - 4.*a.*c)  # square root of the discriminant (doesn't need a separate line anymore)
    # in R, assigning the intermediate variable isn't necessary and can slow down the code
    x <- (-b + dsign * sqrt(b^2 - 4 * a * c)) / (2 * a)
  }
  return(x)
}


#' @title get.computeA
#' \code{get.computeA} Compute the net CO2 assimilation rate using the Farquhar model
#'  Note: even though computeA() is written as a separate function,
#'  the code is, in fact, executed exactly this point in the file (i.e. between the previous if clause and the next section
#'
#' @param Ci is internal CO2 concentration
#' @param Type is photosynthetic pathway type ("C3" or "C4")
#' @param g_m is mesophyll conductance
#' @param Vs_C3 is maximum carboxylation rate for C3 plants
#' @param MM_consts is Michaelis-Menten constants for Rubisco and oxygenase reactions
#' @param Rd is dark respiration rate
#' @param Vcmax  is maximum carboxylation rate
#' @param Gamma_star CO2 compensation point in the absence of mitochondrial respiration
#' @param Je electron transport rate in the chloroplasts
#' @param effcon electron transport to carboxylation efficiency
#' @param atheta curvature parameter for the response of electron transport to irradiance
#' @param kpepcase the relative limitation of electron transport by irradiance
#'
#' @details
#' Even though \code{computeA()}-style logic is written as a separate
#' function here, the calculation is executed exactly at this point in the
#' original SCOPE model's control flow (i.e., between the previous if
#' clause and the next section) — this R port keeps that same ordering.
#'
#' @return A: net CO2 assimilation rate (and related intermediate variables, depending on Type/g_m branch taken).
#' @export
#'
get.computeA <- function(Ci, Type, g_m, Vs_C3, MM_consts, Rd, Vcmax, Gamma_star, Je, effcon, atheta, kpepcase) {
  .scopeinr_state$fcount <- 0

  if (missing(Ci)) {
    .scopeinr_state$fcount <- 0
    return(invisible())
  }

  if (Type == "C3") {
    Vs <- Vs_C3
    if (any(g_m < Inf)) {

      Vc <- sel_root(1 / g_m, -(MM_consts + Ci + (Rd + Vcmax) / g_m), Vcmax * (Ci - Gamma_star + Rd / g_m), -1)
      Ve <- sel_root(1 / g_m, -(Ci + 2 * Gamma_star + (Rd + Je * effcon) / g_m), Je * effcon * (Ci - Gamma_star + Rd / g_m), -1)
      CO2_per_electron <- Ve / Je

    } else {

      Vc <- Vcmax * (Ci - Gamma_star) / (MM_consts + Ci)
      CO2_per_electron <- (Ci - Gamma_star) / (Ci + 2 * Gamma_star) * effcon
      Ve <- Je * CO2_per_electron
    }
  } else {  # C4

    Vc <- Vcmax
    Vs <- kpepcase * Ci
    CO2_per_electron <- effcon
    Ve <- Je * CO2_per_electron
  }


  V <- sel_root(atheta, -(Vc + Ve), Vc * Ve, sign(-Vc))  # i.e., sign(Gamma_star - Ci)
  Ag <- sel_root(0.98, -(V + Vs), V * Vs, -1)
  A <- Ag - Rd

  .scopeinr_state$fcount <- .scopeinr_state$fcount + 1

  if (missing(Ci)) {
    return(invisible())
  } else {
    return(list(A=A, Ag = Ag, Vc = Vc, Vs = Vs, Ve = Ve, CO2_per_electron = CO2_per_electron, fcount = .scopeinr_state$fcount))
  }
}




