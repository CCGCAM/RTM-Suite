#' LIBERTY model
#'
#' Leaf radiative transfer model designed for conifer needles.
#'
#' Note that these coefficients are based on the original values from Dawson et
#' al. (1998). Improved specific absorption coefficients are available from
#' later work by Di Vittorio (2009).
#'
#' @references Dawson, T. P., Curran, P. J., & Plummer, S. E. (1998). LIBERTY—Modeling the
#' Effects of Leaf Biochemical Concentration on Reflectance Spectra. Remote
#' Sensing of Environment, 65(1), 50–60.
#' https://doi.org/10.1016/S0034-4257(98)00007-8
#' @references Di Vittorio, A. V. (2009). Enhancing a leaf radiative transfer
#' model to estimate concentrations and in vivo specific absorption
#' coefficients of total carotenoids and chlorophylls a and b from
#' single-needle reflectance and transmittance. Remote Sensing of Environment,
#' 113(9), 1948–1966. https://doi.org/10.1016/j.rse.2009.05.002
#'

#' @return Length 4 list containing modeled reflectance, transmittance,
#' wavelengths, and R (?).
#' @author Alexey Shiklomanov and modifed by Carlos Camino
#' @export
#' 
liberty <- function(inputLUT)   {
  
  ## parameters for the leaf LIBERTY model
  cell.d=inputLUT[,'cell.d'];  #Cell diameter. Average leaf cell diameter (m-6). Typical value 40 (20-100).
  #Intercellular air space. Determinant for the amount of radiative
  # flux passing between cells. Typical value 0.045 (0.01-0.1)
  inter.c=inputLUT[,'inter.c']; 
  #baseline.abs absorption. Wavelength-independent absorption to
  # compensate for changes in absolute reflectance. Typical value 0.0006 for
  # fresh leaves or or 0.0004 for dry leaves.
  baseline.abs=inputLUT[,'baseline.abs']; 
  #Leaf thickness. Arbitrary value to determine single leaf
  #reflectance and transmittance from infinite reflectance criteria. Typical
  #value 1.6 (1-10).
  leaf.thick=inputLUT[,'leaf.thick'];
  #Albino absorption. Absorption in the visible region due to lignin. Typical value 2 (0-4).
  albino.abs=inputLUT[,'albino.abs'];
  Cab=inputLUT[,'Cab'] *10 ; #Chlorophyll content (mg m-2). Typical value 200 (0-600).
  EWT=inputLUT[,'EWT'] * 10000 ; #Water content (g m-2). Typical value 100 (0-500).
  #Lignin and cellulose content. Combined lignin and cellulose. content (mg m-2). Typical value 40 (10-80).
  lign.cell=inputLUT[,'lign.cell']*10; 
  Nitrogen=inputLUT[,'Nitrogen'] #Nitrogen content (g m-2). Typical value 1  (0.3-2.0).

  
  k.chl <- ToolsRTM::dataspec.liberty[, 1]
  k.water <- ToolsRTM::dataspec.liberty[, 2]
  ke <- ToolsRTM::dataspec.liberty[, 3]
  k.ligcell <- ToolsRTM::dataspec.liberty[, 4]
  k.proteins <- ToolsRTM::dataspec.liberty[, 5]

  nwl <- 420
  transR <- numeric(nwl)
  reflR <- numeric(nwl)
  wavelength <- numeric(nwl)
  RR <- numeric(nwl)

  for (i in seq(1, nwl)) {

    coeff <- cell.d * (baseline.abs +
                  (k.chl[i] * Cab) +
                  (k.water[i] * EWT) +
                  (ke[i] * albino.abs) +
                  (k.ligcell[i] * lign.cell) +
                  (k.proteins[i] * Nitrogen))

    # change of refractive index over wavelength...
    N1 <- 1.4891 - (0.0005 * (i - 1))
    N0 <- 1.0
    in_angle <- 59

    # Index of refraction
    alpha <- in_angle * pi / 180
    beta <- asin((N0 / N1) * sin(alpha))
    # para_r <- (tan(alpha - beta)) / (tan(alpha + beta))
    # vert_r <- -(sin(alpha - beta)) / (sin(alpha + beta))

    me <- 0
    width <- pi / 180
    pp <- seq(1, 90)
    alpha <- pp * pi / 180
    beta <- asin(N0 / N1 * sin(alpha))
    plus <- alpha + beta
    dif <- alpha - beta
    refl <- 0.5 * ((sin(dif))^2 / (sin(plus))^2 + (tan(dif))^2 / (tan(plus))^2)
    me <- 2 * sum(refl * sin(alpha) * cos(alpha) * width)

    mi <- 0
    mint <- 0
    width <- pi / 180
    critical <- asin(N0 / N1) * 180 / pi
    jj <- seq(1, critical)
    alpha <- jj * pi / 180
    beta <- asin(N0 / N1 * sin(alpha))
    plus <- alpha + beta
    dif <- alpha - beta
    refl <- 0.5 * ((sin(dif))^2 / (sin(plus))^2 + (tan(dif))^2 / (tan(plus))^2)
    mint <- sum(refl * sin(alpha) * cos(alpha) * width)
    mi <- (1 - sin(critical * pi / 180)^2) + 2 * mint

    M <- 2 * (1 - (coeff + 1) * exp(-coeff)) / coeff^2
    Tt <- ((1 - mi) * M) / (1 - (mi * M))
    x <- inter.c / (1 - (1 - 2 * inter.c) * Tt)
    a <- me*Tt + x*Tt - me - Tt - x*me*Tt
    b <- 1 + x*me*Tt - 2*x^2*me^2*Tt
    c <- 2*me*x^2*Tt - x*Tt - 2*x*me
    R <- 0.5
    # Initial guess...
    for (iterations in seq(1, 50)) {
      next_R <- -(a * R^2 + c) / b
      R <- next_R
    }

    # The next bit works out transmittance based upon Benford...
    # setting up unchanging parameters...
    rb <- (2 * x * me) + (x * Tt) - (x * Tt * 2 * x * me)
    tb <- sqrt(((R - rb) * (1 - (R * rb))) / R)
    whole <- fix(leaf.thick)
    fraction <- leaf.thick - whole
    #
    # % 	/* The next bit works out the fractional value */
    # % 	/* for the interval between 1 and 2... */
    top <- tb^(1 + fraction) * ((((1 + tb)^2) - rb^2)^(1 - fraction))
    bot1 <- (1 + tb)^(2 * (1 - fraction)) - rb^2
    bot2 <- 1 + (64 / 3) * fraction * (fraction - 0.5) * (fraction - 1) * 0.001
    tif <- top / (bot1 * bot2)
    rif <- (1 + rb^2 - tb^2 - sqrt((1 + rb^2 - tb^2)^2 - 4 * rb^2 * (1 - tif^2))) / (2 * rb)

    # Now to work out for integral thickness greater than 2...
    if (whole >= 2) {
      prev_t <- 1
      prev_r <- 0
      for (step in seq(1, (whole - 1))) {
        cur_t <- (prev_t * tb) / (1 - (prev_r * rb))
        cur_r <- prev_r + (((prev_t * prev_t) * rb) / (1 - (prev_r * rb)))
        prev_t <- cur_t
        prev_r <- cur_r
      }
    } else {
      cur_t <- 1
      cur_r <- 0
    }
    trans <- (cur_t * tif) / (1 - (rif * cur_r))
    refl <- cur_r + ((cur_t^2 * rif) / (1 - (rif * cur_r)))

    transR[i] <- trans
    reflR[i] <- refl
    wavelength[i] <- 400 + ((i - 1) * 5)
    RR[i] <- R
  }

  ### Simulation at 1 nm
  wave_1nm <-c(400:2500)
  reflR_1nm = signal::interp1(wavelength, reflR, wave_1nm, extrap = T)
  transR_1nm = signal::interp1(wavelength, transR, wave_1nm, extrap = T)
  RR_1nm = signal::interp1(wavelength, RR, wave_1nm, extrap = T)
  
  LRT <- list(lambda = wave_1nm,refl = reflR_1nm,tran = transR_1nm,RR = RR_1nm)

  return(LRT)
}

# https://www.mathworks.com/help/matlab/ref/fix.html
fix <- function(x) trunc(x)
