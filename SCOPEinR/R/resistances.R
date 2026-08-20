
#' @title get.resistances
#'
#' @description
#' \code{get.resistances} calculates the aerodynamic and boundary-layer
#' resistances between the soil, the canopy and the reference (measurement)
#' height, following the two-layer resistance scheme of Wallace and Verhoef
#' (2000), including a stability correction for non-neutral atmospheric
#' conditions based on the Monin-Obukhov length.
#'
#' @param data.soil list. Soil properties; must contain \code{rbs} (soil boundary-layer resistance, s m^-1).
#' @param data.canopy list. Canopy structural properties: \code{Cd} (leaf drag coefficient), \code{LAI} (leaf area index, m^2 m^-2), \code{rwc} (within-canopy aerodynamic resistance, s m^-1), \code{zo} (roughness length for momentum, m), \code{d} (zero-plane displacement height, m), \code{hc} (vegetation height, m), \code{leafwidth} (characteristic leaf width, m).
#' @param data.meteo list. Meteorological data: \code{z} (measurement height, m), \code{u} (wind speed at \code{z}, m s^-1), \code{L} (Monin-Obukhov length, m; see \code{\link{get.Monin.Obukhov}}).
#' @return A list \code{resist_out} with (among others): \code{ustar} (friction velocity, m s^-1), \code{uz0} (wind speed at \code{z0m}, m s^-1), \code{Kh} (eddy diffusivity for heat, m^2 s^-1), \code{rai}/\code{rar}/\code{rac} (aerodynamic resistance in the inertial sublayer / roughness sublayer / canopy layer, s m^-1), \code{rws} (aerodynamic resistance within the canopy down to the soil, s m^-1), \code{raa} (total aerodynamic resistance above the canopy, s m^-1), \code{rawc} (total resistance within the canopy air space, s m^-1), \code{raws} (total resistance between the canopy air space and the soil, s m^-1).
#' @export
#'
#' @author  Anne Verhoef, Christiaan van der Tol, Joris Timmermans (Original version in Matlab)
#' @author Carlos Camino (Ported version into R)
#'
#'   Last update: 01 Feb 2008
#'
#' @references
#'  Wallace and Verhoef (2000) 'Modelling interactions in mixed-plant communities: light,
#'  water and carbon dioxide', in: BruceMarshall, Jeremy A. Roberts (ed), 'Leaf Development and Canopy Growth',
#'  Sheffield Academic Press, UK. ISBN 0849397693
#'
#'  ustar:  Tennekes, H. (1973) 'The logaritmic wind profile', J.Atmospheric Science, 30, 234-238
#'
#'  Psih:   Paulson, C.A. (1970), The mathematical representation of wind speed and temperature in the
#'  unstable atmospheric surface layer. J. Applied Meteorol. 9,857-861
#'
#' @examples
#' \dontrun{
#' resist_out <- get.resistances(data.soil, data.canopy, data.meteo)
#' }
get.resistances<-function(data.soil, data.canopy, data.meteo){
  # Input:
  #   resist_in   aerodynamic resistance parameters and wind speed
  #
  # The strucutre resist_in contains the following elements:
  # u         =   windspeed
  # L         =   stability
  # LAI       =   Leaf Area Index

  # rbs       =   Boundary Resistance of soil                         [s m-1]
  # rss       =   Surface resistance of soil for vapour transport     [s m-1]
  # rwc       =   Within canopy Aerodynamic Resistance canopy         [s m-1]

  # z0m       =   Roughness lenght for momentum for the vegetation    [m]
  # d         =   Displacement height (Zero plane)                    [m]
  # z         =   Measurement height                                  [m]
  # h         =   Vegetation height                                   [m]



  # define the output
  # resist_out  aeorodynamic resistances

  resist_out = list()

  ## parameters
  #global constants

  kappa = subset(  SCOPEinR::constants,constant == 'kappa')[[2]]
  Cd = data.canopy[['Cd']]
  LAI = data.canopy[['LAI']]
  rwc = data.canopy[['rwc']]
  z0m = data.canopy[['zo']] # 0.11
  d = data.canopy[['d']] # 1.64
  h = data.canopy[['hc']]
  w = data.canopy[['leafwidth']];
  z = data.meteo[['z']]
  u = max(0.3,data.meteo[['u']]);
  L = data.meteo[['L']];
  rbs = data.soil[['rbs']];
  #rss       =  resist_in.rss;

  # derived parameters
  #zr: top of roughness sublayer, bottom of intertial sublayer
  zr			= 2.5 * h;                   #                            [m]
  #n: dimensionless wind extinction coefficient                       W&V Eq 33
  n			= Cd * LAI /(2 * kappa^2);      #                            []
  ## stability correction for non-neutral conditions


  unstable <- ifelse(L < 0 & L > -500, TRUE, FALSE)
  stable <- ifelse(L > 0 & L < 500, TRUE, FALSE)
  x = (1 - 16 * z / L)^(1 / 4); # only used for unstable

  # stability correction functions, friction velocity and Kh=Km=Kv
  pm_z <- get.psim(z -d, L, unstable, stable,x)
  ph_z <- get.psih(z -d, L,unstable, stable, x)

  pm_h <- get.psim(h -d, L, unstable, stable, x)

  #ph_h <- get.psih(h -d,L,unstable,stable)

  ph_zr <- get.psih(zr-d,L,unstable,stable,x) * (z >= zr) + ph_z * (z < zr)

  phs_zr <- get.phstar(zr,zr,d,L,stable,unstable,x)
  phs_h <- get.phstar(h ,zr,d,L,stable,unstable,x)

  ustar <- max( 0.001,kappa * u /(log(( z - d)/ z0m) - pm_z)) #W&V Eq 30

  Kh <- kappa * ustar * (zr - d);                  #W&V Eq 35


 if (unstable == TRUE) {
   resist_out[['Kh']] <- Kh * (1 - 16 * (h - d)/L) ^ 0.5 # W&V Eq 35
 } else if (stable == TRUE) {
   resist_out[['Kh']] <- Kh * (1 + 5 * (h - d) / L) ^ -1 # W&V Eq 35
 } else {
   resist_out[['Kh']] <-Kh
 }

  ## to avoid NAN from log with negative values

  if (is.na(ustar/kappa * (log((h-d)/z0m) - pm_h)) == TRUE){
    uh1<- 0
  } else {
    uh1<- ustar/kappa * (log((h-d)/z0m) - pm_h)
  }


  ## wind speed at height h and z0m
  uh <- max(uh1, 0.01)
  uz0 <- uh*exp(n*((z0m+d)/h-1))  # W&V Eq 32

  ## resistances
  # The strucutre resist_out contains the following elements:
  # ustar     =   Friction velocity                                   [m s-1]
  # uz0       =   windspeed at z0                                     [m s-1]

  resist_out[['uz0']] <- uz0
  resist_out[['ustar']] <- ustar

  # rai is aerodynamic resistance in inertial sublayer         [s m-1]
  if (z > zr){
    rai = (1 / (kappa * ustar) * (log((z - d) /(zr - d))  - ph_z   + ph_zr)) # W&V Eq 41
  } else {
    rai = 0
  }

  # rar is Aerodynamic resistance in roughness sublayer        [s m-1]
  rar <- 1 / (kappa * ustar) * ((zr - h) / (zr - d)) - phs_zr + phs_h # W&V Eq 39

  # rac is Aerodynamic resistance in canopy layer (above z0+d) [s m-1]
  rac <- h * sinh(n) / (n * Kh) * (log((exp(n) - 1) / (exp(n) + 1)) - log((exp(n * (z0m + d) / h) - 1) / (exp( n * (z0m + d) / h) + 1)))
  # rws is  Aerodynamic resistance within canopy(soil)          [s m-1]
  rws <- h * sinh(n)/ (n * Kh) * (log((exp(n * (z0m + d) / h) - 1) / (exp(n * (z0m + d) / h) + 1)) - log((exp(n * (0.01) / h) - 1)/(exp(n * (0.01) / h) + 1))) # W&V Eq 43

  #rbc = 70/LAI * sqrt(w./uz0);						#		W&V Eq 31, but slightly different

  # The strucutre resist_out contains the following elements:


  # rai       =   Aerodynamic resistance in inertial sublayer         [s m-1]
  # rar       =   Aerodynamic resistance in roughness sublayer        [s m-1]
  # rac       =   Aerodynamic resistance in canopy layer (above z0+d) [s m-1]
  # rws       =   Aerodynamic resistance within canopy(soil)          [s m-1]

  resist_out[['rai']] <- rai
  resist_out[['rar']] <- rar
  resist_out[['rac']] <- rac
  resist_out[['rws']] <- rws

  #resist_out.rbc = rbc;

  raa  = rai + rar + rac
  rawc = rwc # + rbc;
  raws = rws + rbs

  # raa       =   Aerodynamic resistance above the canopy             [s m-1]
  # rawc      =   Total resistance within the canopy (canopy)         [s m-1]
  # raws      =   Total resistance within the canopy (soil)           [s m-1]

  resist_out[['raa']] <- raa  # aerodynamic resistance above the canopy           W&V Figure 8.6
  resist_out[['rawc']] <- rawc			# aerodynamic resistance within the canopy (canopy)
  resist_out[['raws']] <- raws			# aerodynamic resistance within the canopy (soil)

  # rbc       =   Boundary layer resistance (canopy)                  [s m-1]
  # rwc       =   Aerodynamic Resistance within canopy(canopy)(Update)[s m-1]
  # rbs       =   Boundary layer resistance (soil) (Update)           [s m-1]
  # rss       =   Surface resistance vapour transport(soil)(Update)   [s m-1]
  # Kh        =   Diffusivity for heat                                [m2s-1]

  LRT<- resist_out
  return(LRT)


}


##
#' subfunction pm for stability correction (eg. Paulson, 1970)
#'
#' \code{get.psim} computes the stability correction function for momentum
#' transfer, used to correct the logarithmic wind profile for non-neutral
#' atmospheric conditions.
#'
#' @param z numeric. Height above the zero-plane displacement height, in meters.
#' @param L numeric. Monin-Obukhov length, in meters (atmospheric stability scale).
#' @param unstable logical. TRUE if the atmosphere is in an unstable stratification regime.
#' @param stable logical. TRUE if the atmosphere is in a stable stratification regime.
#' @param x numeric. Stability correction intermediate variable, typically \code{(1 - 16*z/L)^0.25}.
#'
#' @author 	Christiaan van der Tol (Original version in Matlab)
#' @author Carlos Camino (Ported version into R)
#'
#' @return numeric. Momentum stability correction function \code{pm}; \code{0} under neutral conditions.
#' @export
#'
#' @examples
#' get.psim(z = 2, L = -50, unstable = TRUE, stable = FALSE, x = 1.2)
get.psim<-function(z,L,unstable,stable,x){
  pm = 0
  if (unstable == TRUE) {
    pm = 2 * log((1 + x) / 2) + log((1 + x^2) / 2) - 2 * atan(x) + pi / 2 #unstable
  } else if (stable == TRUE) {
    pm = -5 * z / L  # stable
  }

return(pm)
}


#' subfunction ph for stability correction (eg. Paulson, 1970)
#'
#' @param z numeric. Height above the surface, in meters.
#' @param L numeric. Monin-Obukhov length, in meters (atmospheric stability scale).
#' @param unstable logical. TRUE if the atmosphere is in an unstable stratification regime.
#' @param stable logical. TRUE if the atmosphere is in a stable stratification regime.
#' @param x numeric. Stability correction intermediate variable, typically (1 - 16*z/L)^0.25.
#'
#' @author 	Christiaan van der Tol (Original version in Matlab)
#' @author Carlos Camino (Ported version into R)
#'
#' @return ph vector
#' @export
#'
#' @examples
#' get.psih(z = 2, L = -50, unstable = TRUE, stable = FALSE, x = 1.2)
get.psih<-function(z,L,unstable,stable,x){

  ph <- 0
  if (unstable == TRUE) {
    ph <- 2 * log(( 1 + x^2) / 2)  # unstable
  } else if (stable == TRUE) {
    ph <- -5 * z / L  # stable
  }
  return(ph)

}

#' subfunction phs for stability correction (eg. Paulson, 1970)
#'
#' @param z numeric. Height above the surface, in meters.
#' @param zR numeric. Reference height, in meters.
#' @param d numeric. Zero-plane displacement height, in meters.
#' @param L numeric. Monin-Obukhov length, in meters (atmospheric stability scale).
#' @param stable logical. TRUE if the atmosphere is in a stable stratification regime.
#' @param unstable logical. TRUE if the atmosphere is in an unstable stratification regime.
#' @param x numeric. Stability correction intermediate variable, typically (1 - 16*z/L)^0.25.
#'
#' @author 	Christiaan van der Tol (Original version in Matlab)
#' @author Carlos Camino (Ported version into R)
#'
#' @return phs vector
#' @export
#'
#' @examples
#' get.phstar(z = 2, zR = 10, d = 0.7, L = -50, stable = FALSE, unstable = TRUE, x = 1.2)
get.phstar<-function(z,zR,d,L,stable,unstable,x){
  phs = 0

  if (unstable == TRUE) {
    phs  = (z - d) / (zR - d) * (x^2 - 1) / (x^2 + 1) #unstable
  } else if (stable == TRUE) {
    phs = -5 * z / L  # stable
  }
  return(phs)
}

