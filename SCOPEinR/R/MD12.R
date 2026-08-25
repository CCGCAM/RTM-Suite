 ##
#' MD12 algorithm for the computation of fluorescence yield
#'
#' \code{MD12} computes the PSII fluorescence quantum yield from the electron
#' transport rate and a set of first-order rate constants describing the
#' competing de-excitation pathways of PSII (photochemistry, fluorescence,
#' constitutive and regulated thermal dissipation). The yield is evaluated
#' both under CO2-limited and under light-limited conditions, and the
#' smaller (limiting) value is returned, following the model described in
#' Van der Tol et al. (2014, Biogeosciences) building on Magnani et al. (2012).
#'
#' @param ps numeric vector. PSII photochemical yield (probability that an absorbed photon leads to charge separation), dimensionless \[E/E\], range 0-1.
#' @param Ja numeric vector. Actual (CO2-limited) electron transport rate, in the same units as \code{Jms} (typically umol electrons m^-2 s^-1).
#' @param Jms numeric vector. Potential (light-limited) electron transport rate, reduced for PSII photodamage, in the same units as \code{Ja}.
#' @param kps numeric vector. Rate constant for photochemistry at PSII (steady state), dimensionless relative rate constant.
#' @param kf numeric vector. Rate constant for fluorescence emission, dimensionless relative rate constant.
#' @param kds numeric vector. Rate constant for constitutive (basal, dark-type) thermal dissipation at PSII, steady state, dimensionless relative rate constant.
#' @param kDs numeric vector. Rate constant for regulated (NPQ-type) thermal dissipation at PSII, steady state, dimensionless relative rate constant.
#'
#' @return numeric vector. PSII fluorescence quantum yield \code{fs}, dimensionless \[E/E\], computed as the minimum of the CO2-limited and light-limited yields.
#' @export
#'
#' @author 	Christiaan van der Tol (Original version in Matlab)
#' @author Carlos Camino (Ported version into R)
#'
#' @examples
#' \dontrun{
#' fs <- MD12(ps, Ja, Jms, kps, kf, kds, kDs)
#' }
MD12<- function(ps,Ja,Jms,kps,kf,kds,kDs){

fs1    = ps * (kf / kps) / (1 - Ja / Jms)         # [E/E]   PSII fluorescence yield under CO2-limited conditions

par1   = kps / (kps - kds)                            # [E/E]   empirical parameter in the relationship under light-limited conditions
par2   = par1 * (kf + kDs + kds)/kf                  # [E/E]   empirical parameter in the relationship under light-limited conditions
fs2    = (par1 - ps) / par2                           # [E/E]   PSII fluorescence yield under light-limited conditions

fs     = min(fs1,fs2)                              # [E/E]   PSII fluorescence yield
return(fs)
}
