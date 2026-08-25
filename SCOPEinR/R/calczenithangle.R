
#' Calculates pi/2-the angle of the sun with the slope of the surface.
#'
#' @param Doy   day of the year
#' @param t  time of the day (hours, GMT)
#' @param Omega_g slope azimuth angle (deg)
#' @param Fi_gm  slope of the surface (deg)
#' @param Long  Longitude (decimal)
#' @param Lat Latitude (decimal)
#'
#' @return A list \code{zenith.angles} intended to contain:
#' \describe{
#'   \item{Omega_s}{Hour angle of the sun (rad).}
#'   \item{Fi_s}{"Classic" solar zenith angle, perpendicular to the horizontal plane (rad).}
#'   \item{Fi_g}{Projected slope of the surface, in the plane through the solar beam and the vertical (rad).}
#'   \item{Fi_gs}{Angle of the sun with the vector perpendicular to the (possibly sloped) surface, i.e. \code{pi/2} minus the surface solar elevation angle (rad).}
#' }
#' @export
#'
#' @details Last updates: Jan 2003; Oct 2008 by Joris Timmermans (corrected equation of time); Oct 2012 (CvdT) comment: input time is GMT, not local time.
#'
#' @author 	Christiaan van der Tol  (Original version in Matlab)
#' @author Carlos Camino (Ported version into R)
#' @examples
#' \dontrun{
#' zenith.angles <- calczenithangle(Doy = 180, t = 12, Omega_g = 210,
#'                                   Fi_gm = 30, Long = 13.75, Lat = 45.5)
#' }
calczenithangle<-function(Doy,t,Omega_g,Fi_gm,Long,Lat){

  zenith.angles<-list()
  #parameters (if not already supplied)
  if (is.null(Long) | is.null(Lat)){
    Long = 13.75;  # longitude
    Lat = 45.5;     # latitude
    if (is.null(Omega_g) | is.null(Fi_gm)){
      Omega_g = 210;  # aspect
      Fi_gm = 30;      # slope angle
    }
  }
  #convert angles into radials
  G = (Doy - 1) / 365 * 2 * pi;           # converts day of year to radials
  Omega_g = Omega_g / 180 * pi;             # converts direction of slope to radials
  Fi_gm = Fi_gm / 180 * pi;               # converts maximum slope to radials
  Lat = Lat / 180 * pi;                 # converts latitude to radials

  #computes the declination of the sun
  d = 0.006918 - 0.399912 * cos(G) + 0.070247 * sin(G) - 0.006758 * cos(2 * G) + 0.000907 * sin(2 * G) - 0.002697 * cos(3 * G) + 0.00148 * sin(3 * G);

  #Equation of time
  Et  = 0.017 + 0.4281 * cos(G) - 7.351 * sin(G) - 3.349 * cos(2 * G) - 9.731 * sin(2 * G);

  #computes the time of the day when the sun reaches its highest angle
  tm = 12+(4*(-Long) - Et ) / 60;      # de Pury and Farquhar (1997), Iqbal (1983)

  # output:
  # Fi_s      'classic' zenith angle: perpendicular to horizontal plane
  # Fi_gs     solar angle perpendicular to surface slope
  # Fi_g      projected slope of the surface in the plane through the solar beam and the vertical
  #
  #computes the hour angle of the sun
  Omega_s = (t - tm) / 12 * pi;
  zenith.angles$Omega_s <- Omega_s
  #computes the zenithangle (equation 3.28 in De Bruin)
  Fi_s = acos(sin(d) * sin(Lat) + cos(d) * cos(Lat) * cos(Omega_s));
  zenith.angles$Fi_s <- Fi_s
  #computes the slope of the surface Fi_g in the same plane as the solar beam
  Fi_g  = atan(tan(Fi_gm) * cos(Omega_s - Omega_g));
  zenith.angles$Fi_g <- Fi_g
  #computes the angle of the sun with the vector perpendicular to the surface
  Fi_gs = Fi_s + Fi_g;
  zenith.angles$Fi_gs <- Fi_gs

  return(zenith.angles)
}
