#' Stern's formula in Lekner & Dorf (1988) gives reflectance for alfa = 90 degrees
#'
#' @param alfa numeric. Maximum incidence solid angle, in degrees.
#' @param nr numeric. Leaf refractive index.
#'
#' @return angles: average transmissivity of a dielectric surface.
#' @export
#'
#' @author 	Christiaan van der Tol (Original version in Matlab)
#' @author Carlos Camino (Ported version into R)
#'
#' @examples
#' tav(59, 1.4)
tav<-function(alfa,nr){

  # Was `nr.^2` -- a MATLAB->R port artifact: MATLAB's `.^` is the
  # element-wise power operator, which R doesn't need/have (`^` already is
  # element-wise); `nr.^2` parses in R as the undefined variable `nr.`
  # (R allows dots in identifiers) raised to the power 2, not `nr^2`.
  n2 =   nr^2;
  np =   n2 + 1;
  nm =   n2 - 1;

  # y1 = (3*n2+2*nr+1)./(3*(nr+1).^2);
  # y2 = 2*nr.^3.*(nr.^2+2*nr-1)./(np.^2.*nm);
  # y3 = n2.*np.*log(nr)./nm.^2;
  # y4 = n2.*nm.^2./np.^3.*log((nr.*(nr+1)./(nr-1)));

  # st = y1-y2+y3-y4;
  a=4
  k=2
  a =   +((nr + 1)^2) / 2;
  k=   -((n2 - 1)^2) / 4;
  # sind() (sine of a degree angle) is a MATLAB built-in with no R
  # equivalent -- was called bare here and never existed in either package.
  sin_a =   sin(alfa * pi / 180);
  #
  if (alfa != 0){
    B2 = sin_a^2 - np / 2;
    if (alfa != 90){
      B1 =  sqrt(B2^2 + k );
    } else{
      B1= 0
    }
    b = B1 - B2;
    b3 =  b^3;
    a3  = a^3;

    ts = (k^2/ (6 * b3) + k / b - b / 2) - (k^2/ (6 * a3) + k / a - a / 2);

    tp1 = -2 * n2 * (b  -  a ) / (np^2);
    tp2 = -2 * n2 * np* (log(b / a)) / (nm^2);
    tp3 = n2 * ( 1 / b - 1 / a ) / 2;

    #     tp4 =   16*n2.^2.* (n2.^2+1) .* ( log(2*np.*b - nm.^2) - log(2*np.*a - nm.^2) ) ./ (np.^3.*nm.^2);
    #     tp5 =   16*n2.^2.* (n2     ) .* ( 1./(2*np.*b - nm.^2) - 1./(2*np.*a - nm.^2)) ./ (np.^3       );

    tp4 =	16 * n2 ^2 * (n2^2 + 1) * (log((2 * np * b - nm^2) /(2 * np * a - nm^2))  ) / (np^3 * nm^2);
    tp5 =  16 * n2^2 * (n2) * ( 1 / (2 * np * b - nm^2) - 1 / (2 * np * a - nm^2)) / (np^3);
    tp =   tp1 + tp2 + tp3 + tp4 + tp5;
    Tav =   (ts + tp)/ (2 * sin_a^2);
  }  else {
    Tav = 4 * nr / ((nr + 1) * (nr + 1));
  }
  return(Tav)
}
