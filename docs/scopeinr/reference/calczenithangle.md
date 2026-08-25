# Calculates pi/2-the angle of the sun with the slope of the surface.

Calculates pi/2-the angle of the sun with the slope of the surface.

## Usage

``` r
calczenithangle(Doy, t, Omega_g, Fi_gm, Long, Lat)
```

## Arguments

- Doy:

  day of the year

- t:

  time of the day (hours, GMT)

- Omega_g:

  slope azimuth angle (deg)

- Fi_gm:

  slope of the surface (deg)

- Long:

  Longitude (decimal)

- Lat:

  Latitude (decimal)

## Value

A list `zenith.angles` intended to contain:

- Omega_s:

  Hour angle of the sun (rad).

- Fi_s:

  "Classic" solar zenith angle, perpendicular to the horizontal plane
  (rad).

- Fi_g:

  Projected slope of the surface, in the plane through the solar beam
  and the vertical (rad).

- Fi_gs:

  Angle of the sun with the vector perpendicular to the (possibly
  sloped) surface, i.e. `pi/2` minus the surface solar elevation angle
  (rad).

## Details

Last updates: Jan 2003; Oct 2008 by Joris Timmermans (corrected equation
of time); Oct 2012 (CvdT) comment: input time is GMT, not local time.

## Author

    Christiaan van der Tol  (Original version in Matlab)

Carlos Camino (Ported version into R)

## Examples

``` r
if (FALSE) { # \dontrun{
zenith.angles <- calczenithangle(Doy = 180, t = 12, Omega_g = 210,
                                  Fi_gm = 30, Long = 13.75, Lat = 45.5)
} # }
```
