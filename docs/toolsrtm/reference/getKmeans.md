# K-means Classification on a RasterStack This function performs a k-means unsupervised classification of a stack of rasters. Number of clusters can be specified, as well as number of iterations and starting sets. An optional geographic weighting system can be turned on that constrains clusters to a geographic area, by including coordinates in the clustering. All variables are normalized before clustering is performed. https://rdrr.io/github/ozjimbob/ecbtools/man/raster.kmeans.html

K-means Classification on a RasterStack This function performs a k-means
unsupervised classification of a stack of rasters. Number of clusters
can be specified, as well as number of iterations and starting sets. An
optional geographic weighting system can be turned on that constrains
clusters to a geographic area, by including coordinates in the
clustering. All variables are normalized before clustering is performed.
https://rdrr.io/github/ozjimbob/ecbtools/man/raster.kmeans.html

## Usage

``` r
getKmeans(stack, k = 12, iter.max = 100, nstart = 10, geo = T, geo.weight = 1)
```

## Arguments

- stack:

  A RasterStack object, or a string pointing to the directory where all
  raster layers are stored. All rasters should have the same extent,
  resolution and coordinate system.

- k:

  Number of clusters to classify to.

- iter.max:

  Maximum number of iterations allowed

- nstart:

  Number of random sets to be chosen.

- geo:

  True/False - should geographic weighting be used?

- geo.weight:

  A weighting multiplier indicating the strength of geographic weighting
  relative to the other variables. A value of 1 gives equal weight.

## Value

Kmenans
