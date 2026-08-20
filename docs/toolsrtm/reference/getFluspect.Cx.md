# Leaf FLUSPECT-B-Cx model

The outputs are refl (reflectance) tran (transmittance) Mb (backward
scattering fluorescence matrix, I for PSI and II for PSII) Mf (forward
scattering fluorescence matrix, I for PSI and II for PSII)

Date: 2007 Update from PROSPECT to FLUSPECT: January 2011 (CvdT)

Nov 2012 (CvdT) Output EF-matrices separately for PSI and PSII 31 Jan
2013 (WV) Adapt to SCOPE v_1.40, using structures for I/O 30 May 2013
(WV) Repair bug in s for non-conservative scattering 24 Nov 2013 (WV)
Simplified doubling routine 25 Nov 2013 (WV) Restored piece of code that
takes final refl and tran outputs as a basis for the doubling routine 03
Dec 2013 (WV) Major upgrade. Border interfaces are removed before the
fluorescence calculation and later added again 23 Dec 2013 (WV) Correct
a problem with N = 1 when calculating k and s; a test on a = Inf was
included 01 Apr 2014 (WV) Add carotenoid concentration (Cca and Kca) 19
Jan 2015 (WV) First beta version for simulation of PRI effect 20 Jan
2021 (CvdT) Include PROSPECT-PRO coefficients

## Usage

``` r
getFluspect.Cx(inputsLeaf, inputsOptipar, version = "Cx")
```

## Arguments

- inputsLeaf:

  a LUT with main plnt traits for Fluspect-B Cx

- inputsOptipar:

  internal parameters

- version:

  Leaf model uses for Kcal: 'Cx'

## Value

a leaf model refl and trans

## Details

`getFluspect.Cx` calculates reflectance and transmittance spectra of a
leaf using FLUSPECT-B, calculates reflectance and transmittance spectra
of a leaf using FLUSPECT-B, plus four excitation-fluorescence matrices

## Author

Wout Verhoef, Christiaan van der Tol, Joris Timmermans, Nastassia Vilfan
(Original version in Matlab)

Carlos Camino (Ported version into R)

## Examples

``` r
inputs = ToolsRTM::inputsRTMs
LUT<-as.data.frame(getLUT(inputs = inputs, nLUT=1, setseed = 1234))
sim <-getFluspect.Cx(inputsLeaf = LUT, inputsOptipar =ToolsRTM::optipar,version='Cx')
```
