### SCOPEinR

**SCOPEinR** is an R package designed for implementing the Soil Canopy
Observation, Photochemistry, and Energy Fluxes (SCOPE) radiative
transfer model. Originally developed in MATLAB, this model allows users
to simulate interactions between soil, canopy, and atmospheric processes
(Van der Tol et al., 2009; Yang et al., 2020).

The **SCOPEinR** package powers the online **RT-Simulator**, providing a
user-friendly interface through a Shiny app that allows users to
simulate canopy reflectance at the Top of Canopy (TOC) level. In
addition to reflectance, it also enables the simulation of chlorophyll
fluorescence emissions using the SCOPE model. For comprehensive
inter-comparison with other key radiative transfer (RT) models, it is
recommended to install the **ToolsRTM** package alongside **SCOPEinR**.

The **SCOPEinR** package executes the Soil Canopy Observation,
Photochemistry, and Energy Fluxes (SCOPE) radiative transfer model,
originally developed in MATLAB by Van der Tol et al. (2009) and further
extended by Yang et al. (2020). This integration allows users to
leverage the capabilities of the SCOPE model within the R environment
for in-depth analyses.

We have seamlessly integrated the **ToolsRTM** and **SCOPEinR** packages
into the online RT-platform, both of which are available under the GNU
license. You can access the online RT simulator at [this
link](https://carlos-camino.shinyapps.io/0-toolsrtm-simulator/).

The repositories for accessing the R packages and platform are as
follows:

- **ToolsRTM package:** <https://gitlab.com/caminoccg/toolsrtm>

- **SCOPEinR package:** <https://gitlab.com/caminoccg/scopeinr>

- **Online RT-platform:**
  <https://gitlab.com/caminoccg/toolsrtm-simulator>

*Note: These modules are currently in the testing phase before public
release. For access, please contact us at caminoccg@gmail.com.*

### Manuals

The manuals for both packages are available through the
[Shinyapp](https://carlos-camino.shinyapps.io/0-toolsrtm-simulator/) You
can also access the documentation directly within the
[ToolsRTM](https://carlos-camino.shinyapps.io/0-toolsrtm-simulator/_w_ef4421a7/Notebooks/R/ToolsRTM/ToolsRTM.html)
and
[SCOPEinR](https://carlos-camino.shinyapps.io/0-toolsrtm-simulator/_w_ef4421a7/Notebooks/R/SCOPEinR/SCOPEinR.html)
packages. See also the **Articles** on the [pkgdown
site](../docs/scopeinr/index.md) — “SCOPE course pipeline: energy
balance, fluorescence, and trait inversion”, including a real
SIF-vs-Vcmax25 experiment.

### Getting started

To install the **SCOPEinR** package, follow these steps in your R
session:

1.  Download the SCOPEinR package as .tar.gz file

&nbsp;

    ## install SCOPEinR
    install.packages('pathWithFile/scopeinr-main.tar.gz',repos = NULL,type = "source")

2.  Check the installed version:

&nbsp;

    # Check the version of ToolsRTM
    packageVersion("SCOPEinR")

2.  Optionally, install the ToolsRTM package for additional features:.

&nbsp;

    ## Install ToolsRTM package.
    install.packages('pathWithFile/toolsrtm-main.tar.gz',repos = NULL,type = "source")

### How to run SCOPE model in R

#### Get options for the SCOPE model

    ## read the setoptions
    table.with.opts<-read.table('input/setoptions.csv',header=T, sep=',')

#### Run the SCOPE model

    ## read the LUT table, here the default SCOPE's LUT ('input_data_default.csv') is used.
    inputLUT=read.table('input/LUT_input.csv',header=T,sep=',')

    ## Execute the SCOPE model, and store the resultant inputs and outputs in a list element.

    db.sim <-SCOPEinR::get.SCOPE(LUT=inputLUT,options.SCOPE=table.with.opts,
                                       optipar=SCOPEinR::optipar2021.Pro.CX,
                                       leaf.model='fluspect-CX',canopy.model='fourSAIL',
                                       get.outputs = 'ALL', get.plots = F)

Note: is needed to have the inputs in same folder as input.

### Runnable pipeline scripts

For the full workflow (simulate N SCOPE runs → convolve to
Sentinel-2/PRISMA → compute indices → invert traits, e.g. Cab, LAI, or
Vcmax25), see
[`Scripts/Pipeline/SCOPE-1-simulate.R`](https://gitlab.com/caminoccg/toolsrtm/-/tree/main/Scripts/Pipeline)
through `SCOPE-3-inversion_deep_learning.R` in the suite repo — the
SCOPEinR counterpart to ToolsRTM’s own pipeline scripts.

The main elements in db.sim are:

|               |               |                  |           |
|:-------------:|:-------------:|:----------------:|:---------:|
|   data.rad    |  data.fluxes  |    data.soil     | data.gap  |
| data.spectral | data.leafbio  |   data.angles    | data.bcu  |
| data.thermal  |  data.canopy  |    data.meteo    | data.bch  |
|   data.opts   | data.profiles | data.directional | iter.ebal |

#### Get SCOPE’s outputs

    ## get main outputs in a folder

    path.outs = 'outs/'

    ## get all the main outs
    get.SCOPE.outputs(data.sim = db.sim, N.sims=N.Samples,
                      LUT=inputLUT,  ## main inputs
                      path.out = path.outs, ## path for outputs
                      get.more.inputs=c('refl','lidf','LIDFb','Ft_Fo','rdo'),
                      get.plots=T)

The `get.SCOPE.outputs` function includes the `get.more.inputs`
parameter, which allows for the addition of extra inputs and outputs in
the output folder. In the previous example, we extracted the
reflectance, LIDF angles used, LIDF-b parameter, Ft-Fo ratio, and TOC
hemispherical-directional reflectance. These additional outputs will be
saved in the `Additional.inputs` folder.

![Fig. 1. Main structure for the output folder. The name of the folder
is taken from time system.](reference/figures/Outputs.png)

**Fig. 1**. Main structure for the output folder. The name of the folder
is taken from time system.

The `get.plots` function allows for the plotting of key outputs from the
SCOPE model, with the default setting set to `FALSE`. This option is
recommended for single simulations. When enabled, the `get.plots`
parameter will generate plots for reflectance, irradiance, radiance,
Sigma-fluorescence, and chlorophyll fluorescence emission.

![Fig. 2. Reflectance.](reference/figures/1-reflectance.png)

**Fig. 2**. Reflectance.

![Fig. 3. Irradiance.](reference/figures/2-Irradiance.png)

**Fig. 3**. Irradiance.

![Fig. 4. Radiance excluding adding
Fluorescence.](reference/figures/3-Radiance_excluding_addingFluorescence_2.png)

**Fig. 4**. Radiance excluding adding Fluorescence.

![Fig. 5. Hemispherically Integrated upwelling
radiance.](reference/figures/5-Hemispherically_Integrated_upwelling_radiance.png)

**Fig. 5**. Hemispherically Integrated upwelling radiance.

![Fig. 6. Chlorophyll
fluorescence.](reference/figures/6-fluorescence.png)

**Fig. 6**. Chlorophyll fluorescence.

### Run the SCOPE model in parallel

This script is designed for parallel simulations to reduce computing
time. The main concept behind parallelization is to divide the data into
smaller chunks and process them simultaneously across multiple cores or
processors. This approach enhances the utilization of available
computational resources, resulting in decreased overall simulation time.

**1 Initialization**

    n.rows <- nrow(LUT)
    print(n.rows)
    chunk_size <- n.rows /5

- `n.rows`: Total number of rows in the data frame `LUT`.

- `chunk_size`: Size of each chunk, calculated by dividing the total
  number of rows by 5 (for parallelization into 5 chunks).

**2 Parallel Simulation Loop**

    for (i in c(1:5)){
      print(i)
      start.row <- ((i - 1) * chunk_size) +1

      end.row <-  min(i * chunk_size, n.rows)

      LUT_ <- LUT[start.row:end.row,]
      # Get random indices to subset your LUT data
      #n.samples =100
      #n.randoms <- sample(1:nrow(LUT), n.samples)

      db.sims <-SCOPEinR::get.SCOPE.parallel(LUT=LUT_,
           options.SCOPE=table.with.opts,optipar=SCOPEinR::optipar2017.ProspectD,
           leaf.model='fluspect-CX',canopy.model='fourSAIL', parallel = T,
           get.outputs='ALL', get.plots = F, get.csv =T n.cores=8)
    }

- `for (i in c(1:5))`: Looping five times, each representing a chunk of
  data.

- `start.row` and `end.row`: Defining the starting and ending rows for
  each chunk.

- `LUT_`: Subset of data for the current iteration.

- [`SCOPEinR::get.SCOPE.parallel`](reference/get.SCOPE.parallel.md):
  Function for running simulations, but parallelization is disabled
  (`parallel = F`).

- `get.outputs='ALL'`, `get.plots = F`, `get.csv = T`: Specifying
  outputs required from the simulations.

**3 Save the outputs in same folder.**

The getCSV function merges all the outputs generated by the SCOPE model
in 5 folders. After generating the simulation, the getCSV will save the
SCOPE’s outputs in a uniqeu folder using the time system for the name.

    ##### Put the simulations in a
    sims <- getCSV(path.out='outs',n.folders=5,  files.names='All')

#### Run the SCOPE model in parallel with more cores

When working on a server with multiple CPU cores, parallelization can
greatly enhance processing speed and overall efficiency. Below is an
explanation of the script that highlights the advantages of parallel
processing in such environments:

    chunk_size <- nrow(LUT) #20,000 simulations
    db.sims <-SCOPEinR::get.SCOPE.parallel(LUT=LUT[1:chunk_size,],options.
         SCOPE=table.with.opts,optipar=SCOPEinR::optipar2017.ProspectD,
         leaf.model='fluspect-CX',canopy.model='fourSAIL', parallel = T,
         get.outputs='ALL', get.plots = F, get.csv =T, n.cores=16)
         

### 1.6 Get some additional plots by main plant trait.

    output.folder = 'path With outputs'
    plant.traits <- c('Vcmax25','EWT','Anth')

    ###  get.plots the options are 'fluorescence'; 'reflectance', and 'radiance'
    get.SCOPE.plots(path.files=output.folder, plant.trait=plant.traits, get.plots='fluorescence')

![Fig. 7. Fluorescence emission for 20 simulations classified by Anth
values.](reference/figures/6-fluorescence_fluorescence_by_Anth.png)

**Fig. 7.** Fluorescence emission for 20 simulations classified by Anth
values.

![Fig. 8. Reflectance (rdo) for 20 simulations classified by EWT
values.](reference/figures/1-reflectance_rdo_by_EWT.png)

**Fig. 8**. Reflectance (rdo) for 20 simulations classified by EWT
values.

![Fig. 9. Reflectance for 20 simulations classified by Vcmax
values.](reference/figures/1-reflectance_refl_by_Vcmax25.png)

**Fig. 9**. Reflectance for 20 simulations classified by Vcmax values.

**Note** Figure 8-10 showed also the effects form other plant traits

### References:

The official SCOPE’s github is available at
<https://github.com/Christiaanvandertol/SCOPE>

Yang, P., E. Prikaziuk, W. Verhoef, and C. van der Tol. 2020. “SCOPE
2.0: A Model to Simulate Vegetated Land Surface Fluxes and Satellite
Signals.” Geoscientific Model Development Discussions 2020: 1–26.
<https://doi.org/10.5194/gmd-2020-251>.

Van der Tol, C., W. Verhoef, J Timmermans, A Verhoef, and Z Su. 2009.
“An Integrated Model of Soil-Canopy Spectral Radiances, Photosynthesis,
Fluorescence, Temperature and Energy Balance.” Biogeosciences 6 (12):
3109–29. <https://doi.org/10.5194/bg-6-3109-2009>.

Other Main References:

G.James Collatz, J.Timothy Ball, Cyril Grivet, and Joseph A Berry.
Physiological and environmental regulation of stomatal conductance,
photosynthesis and transpiration: a model that includes a laminar
boundary layer. Agric. For. Meteorol., 54(2-4):107–136, apr 1991. URL:
<https://www.sciencedirect.com/science/article/pii/0168192391900028>,
<doi:10.1016/0168-1923(91)90002-8>.

GJ Collatz, M Ribas-Carbo, and JA Berry. Coupled Photosynthesis-Stomatal
Conductance Model for Leaves of C sub4/subPlants. Aust. J. Plant
Physiol., 19(5):519, 1992. URL:
<http://www.publish.csiro.au/?paper=PP9920519>, <doi:10.1071/PP9920519>.

Albert Porcar-Castell. A high-resolution portrait of the annual dynamics
of photochemical and non-photochemical quenching in needles of Pinus
sylvestris. Physiol. Plant., 143(2):139–153, oct 2011. URL:
<http://www.ncbi.nlm.nih.gov/pubmed/21615415>
<http://doi.wiley.com/10.1111/j.1399-3054.2011.01488.x>,
<doi:10.1111/j.1399-3054.2011.01488.x>.

G. Schaepman-Strub, M. E. Schaepman, T. H. Painter, S. Dangel, and J. V.
Martonchik. Reflectance quantities in optical remote sensing-definitions
and case studies. Remote Sens. Environ., 103(1):27–42, 2006.
<doi:10.1016/j.rse.2006.03.002>.

Christiaan van der Tol, Micol Rossini, Sergio Cogliati, Wouter Verhoef,
Roberto Colombo, Uwe Rascher, and Gina Mohammed. A model and measurement
comparison of diurnal cycles of sun-induced chlorophyll fluorescence of
crops. Remote Sens. Environ., 186:663–677, dec 2016. URL:
<https://www.sciencedirect.com/science/article/pii/S0034425716303649>,
<doi:10.1016/j.rse.2016.09.021>.

Wout. Verhoef and Nationaal Lucht- en Ruimtevaartlaboratorium
(Netherlands). Theory of radiative transfer models applied in optical
remote sensing of vegetation canopies. \[publisher not identified\],
1998. ISBN 9054858044. URL:
<https://library.wur.nl/WebQuery/wda/945481>.

Wouter Verhoef, Christiaan van der Tol, and Elizabeth M. Middleton.
Hyperspectral radiative transfer modeling to explore the combined
retrieval of biophysical parameters and canopy fluorescence from FLEX –
Sentinel-3 tandem mission multi-sensor data. Remote Sens. Environ.,
204(August 2016):942–963, 2018. URL:
<https://doi.org/10.1016/j.rse.2017.08.006>,
<doi:10.1016/j.rse.2017.08.006>.

Nastassia Vilfan, Christiaan van der Tol, Onno Muller, Uwe Rascher, and
Wouter Verhoef. Fluspect-B: A model for leaf fluorescence, reflectance
and transmittance spectra. Remote Sens. Environ., 186:596–615, 2016.
URL: <http://dx.doi.org/10.1016/j.rse.2016.09.017>,
<doi:10.1016/j.rse.2016.09.017>.

Peiqi Yang, Wout Verhoef, and Christiaan van der Tol. The mSCOPE model:
A simple adaptation to the SCOPE model to describe reflectance,
fluorescence and photosynthesis of vertically heterogeneous canopies.
Remote Sens. Environ., 201:1–11, nov 2017. URL:
<https://www.sciencedirect.com/science/article/pii/S0034425717303954>,
<doi:10.1016/j.rse.2017.08.029>.

XINYOU YIN, JEREMY HARBINSON, and PAUL C. STRUIK. Mathematical review of
literature to assess alternative electron transports and
interphotosystem excitation partitioning of steady-state C3
photosynthesis under limiting light. Plant, Cell Environ.,
29(9):1771–1782, sep 2006. URL:
<http://doi.wiley.com/10.1111/j.1365-3040.2006.01554.x>,
<doi:10.1111/j.1365-3040.2006.01554.x>.

Xinyou Yin and Paul C. Struik. Crop systems biology as an avenue to
bridge applied crop science and fundamental plant biology. In Proc. -
2012 IEEE 4th Int. Symp. Plant Growth Model. Simulation, Vis. Appl. PMA
2012, 15–17. IEEE, oct 2012. URL:
<http://ieeexplore.ieee.org/document/6524806/>,
<doi:10.1109/PMA.2012.6524806>.

Van der Tol, C.V, Berry J. A., Campbell P.K.E., and Rascher U. Models of
fluorescence and photosynthesis for interpreting measurements of
solar-induced chlorophyll fluorescence. J. Geophys. Res. Biogeosciences,
119(12):2312–2327, 2014.

### License

![](https://img.shields.io/badge/License-MIT-yellow.svg)

The **SCOPEinR** package is licensed under the MIT License, allowing for
free use, modification, and distribution. This package is available on
GitLab, and we encourage contributions and collaborations from the
community. For more details, please refer to the LICENSE file in the
repository.
