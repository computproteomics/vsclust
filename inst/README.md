# Overview of additional functionality in the vsclust package

## Command line script to run vsclust analysis
__Folder__: `inst/cmdline`

The two scripts allow direct analysis of a data matrix where the parameters and input files are given by a yaml file like `vsclust.yml` in the folder `inst/extdata`

Call it with `run_clust_app.R vsclust.yml`

## Extra files for testing and command line
__Folder__:  `inst/extdata`
Here you can find the yaml file for parametrizing the command line script and different csv files with data for testing. The file `ArtData.csv` contains artificial data where 50\% of the features are distributed in 5 clusters while the other 50\% of the features just consist of random noise (normally distributed).

## Older scripts and test script
__Folder__: `inst/other`

`runvsclust_algorithmR` can be used as test for vsclust (e.g. in the testing of the conda package)

`FcmClustPEst.R` contains a function to estimate cluster numbers according to the publication "Veit Schwämmle and Ole Nørregaard Jensen. A simple and fast method to determine
the parameters for fuzzy c-means cluster analysis. __Bioinformatics__, _26(22)_:2841–2848,
2010.

## Shiny app 
__Folder__: `inst/shiny`

Here you can find the files and code necessary to run the VSClust Shiny App. In RStudio, you can start it directly from either the `ui.R` or the `server.R` file.

