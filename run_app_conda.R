#!/usr/bin/env Rscript
conda <- Sys.getenv("CONDA_PREFIX")
setwd(paste0(conda,"/share/vsclust"))

  shiny::runApp(port=3838)
