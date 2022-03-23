#!/usr/bin/env Rscript
library(vsclust)
setwd(system.file("shiny",package="vsclust"))
shiny::runApp(port=3838)
