#!/usr/bin/env Rscript

currPath <- getwd()
  initial.options <- commandArgs(trailingOnly = FALSE)
  file.arg.name <- "--file="
  script.name <- sub(file.arg.name, "", initial.options[grep(file.arg.name, initial.options)])
  script.basename <- dirname(script.name)
  cat(paste0("Getting R functions from files located in ", script.basename),"\n")
  setwd(script.basename)
  shiny::runApp(port=3838)
setwd(currPath)


