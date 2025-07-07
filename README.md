# Variance-sensitive clustering of omics data

Feature-based variance-sensitive clustering of omics data. Optimizes cluster assignment by taking into account individual feature variance.

VSClust is available as interactive and user-friendly Shiny webservice and as R package in Bioconductor. For more instructions on both, see below.

VSClust was developed at the

[Protein Research Group](http://www.sdu.dk/en/Om_SDU/Institutter_centre/Bmb_biokemi_og_molekylaer_biologi/Forskning/Forskningsgrupper/Protein.aspx)  
Department of Biochemistry and Molecular Biology  
[University of Southern Denmark](http://www.sdu.dk)  

## Citation
When using VSClust, please cite our paper:  
Veit Schwämmle, Ole N Jensen; VSClust: Feature-based variance-sensitive clustering of omics data, Bioinformatics, bty224, https://doi.org/10.1093/bioinformatics/bty224

## Shiny app

### Web service

You can use the implementation on our web server http://computproteomics.bmb.sdu.dk:  
http://computproteomics.bmb.sdu.dk/Apps/VSClust

Be aware that the tool does allow only one user to run the background R calculations at a time. Therefore the app might become temporarily irresponsive. However, multiple sessions are separated and your data won't be shared between sessions or overwritten. 

### Local implementation

#### Docker
The easiest option is to use the docker image:

```
docker pull veitveit/vsclust
docker run -t -i -p 3838:3838 veitveit/vsclust
```

and access the server through http://localhost:3838

#### Bioconda
Install the package in the command-line
```
conda install -c bioconda vsclust
run_vsclust_app.sh
```
and access the shiny app through http://localhost:3838

#### Manual installation
Install the Bioconductor package `vsclust`
In R:
```
if (!require("BiocManager", quietly = TRUE))
    install.packages("BiocManager")

BiocManager::install("vsclust")
```

You can run the shiny app from the server.R or ui.R files in the "inst/shiny" folder using [Rstudio](http://rstudio.com) or run the app on a [shiny-server](https://www.rstudio.com/products/shiny/shiny-server/)

Be aware that you need to have all files, the R libraries described in Installation *and* the modified e1071 library installed.


### Build and use Docker image
A Dockerfile has been created on top of the rocker/shiny docker image. Copy this repository to a folder and carry out the following command to build the image (takes a while)

`docker build -t veitveit/vsclust .`

You can also just directly download and run the image by

`docker run -t -i -p 3838:3838 veitveit/vsclust`

and access the server through http://localhost:3838

## Command line and R package

All operations but the gene set enrichment can be performed via command line running the R script `runVSClust.R`
or using the functions of the Bioconductor package `vsclust`

### Installation

#### In R: 
```
if (!require("BiocManager", quietly = TRUE))
    install.packages("BiocManager")

BiocManager::install("vsclust")
```

### As conda package
TBD

### Usage

Please take a look at the vignettes and/or in the help packages of the `vsclust` functions

## Contact
For software issues and general questions, please submit an issue.

## License
GPL-2 or higher
