**Welcome to the repository of VSClust**
developed at the

[Protein Research Group](http://www.sdu.dk/en/Om_SDU/Institutter_centre/Bmb_biokemi_og_molekylaer_biologi/Forskning/Forskningsgrupper/Protein.aspx)  
Department of Biochemistry and Molecular Biology  
[University of Southern Denmark](http://www.sdu.dk)  


We provide a shiny app for interactive analysis, a command-line version for running VSClust as script in R and a docker version of the tool to avoid installation issues.

## Shiny app

### Web service

You can use our web server http://computproteomics.bmb.sdu.dk:

http://computproteomics.bmb.sdu.dk/Apps/VSClust

Be aware that the tool does allow only one user to run the background R calculations at a time. Therefore the app might become temporarily irresponsive. However, multiple sessions are separated and your data won't be shared between sessions or overwritten. 

### Implementation on own computer

The easiest option is to use the docker image:

*docker pull veitveit/vsclust

*docker run -t -i -p 3838:3838 veitveit/vsclust

and access the server through http://localhost:3838/VSClust/


You can run the shiny app from the server.R or ui.R files using [Rstudio](http://rstudio.com) or run the app on a [shiny-server](https://www.rstudio.com/products/shiny/shiny-server/)

Be aware that you need to have all files, the R libraries described in Installation *and* the modified e1071 library installed.


### Installation
Install the following R libraries: Mfuzz, matrixStats, limma, qvalue, shiny, clusterProfiler, RDAVIDWebService, parallel, shinyjs and shinythemes.
In R:
*source("https://bioconductor.org/biocLite.R")*   
*biocLite(c("Mfuzz", "matrixStats", "limma", "qvalue", "shiny", "clusterProfiler", "RDAVIDWebService", "parallel", "shinyjs", "shinythemes"))*

Download the files into a folder and install the library *e1071FuzzVec*. You might need to compile the library on your computer:  
Run *R CMD INSTALL e1071FuzzVec_Installation* from the command line. You need to be in the VSClust folder. For Windows users, replace *R* by *InstallationPath/R.exe*.

### Build and use Docker image
A Dockerfile has been created on the basis of an OpenSuse distribution. Copy the repository to a folder and carry out the following command to build the images (takes a while)

*docker build -t veitveit/vsclust .

You can run the image by

*docker run -t -i -p 3838:3838 veitveit/vsclust

and access the server through http://localhost:3838/VSClust/

## Command line 

All operations but the gene set enrichment can be performed via command line running the R script runVSClust.R

### Usage

Given that you installed the *e1071FuzzVec* library, open the R script and change the relevant file names and parameters. The different operations are described in the R file. 

### Installation
The *e1071FuzzVec* library needs to be compiled and installed, see above.

## Contact
For software issues and general questions, please submit an issue.

## License
GPL-2 or higher
