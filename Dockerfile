FROM rocker/shiny:4.1.2
LABEL maintainer="Veit Schwaemmle <veits@bmb.sdu.dk>"
LABEL description="Docker image of VSClust implementation on top of shiny-server. The number of to-be-installed R packages requires patience when building this image."

RUN apt-get update && apt-get install -y libssl-dev liblzma-dev libbz2-dev libicu-dev libxml2-dev openjdk-8-jdk tk tk-dev  libglpk-dev  && apt-get clean 

RUN R -e "install.packages('BiocManager', repos='http://cran.us.r-project.org'); \
  update.packages(ask=F); \
  BiocManager::install(c('BiocManager', 'devtools', 'RJDBC', 'dplyr', 'plotly', 'RPostgreSQL','rJava', 'lubridate', 'DT'),ask=F)"
RUN R CMD javareconf

RUN bash

RUN R -e "library(BiocManager); BiocManager::install(c('genefilter', 'Rcpp',  'clusterProfiler','qvalue','limma','matrixStats','yaml','shinyjs','shinythemes','graph', 'GOstats', 'Category', 'RBGL'),ask=F, force=T)"

RUN rm -rf /srv/shiny-server
RUN mkdir /srv/shiny-server
COPY inst/shiny/*  /srv/shiny-server/
COPY inst/shiny/www /srv/shiny-server/www
COPY inst/other/bioconductor-rdavidwebservice_1.28.0_src_all.tar.gz .
RUN tar -xzf bioconductor-rdavidwebservice_1.28.0_src_all.tar.gz
RUN R CMD INSTALL RDAVIDWebService

# installing customized library
COPY .  /srv/shiny-server/vsclust
RUN R CMD INSTALL /srv/shiny-server/vsclust
