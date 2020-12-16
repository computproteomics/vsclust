FROM rocker/shiny
LABEL maintainer="Veit Schwaemmle <veits@bmb.sdu.dk>"
LABEL description="Docker image of VSClust implementation on top of shiny-server. The number of to-be-installed R packages requires patience when building this image."

      RUN apt-get update && apt-get install -y libssl-dev liblzma-dev libbz2-dev libicu-dev libxml2-dev openjdk-8-jdk tk tk-dev  libglpk-dev  && apt-get clean 

RUN R -e "install.packages('BiocManager', repos='http://cran.us.r-project.org'); \
  update.packages(ask=F); \
  BiocManager::install(c('BiocManager', 'devtools', 'RJDBC', 'dplyr', 'plotly', 'RPostgreSQL','rJava', 'lubridate', 'DT'),ask=F)"
RUN R CMD javareconf

RUN bash

  
RUN R -e "library(BiocManager); BiocManager::install(c('geneilter', 'clusterProfiler','qvalue','limma','matrixStats','yaml','shinyjs','shinythemes','RDAVIDWebService','Mfuzz'),ask=F)"




RUN rm -rf /srv/shiny-server
RUN mkdir /srv/shiny-server
COPY *R  /srv/shiny-server/
COPY *csv  /srv/shiny-server/
RUN mkdir /srv/shiny-server/www
COPY www/* /srv/shiny-server/www/

# installing customized library
COPY e1071FuzzVec_Installation/ /srv/shiny-server/e1071FuzzVec_Installation/
RUN chmod a+x /srv/shiny-server/e1071FuzzVec_Installation/configure
RUN R CMD INSTALL /srv/shiny-server/e1071FuzzVec_Installation
