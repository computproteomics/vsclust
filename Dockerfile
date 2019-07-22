FROM opensuse/leap:15.1
LABEL maintainer="Veit Schwaemmle <veits@bmb.sdu.dk>"
LABEL description="Docker image of VSClust implementation on top of shiny-server and OpenSuse tumblewee. The number of to-be-installed R packages requires patience when building this image."

RUN zypper --non-interactive  update

RUN zypper --non-interactive install java-1_8_0-openjdk java-1_8_0-openjdk-devel tar wget

# R devtools pre-requisites:
RUN zypper --non-interactive install  git libxml2 libxml2-devel curl curl-devel openssl-devel pandoc postgresql


RUN zypper --non-interactive install  libcurl-devel R-base R-base-devel

RUN zypper --non-interactive install gcc gcc-c++ curl pcre-devel libbz2-1 libbz2-devel pcre icu libicu-devel gcc-fortran

RUN R CMD javareconf

RUN bash

RUN R -e "install.packages('rJava', repos='http://cran.r-project.org', INSTALL_opts='--no-html')"

RUN R -e "install.packages(c( 'BiocManager', 'shiny', 'rmarkdown', 'devtools', 'RJDBC', 'dplyr', 'plotly', 'RPostgreSQL', 'lubridate', 'DT'), repos='http://cran.r-project.org', INSTALL_opts='--no-html')"

RUN R -e "library(BiocManager); install(); install(c('geneFilter', 'clusterProfiler','qvalue','limma','matrixStats','shinyjs','shinythemes','RDAVIDWebService','Mfuzz'))"
RUN wget https://download3.rstudio.org/centos6.3/x86_64/shiny-server-1.5.9.923-x86_64.rpm 

#RUN zypper addrepo http://download.opensuse.org/repositories/server:monitoring/openSUSE_Tumbleweed/server:monitoring.repo 
#RUN zypper --non-interactive --no-gpg-checks refresh
RUN zypper --non-interactive install libffi-devel libffi 

RUN wget http://download.opensuse.org/repositories/science/openSUSE_Leap_15.1/x86_64/udunits2-2.2.26-lp151.1.1.x86_64.rpm
RUN wget http://download.opensuse.org/repositories/science/openSUSE_Leap_15.1/x86_64/udunits2-devel-2.2.26-lp151.1.1.x86_64.rpm
RUN rpm -ivh --nodeps  shiny-server*.rpm udunits2*.rpm

RUN R -e "library(BiocManager); install(); install(c('genefilter', 'clusterProfiler','qvalue','limma','matrixStats','shinyjs','shinythemes','RDAVIDWebService','Mfuzz'))"

RUN mkdir -p /var/log/shiny-server
RUN mkdir -p /home/shiny
RUN chown shiny.shiny /var/log/shiny-server
RUN chown -R shiny:shiny /srv/shiny-server
RUN chown -R shiny:shiny /var/lib/shiny-server

RUN mkdir /srv/shiny-server/VSClust

ADD  *.R /srv/shiny-server/VSClust/
ADD www /srv/shiny-server/VSClust/
ADD ArtData.csv /srv/shiny-server/VSClust/
RUN mkdir /srv/shiny-server/VSClust/e1071FuzzVec_Installation
COPY e1071FuzzVec_Installation/ /srv/shiny-server/VSClust/e1071FuzzVec_Installation/


RUN chmod a+x /srv/shiny-server/VSClust/e1071FuzzVec_Installation/configure

RUN R CMD INSTALL /srv/shiny-server/VSClust/e1071FuzzVec_Installation


RUN R CMD javareconf

EXPOSE 8787 3838


USER shiny
CMD shiny-server




