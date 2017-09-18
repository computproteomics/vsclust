## TODO Clean accession numbers better and separate ()?
library(matrixStats)
library(Mfuzz)
library(limma)
library(parallel)
library(qvalue)
library(e1071FuzzVec)
library(shinyjs)
library(clusterProfiler)
library(RDAVIDWebService)
source("FcmClustPEst.R")
source("mfuzz.plotpdf.R")
source("HelperFuncs.R")

options(shiny.maxRequestSize=2000*1024^2,shiny.trace=T,shiny.sanitize.errors=F) 

options(java.parameters="-Xss2560k")

shinyServer(function(input, output,clientData,session) {
  
  dat <- NULL
  pars <- NULL
  pars$m <- NULL
  maxClust <- 25
  NumCond <- NULL
  cores <- 2
  
  # # Resetting views
  # reset_memory <- 0
  # observeEvent(input$reset,{
  # if (input$reset - reset_memory > 0) {
  #   print("Hello")
  #       pars <<- NULL
  #   dat <<- NULL
  #   reset_memory <<- input$reset
  # }
  #   
  # })
  # 
  # setting the right tab
  observe({
    input$NumReps
    input$NumCond
    input$isPaired
    input$isStat
    input$row.names
    input$infile
    input$protnames
    input$examplefile
    input$is_header
    if (!is.null(v$dat) & is.numeric(input$NumReps) & is.numeric(input$NumCond) & input$isStat) {
      if (is.numeric(ifelse(input$protnames,as.matrix(v$dat[,2:ncol(v$dat)]),as.matrix(v$dat)))) {
        if (ifelse(input$protnames,ncol(v$dat)-1,ncol(v$dat))==input$NumReps*input$NumCond) {
          print("go to stats")
          updateTabsetPanel(session, inputId ="tabset", selected = "stat")
        }
      }
      maxSamp <- ifelse(input$protnames,ncol(v$dat)-1,ncol(v$dat))
      updateNumericInput(session,"NumCond",max=maxSamp)
      updateNumericInput(session,"NumReps",max=maxSamp)
      output$RepsCond <- renderText(paste("Total number of samples (replicates * conditions):", maxSamp))
    } else {
      if (!input$isStat & !is.null(v$dat)) {
        updateTabsetPanel(session, inputId ="tabset", selected = "stat")
      } else  {
        updateTabsetPanel(session, inputId ="tabset", selected = "fin") 
      }
    }
  })
  
  output$tparame <- renderText("Parameter estimation")
  output$tresults <- renderText("Clustering results")
  output$tdescr <- renderText("Description")
  output$ui <- renderUI({
    if (is.null(input$isStat) | input$isStat) {
      p(
        checkboxInput(inputId="isPaired", label="Paired tests",value=T),
        textOutput("RepsCond"),
        numericInput("NumReps",min=2,max=20,value=2,label="Number of replicates",step=1),
        numericInput("NumCond",min=2,max=20,value=3,label="Number of conditions",step=1))
    }
  })
  v <- reactiveValues(dat = NULL, example = F, clustOut = NULL)
  observe({
    input$reset
    session$sendCustomMessage(type = "resetFileInputHandler", "in_file")  
    v <- NULL
    reset("in_file")
  })
  observeEvent(input$examplefile,{
    dat <- read.csv("ArtData.csv",row.names=1)
    dat <- dat[1:nrow(dat),]
    print("example")
    updateNumericInput(session,"NumCond",value=10,max=10)
    updateNumericInput(session,"NumReps",value=5,max=5)
    updateCheckboxInput(session,"isPaired",value=F)
    updateCheckboxInput(session,"isStat",value=T)
    v$dat <- dat
    v$example <- T
  })
  observeEvent(input$in_file,{
    ## test for right replicate and condition numbers, min 50 features, ...
    dat <- NULL
    withProgress(message="Reading file ...", min=0,max=1, value=0.5,  {
      try(dat <- read.csv(input$in_file$datapath,row.names=1,header=input$is_header))
      dat <- dat[rownames(dat)!="",]
      v$example <- F
      v$dat <- dat
      updateNumericInput(session,"NumCond",max=ifelse(input$protnames,ncol(dat)-1,ncol(dat)))
      updateNumericInput(session,"NumReps",max=ifelse(input$protnames,ncol(dat)-1,ncol(dat)))
    })
  })
  output$plot0 <- renderPlot({
    
    proteins <- NULL
    dat <- v$dat
    if(input$protnames & !v$example) {
      proteins <- dat[,1]
      dat <- dat[,2:ncol(dat)]
      names(proteins) <- rownames(dat)
    }
    
    print(dim(dat))
    print(input$in_file)
    validate(need(!is.null(dat),"no data"))
    validate(need(is.numeric(as.matrix(dat)), "The data table contains non-numerical  values!"))
    if(!is.null(dat)) {
      withProgress(message="Statistics + PCA ...", min=0,max=1, value=0.5,  {
        NumReps <- input$NumReps
        NumCond <- input$NumCond
        print(NumReps)
        print(NumCond)
        dat[!is.finite(as.matrix(dat))] <- NA
        num_miss <- sum(is.na(dat))
        if (input$isStat) {
          validate(need(ncol(dat)==NumReps*NumCond, "Number of data columns must correspond to product of conditions and replicates!"))
        }
        validate(need(try(statOut <- statWrapper(dat, NumReps, NumCond, input$isPaired, input$isStat)), 
                      "Please remove the following items from your input file:\na) empty columns or rows\nb) non-numerical or infinite values\nc) commenting characters (e.g. #)"))
        
        dat <- statOut$dat
        Sds <- dat[,ncol(dat)]
        output$data_summ <- renderUI({HTML(paste("Features:",nrow(dat),"<br/>Missing values:",
                                                 num_miss,"<br/>Median standard deviations:",
                                                 round(median(Sds,na.rm=T),digits=3)))})
        
        pars$dat <<- dat 
        pars$proteins <<- proteins
        
        # file to download q-values
        output$downloadDataLimma <- downloadHandler(
          filename = function() {
            paste("LimmaResults", Sys.Date(), ".csv", sep="");
          },
          content = function(file) {
            write.csv(statOut$statFileOut, file)
          })
      })
    }
  })
  
  output$plot1 <- renderPlot({
    if(input$clButton1 == 0) 
      return()
    isolate({
      if((input$clButton1)) {
        dat <- pars$dat
        clustNumOut <- NULL
        withProgress(message="Calculating ...", min=3,max=maxClust, value=2,  {
          clustNumOut <- estimClustNum(dat, input$maxclust, cores)
        })
        
        output$downloadParamEst <- downloadHandler(
          filename = function() {
            paste("EstimatedClustNumber", Sys.Date(), ".pdf", sep="");
          },
          content = function(file) {
            pdf(file,height=6,width=15)
            print(clustNumOut$p)
            dev.off()  
          })
      }
    })},height=600)
  
  
  output$plot2 <- renderPlot({
    if(input$clButton2 == 0) 
      return()
    isolate({
      if((input$clButton2)) {
        if (!is.null(pars$dat)) {
          withProgress(message="Calculating ...", min=0,max=1, {
            print(input$nclust1)
            Nclust <- input$nclust1
            # pars$dat <- pars$dat[rowSums(is.na(pars$dat))==0,]
            ClustOut <- runClustWrapper(pars$dat, Nclust, pars$proteins, VSClust=T, cores)
            output$downloadFigure <- downloadHandler(
              filename = function() {
                paste("FCMVarMResults", Sys.Date(), ".pdf", sep="");
              },
              content = function(file) {
                pdf(file,height=5*round(sqrt(Nclust)),width=5*ceiling(sqrt(Nclust)))
                replayPlot(ClustOut$p)
                dev.off()
              })
            
            Bestcl <- ClustOut$Bestcl
            pars$Bestcl2 <<- Bestcl  
            
            output$downloadData2 <- downloadHandler(
              filename = function() {
                paste("FCMVarMResults", Sys.Date(), ".csv", sep="");
              },
              content = function(file) {
                write.csv(data.frame(cluster=Bestcl$cluster,ClustOut$outFileClust,isClusterMember=rowMaxs(Bestcl$membership)>0.5,maxMembership=rowMaxs(Bestcl$membership),
                                     Bestcl$membership), file)
              })
            output$downloadCentroid2 <- downloadHandler(
              filename = function() {
                paste("FCMVarMResultsCentroids", Sys.Date(), ".csv", sep="");
              },
              content = function(file) {
                write.csv(Bestcl$centers, file)
              })
            output$clustinf1 <- renderDataTable(as.data.frame(ClustOut$ClustInd))
          })
        }}
    })
  })
  
  
  output$plot3 <- renderPlot({
    if(input$clButton3 == 0) 
      return()
    isolate({
      if((input$clButton3)) {
        if (!is.null(pars$dat)) {
          withProgress(message="Calculating ...", min=0,max=1, {
            Nclust <- input$nclust2
            # pars$dat <- pars$dat[rowSums(is.na(pars$dat))==0,]
            
            ClustOut <- runClustWrapper(pars$dat, Nclust, pars$proteins, VSClust=F, cores)
            output$downloadFigure2 <- downloadHandler(
              filename = function() {
                paste("FCMResults", Sys.Date(), ".pdf", sep="");
              },
              content = function(file) {
                pdf(file,height=5*round(sqrt(Nclust)),width=5*ceiling(sqrt(Nclust)))
                replayPlot(ClustOut$p)
                dev.off()
              })
            
            Bestcl <- ClustOut$Bestcl
            pars$Bestcl1 <<- Bestcl  
            
            output$downloadData3 <- downloadHandler(
              filename = function() {
                paste("FCMResults", Sys.Date(), ".csv", sep="");
              },
              content = function(file) {
                write.csv(data.frame(cluster=Bestcl$cluster,ClustOut$outFileClust,isClusterMember=rowMaxs(Bestcl$membership)>0.5,maxMembership=rowMaxs(Bestcl$membership),
                                     Bestcl$membership), file)
              })
            output$downloadCentroid3 <- downloadHandler(
              filename = function() {
                paste("FCMResultsCentroids", Sys.Date(), ".csv", sep="");
              },
              content = function(file) {
                write.csv(Bestcl$centers, file)
              })
            output$clustinf2 <- renderDataTable(as.data.frame(ClustOut$ClustInd))
          })
        }}
    })
  })
  
  output$plot4 <- renderPlot({
    if(input$goButton == 0)
      return()
    validate(need(!is.null(pars$Bestcl2),"Run clustering first"))
    if (!is.null(pars$Bestcl2)) {
      isolate({
        if(input$goButton) {
          Accs <- list()
          cl <- pars$Bestcl2
          dat <- pars$dat[,1:(ncol(pars$dat)-1)]
          for (c in 1:max(cl$cluster)) {
            Accs[[c]] <- names(which(cl$cluster==c & rowMaxs(cl$membership)>0.5))
            Accs[[c]] <- Accs[[c]][Accs[[c]]!=""]
            if (length(Accs[[c]])>0) {
              if (length(Accs[[c]])>1) {
                tdat <- (dat[Accs[[c]],])
              } else {
                tdat <- t(dat[Accs[[c]],])
              }
              Accs[[c]] <- sub("-[0-9]","",Accs[[c]])
            } else {
              tdat <- t(rep(NA,ncol(dat)+1))
            }
            if (!is.null(pars$proteins)) {
              Accs[[c]] <- as.character(pars$proteins[Accs[[c]]])
            }
          }
          names(Accs) <- paste("Cluster",1:length(Accs))
          Accs <- lapply(Accs,function(x) unique(ifelse(is.na(x),"B3",x)))
          Accs <- Accs[lapply(Accs,length)>0]
          print(lapply(Accs,length))
          withProgress(message="Waiting for data (1/2)...", min=0,max=1, {
            x <- NULL
            try(x <- compareCluster(Accs, fun="enrichDAVID", annotation=input$infosource,
                                idType=input$idtype,
                                listType="Gene", david.user = "veits@bmb.sdu.dk"))
            validate(need(!is.null(x),"No result. Wrong ID type?"))
            incProgress(0.7, detail = "received")
            x@compareClusterResult <- cbind(x@compareClusterResult,log10padval=log10(x@compareClusterResult$p.adjust))
            incProgress(0.8, detail = "plotting")
            validate(need(nrow(x@compareClusterResult)>1,"No significant results"))
            output$downloadGOData1 <- downloadHandler(
              filename = function() {
                paste("DAVIDResultsVarM", Sys.Date(), ".csv", sep="");
              },
              content = function(file) {
                write.csv(as.data.frame(x@compareClusterResult), file)
              })
            incProgress(0.9, detail = "calculating BHI")
            BHI <- calcBHI(Accs,x)
              y <- new("compareClusterResult",compareClusterResult=x@compareClusterResult)
              if (length(unique(y@compareClusterResult$ID)) > 20) {
                print("Reducing number of DAVID results")
                y@compareClusterResult <- y@compareClusterResult[
                  order(y@compareClusterResult$p.adjust)[1:20],]
                # print(x@compareClusterResult)
              }
              plot(y,title=paste("BHI:",BHI),showCategory=1000,colorBy="log10padval",font.size=10)
          })
        }})}})
  
  output$plot5 <- renderPlot({
    if(input$goButton == 0)
      return()
    validate(need(!is.null(pars$Bestcl1),"Run clustering first"))
    if (!is.null(pars$Bestcl1)) {
      isolate({
        if(input$goButton) {
          Accs <- list()
          cl <- pars$Bestcl1
          dat <- pars$dat[,1:(ncol(pars$dat)-1)]
          for (c in 1:max(cl$cluster)) {
            Accs[[c]] <- names(which(cl$cluster==c & rowMaxs(cl$membership)>0.5))
            Accs[[c]] <- Accs[[c]][Accs[[c]]!=""]
            if (length(Accs[[c]])>0) {
              if (length(Accs[[c]])>1) {
                tdat <- (dat[Accs[[c]],])
              } else {
                tdat <- t(dat[Accs[[c]],])
              }
              Accs[[c]] <- sub("-[0-9]","",Accs[[c]])
            } else {
              tdat <- t(rep(NA,ncol(dat)+1))
            }
            if (!is.null(pars$proteins)) {
              Accs[[c]] <- as.character(pars$proteins[Accs[[c]]])
            }
          }
          names(Accs) <- paste("Cluster",1:length(Accs))
          Accs <- lapply(Accs,function(x) unique(ifelse(is.na(x),"B3",x)))
          Accs <- Accs[lapply(Accs,length)>0]
          print(lapply(Accs,length))
          withProgress(message="Waiting for data (2/2) ...", min=0,max=1, {
            x <- NULL
            try(x <- compareCluster(Accs, fun="enrichDAVID", annotation=input$infosource,
                                    idType=input$idtype,
                                    listType="Gene", david.user = "veits@bmb.sdu.dk"))
            validate(need(!is.null(x),"No result. Wrong ID type?"))
            incProgress(0.7, detail = "received")
            x@compareClusterResult <- cbind(x@compareClusterResult,log10padval=log10(x@compareClusterResult$p.adjust))
            incProgress(0.8, detail = "plotting")
            validate(need(nrow(x@compareClusterResult)>1,"No significant results"))
            output$downloadGOData2 <- downloadHandler(
              filename = function() {
                paste("DAVIDResultsStand", Sys.Date(), ".csv", sep="");
              },
              content = function(file) {
                write.csv(as.data.frame(x@compareClusterResult), file)
              }) 
            incProgress(0.9, detail = "calculating BHI")
            BHI <- calcBHI(Accs,x)
              input$goButton
              y <- new("compareClusterResult",compareClusterResult=x@compareClusterResult)
              if (length(unique(y@compareClusterResult$ID)) > 20) {
                print("Reducing number of DAVID results")
                y@compareClusterResult <- y@compareClusterResult[
                  order(y@compareClusterResult$p.adjust)[1:20],]
              }
              plot(y,title=paste("BHI:",BHI),showCategory=1000,colorBy="log10padval",font.size=10)
          })
        }})}
  })
  
  output$intro <- renderUI({HTML("We present a new method to apply fuzzy c-means clustering to data sets that exhibit 
                                  non-constant variance of their features. The individual variation is incorporated into 
the estimation of the fuzzifier (parameter m). For details see<br/> Schwämmle, V. and Jensen, O. N., &quotVariance sensitive fuzzy clustering&quot (<i>in preparation</i>)<br/> 
                       The algorithm for the estimation of the parameters fuzzifier and cluster numbers is furthermore based on
<br/>Schwämmle, V. and <a href='http://www.ncbi.nlm.nih.gov/pubmed/20880957'>Jensen, O. N. &quotA simple and fast method to determine the parameters for fuzzy c-means cluster analysis&quot. <i>Bioinformatics</i>, 2010, <b>26</b>, 2841-2848<br/></a>
                                 Please note that this method reveals its power for 8 or more different conditions (dimensions). Lower numbers yields results nearly identical to standard fuzzy c-means clustering.<br/>
                                 <b><i>Example data set:</i></b> You can test the method with an artificial data set by clicking on <i>Load example</i>. The data set contains 500 features where half of them 
                                 make part of five predefined clusters while the rest of the data is randomly distributed.")})
  output$finput <- renderUI({HTML("Input format is restricted to comma-delimited files (.csv). Files are required to contain only the numerical data that will
be analyzed in addition to the following contents: A one-row header containing the 
                                  column names and additionally defined feature names for the rows (e.g. gene names, transcripts, peptides and proteins). The latter
                                  need to be unique, duplicated feature names are not accepted. Additionally, in the case of data containing multiple measurements of the same gene/protein, a second column can contain this information
that then will be used for the GO term enrichment.<br/>
                                  The order of the columns is crucial to distinguish conditions from replicated values. With numbers denoting conditions (1-4) and letters
corresponding to replicates (A-C), the columns are arranged: A1, A2, A3, A4, B1, B2, B3, B4, C1, C2, C3, C4. <br/>
                                  In the case of an input file that already contains estimated standard deviations, all columns but the last ones are considered 
                                  different conditions while the last column contains the standard deviations.")})
  output$stat <- renderUI({HTML("We offer to calculate feature standard deviations using the <a href='https://bioconductor.org/packages/limma'>limma package</a>, well-performing for microarray, mass spectrometry and RNA-seq data. 
                                The method corrects feature variation by ensuring a minimal, adapted variation that is estimated from the entire data set.  
                                Statistics and estimation of feature variation can be carried out using paired and unpaired tests. <br/>
                                The plot on the &quotStatistics and variance&quot tab shows a simplified visualization of the data projected on the first principal components (standard PCA). Point sizes correspond to 
                                estimated feature standard deviation.")})
  output$pest <- renderUI({HTML("Estimation of an optimal cluster number can become tricky and many validation methods have been proposed
                                to pick out the right number. The tab &quotEstimation of cluster number&quot shows the values of the minimum 
                                centroid distance and the Xie-Beni index over the user-defined range of cluster numbers. Suitable cluster numbers are given by large jumps
                                (minimal centroid distance) and minima (Xie-Beni index). The best estimates for the variance-sensitive method are marked by black squares. Additionally, the 3rd figures depticts the number of assigned features (features with membership values above 0.5)")})
  output$fclust <- renderUI({HTML("Clustering results are available in 2 tabs, showing the results for our variance-sensitive algorithm and 
                                  standard c-means clustering. Users are required to set the number of clusters. Expression profiles are shown in separate panels where the color corresponds to the membership values.
                                  The plot shows only features confidently assigned to a cluster requiring a minimum membership value of 0.5. Users can download 
                                  pdf-files of the figure, data with membership values and principal expression profiles of the clusters (cluster centroids). Furthermore, the number
                                  of cluster members is shown in a table. This can be useful when comparing the two methods, where the variance-based algorithm should yield 
                                  lower numbers of cluster members as features with large standard deviations often get discarded. The clustering is implemented using <a href='https://bioconductor.org/packages/release/bioc/html/Mfuzz.html'>Mfuzz</a> and a modified version of <a href='https://cran.r-project.org/web/packages/e1071/index.html'>e1071</a>.")})
  output$goterms <- renderUI({HTML("Simple viewer for GO term and pathway enrichment based on <a  target='_blank'' href='https://bioconductor.org/packages/release/bioc/html/clusterProfiler.html'>clusterProfiler</a> 
                                    configured to access the <a href='https://david.ncifcrf.gov/content.jsp?file=WS.html'>DAVID</a> web service. <i>This can take a while (response times of DAVID server)</i>. We furthermore offer a rough estimate of enrichment efficiency, 
                                   given by the biological homogeneity index (BHI, see also <a href='http://www.ncbi.nlm.nih.gov/pmc/articles/PMC1590054/'>paper</a>)<br/>
#<b>Important: </b>In the case of more than 2000 features in cluster, they will be limited to the 2000 with the highest membership values.
    ")})
  output$reading <- renderUI({HTML("The source code of VSClust is available at <a href='https://bitbucket.org/veitveit/vsclust'>bitbucket</a>. VSClust was developed and implemented 
                                   at the <a href='http://www.sdu.dk/en/Om_SDU/Institutter_centre/Bmb_biokemi_og_molekylaer_biologi/Forskning/Forskningsgrupper/Protein.aspx'>Protein Research Group</a> of the University of Southern Denmark. See also <a href='computproteomics.bmb.sdu.dk'>computproteomics.bmb.sdu.dk</a> for more information.")})
  
})
