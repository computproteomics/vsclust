## TODO Clean accession numbers better and separate ()?
library(vsclust)
library(shinyjs)
library(clusterProfiler)
library(jsonlite)
library(matrixStats)
library(DT)

validate <- shiny::validate

options(shiny.maxRequestSize=2000*1024^2,shiny.sanitize.errors=F) 

options(java.parameters="-Xss2560k")
shinyServer(function(input, output,clientData,session) {
  
  dat <- NULL
  pars <- NULL
  pars$m <- NULL
  maxClust <- 25
  cores <- 4
  shiny_threads <- as.numeric(Sys.getenv("SHINY_THREADS"))                                                       
  if (!is.na(shiny_threads)) {                                                                                   
    cores <- shiny_threads                                                                                 
    print(paste("Set number of threads to",cores))                                                          
  }                                                                                                              
  
  
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
  
  #global variables
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
  
  # ## Reading in external data from e.g. a call of the window + message
  observe({
    if (!is.null(input$extdata)) {
      isolate({
        #        withProgress(message="Gathering external data ...", value="please wait",  {
        jsonmessage <- fromJSON(input$extdata)
        # Loading parameters
        NumCond <- jsonmessage[["numcond"]]
        NumReps <- jsonmessage[["numrep"]]
        isPaired <- jsonmessage[["paired"]]
        isGrouped <- jsonmessage[["grouped"]]
        hasProtnames <- jsonmessage[["modsandprots"]]
        # reading data matrix
        expr_matr <- jsonmessage[["expr_matrix"]]
        first_col <- ifelse(hasProtnames,3,2)
        output$fileInText <- renderText({
          validate(need(!is.null(expr_matr), "Uploaded data empty"))
          validate(need(length(expr_matr)>1, "Uploaded data does not contain multiple samples"))
          validate(need(sum(duplicated(expr_matr[[1]]),na.rm=T)==0,"Duplicated feature names in first column!"))
          tdat <- matrix(NA,nrow=length(expr_matr[[1]]),ncol=length(expr_matr)-first_col+1,dimnames=list(rows=expr_matr[[1]], cols=names(expr_matr)[first_col:length(expr_matr)]))
          for (i in first_col:length(expr_matr)) {
            validate(need(length(expr_matr[[i]]) == nrow(tdat),
                          paste("Wrong array length of sample", names(expr_matr)[i])))
            tdat[,i-first_col+1] <- as.numeric(expr_matr[[i]])
            
          }
          validate(need(is.numeric(tdat),"The uploaded data table contains non-numerical values!"))
          if (hasProtnames) {
            tdat <- data.frame(expr_matr[[2]],tdat)
          }
          v$example <- F
          v$dat <- tdat
          updateNumericInput(session,"NumCond",value=NumCond)
          updateNumericInput(session,"NumReps",value=NumReps)
          updateCheckboxInput(session, "isPaired", value=isPaired)
          updateCheckboxInput(session, "qcol_order", value=isGrouped)
          updateCheckboxInput(session, "protnames", value=hasProtnames)
          return(paste("Loaded external data"))
        })
        
        #       })
      })
    }
  })
  
  
  observeEvent(input$in_file,{
    ## test for right replicate and condition numbers, min 50 features, ...
    dat <- NULL
    withProgress(message="Reading file ...", min=0,max=1, value=0.5,  {
      try(dat <- read.csv(input$in_file$datapath,header=input$is_header))
      output$fileInText <- renderText({
        validate(need(sum(duplicated(dat[,1]),na.rm=T)==0,"Duplicated feature names in first column!"))
        rownames(dat) <- dat[,1]
        dat <- dat[,2:ncol(dat)]
        dat <- dat[rownames(dat)!="",]
        if(input$protnames) {
          validate(need(is.numeric(as.matrix(dat[,2:ncol(dat)])), "The data table contains non-numerical  values!"))
        } else {
          validate(need(is.numeric(as.matrix(dat)), "The data table contains non-numerical  values!"))
        }
        v$example <- F
        v$dat <- dat
        updateNumericInput(session,"NumCond",max=ifelse(input$protnames,ncol(dat)-1,ncol(dat)))
        updateNumericInput(session,"NumReps",max=ifelse(input$protnames,ncol(dat)-1,ncol(dat)))
        return("")
        
      })
    })
  })
  
  
  
  output$plot0 <- renderPlot({
    
    proteins <- NULL
    dat <- v$dat
    if(input$protnames & !v$example) {
      print(head(dat))
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
          if(!input$qcol_order) {
            dat <- dat[,rep(0:(NumCond-1),NumReps)*NumReps+rep(1:(NumReps), each=NumCond)]
          }
        }
        fulldat <- dat
        validate(need(try(statOut <- PrepareForVSClust(dat, NumReps, NumCond, input$isPaired, input$isStat)), 
                      "Please remove the following items from your input file:\na) empty columns or rows\nb) non-numerical or infinite values\nc) commenting characters (e.g. #)"))
        
        dat <- statOut$dat
        Sds <- dat[,ncol(dat)]
        output$data_summ <- renderUI({HTML(paste("Features:",nrow(dat),"<br/>Missing values:",
                                                 num_miss,"<br/>Median standard deviations:",
                                                 round(median(Sds,na.rm=T),digits=3)),"<br/>",
                                           paste("<i>Condition ",1:NumCond,":</i>", sapply(1:NumCond, function(x) 
                                             paste(colnames(fulldat)[(0:(NumReps-1))*NumCond+x],collapse=", ")),"<br/>",collapse=""))})
        
        
        
        pars$dat <<- dat 
        pars$proteins <<- proteins
        
        # file to download q-values
        output$downloadDataLimma <- downloadHandler(
          filename = function() {
            paste("LimmaResults", Sys.Date(), ifelse(input$limma_tsv,".tsv",".csv"), sep="");
          },
          content = function(file) {
            if (input$limma_tsv) {
              write.table(statOut$statFileOut, file, sep="\t", quote=F) 
            } else {
              write.csv(statOut$statFileOut, file)
            }
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
        clustInd <- NULL
        withProgress(message="Calculating ...", min=3,max=maxClust, value=2,  {
          clustInd <- estimClustNum(dat, input$maxclust, cores)
          estimClust.plot(clustInd)
          updateSliderInput(session,"nclust1",value=optimalClustNum(clustInd))
          updateSliderInput(session,"nclust2",value=optimalClustNum(clustInd, method="FCM"))
        })
        
        output$downloadParamEst <- downloadHandler(
          filename = function() {
            paste("EstimatedClustNumber", Sys.Date(), ".pdf", sep="");
          },
          content = function(file) {
            pdf(file,height=6,width=15)
            estimClust.plot(clustInd)
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
                paste("FCMVarMResults", Sys.Date(), ifelse(input$clustvar_tsv,".tsv",".csv"), sep="");
              },
              content = function(file) {
                if (input$clustvar_tsv) {
                  write.table(data.frame(cluster=Bestcl$cluster,ClustOut$outFileClust,isClusterMember=rowMaxs(Bestcl$membership)>0.5,maxMembership=rowMaxs(Bestcl$membership),
                                         Bestcl$membership), file, quote=F, sep="\t")
                } else {
                  write.csv(data.frame(cluster=Bestcl$cluster,ClustOut$outFileClust,isClusterMember=rowMaxs(Bestcl$membership)>0.5,maxMembership=rowMaxs(Bestcl$membership),
                                       Bestcl$membership), file)
                }
              })
            output$downloadCentroid2 <- downloadHandler(
              filename = function() {
                paste("FCMVarMResultsCentroids", Sys.Date(), ifelse(input$clustvar_tsv,".tsv",".csv"), sep="");
              },
              content = function(file) {
                if (input$clustvar_tsv) {
                  write.table(Bestcl$centers, file, quote=F, sep="\t")
                } else {
                  write.csv(Bestcl$centers, file)
                }
              })
            output$clustinf1 <- DT::renderDataTable(as.data.frame(ClustOut$ClustInd))
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
                paste("FCMResults", Sys.Date(), ifelse(input$cluststd_tsv,".tsv",".csv"), sep="");
              },
              content = function(file) {
                if (input$cluststd_tsv) {
                  write.table(data.frame(cluster=Bestcl$cluster,ClustOut$outFileClust,isClusterMember=rowMaxs(Bestcl$membership)>0.5,maxMembership=rowMaxs(Bestcl$membership),
                                         Bestcl$membership), file, quote=F, sep="\t")
                } else {
                  write.csv(data.frame(cluster=Bestcl$cluster,ClustOut$outFileClust,isClusterMember=rowMaxs(Bestcl$membership)>0.5,maxMembership=rowMaxs(Bestcl$membership),
                                       Bestcl$membership), file)
                  
                }
              })
            output$downloadCentroid3 <- downloadHandler(
              filename = function() {
                paste("FCMResultsCentroids", Sys.Date(), ifelse(input$cluststd_tsv,".tsv",".csv"), sep="");
              },
              content = function(file) {
                if (input$cluststd_tsv) {
                  write.table(Bestcl$centers, file, quote=F, sep="\t")
                } else {
                  write.csv(Bestcl$centers, file)
                }
              })
            output$clustinf2 <- DT::renderDataTable(as.data.frame(ClustOut$ClustInd))
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
          cl <- pars$Bestcl2
          dat <- pars$dat[,1:(ncol(pars$dat)-1)]
          withProgress(message="Waiting for data (1/2)...", min=0,max=1, {
            
            enriched <- runFuncEnrich(cl, dat, pars$proteins, input$idtype, input$infosource)
            x <- enriched$fullFuncs
            y <- enriched$redFuncs
            BHI <- enriched$BHI
            validate(need(!is.null(x),"No result. Wrong ID type?"))
            incProgress(0.7, detail = "received")
            incProgress(0.8, detail = "plotting")
            validate(need(nrow(x@compareClusterResult)>1,"No significant results"))
            output$downloadGOData1 <- downloadHandler(
              filename = function() {
                paste("DAVIDResultsVarM", Sys.Date(), ".csv", sep="");
              },
              content = function(file) {
                write.csv(as.data.frame(x@compareClusterResult), file)
              })
            dotplot(y,title=paste("BHI:",round(BHI,digits=3)),showCategory=20,font.size=10)
          })
        }})}})
  
  output$plot5 <- renderPlot({
    if(input$goButton == 0)
      return()
    validate(need(!is.null(pars$Bestcl1),"Run clustering first"))
    if (!is.null(pars$Bestcl1)) {
      isolate({
        if(input$goButton) {
          cl <- pars$Bestcl1
          dat <- pars$dat[,1:(ncol(pars$dat)-1)]
          withProgress(message="Waiting for data (2/2)...", min=0,max=1, {
            
            enriched <- runFuncEnrich(cl, dat, pars$proteins, input$idtype, input$infosource)
            x <- enriched$fullFuncs
            y <- enriched$redFuncs
            BHI <- enriched$BHI
            validate(need(!is.null(x),"No result. Wrong ID type?"))
            incProgress(0.7, detail = "received")
            incProgress(0.8, detail = "plotting")
            validate(need(nrow(x@compareClusterResult)>1,"No significant results"))
            output$downloadGOData2 <- downloadHandler(
              filename = function() {
                paste("DAVIDResultsStand", Sys.Date(), ".csv", sep="");
              },
              content = function(file) {
                write.csv(as.data.frame(x@compareClusterResult), file)
              }) 
            #print(y@compareClusterResult)
            dotplot(y,title=paste("BHI:",round(BHI,digits=3)),showCategory=20,font.size=10)
          })
        }})}
  })
  
  output$intro <- renderUI({HTML("We present a new method to apply fuzzy c-means clustering to data sets that exhibit 
                                  non-constant variance of their features. The individual variation is incorporated into 
the estimation of the fuzzifier (parameter m). For details see<br/><a href='https://doi.org/10.1093/bioinformatics/bty224'>Veit Schw&auml;mmle, Ole N Jensen; VSClust: Feature-based variance-sensitive clustering of omics data, Bioinformatics, 2018, bty224, https://doi.org/10.1093/bioinformatics/bty224</a><br/> 
                       The algorithm for the estimation of the parameters fuzzifier and cluster numbers is furthermore based on
<br/><a href='http://www.ncbi.nlm.nih.gov/pubmed/20880957'>Schw&auml;mmle, V. and Jensen, O. N. &quotA simple and fast method to determine the parameters for fuzzy c-means cluster analysis&quot. <i>Bioinformatics</i>, 2010, <b>26</b>, 2841-2848<br/></a>
                                 Please note that this method reveals its power for 8 or more different conditions (dimensions). Lower numbers yields results nearly identical to standard fuzzy c-means clustering.<br/>
                                 <b><i>Example data set:</i></b> You can test the method with an artificial data set by clicking on <i>Load example</i>. The data set contains 500 features where half of them 
                                 make part of five predefined clusters while the rest of the data is randomly distributed. The different features 
                                 do not contain identifiers for subsequent gene enrichement analysis (long waiting times for response from the DAVID web service would compromise VSClust performance by multiple users).")})
  output$finput <- renderUI({HTML("Input format is restricted to comma-delimited files (.csv). Files are required to contain only the numerical data that will
be analyzed in addition to the following contents: A one-row header containing the 
                                  column names and additionally defined feature names for the rows (e.g. gene names, transcripts, peptides and proteins). The latter
                                  need to be unique, duplicated feature names are not accepted. Additionally, in the case of data containing multiple measurements of the same gene/protein, a second column can contain this information
that then will be used for the GO term enrichment.<br/>
                                  The order of the columns is crucial to distinguish conditions from replicated values. With numbers denoting conditions (1-4) and letters
corresponding to replicates (A-C), grouped columns are arranged A1, A2, A3, A4, B1, B2, B3, B4, C1, C2, C3, C4 and ungrouped columns as A1, B1, C1,..., A4, B4, C4,...  <br/>
                                  In the case of an input file that already contains estimated standard deviations, all columns but the last ones are considered 
                                  different conditions while the last column contains the standard deviations. In this case, untick the 'Estimate variance levels from replicated quantifications' checkbox.")})
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
                                  lower numbers of cluster members as features with large standard deviations often get discarded.")})
  output$goterms <- renderUI({HTML("Simple viewer for GO term and pathway enrichment based on <a  target='_blank'' href='https://bioconductor.org/packages/release/bioc/html/clusterProfiler.html'>clusterProfiler</a> 
                                    configured to access the <a href='https://david.ncifcrf.gov/content.jsp?file=WS.html'>DAVID</a> web service. <i>This can take a while (response times of DAVID server)</i>. We furthermore offer a rough estimate of enrichment efficiency, 
                                   given by the biological homogeneity index (BHI, see also <a href='http://www.ncbi.nlm.nih.gov/pmc/articles/PMC1590054/'>paper</a>)<br/>
<b>Important: </b>In the case of more than 2000 features in cluster, they will be limited to the 2000 with the highest membership values.
    ")})
  output$reading <- renderUI({HTML("The source code of VSClust is available at <a href='https://bitbucket.org/veitveit/vsclust'>bitbucket</a>. VSClust was developed and implemented 
                                   at the <a href='http://www.sdu.dk/en/Om_SDU/Institutter_centre/Bmb_biokemi_og_molekylaer_biologi/Forskning/Forskningsgrupper/Protein.aspx'>Protein Research Group</a> of the University of Southern Denmark. See also <a href='computproteomics.bmb.sdu.dk'>computproteomics.bmb.sdu.dk</a> for more information.")})
  output$DownloadExample <- downloadHandler(
    filename = function() {
      paste("ExampleFile.csv")
    },
    content = function(file) {
      file.copy("ArtData.csv",file)
    },
    contentType = "application/csv"
  )
  
})
