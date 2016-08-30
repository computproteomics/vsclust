library(matrixStats)
library(Mfuzz)
library(limma)
library(qvalue)
library(e1071)
library(clusterProfiler)
library(RDAVIDWebService)
source("FcmClustPEst.R")
source("mfuzz.plotpdf.R")
source("HelperFuncs.R")

options(shiny.maxRequestSize=500*1024^2) 

dat <- NULL
pars <- NULL
pars$m <- NULL
maxClust <- 25
NumCond <- NULL

shinyServer(function(input, output,clientData,session) {
  
  
  
  # setting the right tab
  observe({
    input$NumReps
    input$NumCond
    input$isPaired
    input$isStat
    input$row.names
    input$infile
    input$is_header
    updateTabsetPanel(session, inputId ="tabset", selected = "stat")
  })
  
  output$description <- renderText("A new method for fuzzy c-means clustering adapted to data with non-constant feature variance. Manuscript in preparation, the method is based on a previous work that estimates the parameter values (fuzzifier and number of clusters) for standard fuzzy c-means clustering.")
  output$tparame <- renderText("Parameter estimation")
  output$tresults <- renderText("Clustering results")
  output$tdescr <- renderText("Description")
  output$ui <- renderUI({
    if (is.null(input$isStat) | input$isStat) {
      p(
        checkboxInput(inputId="isPaired", label="Paired tests",value=T),
        sliderInput("NumReps",min=2,max=20,value=2,label="Number of replicates",step=1),
        sliderInput("NumCond",min=2,max=20,value=3,label="Number of conditions",step=1))
    }
  })
  v <- reactiveValues(dat = NULL, example = F)
  observeEvent(input$examplefile,{
    dat <- read.csv("ArtData.csv",row.names=1)
    dat <- dat[1:nrow(dat),]
    print("example")
    updateNumericInput(session,"NumCond",value=10,max=10)
    updateNumericInput(session,"NumReps",value=5,max=5)
    updateCheckboxInput(session,"isPaired",value=F)
    v$dat <- dat
    v$example <- T
  })
  observeEvent(input$in_file,{
    ## test for right replicate and condition numbers, min 50 features, ...
    dat <- NULL
    try(dat <- read.csv(input$in_file$datapath,row.names=1,header=input$is_header))
    dat <- dat[rownames(dat)!="",]
    v$example <- F
    v$dat <- dat
    updateNumericInput(session,"NumCond",max=ifelse(input$protnames,ncol(dat)-1,ncol(dat)))
    updateNumericInput(session,"NumReps",max=ifelse(input$protnames,ncol(dat)-1,ncol(dat)))
  })
  output$plot0 <- renderPlot({
    dat <- NULL
    proteins <- NULL
    dat <- v$dat
    if(input$protnames & !v$example) {
      proteins <- dat[,1]
      dat <- dat[,2:ncol(dat)]
      names(proteins) <- rownames(dat)
    }
    fdat <- dat
    maxClust <- input$maxclust
    print(dim(dat))
    print(input$in_file)
    validate(need(!is.null(dat),"no data"))
    validate(need(is.numeric(as.matrix(dat)), "The data table contains non-numerical  values!"))
    if(!is.null(dat)) {
      NumReps <- input$NumReps
      NumCond <- input$NumCond
      print(NumReps)
      print(NumCond)
      qvals <- NULL
      if (input$isStat) {
        validate(need(ncol(dat)==NumReps*NumCond, "Number of data columns must correspond to product of conditions and replicates!"))
        if (input$isPaired) {
          ttt <- SignAnalPaired(dat,NumCond, NumReps)
        } else {
          ttt <- SignAnal(dat,NumCond, NumReps)
        }
        tdat<-rowMeans(dat[,seq(1,NumReps*NumCond,NumCond)],na.rm=T)
        Sds <- ttt$Sds
        qvals <- ttt$plvalues
        for (i in 2:NumCond) {
          tdat<-cbind(tdat,rowMeans(dat[,seq(i,NumReps*NumCond,NumCond)],na.rm=T))
        }
        colnames(tdat)<-paste("Mean of log T",0:(NumCond-1),sep="")
        dat <- cbind(tdat,Sds=Sds)
        #      dat <- dat[rowSums(is.na(ttt$plvalues))==0,]
        
      } else {
        Sds <- dat[,ncol(dat)]
        tdat <- dat[,1:(ncol(dat)-1)]
      }
      print(colnames(dat))
      output$data_summ <- renderUI({HTML(paste("Features:",nrow(dat),"<br/>Missing values:",
                                               sum(is.na(v$dat)),"<br/>Median standard deviations:",
                                               round(median(Sds,na.rm=T),digits=3)))})
      #    updateTabsetPanel(session, inputId ="tabset", selected = "pest")
      dat[,ncol(dat)] <- dat[,ncol(dat)]/ rowSds(as.matrix(tdat),na.rm=T)
      pars$dat <<- dat 
      pars$proteins <<- proteins
      tdat <- dat[,1:(ncol(dat)-1)]
      Sds <- dat[,ncol(dat)]
      
      # file to download q-values
      output$downloadDataLimma <- downloadHandler(
        filename = function() {
          paste("LimmaResults", Sys.Date(), ".csv", sep="");
        },
        content = function(file) {
          dim(qvals)
          colnames(qvals) <- paste("qvalue T",1:(NumCond-1),"vsT0",sep="")
          outdat <- cbind(dat,qvals)
          write.csv(outdat, file)
        })
      
      pcaDat <- t(fdat[complete.cases(fdat),])
      validate(need(length(pcaDat)> 0, "Principal component analysis not shown as too many missing values"))      
      validate(need(nrow(pcaDat)> 10, "Principal component analysis not shown as too many missing values"))      
      pca<-prcomp(pcaDat,scale=T,retx=T)
      scores <- pca$x
      loadings <- pca$rotation
      par(mfrow=c(1,2))
      plot(loadings,cex=Sds/5,pch=19)
      title(main="Principal component analysis of data set (loadings)",sub="The point size corresponds to the estimated standard deviation")
      plot(scores,pch=19,col=rainbow(NumCond)[rep(1:NumCond,NumReps)])
      title(main="Principal component analysis of data set (scores)",sub="Colors denote different conditions")
      legend("topright",paste("Condition",1:NumCond),pch=rep(19,NumCond),col=rainbow(NumCond)[1:NumCond])
    }
  })
  
  output$plot1 <- renderPlot({
    if(input$clButton1 == 0) 
      return()
    isolate({
      if((input$clButton1)) {
        dat <- pars$dat
        maxClust <- input$maxclust
        ClustInd<-matrix(NA,nrow=maxClust,ncol=4)
        withProgress(message="Calculating ...", min=3,max=maxClust, value=2,  {
          for (NClust in 3:maxClust) {
            print(NClust)
            incProgress(1, detail = paste("Running cluster number",NClust))
            clustout <- ClustComp(dat[,1:(ncol(dat)-1)],NClust=NClust,Sds=dat[,ncol(dat)],NSs=32)
            ClustInd[NClust,]<-clustout$indices
          }
        })
        dmindist <- c(which.max(ClustInd[3:(maxClust-2),1]-ClustInd[4:(maxClust-1),1])+2,
                      which.max(ClustInd[3:(maxClust-2),3]-ClustInd[4:(maxClust-1),3])+2)
        dxiebeni <- c(which.min(ClustInd[3:(maxClust-1),2])+2,
                      which.min(ClustInd[3:(maxClust-1),4])+2)
        print(dmindist)
        par(mfrow=c(2,1))
        plot(3:maxClust,ClustInd[3:(maxClust),1],col=2 , type="b", 
             main="Min. centroid distance\n(Highest jump is best)",
             xlab="Number of clusters", ylab="Index",ylim=c(min(ClustInd[,c(1,3)],na.rm=T),max(ClustInd[,c(1,3)],na.rm=T)))
        lines(3:maxClust,ClustInd[3:(maxClust),3],col=3,type="b")
        points(dmindist[1],ClustInd[dmindist[1],1],pch=15,col=1,cex=2)
        legend("topright",legend = c("Variance-based method","Standard"), lty=c(1,1),col=2:3)
        plot(3:maxClust,ClustInd[3:(maxClust),2], col=2, type="b", main="Xie-Beni index\n(Lowest is best)",
             xlab="Number of clusters", ylab="Index",ylim=c(min(ClustInd[,c(2,4)],na.rm=T),max(ClustInd[,c(2,4)],na.rm=T)))
        lines(3:maxClust,ClustInd[3:(maxClust),4],type="b",col=3)
        points(dxiebeni[1],ClustInd[dxiebeni[1],2],pch=15,col=1,cex=2)
        legend("topright",legend = c("Variance-based method","Standard"), lty=c(1,1),col=2:3)
        
        output$downloadParamEst <- downloadHandler(
          filename = function() {
            paste("EstimatedClustNumber", Sys.Date(), ".pdf", sep="");
          },
          content = function(file) {
            pdf(file,height=12)
            par(mfrow=c(2,1))
            plot(3:maxClust,ClustInd[3:(maxClust),1],col=2 , type="b", 
                 main="Min. centroid distance\n(Highest jump is best)",
                 xlab="Number of clusters", ylab="Index",ylim=c(min(ClustInd[,c(1,3)],na.rm=T),max(ClustInd[,c(1,3)],na.rm=T)))
            lines(3:maxClust,ClustInd[3:(maxClust),3],col=3,type="b")
            points(dmindist[1],ClustInd[dmindist[1],1],pch=15,col=1,cex=2)
            legend("topright",legend = c("Variance-based method","Standard"), lty=c(1,1),col=2:3)
            plot(3:maxClust,ClustInd[3:(maxClust),2], col=2, type="b", main="Xie-Beni index\n(Lowest is best)",
                 xlab="Number of clusters", ylab="Index",ylim=c(min(ClustInd[,c(2,4)],na.rm=T),max(ClustInd[,c(2,4)],na.rm=T)))
            lines(3:maxClust,ClustInd[3:(maxClust),4],type="b",col=3)
            points(dxiebeni[1],ClustInd[dxiebeni[1],2],pch=15,col=1,cex=2)
            legend("topright",legend = c("Variance-based method","Standard"), lty=c(1,1),col=2:3)
            dev.off()            
          })
        
        #         updateNumericInput(session,"nclust1",max=maxClust,value=dmindist[1])
        #         updateNumericInput(session,"nclust2",max=maxClust,value=dmindist[1])
        #updateTabsetPanel(session, inputId ="tabset", selected = "clust1")
      }
    })},height=600)
  observeEvent(input$clButton2, {
    par(mar=c(2,2,1,3),mgp=c(2,1,0))
    if (!is.null(pars$dat)) {
      print(input$nclust1)
      Nclust <- input$nclust1
      pars$dat <- pars$dat[rowSums(is.na(pars$dat))==0,]
      
      dat <- new("ExpressionSet",expr=as.matrix(pars$dat[,1:(ncol(pars$dat)-1)]))
      PExpr <- new("ExpressionSet",expr=as.matrix(dat))
      PExpr.r                                                                                                                                                    <- filter.NA(PExpr, thres = 0.25)
      PExpr <- fill.NA(PExpr.r,mode = "mean")
      tmp <- filter.std(PExpr,min.std=0,visu=F)
      dat <- standardise(PExpr)
      
      withProgress(message="Calculating ...", min=0,max=1, {
        
        clustout <- ClustComp(exprs(dat),NClust=Nclust,Sds=pars$dat[,ncol(pars$dat)],NSs=32)
        Bestcl <- clustout$Bestcl
        Bestcl <- SwitchOrder(Bestcl,Nclust)
        pars$Bestcl2 <<- Bestcl  
        
        # sorting for membership values (globally)
        Bestcl$cluster <- Bestcl$cluster[order(rowMaxs(Bestcl$membership,na.rm=T))]
        Bestcl$membership <- Bestcl$membership[order(rowMaxs(Bestcl$membership,na.rm=T)),]
        dat <- dat[names(Bestcl$cluster),]
        
        incProgress(0.7, detail = paste("Plotting",Nclust))
        output$plot2 <- renderPlot({
          par(lwd=0.25)
          oldmar <- par("mar")
          par(mar=par("mar")/max(1,Nclust/20))
          mfuzz.plot2(dat,cl=Bestcl,mfrow=c(round(sqrt(Nclust)),ceiling(sqrt(Nclust))),min.mem=0.5,x11=F,colo="fancy")
          par(lwd=1,mar=oldmar)
        })
        output$downloadData2 <- downloadHandler(
          filename = function() {
            paste("FCMVarMResults", Sys.Date(), ".csv", sep="");
          },
          content = function(file) {
            colnames(Bestcl$membership) <- paste("membership of cluster",colnames(Bestcl$membership))
            outdat <- exprs(dat)
            if (!is.null(pars$proteins)) {
              outdat <- cbind(outdat,names=as.character(pars$proteins[rownames(outdat)]))
            }
            write.csv(data.frame(cluster=Bestcl$cluster,outdat,isClusterMember=rowMaxs(Bestcl$membership)>0.5,maxMembership=rowMaxs(Bestcl$membership),
                                 Bestcl$membership), file)
          })
        output$downloadCentroid2 <- downloadHandler(
          filename = function() {
            paste("FCMVarMResultsCentroids", Sys.Date(), ".csv", sep="");
          },
          content = function(file) {
            rownames(Bestcl$centers) <- paste("Cluster",rownames(Bestcl$centers))
            write.csv(Bestcl$centers, file)
          })
        output$downloadFigure <- downloadHandler(
          filename = function() {
            paste("FCMVarMResults", Sys.Date(), ".pdf", sep="");
          },
          content = function(file) {
            pdf(file,height=5*round(sqrt(Nclust)),width=5*ceiling(sqrt(Nclust)))
            par(lwd=0.25)
            oldmar <- par("mar")
            par(mar=par("mar")/max(1,Nclust/20))
            mfuzz.plot2(dat,cl=Bestcl,mfrow=c(round(sqrt(Nclust)),ceiling(sqrt(Nclust))),min.mem=0.5,x11=F,colo="fancy")
            par(lwd=1,mar=oldmar)
            dev.off()
            # mfuzz.plotpdf(dat,cl=Bestcl,mfrow=c(round(sqrt(Nclust)),ceiling(sqrt(Nclust))),min.mem=0.5,filename=file)
          })
        ClustInd <- as.data.frame(table(Bestcl$cluster[rowMaxs(Bestcl$membership)>0.5]))
        colnames(ClustInd) <- c("Cluster","Members")
        output$clustinf1 <- renderDataTable(ClustInd)
      })
    }})
  observeEvent(input$clButton3, {
    if (!is.null(pars$dat)) {
      par(mar=c(2,2,1,3),mgp=c(2,1,0))
      Nclust <- input$nclust2
      pars$dat <- pars$dat[rowSums(is.na(pars$dat))==0,]
      
      dat <- new("ExpressionSet",expr=as.matrix(pars$dat[,1:(ncol(pars$dat)-1)]))
      PExpr <- new("ExpressionSet",expr=as.matrix(dat))
      PExpr.r <- filter.NA(PExpr, thres = 0.25)
      PExpr <- fill.NA(PExpr.r,mode = "mean")
      tmp <- filter.std(PExpr,min.std=0,visu=F)
      dat <- standardise(PExpr)
      
      withProgress(message=paste("Calculating ..."), min=0,max=1, {
        clustout <- ClustComp(exprs(dat),NClust=Nclust,Sds=pars$dat[,ncol(pars$dat)],NSs=32)
        Bestcl <- clustout$Bestcl2
        Bestcl <- SwitchOrder(Bestcl,Nclust)
        pars$Bestcl1 <<- Bestcl  
        
        # sorting for membership values (globally)
        Bestcl$cluster <- Bestcl$cluster[order(rowMaxs(Bestcl$membership,na.rm=T))]
        Bestcl$membership <- Bestcl$membership[order(rowMaxs(Bestcl$membership,na.rm=T)),]
        dat <- dat[names(Bestcl$cluster),]
        
        incProgress(0.7, detail = paste("Plotting",Nclust))
        output$plot3 <- renderPlot({
          par(lwd=0.25)
          mfuzz.plot2(dat,cl=Bestcl,mfrow=c(round(sqrt(Nclust)),ceiling(sqrt(Nclust))),min.mem=0.5,x11=F,colo="fancy")
          par(lwd=1)
        })
        output$downloadData3 <- downloadHandler(
          filename = function() {
            paste("FCMResults", Sys.Date(), ".csv", sep="");
          },
          content = function(file) {
            colnames(Bestcl$membership) <- paste("membership of cluster",colnames(Bestcl$membership))
            outdat <- exprs(dat)
            if (!is.null(pars$proteins)) {
              outdat <- cbind(outdat,names=as.character(pars$proteins[rownames(outdat)]))
            }
            write.csv(data.frame(cluster=Bestcl$cluster,outdat,isClusterMember=rowMaxs(Bestcl$membership)>0.5,maxMembership=rowMaxs(Bestcl$membership),Bestcl$membership), file)
          })
        output$downloadCentroid3 <- downloadHandler(
          filename = function() {
            paste("FCMResultsCentroids", Sys.Date(), ".csv", sep="");
          },
          content = function(file) {
            rownames(Bestcl$centers) <- paste("Cluster",rownames(Bestcl$centers))
            write.csv(Bestcl$centers, file)
          })        
        output$downloadFigure2 <- downloadHandler(
          filename = function() {
            paste("FCMResults", Sys.Date(), ".pdf", sep="");
          },
          content = function(file) {
            pdf(file,height=5*round(sqrt(Nclust)),width=5*ceiling(sqrt(Nclust)))
            par(lwd=0.25)
            mfuzz.plot2(dat,cl=Bestcl,mfrow=c(round(sqrt(Nclust)),ceiling(sqrt(Nclust))),min.mem=0.5,x11=F,colo="fancy")
            par(lwd=1)
            dev.off()
            # mfuzz.plotpdf(dat,cl=Bestcl,mfrow=c(round(sqrt(Nclust)),ceiling(sqrt(Nclust))),min.mem=0.5,filename=file)
          })
        ClustInd <- as.data.frame(table(Bestcl$cluster[rowMaxs(Bestcl$membership)>0.5]))
        colnames(ClustInd) <- c("Cluster","Members")
        output$clustinf2 <- renderDataTable(ClustInd) 
      })
    }})
  
  
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
          # print(Accs)
          withProgress(message="Waiting for data (1/2)...", min=0,max=1, {
            
            x <- compareCluster(Accs, fun="enrichDAVID2", annotation=input$infosource,
                                idType=input$idtype,
                                listType="Gene", david.user = "veits@bmb.sdu.dk")
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
          # print(Accs)
          withProgress(message="Waiting for data (2/2) ...", min=0,max=1, {
            x <- compareCluster(Accs, fun="enrichDAVID2", annotation=input$infosource,
                                idType=input$idtype,
                                listType="Gene", david.user = "veits@bmb.sdu.dk")
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
<br/>Schwämmle, V. and Jensen, O. N. &quotA simple and fast method to determine the parameters for fuzzy c-means cluster analysis&quot. <i>Bioinformatics</i>, 2010, <b>26</b>, 2841-2848<br/>
                                 Please note that this method reveals its power for 8 or more different conditions (dimensions). Lower numbers yields results nearly identical to standard fuzzy c-means clustering.<br/>
                                 <i>Example data set:</i> You can test the method with an artificial data set by clicking on <i>Load example</i>. The data set contains 500 features where half of them 
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
  output$stat <- renderUI({HTML("We offer to calculate feature standard deviations using the limma package, well-performing for micro-array, mass spectrometry and RNA-seq data. 
                                The method corrects feature variation by ensuring a minimal, adapted variation that is estimated from the entire data set.  
                                Statistics and estimation of feature variation can be carried out using paired and unpaired tests. <br/>
                                The plot on the &quotStatistics and variance&quot tab shows a simplified visualization of the data projected on the first principal components (standard PCA). Point sizes correspond to 
                                estimated feature standard deviation.")})
  output$pest <- renderUI({HTML("Estimation of an optimal cluster number can become tricky and many validation methods have been proposed
                                to pick out the right number. The tab &quotEstimation of cluster number&quot shows the values of the minimum 
                                centroid distance and the Xie-Beni index over the user-defined range of cluster numbers. Suitable cluster numbers are given by large jumps
(minimal centroid distance) and minima (Xie-Beni index). The best estimates for the variance-sensitive method are marked by black squares.")})
  output$fclust <- renderUI({HTML("Clustering results are available in 2 tabs, showing the results for our variance-sensitive algorithm and 
                                  standard c-means clustering. Users are required to set the number of clusters. Expression profiles are shown in separate panels where the color corresponds to the membership values.
                                  The plot shows only features confidently assigned to a cluster requiring a minimum membership value of 0.5. Users can download 
                                  pdf-files of the figure, data with membership values and principal expression profiles of the clusters (cluster centroids). Furthermore, the number
                                  of cluster members is shown in a table. This can be useful when comparing the two methods, where the variance-based algorithm should yield 
                                  lower numbers of cluster members as features with large standard deviations often get discarded.")})
  output$goterms <- renderUI({HTML("Simple viewer for GO term enrichment based on clusterProfiler configured to access DAVID .... We furthermore offer a rough estimate of enrichment efficiency, 
                                   given by the biological homogeneity index (BHI)")})
  output$reading <- renderUI({HTML("Used libraries, further reading, my web page, PRG, paper ..")})
  
})
