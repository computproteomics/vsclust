# Function to calculate the parameters for the fuzzy c-means clustering algorithm. 
# For details, see V. Schwämmle and O. N. Jensen, Bioinformatics, Vol. 26, Nr. 22, 2010, 2841. 
# The Mfuzz library is required
# This package estimates the fuzzifier using the a specific formula the works well for >=4 dimensions and more than 1000 data points
# The cluster number is estimate using the well-known Xie-Beni validation index and an own measure here called minimum cluster distance
# In the case of now agreement of both measures, the latter is taken.
# Input values are: 
# eset: ExpressionSet of your data set
# maxc: maximal number of clusters to be taken into account its estimation
# samples: number of samples taken for calculating a best estimation of optimal cluster results
FcmClustPEst <- function (eset, maxc=20, samples=10) 
{

library(Mfuzz)  

cvalidate.xiebeni <- function(clres,m) {                                                
        xrows <- dim(clres$me)[1]                                                
        minimum <- -1                                                            
        error <- clres$within                                                    
        ncenters <- dim(clres$centers)[1]                                        
        for (i in 1:(ncenters - 1)) {                                            
            for (j in (i + 1):ncenters) {                                        
                diff <- clres$ce[i, ] - clres$ce[j, ]                            
                diffdist <- t(diff) %*% t(t(diff))                               
                if (minimum == -1)                                               
                  minimum <- diffdist                                            
                if (diffdist < minimum)                                          
                  minimum <- diffdist                                            
            }                                                                    
        }                                                                        
        xiebeni <- error/(xrows * minimum)                                       
        return(xiebeni)                                                          
}                                                                            


dims<-dim(exprs(eset))
MaxClust <- maxc
# calculate from magic formula
m <- 1 + (1418/dims[1] + 22.05)* dims[2]^(-2) + (12.33/dims[1] + 0.243)*dims[2]^(-0.0406*log(dims[1])-0.1134)
print(paste("Estimated m=",m))
NumS <- samples
MaxC <- vector(,1)
Wk <- vector(,MaxClust)
Bestcl<-NULL
for (NClust in 2:MaxClust) {
      MaxC[m] <- NClust-1
      clBest <- 100000
      print(paste("fcm-clustering for",NClust,"clusters"))
      for(NS in 1:NumS) {
        cl <- mfuzz(eset,c = NClust, m=m)
        if (cl$withinerror < clBest) {
          Bestcl[[NClust]] <- cl
          clBest <- cl$withinerror
        }
        dimsc <- dim(cl$membership)
      }
}


xb<-vector(,MaxClust)
 mindistc <- vector(,MaxClust)

for (NClust in 2:MaxClust) {
  xb[NClust] <- cvalidate.xiebeni(Bestcl[[NClust]],m)
  mindistc[NClust] <- min(dist(Bestcl[[NClust]]$centers))
}

# X11()
par(mfrow=c(2,1),ps=15)
plot(2:NClust,mindistc[2:NClust],main="minimum centroid distance",xlab="c",ylab="")
plot(2:NClust,log(xb[2:NClust]), main="log(xie beni index)",xlab="c",ylab="")

dmindist<-vector(,(NClust-1))
dmindist<-mindistc[3:(NClust-2)]-mindistc[4:(NClust-1)]
c1<-which.max(dmindist)+2
c2<-which.min(xb[2:NClust])+1

print(paste("best c from mindist: ",c1))
print(paste("best c from xie beni index: ",c2))
ccc<-NULL
ccc<-c(m,c1)
return(ccc)
}