###############################################################
###    SIMULATIONS: TYPE I, TYPE II AND MEAN ERROR RATES    ###
###############################################################


## PRELIMINARY

#HPC
args <- commandArgs(TRUE)
conditionline <- as.numeric(as.character(args))[1]

#Paths
scripts <- "path where your scripts are located"
datastor <- "path where maps and results will be stored"


#Libraries
library(neuRosim)
library(oro.nifti)
library(AnalyzeFMRI)
library(lattice)
library(MASS)
library(smoothie)


#Source files
setwd(scripts)
source("functions_GLM_2D.R")

#Define parameter values
conditions <- read.csv(paste(scripts,"conditions_hpc.txt",sep=""))
delta <- conditions[conditionline, 1]*100
nscan <- conditions[conditionline, 2]
perc <- conditions[conditionline, 3]
es <- 0.02
cnr <- 0.46
sigmanoise <- es/cnr 
nsim <- 1000
base <- 100
alpha <- 0.05
FDRqval <- 0.05
k <- 1.5 ##Should be adjusted to the treshold of your choice. 
#For this paper: k=8 for LR methods and k=0.88,1.5 or 3.68




#Image characteristics
dim1 <- 32
dim2 <- 32 
imdim <- c(dim1,dim2) 
voxdim <- c(3.5,3.5)
TR <- 2
nregio <- 2


#CREATING GROUND TRUTH

truth1 <- array(0,dim=imdim)
truth2 <- array(0,dim=imdim)
truth1[c(5:9),c(4:10)] <- 1 
truth2[c(15:25),c(20:25)] <- 1 
truth <- truth1 + truth2 
actvox <- sum(truth)
inactvox <- prod(imdim)-actvox
inact <- ifelse(truth == 1, 0, 1)

#Storage of results
counts <- c()
  
tmp.results <- array(NA,dim=c(nsim,15))
colnames(tmp.results) <- c("FP.bon","FN.bon","AVG.bon","FP.LR","FN.LR","AVG.LR","FP.fdr","FN.fdr","AVG.fdr","FP.mlr","FN.mlr","AVG.mlr","FP.lrdel","FN.lrdel","AVG.lrdel")
  

set.seed(171018)

  for(i in 1:nsim) {
    
    print(i)
   
    ## CREATING DESIGN AND SIMULATED fMRI TIME SERIES 
    
    #putting in the temporal signal: block design 20s ON/OFF
    total.time <- nscan * TR
    dur <- total.time/8
    onsets <- c(1, total.time/4+1, 2*total.time/4+1, 3*total.time/4+1)
    X <- simprepTemporal(total.time,1,onsets=onsets,effectsize = 100,durations=dur,TR=TR,acc=0.1, hrf="double-gamma") 
    pred <- simTSfmri(design=X, base=base, SNR=1, noise="none", verbose=FALSE) 
    
    #Active voxels
    design.act <- es * (pred-base) + base

    #Inactive voxels 
    design.inact <- rep(base,times=length(design.act)) 
    
    #Smoothed signal and noise
    rawsignal1 <- inact %o% design.inact
    rawsignal2 <- truth %o% design.act
    signal <- rawsignal1 + rawsignal2
    
    noisim <- array(rnorm(prod(imdim)*nscan,0,sigmanoise),dim=c(imdim,nscan))
    data <- noisim + signal
    
    smdata <- array(NA, dim=dim(data))
    
    for(j in 1:nscan) {
      
      smdata[,,j] <- kernel2dsmooth(data[,,j], kernel.type="gauss", sigma=1.5, nx=32, ny=32)
      
    }
    
    
    ## CALCULATING B1, S(B1), T and P FOR EACH VOXEL
    b1 <- tbeta(data,pred)
    sb1 <- tsebeta(b1,data,pred)
    tmap <- b1/sb1
    pmap.null <- pt(tmap, nscan-2,lower.tail=FALSE)
    
    

    ##ANALYZING DATA USING BONFERRONI CORRECTION AND LR+95th PERCENTILE
    
    #Bonferroni correction
    active.bonf <- ifelse(pmap.null <= alpha/prod(dim(tmap)), 1, 0)
    
    #LR + xth percentile (Kang et al., 2015)
    LR <- array(0, dim=dim(b1))
    teller<-dnorm(b1, mean=quantile(b1,perc,na.rm=TRUE), sd=sb1)
    noemer<-dnorm(b1, mean=0, sd=sb1)
    LR <- teller/noemer
    active.LR <- ifelse(LR >= k , 1, 0)
    
    #FDR corrected
    fdr.thresh <- Threshold.FDR(c(tmap), cV.type=2, type="t",df1=nscan-2,q=0.05)
    active.fdr <- ifelse(tmap >= fdr.thresh, 1, 0)
    
    #mLR
    b1stats <- c(b1)
    sb1stats <- c(sb1)
    
    teller_mLR <- c()
    noemer_mLR <- c()
    
    for(l in 1:length(b1stats)) {
      
      if (b1stats[l] > delta) {teller_mLR[l] <- dnorm(b1stats[l], mean=b1stats[l], sd=sb1stats[l])}
      if (b1stats[l] <= delta) {teller_mLR[l] <- dnorm(b1stats[l], delta, sd=sb1stats[l])}
      if (b1stats[l] < delta) {noemer_mLR[l]<-dnorm(b1stats[l], mean=b1stats[l], sd=sb1stats[l])}
      if (b1stats[l] >= delta){noemer_mLR[l]<-dnorm(b1stats[l], delta, sd=sb1stats[l])}
      
    }
    
    mLR <-teller_mLR/noemer_mLR
    mLR <- array(mLR, dim=dim(tmap))
    active.mLR <- ifelse(mLR >= k , 1, 0)
    active.mLR[is.na(active.mLR)] <- 0
    
    
    #LR + delta
    LRdel <- array(0, dim=dim(b1))
    tellerdel<-dnorm(b1, mean=delta, sd=sb1)
    noemerdel<-dnorm(b1, mean=0, sd=sb1)
    LRdel <- tellerdel/noemerdel
    active.LRdel <- ifelse(LRdel >= k , 1, 0)
    
    

    ##COMPUTE ERRORS
    
    emap.bonf <- array(NA, dim=dim(tmap)) 
    emap.LR <- array(NA, dim=dim(tmap)) 
    emap.fdr <- array(NA, dim=dim(tmap)) 
    emap.mLR <- array(NA, dim=dim(tmap)) 
    emap.LRdel <- array(NA, dim=dim(tmap))
    
    emap.bonf <- ifelse(active.bonf==1 & truth==1, 1, NA) #TP
    emap.bonf[active.bonf==1 & truth==0] <- 2 #FP  
    emap.bonf[active.bonf==0 & truth==1] <- 3 #FN
    emap.bonf[active.bonf==0 & truth==0] <- 4 #TN
    
    emap.LR <- ifelse(active.LR==1 & truth==1, 1, NA) #TP
    emap.LR[active.LR==1 & truth==0] <- 2 #FP  
    emap.LR[active.LR==0 & truth==1] <- 3 #FN
    emap.LR[active.LR==0 & truth==0] <- 4 #TN    
    
    emap.mLR <- ifelse(active.mLR==1 & truth==1, 1, NA) #TP
    emap.mLR[active.mLR==1 & truth==0] <- 2 #FP  
    emap.mLR[active.mLR==0 & truth==1] <- 3 #FN
    emap.mLR[active.mLR==0 & truth==0] <- 4 #TN 
    
    emap.fdr <- ifelse(active.fdr==1 & truth==1, 1, NA) #TP
    emap.fdr[active.fdr==1 & truth==0] <- 2 #FP  
    emap.fdr[active.fdr==0 & truth==1] <- 3 #FN
    emap.fdr[active.fdr==0 & truth==0] <- 4 #TN  
    
    emap.LRdel <- ifelse(active.LRdel==1 & truth==1, 1, NA) #TP
    emap.LRdel[active.LRdel==1 & truth==0] <- 2 #FP  
    emap.LRdel[active.LRdel==0 & truth==1] <- 3 #FN
    emap.LRdel[active.LRdel==0 & truth==0] <- 4 #TN    
    
    tmp.results[i,1] <- sum(emap.bonf[!is.na(emap.bonf)]==2)/inactvox
    tmp.results[i,2] <- sum(emap.bonf[!is.na(emap.bonf)]==3)/actvox    
    tmp.results[i,3] <- (tmp.results[i,1] + tmp.results[i,2])/2  
    tmp.results[i,4] <- sum(emap.LR[!is.na(emap.LR)]==2)/inactvox
    tmp.results[i,5] <- sum(emap.LR[!is.na(emap.LR)]==3)/actvox    
    tmp.results[i,6] <- (tmp.results[i,4] + tmp.results[i,5])/2
    tmp.results[i,7] <- sum(emap.fdr[!is.na(emap.fdr)]==2)/inactvox
    tmp.results[i,8] <- sum(emap.fdr[!is.na(emap.fdr)]==3)/actvox    
    tmp.results[i,9] <- (tmp.results[i,7] + tmp.results[i,8])/2
    tmp.results[i,10] <- sum(emap.mLR[!is.na(emap.mLR)]==2)/inactvox
    tmp.results[i,11] <- sum(emap.mLR[!is.na(emap.mLR)]==3)/actvox    
    tmp.results[i,12] <- (tmp.results[i,10] + tmp.results[i,11])/2
    tmp.results[i,13] <- sum(emap.LRdel[!is.na(emap.LRdel)]==2)/inactvox
    tmp.results[i,14] <- sum(emap.LRdel[!is.na(emap.LRdel)]==3)/actvox    
    tmp.results[i,15] <- (tmp.results[i,13] + tmp.results[i,14])/2
  }
  
  counts <- c(delta, nscan, perc, mean(tmp.results[,1],na.rm=TRUE), mean(tmp.results[,2],na.rm=TRUE),
   mean(tmp.results[,3],na.rm=TRUE), mean(tmp.results[,4],na.rm=TRUE), mean(tmp.results[,5],na.rm=TRUE),
   mean(tmp.results[,6],na.rm=TRUE), mean(tmp.results[,7],na.rm=TRUE), mean(tmp.results[,8],na.rm=TRUE),
   mean(tmp.results[,9],na.rm=TRUE), mean(tmp.results[,10],na.rm=TRUE), mean(tmp.results[,11],na.rm=TRUE),
   mean(tmp.results[,12],na.rm=TRUE), mean(tmp.results[,13],na.rm=TRUE), mean(tmp.results[,14],na.rm=TRUE),
   mean(tmp.results[,15],na.rm=TRUE))

  
  resultlabels <- c("delta","nscan","perc","FP.bon","FN.bon","AVG.bon","FP.LR","FN.LR","AVG.LR","FP.fdr","FN.fdr","AVG.fdr","FP.mlr","FN.mlr","AVG.mlr","FP.lrdel","FN.lrdel","AVG.lrdel")
  results <- array(NA, dim=c(1,length(resultlabels)))
  results[1,] <- counts
  colnames(results) <- resultlabels

  write.csv(results, file=paste(datastor,conditionline,".csv",sep=""))
