
################################################################
##                ANALYSIS WITH LR and mLR                    ##
################################################################



## PRELIMINARY
args <- commandArgs(TRUE)
deltaline <- as.numeric(as.character(args))[1]
stor <- as.numeric(as.character(args))[2] #To use the copes and varcopes of 2,3,4 or 5 blocks combined
#Relevant directories

data <- "pathtocopesandvarcopes"
datastor <- "pathtostoreanalyzedSPMs"

#Libraries
library(neuRosim)
library(oro.nifti)
library(AnalyzeFMRI)
library(lattice)


##COMPUTING LR AND MLR

for(i in 1:10) {

  print(i)

  #Delta calculation
  EStruth <- readNIfTI(paste("pathtogroundtruthsofeffecsizes/truth_ES_",i,".nii.gz",sep=""))
  mask <- readNIfTI(paste("pathtosubjectmask/subj_1_OVERALLMASK.nii",sep=""))
  mask <- ifelse(mask==0, NA, mask)

  truth <- EStruth * mask

  es1 <- quantile(truth, 0.60, na.rm=TRUE)
  es2 <- quantile(truth, 0.75,na.rm=TRUE)
  es3 <- quantile(truth, 0.90,na.rm=TRUE)


  deltavec <- c(es1,es2,es3)
  deltaperc <- c(600, 750, 900)


  ## READING IN B1, S(B1), TRUTH
  b1 <- readNIfTI(paste(data,"cope_",i,".nii.gz",sep=""),reorient=FALSE)
  var1 <- readNIfTI(paste(data,"varcope_",i,".nii.gz",sep=""),reorient=FALSE)
  sb1 <- sqrt(var1)
  mask <- ifelse(b1 ==0, NA, 1)
  b1mask <- b1 * mask
  sb1mask <- sb1 * mask


  ## LR AND MLR

  delta <- deltavec[deltaline]
  
  #LR + delta
  LR <- array(0, dim=dim(b1))

  teller<-dnorm(b1mask, mean=delta, sd=sb1)
  noemer<-dnorm(b1mask, mean=0, sd=sb1)

  LR <- teller/noemer
  LR[is.na(LR)] <- 0


  #Maximized LR method + delta
  b1stats <- c(b1)
  sb1stats <- c(sb1)

  teller_mLRdelta <- c()
  noemer_mLRdelta <- c()

  for(l in 1:length(b1stats)) {
  
    if (b1stats[l] > delta) {teller_mLRdelta[l] <- dnorm(b1stats[l], mean=b1stats[l], sd=sb1stats[l])}
    if (b1stats[l] <= delta) {teller_mLRdelta[l] <- dnorm(b1stats[l], delta, sd=sb1stats[l])}
    if (b1stats[l] < delta) {noemer_mLRdelta[l]<-dnorm(b1stats[l], mean=b1stats[l], sd=sb1stats[l])}
    if (b1stats[l] >= delta){noemer_mLRdelta[l]<-dnorm(b1stats[l], delta, sd=sb1stats[l])}

  }

  mLRdelta <-teller_mLRdelta/noemer_mLRdelta
  mLRdelta <- array(mLRdelta, dim=dim(b1))
 

  ## WRITING AWAY THE DATA

  writeNIfTI(LR, paste(datastor,"LR_delta_",deltaperc[deltaline],"_",i,sep=""),gzipped=TRUE)
  writeNIfTI(mLRdelta, paste(datastor,"mLR_delta_",deltaperc[deltaline],"_",i,sep=""),gzipped=TRUE)

}

