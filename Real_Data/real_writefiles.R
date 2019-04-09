##############################################################
##    READ IN DATA WRITE AWAY OVERALL MASK FOR SUBJECT      ##
##############################################################

#Working directories
data <- "PATH1" #path where copes and varcopes are stored

#Libraries
library(oro.nifti)
library(AnalyzeFMRI)
library(lattice)
library(qvalue)
library(locfdr)

## OVERALL BRAIN MASK

setwd(data)
cope <- readNIfTI("cope_1.nii",reorient=FALSE)
mask <- ifelse(cope!=0, 1, 0)

for (i in 2:10) {

	cope <- readNIfTI(paste("cope_",i,".nii",sep=""),reorient=FALSE)
	cope <- ifelse(cope!=0, 1, 0)
	mask <- mask + cope
}

mask <- ifelse(mask==10, 1, 0)
f.write.nifti(mask,"real_OVERALLMASK",size="float", nii=TRUE)
