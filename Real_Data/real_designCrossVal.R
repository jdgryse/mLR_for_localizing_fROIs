####################
#### TITLE:     Create the design.mat, design.con and design.grp files
#### Contents: 	
#### 
#### Source Files: //ABT_for_localizing_fROIs
#### First Modified: 2016
#### Notes: 
#################

##
###############
### Notes
###############
##

# We also need a mask, so we read in a contrast file and put all non zero elements to 1

##
###############
### Preparation
###############
##

# Library
library(oro.nifti)
	# Print which libraries are needed
	print("Need package: oro.nifti")

setwd("PATH1") #path where your copes etc are stored

##
###############
### Write design.mat
###############
##

# File connection
fileCon <- "design.mat"
# Text to be written to the file
cat('/NumWaves\t1
/NumPoints\t',paste(9,sep=''),'
/PPheights\t\t1.000000e+00

/Matrix
',rep("1.000000e+00\n",9),file=fileCon)


##
###############
### Write design.con
###############
##

fileCon <- file("design.con")
	writeLines('/ContrastName1	Group Average 
/NumWaves	1
/NumContrasts	1
/PPheights		1.000000e+00
/RequiredEffect		5.034

/Matrix
1.000000e+00 
',fileCon)
close(fileCon)


##
###############
### Write design.grp
###############
##

# File connection
fileCon <- "design.grp"
# Text to be written to the file
cat('/NumWaves\t1
/NumPoints\t',paste(9,sep=''),'

/Matrix
',rep("1\n",9),file=fileCon)


##
###############
### Create masks
###############
##

#Read in the cope image: all copes from each subject are read in.
#	# Then we take the union of the voxels being activated in all subjects. 
# 		# First process the cope images such that:
# 			# Voxels with any value get 1.
# 			# Voxels with NaN get 0.
 	cope[which(cope!=0)] <- 1
 	cope[which(is.nan(cope)==TRUE)] <- 0










