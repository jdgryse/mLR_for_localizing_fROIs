# Welcome

1)  This repository holds both code and data files that can be used across platforms (linux, Windows, macOS) to perform all methods mentioned in the original paper. The results of the paper are obtained with the following operating system: macOS Mojave, version 10.14.4.
2)  - All .R files contain code that should be run with the statistal program R. In order to run this code, the user should install R or RStudio. The version of R used for the results in the paper is R 3.5.3. Once R is installed, the code can readily be run in the console, with the exception for the paths where data are stored locally, which should be adjusted accordingly,  and additional R packages that may not already be installed locally.
    - All .sh files are shell scripts (Unix language), with code that can be run in a Bash shell. In order for most
       commands to work, FSL should be installed on your computer. The results in this paper were obtained with FSL 5.0.9.

# Simulations for the maximized likelihood ratio paradigm for the definition of fROIs

## simulations.R

In this file, the subject maps are simulated analyzed and thresholded using NHST, the likelihood ratio (LR) and the  maximized likelihood ratio (mLR) method. Their Type I, Type II and average error rates are computed. This is described in Section X of the paper. The resulting error rates are used to produce Figure X in the paper.

## make_conditions.R

This file makes a table with all parameter configurations needed to reproduce the simulation study of the paper in R. All possible values for the parameters that were varied in the simulation study (the value for Delta for the mLR method, the number of scans and the specific quantile used to estimate the effect size define in the alternative for the LR method) are combined. The resulting table is conditions.txt. This table is used in simulations.R.


## functions_GLM_2D.R

This file contains the function to estimate the parameters of the GLM, their standard error and the tmap needed to perform the analyses done in simulations.R. No adjustments are required for this file. 



# Real data example and evaluation of the maximized likelihood ratio method for the definition of fROIs


## Data

A toy data set in order to perform the cross-validation performed in the original study, as well as an overall brain mask for this “simulated individual”, since the data described in Section X of the original paper are not freely available. These are all NIFTI files that can be read into R or FSL or the fMRI data analysis program of your choice. 

## design.

These files contain information concerning the design of the fixed effects analysis in FSL, needed for the flameo command in real_es_groundtruth.sh and do not need to be adjusted in order to work for the toy data set.

## real_analysis_LR_mLR.R

In this file, the LR and mLR method are performed on the left out run under different parameter configurations, which is described in Section X of the paper. The result of this file are SPMs with either LRs or mLRs. The resulting SPMs were then used to produce Figure X in the paper.

## real_designCrossVal.R

Here, the design. files mentioned above are made.

## real_es_groundtruth.sh

Here, the ground truth of effect sizes for each step of the cross-validation is constructed. The ground truth of effect sizes is then used in real_analysis_LR_mLR.R.

## real_writefiles.R

In this file, the overall brain mask for the “simulated individual” is constructed. This is then used in  real_analysis_LR_mLR.R.
