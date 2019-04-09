# Simulations for the maximized likelihood ratio paradigm for the definition of fROIs

## simulations.R

In this file, the subject maps are simulated and analyzed using NHST, the LR and the mLR method. Their Type I, Type II and average error rates are computed.

## sims_generation.sh

This file is used to submit to the linux-based supercomputer in order to simulate the brain maps.


## make_conditions.R

This file makes a table with all parameter configurations needed in simulations.R. The resulting table is conditions.txt


## functions_GLM_2D.R

This file contains the function to estimate the parameters of the GLM, their standard error and the tmap



# Real data example and evaluation of the maximized likelihood ratio method for the definition of fROIs

## Data

A toy data set in order to perform the cross-validation performed in the original study, as well as an overall brain mask for this “simulated individual”.

## design.

These files contain information concerning the design of the fixed effects analysis in FSL, needed for the flameo command in real_es_groundtruth.sh

## real_analysis_LR_mLR.R

In this file, the LR and mLR method are performed on the left out run under different parameter configurations.

## real_designCrossVal.R

Here, the design files mentioned above are made.

## real_es_groundtruth.sh

Here, the ground truth of effect sizes for each step of the cross-validation is constructed.

## real_writefiles.R

In this file, the overall brain mask for the “simulated individual” is constructed.
