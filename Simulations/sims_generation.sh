####################################################################################
###### Script to start R script with information about which line			  ######
####################################################################################

#!/bin/sh

#PBS -o outputSims/outputSims.file

#PBS -e errorSims/errorSims.file

#PBS -l walltime=71:00:00

#PBS -l nodes=1:ppn=1

#PBS -l vmem=10GB

module load R

Rscript pathtoscripts/simulations.R $PBS_ARRAYID

echo "job finished"

