#!/bin/bash
#PBS -j oe
#PBS -N train_model
#PBS -l walltime=15:00:00
#PBS -l select=1:ncpus=40:mem=120gb

cd $PBS_O_WORKDIR
shopt -s expand_aliases
source /etc/profile.d/modules.sh
echo "Job identifier is $PBS_JOBID"
echo "Working directory is $PBS_O_WORKDIR"

module load R/3.6.1
Rscript train_model.R outfile=Mus_musculus_model.rds train=Mus_musculus.rds ncores=70