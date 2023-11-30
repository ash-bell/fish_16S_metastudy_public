#!/bin/bash
#SBATCH --job-name=BIOPROJECT_DADA2
#SBATCH --export=ALL
#SBATCH --chdir=.
#SBATCH --partition=pq # submit to the parallel queue
#SBATCH --time=168:00:00 # maximum walltime for the job
#SBATCH --account=Research_Project-###### # research project to submit under
#SBATCH --nodes=1 # specify number of nodes
#SBATCH --ntasks-per-node=16 # specify number of processors per node
#SBATCH -o BIOPROJECT_trim.ASV_assign.taxa_prevelence.filter.log

source ~/.bashrc

conda activate dada2

cd WKDIR/BIOPROJECT

Rscript BIOPROJECT_trim.ASV_assign.taxa_prevelence.filter.R
