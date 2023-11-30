#!/bin/bash
#SBATCH --job-name=merge_PS
#SBATCH --export=ALL
#SBATCH --chdir=.
#SBATCH --partition=pq # submit to the parallel queue
#SBATCH --time=168:00:00 # maximum walltime for the job
#SBATCH --account=Research_Project-###### # research project to submit under
#SBATCH --nodes=1 # specify number of nodes
#SBATCH --ntasks-per-node=16 # specify number of processors per node
#SBATCH -o merge_phyloseq.log
#SBATCH -p phi

source ~/.bashrc

conda activate dada2

cd WKDIR

Rscript merge_phyloseq.R
