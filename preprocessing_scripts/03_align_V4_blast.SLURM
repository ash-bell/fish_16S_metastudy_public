#!/bin/bash
#SBATCH --job-name=BIOPROJECT_blast
#SBATCH --export=ALL
#SBATCH --chdir=.
#SBATCH --partition=pq # submit to the parallel queue
#SBATCH --time=168:00:00 # maximum walltime for the job
#SBATCH --account=Research_Project-###### # research project to submit under
#SBATCH --nodes=1 # specify number of nodes
#SBATCH --ntasks-per-node=16 # specify number of processors per node
#SBATCH -o BIOPROJECT_V4_blastn.log

source ~/.bashrc

conda activate dada2

cd WKDIR

# If the blastdb file doesn't exist make it, else don't make the blastdb file
# Always make V4 blast database from PRJEB46984 as it is confirmed to only include the V4 region
#[ ! -f blastDB_V4.ndb ] && makeblastdb -in PRJEB46984/PRJEB46984_ASVs.fna -dbtype nucl -out blastDB_V4
blastn -query BIOPROJECT/BIOPROJECT_ASVs.fna -db blastDB_V4 -outfmt '6 qseqid sseqid qstart qend sstart send' -num_threads 16 -out BIOPROJECT/BIOPROJECT_blastout.txt -word_size 4 -max_target_seqs 50
