#!/bin/bash
#SBATCH --job-name=BIOPROJECT_ACCESSION_NUMBER_DL
#SBATCH --export=ALL
#SBATCH --chdir=.
#SBATCH --partition=pq # submit to the parallel queue
#SBATCH --time=001:00:00 # maximum walltime for the job
#SBATCH --account=Research_Project-###### # research project to submit under
#SBATCH --nodes=1 # specify number of nodes
#SBATCH --ntasks-per-node=1 # specify number of processors per node
#SBATCH -o logfiles/BIOPROJECT_ACCESSION_NUMBER_DL.log

source ~/.bashrc

conda activate bbmap

mkdir WKDIR/BIOPROJECT
mkdir WKDIR/BIOPROJECT/reads

cd WKDIR/BIOPROJECT/reads

prefetch ACCESSION_NUMBER
vdb-validate ACCESSION_NUMBER
fastq-dump --split-files ACCESSION_NUMBER

#Trim adapters, synthetic sequencing artifacts and phix spike-ins by kmer-matching
if [ ! -f ACCESSION_NUMBER_2.fastq ]; 
then 
    bbduk.sh in=ACCESSION_NUMBER_1.fastq out=BIOPROJECT_ACCESSION_NUMBER_#.trm.fq.gz ktrim=r k=23 mink=11 hdist=1 tbo=t tpe=t minlen=70 ref=adapters,artifacts,phix threads=1 maxns=0 overwrite=t
else 
    bbduk.sh in=ACCESSION_NUMBER_#.fastq out=BIOPROJECT_ACCESSION_NUMBER_#.trm.fq.gz ktrim=r k=23 mink=11 hdist=1 tbo=t tpe=t minlen=70 ref=adapters,artifacts,phix threads=1 maxns=0 overwrite=t
fi

#Quaility trimming
bbduk.sh in=BIOPROJECT_ACCESSION_NUMBER_#.trm.fq.gz out=BIOPROJECT_ACCESSION_NUMBER_#.QC.trm.fq.gz qtrim=rl trimq=10 maq=10 minlen=70 threads=1

rm -r ACCESSION_NUMBER_1.fastq ACCESSION_NUMBER_2.fastq ACCESSION_NUMBER/ BIOPROJECT_ACCESSION_NUMBER_1.trm.fq.gz BIOPROJECT_ACCESSION_NUMBER_2.trm.fq.gz
