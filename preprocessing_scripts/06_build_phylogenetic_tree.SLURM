#!/bin/bash
#SBATCH --job-name=RAxML
#SBATCH --export=ALL
#SBATCH --chdir=.
#SBATCH --partition=pq # submit to the parallel queue
#SBATCH --time=168:00:00 # maximum walltime for the job
#SBATCH --account=Research_Project-###### # research project to submit under
#SBATCH --nodes=1 # specify number of nodes
#SBATCH --ntasks-per-node=16 # specify number of processors per node
#SBATCH -o build_phylo_tree.log

source ~/.bashrc

# originally use raxml-ng but took too long so switched to IQ-TREE
conda activate raxml-ng

cd /lustre/projects/Research_Project-172179/metag/ashley/fish_16S/

mafft --auto --thread 16 all_studies_ASVs.fna > all_studies_ASVs.msa
modeltest-ng -i all_studies_ASVs.msa -o model_test_output.txt -p 16
iqtree -s all_studies_ASVs.msa -nt AUTO -m TPM3uf+I+G4 -pre TPM3uf_IQTREE
