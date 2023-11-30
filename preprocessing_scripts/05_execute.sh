sed "s|WKDIR|/lustre/home/bt273/BIOSSCOPE/metag/ashley/fish_16S|g" 05_merge_phyloseq.R > merge_phyloseq.R;
sed "s|WKDIR|/lustre/home/bt273/BIOSSCOPE/metag/ashley/fish_16S|g" 05_SLURM_script.SLURM > merge_phyloseq.SLURM;
sbatch merge_phyloseq.SLURM
