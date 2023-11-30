for j in $(cat BIOPROJECT_list.txt); do
sed "s/BIOPROJECT/$j/g;s|WKDIR|/lustre/home/bt273/BIOSSCOPE/metag/ashley/fish_16S|g" 04_trim.ASV_assign.taxa_prevelence.filter.R > $j/${j}_trim.ASV_assign.taxa_prevelence.filter.R;
sed "s/BIOPROJECT/$j/g;s|WKDIR|/lustre/home/bt273/BIOSSCOPE/metag/ashley/fish_16S|g" 04_SLURM_script.SLURM > $j/tmp.SLURM;
sbatch $j/tmp.SLURM;
done
