for j in $(cat BIOPROJECT_list.txt); do
sed "s/BIOPROJECT/$j/g;s|WKDIR|/lustre/home/bt273/BIOSSCOPE/metag/ashley/fish_16S|g" 02_dada2_filter.R > $j/${j}_dada2_process.R;
sed "s/BIOPROJECT/$j/g;s|WKDIR|/lustre/home/bt273/BIOSSCOPE/metag/ashley/fish_16S|g" 02_SLURM_script.SLURM > $j/tmp.SLURM;
sbatch $j/tmp.SLURM;
done

