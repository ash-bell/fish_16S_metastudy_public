for j in $(cat BIOPROJECT_list.txt); do
sed "s/BIOPROJECT/$j/g;s|WKDIR|/lustre/home/bt273/BIOSSCOPE/metag/ashley/fish_16S|g" 03_align_V4_blast.SLURM > $j/tmp.SLURM;
sbatch $j/tmp.SLURM;
done
