for i in $(cat BIOPROJECT_list.txt); do mkdir -p $i; done

for j in $(cat BIOPROJECT_list.txt); do
for i in $(cat $j/${j}_SRA_Acc_list.txt); do
sed "s/BIOPROJECT/$j/g;s/ACCESSION_NUMBER/$i/g;s|WKDIR|/lustre/home/bt273/BIOSSCOPE/metag/ashley/fish_16S|g" 01_get_seq_NCBI.SLURM > $j/tmp.SLURM;
sbatch $j/tmp.SLURM;
done;
done
for i in PRJ*metadata*; do sed -i '1s/^\xEF\xBB\xBF//' $i ; done # removes byte order mark from files
