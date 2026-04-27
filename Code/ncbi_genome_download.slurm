#!/bin/bash --login

#SBATCH --job-name ncbi_datasets
#SBATCH --partition=cpu
#SBATCH --array=1-169%20
#SBATCH --nodes=1
#SBATCH --ntasks=1
#SBATCH --cpus-per-task=2
#SBTACH --mem-per-cpu=2G
#SBATCH --output=ncbi_datasets.%A_%a.out
#SBATCH --error=ncbi_datasets.%A_%a.out
#SBATCH --time 0-24:00:00

module load gcc miniconda3 mafft/7.481 trimal/1.4.1 emboss/6.6.0

conda activate ~/myproject_envs/Pigmentation_env/

#Array indexing
acc=$(awk '{print $1}' sp_acc_list.tsv | sed -n "$SLURM_ARRAY_TASK_ID"p)
sp=$(awk '{print $2}' sp_acc_list.tsv | sed -n "$SLURM_ARRAY_TASK_ID"p)

just_dat_dir="/work/FAC/FBM/DEE/mrobinso/evolseq/Agneesh/Projects/Metamorphosis/00.Orthfinder/updated_sp_list"
prot_dir="./updated_sp_list/Protein_dat"

mkdir $prot_dir/${sp}

#### download protein data from ncbi #####
./datasets download genome accession ${acc} --include protein --filename $prot_dir/${sp}/${sp}".zip"
unzip $prot_dir/${sp}/${sp}."zip" -d $prot_dir/${sp}
mv $prot_dir/${sp}/ncbi_dataset/data/${acc}/protein.faa $prot_dir/${sp}/${sp}"_protein.fa"
cp $prot_dir/${sp}/${sp}"_protein.fa" $just_dat_dir
rm $prot_dir/${sp}/README.md
rm -r $prot_dir/${sp}/ncbi_dataset
