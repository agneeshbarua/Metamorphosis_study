#!/bin/bash --login

#SBATCH --job-name kal_index
#SBATCH --partition=cpu
#SBATCH --array=1-1%2
#SBATCH --nodes=1
#SBATCH --ntasks=1
#SBATCH --cpus-per-task=4
#SBTACH --mem-per-cpu=2G
#SBATCH --output=kal.%A_%a.out
#SBATCH --error=kal.%A_%a.err
#SBATCH --time 0-24:00:00

module load gcc kallisto/0.48.0


#Array indexing
sp=$(sed -n "$SLURM_ARRAY_TASK_ID"p sp_names.txt)

#kallisto index -i ${sp}".idx" /work/FAC/FBM/DEE/mrobinso/default/Agneesh/Projects/Metamorphosis_project/00.Download_data/${sp}/${sp}"_cds.fa" --make-unique

kallisto index -i clown.idx  ./Amphiprion_ocellaris.ASM2253959v1.cds.all.fa --make-unique
