#!/bin/bash --login

#SBATCH --job-name kal_quant
#SBATCH --partition=cpu
#SBATCH --array=1-1%1
#SBATCH --nodes=1
#SBATCH --ntasks=1
#SBATCH --cpus-per-task=4
#SBTACH --mem-per-cpu=2G
#SBATCH --output=kal_quant.%A_%a.out
#SBATCH --error=kal_quant.%A_%a.err
#SBATCH --time 0-48:00:00

module load gcc  kallisto/0.48.0
#Array indexing
file=$(sed -n "$SLURM_ARRAY_TASK_ID"p namelist_sra.txt)


#cd /work/FAC/FBM/DEE/mrobinso/default/Agneesh/Projects/Metamorphosis_project/00.Download_data/${sp}
#cd /work/FAC/FBM/DEE/mrobinso/default/Agneesh/Projects/Metamorphosis_project/00.Download_data/Goby/RhinoBB/gigafile-1015-da782e73a601cb0629ffa1d84c239dc3


#for F1 in *; do
    #F2=$(echo $F1 | sed 's/_1.fastq.gz/_2.fastq.gz/')
    #ACC=$(echo $F1 | sed 's/_1.fastq.gz//')


for F1 in *_1.fastq; do    
       	F2=$(echo $F1 | sed 's/_1.fastq/_2.fastq/') 
	ACC=$(echo $F1 | sed 's/_1.fastq//') ;     
       	
#	echo ${F1};     
#	echo ${F2}; 
#	echo ${ACC};  


   kallisto quant -t 4 -i clown.idx -o ${ACC}"_kal_quant" -b 100 ${F1} ${F2}
   
   mv ${ACC}"_kal_quant" clown_res;

done
