#!/bin/bash --login

#SBATCH --job-name HyPhy_FUBAR
#SBATCH --partition=cpu
#SBATCH --array=1-1%1
#SBATCH --nodes=1
#SBATCH --ntasks=1
#SBATCH --cpus-per-task=2
#SBTACH --mem-per-cpu=6G
#SBATCH --output=%A_%a.out
#SBATCH --error=%A_%a.err
#SBATCH --time 72:00:00

module load gcc

conda activate ~/myproject_envs/hyphy-2.5.6x/

#Array indexing
og=$(sed -n "$SLURM_ARRAY_TASK_ID"p non_homologous_prefix_06)

#make process dir
rm -r ${og}
mkdir ${og}

#copy alignment
cp /work/FAC/FBM/DEE/mrobinso/evolseq/Agneesh/Projects/Metamorphosis/04.NCBI_datasets/Seq_downloads/Procesed_dat/${og}/${og}"_cds_treeinput.fa" ${og}
#copy tree
cp /work/FAC/FBM/DEE/mrobinso/evolseq/Agneesh/Projects/Metamorphosis/04.NCBI_datasets/Seq_downloads/Procesed_dat/${og}/${og}"_generax_complete.nwk" ${og}

#analysis code
hyphy fubar \
  --alignment ${og}/${og}"_cds_treeinput.fa" \
  --tree ${og}/${og}"_generax_complete.nwk" 

#move output
mv ${og}/${og}"_cds_treeinput.fa.FUBAR.json" ../hyphy_fubar_outputs

rm -r ${og}
