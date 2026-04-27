#!/bin/bash --login

#SBATCH --job-name fasterq_dump
#SBATCH --partition=cpu
#SBATCH --array=1-21%3
#SBATCH --nodes=1
#SBATCH --ntasks=1
#SBATCH --cpus-per-task=4
#SBTACH --mem-per-cpu=2G
#SBATCH --output=kal.%A_%a.out
#SBATCH --error=kal.%A_%a.err
#SBATCH --time 0-24:00:00

module load gcc sratoolkit/3.0.0

#Bulk download (adjust bash array job parameters for bulk)

name=$(sed -n "$SLURM_ARRAY_TASK_ID"p namelist_sra.txt)

prefetch $name
fasterq-dump $name
