#!/bin/bash --login

#SBATCH --job-name=CAFE5
#SBATCH --partition=cpu
#SBATCH --array=1-1%1
#SBATCH --nodes=1
#SBATCH --ntasks=1
#SBTACH --mem-per-cpu=2G
#SBATCH --cpus-per-task=4
#SBATCH --time=72:00:00
#SBATCH --output=cafeout_%A_%a.out
#SBATCH --error=caferror_%A_%a.err


module load gcc 

conda activate ./cafe5_env/
#Array indexing
file=$(sed -n "$SLURM_ARRAY_TASK_ID"p ./sample_list.txt) #Change array index for the number of lines here

res="results_cafe_"$file

#Make all the required directories
mkdir ./cafeOut/$res
mkdir ./cafeOut/$res/Base_model
mkdir ./cafeOut/$res/Multi_lambda_model_1
mkdir ./cafeOut/$res/Multi_lambda_model_2
mkdir ./cafeOut/$res/Multi_lambda_model_3

#base model and error estimation
cafe5 -i ./new_dat/"cafe_input_"$file".tsv" -t ./new_dat/cafe_no_WGD_tree.nwk -p -e -k 1 -o ./cafeOut/$res/Base_model

#multi labmda model; tests for difference in rates in acanthomorpha
cafe5 -i ./new_dat/"cafe_input_"$file".tsv" -t ./new_dat/cafe_no_WGD_tree.nwk -y ./new_dat/cafe_no_WGD_rates_tree.nwk -p -e./cafeOut/$res/Base_model/Base_error_model.txt -o ./cafeOut/$res/Multi_lambda_model_1

cafe5 -i ./new_dat/"cafe_input_"$file".tsv" -t ./new_dat/cafe_no_WGD_tree.nwk -y ./new_dat/cafe_no_WGD_rates_tree.nwk -p -e./cafeOut/$res/Base_model/Base_error_model.txt -o ./cafeOut/$res/Multi_lambda_model_2

cafe5 -i ./new_dat/"cafe_input_"$file".tsv" -t ./new_dat/cafe_no_WGD_tree.nwk -y ./new_dat/cafe_no_WGD_rates_tree.nwk -p -e./cafeOut/$res/Base_model/Base_error_model.txt -o ./cafeOut/$res/Multi_lambda_model_3
