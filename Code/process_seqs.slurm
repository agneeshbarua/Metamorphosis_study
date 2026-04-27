#!/bin/bash --login

#SBATCH --job-name csubst_pipeline_x02
#SBATCH --partition=cpu
#SBATCH --array=1-1%1
#SBATCH --nodes=1
#SBATCH --ntasks=1
#SBATCH --cpus-per-task=5
#SBTACH --mem-per-cpu=2G
#SBATCH --output=csubst.%A_%a.out
#SBATCH --error=csubst.%A_%a.out
#SBATCH --time 0-72:00:00

#Environment setup
module load gcc mafft/7.505 trimal/1.4.1 
conda activate ~/myproject_envs/Pigmentation_env/

#Array indexing
name=$(sed -n "$SLURM_ARRAY_TASK_ID"p file_list.txt)

rm -r Procesed_dat/${name}
mkdir Procesed_dat/${name}

#make nucleotide sequences in-frame by head and tail paddings
~/seqkit sort -n OG_SEQS/cds_seqs/${name}"_cds.fa" | cdskit pad --codontable 1 > Procesed_dat/${name}/${name}"_cds_pad.fa"

#translate cds into protein sequences
transeq -sequence Procesed_dat/${name}/${name}"_cds_pad.fa" -outseq Procesed_dat/${name}/${name}"_cds_trans.fa" -frame 1 -table 0

#remove the '_1' from header made by transeq
sed -i 's/_1$//' Procesed_dat/${name}/${name}"_cds_trans.fa"

#align the protein sequences
mafft --auto --amino Procesed_dat/${name}/${name}"_cds_trans.fa" > Procesed_dat/${name}/${name}"_cds_trans_pep_aligned.fa"

#translate align the protein alignment with unaigned cds
tranalign -table 0 -asequence Procesed_dat/${name}/${name}"_cds_pad.fa" -bsequence Procesed_dat/${name}/${name}"_cds_trans_pep_aligned.fa" -outseq Procesed_dat/${name}/${name}"_cds_tranaligned.fa"

#trim poorly aligned codon sites
trimal -in Procesed_dat/${name}/${name}"_cds_tranaligned.fa" -backtrans Procesed_dat/${name}/${name}"_cds_pad.fa" -out Procesed_dat/${name}/${name}"_cds_trimal_aligned.fa" -ignorestopcodon -automated1

#mask ambiguous and/or stop codons
cdskit mask --seqfile Procesed_dat/${name}/${name}"_cds_trimal_aligned.fa" --outfile Procesed_dat/${name}/${name}"_cds_mask_tranalign.fa"

#remove uninformative sites
cdskit hammer --seqfile Procesed_dat/${name}/${name}"_cds_mask_tranalign.fa" --outfile Procesed_dat/${name}/${name}"_cds_hammer.fa" --nail 4 --codontable 1

#remove problematic sequences
cdskit rmseq -s Procesed_dat/${name}/${name}"_cds_hammer.fa" --problematic_percent 60 -o Procesed_dat/${name}/${name}"_cds_treeinput.fa"

#iqtree to construct gene trees
iqtree -s Procesed_dat/${name}/${name}"_cds_treeinput.fa" -m MFP -bb 1000 -nt AUTO -bnni --prefix Procesed_dat/${name}/${name}"_cds_treeinput.fa" -wsr -asr

#run generax
bash ./run_generax.sh Procesed_dat/${name}/$name"_cds_treeinput.fa"

#convert nhx format to newick format
nwkit nhx2nwk --infile Procesed_dat/${name}/$name"_generax"/reconciliations/family_1_reconciliated.nhx --outfile Procesed_dat/${name}/$name"_generax.nwk" 
nwkit drop --infile Procesed_dat/${name}/$name"_generax.nwk" --target intnode --name yes --outfile Procesed_dat/${name}/$name"_generax_complete.nwk"

#Remove intermediate files
rm Procesed_dat/${name}/${name}"_cds_pad.fa"
rm Procesed_dat/${name}/${name}"_cds_trans.fa"
rm Procesed_dat/${name}/${name}"_cds_trans_pep_aligned.fa"
rm Procesed_dat/${name}/${name}"_cds_tranaligned.fa"
rm Procesed_dat/${name}/${name}"_cds_trimal_aligned.fa"
rm Procesed_dat/${name}/${name}"_cds_mask_tranalign.fa"
rm Procesed_dat/${name}/${name}"_cds_hammer.fa"
rm Procesed_dat/${name}/$name"_generax.nwk" 
