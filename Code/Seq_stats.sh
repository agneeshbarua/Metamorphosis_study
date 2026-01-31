for f in $ll ../04.NCBI_datasets/Seq_downloads/Procesed_dat/*/*.fa; do ~/seqkit stat $f >> seq_stats.txt; done

sed '/^file[[:space:]]\+/d' seq_stats.txt | awk 'BEGIN{OFS="\t"} NF>=8 {
  n = split($1,a,"/")
  id = a[n-1]
  printf "%s", id
  for(i=2;i<=NF;i++) printf "%s%s", OFS, $i
  printf "\n"
}' | sed 's/,//g' > seq_dat_clean.tsv