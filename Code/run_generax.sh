generate_generax_mapfile () {
		# https://github.com/BenoitMorel/GeneRax/wiki/Gene-to-species-mapping
		my_aln_file=$1
		name=$(echo ${my_aln_file} | sed 's/\/OG[0-9]*_cds_hammer.fa//')
		cat ${my_aln_file} | grep "^>" | sed -e "s/>//" > $name/$name"_tmp.gene_names.txt"
		cat $name/$name"_tmp.gene_names.txt" | sed -e "s/_/|/" -e "s/_.*//" -e "s/|/_/" > $name/$name"_tmp.species_names.txt"
		paste $name/$name"_tmp.gene_names.txt" $name/$name"_tmp.species_names.txt" > $name/$name"_generax_map.txt"
		rm $name/$name"_tmp.gene_names.txt" $name/$name"_tmp.species_names.txt"
	}
	generate_generax_mapfile $1


echo """
	[FAMILIES]
	- family_1
	starting_gene_tree = $name/$name"_cds_hammer.fa.treefile"
	alignment = $name/$name"_cds_hammer.fa"
	mapping = $name/$name"_generax_map.txt"
	subst_model = GTR+G4
	""" | sed -e "s/^[[:space:]]*//" | grep -v "^$" > $name/$name"_generax_families.txt"


mpiexec -n 20 --oversubscribe generax \
--species-tree tree.newick \
--families $name/$name"_generax_families.txt" \
--strategy SPR \
--reconcile \
--rec-model UndatedDL \
--prefix $name/$name"_generax" \
--per-family-rates \
--seed 12345
