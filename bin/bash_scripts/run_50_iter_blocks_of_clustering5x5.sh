
#!/bin/bash
for i in `seq 1 25`;
do

python bin/python_scripts/som_to_kmeans_clustering.py --som_x 5 --som_y 5 --out_dir data/2014-1-5.1_clusters_from_5x5-som_to_kmeans/ --input_file data/2013-12-10.4_vectors_cleaned_laccases_vec/lac_pt.faa.pruneBJXZUO.fasta.SequenceOrderCouplingNumberTotal.EL_DESC_VAL.vec
python bin/python_scripts/som_to_kmeans_clustering.py --som_x 5 --som_y 5 --out_dir data/2014-1-5.1_clusters_from_5x5-som_to_kmeans/ --input_file data/2013-12-10.4_vectors_cleaned_laccases_vec/lac_pt.faa.pruneBJXZUO.fasta.shenC-shenT-shenD.EL_DESC_VAL.vec
python bin/python_scripts/som_to_kmeans_clustering.py --som_x 5 --som_y 5 --out_dir data/2014-1-5.1_clusters_from_5x5-som_to_kmeans/ --input_file data/2013-12-10.4_vectors_cleaned_laccases_vec/lac_pt.faa.pruneBJXZUO.fasta.SPLITAA.EL_DESC_VAL.vec
python bin/python_scripts/som_to_kmeans_clustering.py --som_x 5 --som_y 5 --out_dir data/2014-1-5.1_clusters_from_5x5-som_to_kmeans/ --input_file data/2013-12-10.4_vectors_cleaned_laccases_vec/lac_pt.faa.pruneBJXZUO.fasta.T.EL_DESC_VAL.vec

done
