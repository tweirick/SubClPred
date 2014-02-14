
python3 bin/python_scripts/feature_generator/feature_generator.py --file_set data/2013-01-20.1_test_feruloyl_esterases_seqs_faa/feruloyl_esterases.faa.pruneBJXZU.fasta.greaterthan61chars.faa  --vectors_to_generate "AACOMP,DIPEP,TRIPEP,SPLITAA,MoreauBroto,Moran,Geary,C,T,D,QuasiSequenceOrder,SequenceOrderCouplingNumberTotal,PLPredPhysChem22_nolog,Conjoint_Triad,shenC-shenD-shenT" --output_format EL_DESC_VAL  --line_name_format SEQ_ID

mv data/2013-01-20.1_test_feruloyl_esterases_seqs_faa/*.vec data/2013-01-20.2_test_feruloyl_esterase_vec/

