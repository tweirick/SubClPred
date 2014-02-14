

python3 bin/python_scripts/feature_generator/feature_generator.py --file_set data/2013-01-23.1_laccases_seqs_faa/LO1_tp_combined.faa.pruneBJXZU.fasta.greaterthan61chars.faa --vectors_to_generate "AACOMP,DIPEP,TRIPEP,SPLITAA,MoreauBroto,Moran,Geary,C,T,D,QuasiSequenceOrder,SequenceOrderCouplingNumberTotal,PLPredPhysChem22_nolog,Conjoint_Triad,shenC-shenD-shenT" --output_format EL_DESC_VAL  --line_name_format SEQ_ID

mv data/2013-01-23.1_laccases_seqs_faa/*.vec data/2013-01-24.1_laccases_vec/
