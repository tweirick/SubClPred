#!/bin/bash
#This will vectorize the sequences in the cleaned laccase file. 

mkdir data/2013-12-10.4_vectors_cleaned_laccases_vec/ 

python3 bin/python_scripts/feature_generator/feature_generator.py --vectors_to_generate "AACOMP,DIPEP,TRIPEP,SPLITAA,MoreauBroto,Moran,Geary,C,T,D,QuasiSequenceOrder,SequenceOrderCouplingNumberTotal,PLPredPhysChem22_nolog,Conjoint_Triad,shenC-shenT-shenD" --output_format EL_DESC_VAL --line_name_format SEQ_ID --file_set data/2013-12-10.3_cleaned_laccases_faa/lac_pt.faa.pruneBJXZUO.fasta --out_dir data/2013-12-10.4_vectors_cleaned_laccases_vec/
