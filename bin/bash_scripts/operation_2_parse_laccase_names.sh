#!/bin/bash
#Want to analize the names of sequences obtained just to be sure.
mkdir data/2013-12-10.2_laccase_names_tsv 
python bin/python_scripts/tabulate_uniprot_fasta_prot_descs.py --file_set data/2013-12-10.1_laccase_faa/lac_pt.faa > data/2013-12-10.2_laccase_names_tsv/lac_pt.faa.names.tsv
