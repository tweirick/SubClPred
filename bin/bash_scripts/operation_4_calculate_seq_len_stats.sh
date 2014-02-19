#!/bin/bash
#This will calculate the statistics about the lengths of sequences present in the data set. 
python bin/python_scripts/count_lengths_of_fasta_lines.py --file_set data/2013-12-10.1_laccase_faa/lac_pt.faa --out_dir data/2013-12-10.3_cleaned_laccases_faa/ --title_text Lengths_of_Laccase_Fasta_Entries  --graph_out_dir docs/
