
python3 bin/python_scripts/repair_fasta.py --file_set data/2013-01-23.1_laccases_seqs_faa/LO1_tp_combined.faa  --remove_over_x_percent 0.0

python3 bin/python_scripts/remove_fastas_less_than_x_chars.py --file_set data/2013-01-23.1_laccases_seqs_faa/LO1_tp_combined.faa.pruneBJXZU.fasta --min_chars 61
