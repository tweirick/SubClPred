'''
This program will make a graphical and text based confusion matrix. 
'''
import argparse
from glob import glob
parser = argparse.ArgumentParser(description=__doc__)
parser.add_argument('--pred_files',
    required=True,
    help='A set of fastas, fastas will be removed from this file.')

args       = parser.parse_args()
pred_files = args.pred_files

for file_name in glob(pred_files):
    class_name = file_name.split(".faa.0.9.fasta.")[0]
    
    print(class_name)
    for line in open(file_name,'r'):
        print( line.split(  )[1] )
 
