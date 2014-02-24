'''
This program will calculate averages of the SOM which occur in this format.
lac_pt.faa.pruneBJXZUO.fasta.AACOMP.EL_DESC_VAL.vec     number_of_clusters      11       DBINDEX        0.910282890954   inter_cluster  0.000107045941178        intra_cluster  0.000285169629542
 
'''

import argparse
from glob import glob
import numpy as np

#Get arguments from command line, output will be printed.  
parser = argparse.ArgumentParser(description=__doc__)
parser.add_argument('--file_set',
                    required=True,
                    help='''A file or regex must be in format show in doc string.''' ) 
args      = parser.parse_args()
file_glob = glob(args.file_set)

#Gather all input data 
stat_dict = {}
for file_name in file_glob: 
    for line in open(file_name,'r'):
        sp_line = line.strip().split("\t")    

        in_file_name  = sp_line[0]
        numb_of_cltrs = float(sp_line[2] )
        dbindex       = float(sp_line[4])
        inter_c       = float(sp_line[6])
        intra_c       = float(sp_line[8])

        out_arr = [numb_of_cltrs,dbindex,inter_c,intra_c,inter_c-intra_c]
        
        if in_file_name in stat_dict:
             stat_dict[in_file_name].append( out_arr  )
        else: 
            stat_dict[in_file_name]=[ out_arr ]


for vec_type in sorted(stat_dict):
   data_set = np.array( stat_dict[vec_type] )
   
   avg_inter_intra_cluster_diff = np.average(data_set[:,-1],axis=0)
   avg_cluster_numb             = np.average(data_set[:,0],axis=0) 
   print(vec_type+"\t"+str(avg_inter_intra_cluster_diff)+"\t"+str(avg_cluster_numb))

