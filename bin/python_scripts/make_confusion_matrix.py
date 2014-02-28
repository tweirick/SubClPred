'''
This program will make a graphical and text based confusion matrix. 
'''
import argparse
from glob import glob
import numpy as np
import re
import pylab as pl


def natural_sort(l): 
    convert = lambda text: int(text) if text.isdigit() else text.lower() 
    alphanum_key = lambda key: [ convert(c) for c in re.split('([0-9]+)', key) ] 
    return sorted(l, key = alphanum_key)


parser = argparse.ArgumentParser(description=__doc__)
parser.add_argument('--pred_files',
    required=True,
    help='A set of fastas, fastas will be removed from this file.')

args       = parser.parse_args()
pred_files = args.pred_files

check_for_sim   = {}
class_order     = {}
class_list = []


seq_dict            = {}
seq_true_class_dict = {}

for file_name in sorted(glob(pred_files)):
    class_name = file_name.split(".faa.0.9.fasta.")[0]
    
    preded_dict       = {}
    true_class_dict   = {}
    seq_in_class_dict = {} 
   
    for line in open(file_name,'r'):
                     
        spl        = line.split()
        seq_class  = spl[0].split(".cluster-members.txt.")[-1].split( "faa.0.1.fasta."  )[0]
        seq_id     = spl[1]
        true_class = spl[2]
        pred_val   = spl[3]
        
        if seq_id in seq_dict:
            seq_dict[seq_id].update( { class_name: pred_val  }  )   
            seq_true_class_dict[seq_id].update( { class_name: true_class  }  )
        else:
            seq_dict.update({seq_id:{class_name:pred_val}})
            seq_true_class_dict.update({seq_id : { class_name: true_class }})
            
    class_list.append( [class_name,preded_dict,true_class_dict,seq_in_class_dict] )
        

#Sort data
tps=0
fps=0
fns=0
tns=0

confusion_dict = {}
for seq in seq_dict: 

    sorted_col = sorted(seq_dict[seq], key=seq_dict[seq].get,reverse=True)
    first_is_true = "1" == seq_true_class_dict[seq][sorted_col[0]]
    max_class = sorted_col[0]
    for class_name in natural_sort(seq_dict[seq]):
        if seq_true_class_dict[seq][class_name] == "1":
            true_class = class_name
            if first_is_true:
                tps+=1 
            else:
                fns+=1
        else: #seq_true_class_dict[seq][class_name] == "-1":
            if first_is_true: 
                tns+=1
            else:
                if class_name == max_class:
                    fps+=1 
                else:
                    tns+=1 
    if not true_class in confusion_dict:
        confusion_dict.update({true_class:{}})
    if not max_class in confusion_dict[true_class]:
        confusion_dict.update({true_class:{max_class:0}})
    confusion_dict[true_class][max_class]+=1
                      

out_matrix = []
matrix = []
for true_el in sorted(confusion_dict, key=confusion_dict.get,reverse=True):
   matrix_col     = []
   matrix_txt_col = []
   for pred_el in sorted(confusion_dict, key=confusion_dict.get,reverse=True):
       if pred_el in confusion_dict[true_el]:
           matrix_col.append( confusion_dict[true_el][pred_el] )
           matrix_txt_col.append( str(confusion_dict[true_el][pred_el] ))
       else: 
            matrix_col.append( 0 )
            matrix_txt_col.append( "0" )

   matrix.append(matrix_col)
   out_matrix.append("\t".join(matrix_txt_col))


print("\n".join(out_matrix))
print(len(confusion_dict)) 
print(confusion_dict)    
print(fns,fps,tns,tps)


# Show confusion matrix in a separate window
pl.matshow(matrix)
pl.title('Confusion matrix')
pl.colorbar()
pl.ylabel('True label')
pl.xlabel('Predicted label')
pl.savefig('foo.png', bbox_inches='tight')
