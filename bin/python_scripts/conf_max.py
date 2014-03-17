'''
This program will make a graphical and text based confusion matrix. 
'''
import argparse
from   glob import glob
import numpy as np
import re
import pylab as pl
from math import sqrt



class PerformanceCalculation():

    def __init__(self,FN=0,FP=0,TN=0,TP=0):
        self.FN = float(FN)
        self.FP = float(FP)
        self.TP = float(TP)
        self.TN = float(TN)
        self.run_params = ""

    def getaccuracy(self):
        numerator = float(100*(self.TP+self.TN))
        denominator = float((self.TP+self.TN+self.FP+self.FN))
        if denominator == 0:
            return 0
        else:
            return numerator/denominator
            
    def getprecision(self):
        numerator = float(100*self.TP)
        denominator = float(self.TP+self.FP)
        if denominator == 0:
            return 0
        else:
            return numerator/denominator
        
    def getsensitivity(self):
        #{"Sensitivity":100*true_pos/(true_pos+false_neg)})
        numerator = float(100*self.TP)
        denominator = float(self.TP+self.FN)
        if denominator == 0:
            return 0
        else:
            return numerator/denominator
        
    def getspecificity(self):
        #{"Specificity":100*true_neg/(true_neg+false_pos)})
        numerator = float(100*self.TN)
        denominator = float(self.TN+self.FP)
        if denominator == 0:
            return 0
        else:
            return numerator/denominator

    def ifzeronone(self,x):
        if x == 0.0:
            return 1.0
        else:
            return x

    def getMCC(self):
        #numerator = (true_pos*true_neg)-(false_pos*false_neg)
        #denominator = (true_pos+false_pos)*(true_pos+false_neg)*(true_neg+false_pos)*(true_neg+false_neg)
        numerator = float(self.TP*self.TN-self.FP*self.FN)

        a=self.ifzeronone(self.TP+self.FP)
        b=self.ifzeronone(self.TP+self.FN)
        c=self.ifzeronone(self.TN+self.FP)
        d=self.ifzeronone(self.TN+self.FN)

        denominator = float(sqrt(a*b*c*d))

        if denominator == 0:
            return 0
        else:
            return numerator/denominator
    
    def geterror(self):
        #({"Error": 100*numerator/denominator }) numerator = 100*TP
        #numerator = (false_pos+false_neg)
        #denominator = (true_pos+true_neg+false_pos+false_neg)
        numerator = float(100*self.FP+self.FN)
        denominator = float(self.TP+self.TN+self.FP+self.FN)
        if denominator == 0:
            denominator=-1
        return numerator/denominator

    def getperformance(self):
        rtrn_str = (
        self.run_params+"\t"+
        'accuracy: {number:.{digits}f}'.format(number=(self.getaccuracy()),digits=2)+"\t"+
        'error: {number:.{digits}f}'.format(number=(self.geterror()),digits=2)+"\t"+
        'getMCC: {number:.{digits}f}'.format(number=(self.getMCC()),digits=5)+"\t"+
        'precision: {number:.{digits}f}'.format(number=(self.getprecision()),digits=2)+"\t"+
        'sensitivity: {number:.{digits}f}'.format(number=(self.getsensitivity()),digits=2)+"\t"+
        'specificity: {number:.{digits}f}'.format(number=(self.getspecificity()),digits=2)+"\t"+
        'FN: '+str(self.FN)+"\t"+
        'FP: '+str(self.FP)+"\t"+
        'TP: '+str(self.TP)+"\t"+
        'TN: '+str(self.TN)+"\t"
         )
        return rtrn_str
    
    def getperformanceasdict(self):
        rtrn_str = {
        "accuracy" :self.getaccuracy(),
        "error" :self.geterror(),
        "mcc" :self.getMCC(),
        "precision" :self.getprecision(),
        "sensitivity":self.getsensitivity(),
        "specificity":self.getspecificity(),
        "FN" :self.FN,
        "FP" :self.FP,
        "TP":self.TP,
        "TN":self.TN
         }
        return rtrn_str

    def gettitle(self):
        rtrn_str = (
        'accuracy \t'+
        'error \t'+
        'MCC \t'+
        'precision \t'+
        'sensitivity\t'+
        'specificity')
        
        return rtrn_str



def natural_sort(l): 
    convert = lambda text: int(text) if text.isdigit() else text.lower() 
    alphanum_key = lambda key: [ convert(c) for c in re.split('([0-9]+)', key) ] 
    return sorted(l, key = alphanum_key)


parser = argparse.ArgumentParser(description=__doc__)
parser.add_argument('--pred_files',
    required=True,
    help='A set of fastas, fastas will be removed from this file.')

parser.add_argument('--title',
    required=False,
    default="",
    help='add extra to title of graph')

args       = parser.parse_args()
pred_files = args.pred_files
title      = args.title

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
        seq_class  = spl[0].split(".cluster-members.txt.")[-1].split( "faa.0."  )[0]
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
    pred_val  = float(seq_dict[seq][ sorted_col[0]]) 
    max_class = sorted_col[0]
    for class_name in natural_sort(seq_dict[seq]):
        if seq_true_class_dict[seq][class_name] == "1":
            true_class = class_name
            if first_is_true: #and pred_val >= 0:
                tps+=1
            else:
                fns+=1
        else:
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
        confusion_dict[true_class].update({max_class:0})
    print(true_class,max_class)
    confusion_dict[true_class][max_class]+=1
                      


print(confusion_dict)



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
pref_calc = PerformanceCalculation(FN=fns,FP=fps,TN=tns,TP=tps)
print(pref_calc.getperformance())

# Show confusion matrix in a separate window
pl.matshow(matrix)
pl.title('Confusion matrix '+title)
pl.colorbar()
pl.ylabel('True label')
pl.xlabel('Predicted label')
pl.savefig(title.replace(" ","_")+'.png', bbox_inches='tight')



