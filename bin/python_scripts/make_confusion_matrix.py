'''
This program will make a graphical and text based confusion matrix and ROC 
curves.  
'''
import argparse
from glob import glob
import numpy as np
import re
import pylab as pl
from   sklearn.metrics import roc_curve, auc
from   sklearn.cross_validation import train_test_split
from   sklearn.preprocessing import label_binarize
from   sklearn.multiclass import OneVsRestClassifier
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
    '''Sort strings containing numbers so that numbers are sorted by size. 
    '''
    convert = lambda text: int(text) if text.isdigit() else text.lower() 
    alphanum_key = lambda key:[convert(c) for c in re.split('([0-9]+)',key)] 
    return sorted(l, key = alphanum_key)

parser = argparse.ArgumentParser(description=__doc__)
parser.add_argument('--pred_files',required=True,
    help='A set of pred files of the format file_name\tAC\ttrue_pred\tscore')

parser.add_argument('--graph_title',required=False,
    help='Text here will be appended to the graph title.')

args           = parser.parse_args()
pred_file_glob =  natural_sort(glob(args.pred_files))
title_text     = args.graph_title

check_for_sim   = {}
class_order     = {}
class_list      = []
#Contains 
seq_dict             = {}
seq_true_class_dict  = {}
score_list           = []
true_class_list      = []
true_class_name_list = []

#Build data structures from pred files
#Order in pred files is not preseved, so we must re-order the preds

class_len = len(pred_file_glob)
f         = open(pred_file_glob[0],'r')
line_len  = len( f.readlines() ) 
f.close()


arr_dimens = (line_len,class_len)
score_array   = np.zeros( arr_dimens  )
true_pred_arr = np.zeros( arr_dimens  )

f_cnt = 0
for file_name in pred_file_glob:
    class_name = file_name.split(".faa.0.")[0]
    preded_dict       = {}
    true_class_dict   = {}
    seq_in_class_dict = {}    
    
    for line in open(file_name,'r'):
        spl        = line.split()
        seq_class  = spl[0].split(".cluster-members.txt.")[-1].split( "faa.0." )[0]
        seq_id     = spl[1]
        true_class = spl[2]
        pred_val   = spl[3]
        
        if seq_id in seq_dict:
            seq_dict[seq_id].update( { class_name: pred_val  }  )   
            seq_true_class_dict[seq_id].update( { class_name: true_class  }  )
        else:
            seq_dict.update({seq_id:{class_name:pred_val}})
            seq_true_class_dict.update({seq_id : { class_name: true_class }} )
    
    line_cnt = 0
    for seq_id in sorted(seq_dict):
        
        score_array[line_cnt][f_cnt] = float(seq_dict[seq_id][class_name] )
        tcn = seq_true_class_dict[seq_id][class_name]
        if tcn   == "1":
            true_pred_arr[line_cnt][f_cnt] = 1.0
        elif tcn == "-1":
            true_pred_arr[line_cnt][f_cnt] = 0.0
        else: 
            print("ERROR",seq_true_class_dict[seq_id] )     
        line_cnt+=1
    f_cnt+=1
    
    class_list.append( [class_name,preded_dict,true_class_dict,seq_in_class_dict] )


# Compute ROC curve and ROC area for each class
fpr     = dict()
tpr     = dict()
roc_auc = dict()
for i in range(0,len(pred_file_glob)):
    fpr[i], tpr[i], _ = roc_curve(true_pred_arr[:, i], score_array[:, i])
    roc_auc[i] = auc(fpr[i], tpr[i])

# Compute micro-average ROC curve and ROC area
fpr["micro"], tpr["micro"], _ = roc_curve(true_pred_arr.ravel(), score_array.ravel())
roc_auc["micro"] = auc(fpr["micro"], tpr["micro"])

# Plot ROC curve
pl.clf()
pl.plot(fpr["micro"], tpr["micro"],
        label='micro-average ROC curve (area = {0:0.2f})'
              ''.format(roc_auc["micro"]))

for i in range( 0,len(pred_file_glob)  ):
    #print( fpr[i],tpr[i],  roc_auc[i] )
    pl.plot(fpr[i], tpr[i], label='class {0} (area = {1:0.2f})'
                                  ''.format(i, roc_auc[i]))
pl.plot([0, 1], [0, 1], 'k--')
pl.xlim([0.0, 1.0])
pl.ylim([0.0, 1.05])
pl.xlabel('False Positive Rate')
pl.ylabel('True Positive Rate')
pl.title('Receiver Operating Characteristics '+title_text)
pl.legend(loc="lower right")
pl.savefig('foo_roc.png', bbox_inches='tight')

tps,fps,fns,tns=0,0,0,0
confusion_dict = {}

true_class = None
for seq in sorted(seq_dict): 
    sorted_col = sorted(seq_dict[seq], key=seq_dict[seq].get,reverse=True)
    first_is_true = "1" == seq_true_class_dict[seq][sorted_col[0]]
    max_class = sorted_col[0]
    
    for class_name in sorted(seq_dict[seq]):
        if seq_true_class_dict[seq][class_name] == "1":
            true_class = class_name
            if first_is_true:
                tps+=1 
            else:
                fns+=1
        else: 
            if not first_is_true and class_name == max_class:
                fps+=1 
            else:
                tns+=1
 
    if true_class == None:
         print("ERROR: true class not found for this sequence. Exiting.")
         exit()
    if not true_class in confusion_dict:
        confusion_dict.update({true_class:{}})
    
    if not max_class in confusion_dict[true_class]:
        confusion_dict[true_class].update({max_class:0})

    confusion_dict[true_class][max_class]+=1
                      

out_matrix = []
matrix     = []

for true_el in natural_sort(sorted(confusion_dict, key=confusion_dict.get,reverse=True)):
   
   matrix_col     = []
   matrix_txt_col = []
   for pred_el in natural_sort(sorted(confusion_dict, key=confusion_dict.get,reverse=True)):
       
       if pred_el in confusion_dict[true_el]:
           matrix_col.append( confusion_dict[true_el][pred_el] )
           matrix_txt_col.append( str(confusion_dict[true_el][pred_el] ))
       else: 
            matrix_col.append( 0 )
            matrix_txt_col.append( "0" )

   matrix.append(matrix_col)
   out_matrix.append("\t".join(matrix_txt_col))


print(matrix)
print("\n".join(out_matrix))
print(len(confusion_dict)) 
print(confusion_dict)    
#print(fns,fps,tns,tps)
#FN=0,FP=0,TN=0,TP=0
pref = PerformanceCalculation(FN=fns,FP=fps,TN=tns,TP=tps)
print( pref.getperformance() )

# Show confusion matrix in a separate window
pl.matshow(matrix)
pl.title('Confusion matrix')
pl.colorbar()
pl.ylabel('True label')
pl.xlabel('Predicted label')
pl.savefig('foo.png', bbox_inches='tight')
