"""
@author: Tyler Weirick
@date: 2013-5-20
This script used the pysvmlight library to preforme n-fold testing.

"""
from time import time
from math      import sqrt
from sys       import hexversion
from glob      import glob
#SVM light library downloaded from.
#https://bitbucket.org/wcauchois/pysvmlight
import svmlight
import argparse



class PerformanceCalculation():

    def __init__(self,FN=0,FP=0,TN=0,TP=0):
        self.FN = float(FN)
        self.FP = float(FP)
        self.TP = float(TP)
        self.TN = float(TN)
        self.run_params = ""        
    def getaccuracy(self):
        numerator   = float(100*(self.TP+self.TN))
        denominator = float((self.TP+self.TN+self.FP+self.FN))
        if denominator == 0:
            return 0
        else:
            return numerator/denominator
            
    def getprecision(self):
        numerator   = float(100*self.TP)
        denominator = float(self.TP+self.FP)
        if denominator == 0:
            return 0
        else:
            return numerator/denominator
        
    def getsensitivity(self):
        #{"Sensitivity":100*true_pos/(true_pos+false_neg)})
        numerator   = float(100*self.TP)
        denominator = float(self.TP+self.FN)
        if denominator == 0:
            return 0
        else:
            return numerator/denominator
        
    def getspecificity(self):
        #{"Specificity":100*true_neg/(true_neg+false_pos)})
        numerator   = float(100*self.TN)
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
        numerator   = float(self.TP*self.TN-self.FP*self.FN)

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
        #({"Error": 100*numerator/denominator })  numerator   = 100*TP
        #numerator   = (false_pos+false_neg)
        #denominator = (true_pos+true_neg+false_pos+false_neg)
        numerator   = float(100*self.FP+self.FN)
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
        "accuracy"   :self.getaccuracy(),
        "error"      :self.geterror(),
        "mcc"        :self.getMCC(),
        "precision"  :self.getprecision(),
        "sensitivity":self.getsensitivity(),
        "specificity":self.getspecificity(),
        "FN"        :self.FN,
        "FP"  :self.FP,
        "TP":self.TP,
        "TN":self.TN
         }
        return rtrn_str

    def gettitle(self):
        rtrn_str = (
        'accuracy   \t'+
        'error      \t'+
        'MCC        \t'+
        'precision  \t'+
        'sensitivity\t'+
        'specificity')
        
        return rtrn_str



#=============================================================================
#                           Functions 
#=============================================================================

from math import ceil

def splitlistequally(l, n):
    return [ l[i::n] for i in xrange(n) ]

def parsevectorfiletolist(file_name,label_val,vector_list=[]):
    """
    (<label>, [(<feature>, <value>), ...], <queryid>)
    """ 
    vector_list=[]
    vect_dict = {}
    cnt = 0 
    for line in open(file_name,'r'):
        sp_line = line.strip().split()
        vec_id = sp_line[0]        
        tmp_vec_list = []
        cnt+=1
        #Count from one to account for missing first el
        for i in range(1,len(sp_line[1:])):
            vec_el = sp_line[i].split(":")

            #There should be at most one colon in the elment.
            assert len(vec_el) <= 2
            #Getting the last value allows for the case of no colons to 
            #be handled as well
            vec_val = float(vec_el[-1])
            tmp_vec_list.append( (i,vec_val) )
        vector_list.append( (label_val,tmp_vec_list, hash(file_name.split("/")[-1]+"\t"+vec_id))  )
        vect_dict.update({ hash(file_name.split("/")[-1]+"\t"+vec_id) : file_name.split("/")[-1]+"\t"+vec_id } )
        #vector_list.append( (label_val,tmp_vec_list)  )
    #print(cnt)
    #print(len(vector_list))
    return vector_list,vect_dict


def drange(start, stop, step):
    r = start
    while r < stop:
        yield r
        r += step

def getpredtype(pred,known_val,fold_fn,fold_fp,fold_tn,fold_tp):
    #Use list index fold_number to interperet.
    if pred >= 0 and known_val > 0:
        #True Positive
        fold_tp+=1
    elif pred < 0 and known_val > 0:
        #False Negative
        fold_fn+=1
    elif pred >= 0 and known_val < 0:
        #False Positive
        fold_fp+=1
    elif pred < 0 and known_val < 0:
        # True Negative
        fold_tn+=1
    else:
        #Error
        print(pred,known_val)
        print("Error this position should not be reached.")
        exit()
    return fold_fn,fold_fp,fold_tn,fold_tp


parser = argparse.ArgumentParser(description='Do training.')
parser.add_argument(
            '--pos_vecs_file',
            required=True, 
            type=str, 
            help='')

parser.add_argument(
            '--neg_vecs_file_set',
            required=True,
            type=str,
            help='')

parser.add_argument(
           '--out_base',
            type=str,
            help='')

parser.add_argument(
            '--k_folds', 
            type=int, 
            default=str,
            help='')

parser.add_argument(
            '--kernal',type=str,default="RBF",
             help='')

parser.add_argument(
            '--verbose',type=str,default=False,
             help='')


args = parser.parse_args()
pos_vecs_file     = args.pos_vecs_file
neg_vecs_file_set = sorted(glob(args.neg_vecs_file_set))
verbose           = args.verbose

if pos_vecs_file in neg_vecs_file_set: 
    #print(len(neg_vecs_file_set))
    neg_vecs_file_set.remove(pos_vecs_file)
    #print(len(neg_vecs_file_set))
    if verbose: print('Removing the file"'+pos_vecs_file+'" from the negative set.')

out_base          = args.out_base
kernal_type       = args.kernal
k_folds           = args.k_folds

#Make positive folds.
pos_vecs_list,pos_vec_dict = parsevectorfiletolist(pos_vecs_file,1)
kfold_pos_vecs_list = splitlistequally(pos_vecs_list, k_folds)

#Make negative folds.
neg_vecs_list_kfold_list = [[] for i in range( k_folds)]
neg_vec_dict = dict()
for file_name in neg_vecs_file_set:
    #print(file_name)
    neg_vecs_list = []

    neg_vecs_list,neg_tmp_dict = parsevectorfiletolist(file_name,-1)
    neg_vec_dict.update(neg_tmp_dict)

    #print("len neg vec",len(neg_vecs_list),k_folds)
    #print(neg_vecs_list[0])
    tmp_folds = splitlistequally(neg_vecs_list, k_folds)   
 
    assert len(tmp_folds) == k_folds, tmp_folds[-1]
    for fold_i in range(len(tmp_folds)):
        neg_vecs_list_kfold_list[fold_i]+=tmp_folds[fold_i]
        
#Add combine postive and negative folds. 
main_fold_list = [[] for i in range( k_folds)]
for i in range(k_folds):
    main_fold_list[i] = kfold_pos_vecs_list[i] + neg_vecs_list_kfold_list[i] 
    

pos_vec_dict.update(neg_vec_dict)
main_hash_dict = pos_vec_dict

def maketrainandmodelfold(train_fold_i,k_folds,main_fold_list):
    model_list,test_list = [] , []
    for fold_i in range(k_folds):
        if fold_i == train_fold_i:
            test_list+=main_fold_list[fold_i]
        else:
            model_list+=main_fold_list[fold_i]
    return model_list,test_list


best_stats = None
if "LIN" == kernal_type or  "linear" == kernal_type:
    print("Comming soon.")
elif "POLY" == kernal_type or  "polynomial" == kernal_type:
    print("Comming soon.")
elif "RBF" == kernal_type or  "radial" == kernal_type:
    if True:  
        gamma_iter = drange(0.01,450.0,30.0)#  # 1.0,350.0,5.0   #1.0,600.0,25.0
        j_iter     = drange(0.01,15.0,1.5)  # # 1.0,15.0,1.5
        c_iter     = drange(0.01,500.0,35.0)# # 1.0,510.0,30.0  #1.0,800.0,30.0
    else: 
        print("to do, add file input.")
        gamma_iter = []
        j_iter     = []
        c_iter     = []
    pred_list = []
    training_history = [] 
    gamma_iter = list(gamma_iter)
    j_iter     = list(j_iter)
    c_iter     = list(c_iter)

    for g in gamma_iter:
        for j in j_iter:
            for c in c_iter:
                #For overall stats
                total_fn,total_fp,total_tn,total_tp = 0,0,0,0
                #A list of all sequences considered and the values for their pred. 
                pred_results_list = []
                fold_stats = []
                for fold_i in range(k_folds):
                    #For fold specific stats
                    f_fn,f_fp,f_tn,f_tp = 0,0,0,0
                     
                    #Get the training and testing sets. 
                    model_list,test_list = maketrainandmodelfold(
                    fold_i,k_folds,main_fold_list)
                      
                    #Make a model.
                    fold_model = svmlight.learn(
                        model_list,
                        type='classification',kernel='rbf',
                        C=c,rbf_gamma=g,costratio=j)
                     
                    #Get predicitons from test data.
                    #It is a bit faster to do all calcs at once but only about
                    #~0.005s faster  
                    for vector_to_pred in test_list:
                        true_value = vector_to_pred[0] 
                        raw_vector = vector_to_pred[1]
                        seq_id     = main_hash_dict[vector_to_pred[2]]
                        pred_list = svmlight.classify(
                            fold_model,
                            [(0,raw_vector)],
                        )
                        pred = pred_list[0]
                        #Interpret the value from the prediction.
                        f_fn,f_fp,f_tn,f_tp = getpredtype(pred,true_value,f_fn,f_fp,f_tn,f_tp)  
                        pred_results_list.append(seq_id+"\t"+str(true_value)+"\t"+str(pred))                   
                     
                    #Calculate fold stats. 
                    pref_calc = PerformanceCalculation(f_fn,f_fp,f_tn,f_tp)
                    pref_calc.run_params = "fold:"+str(fold_i)+"\t--g:"+str(g)+"\t--j:"+str(j)+"\t--c:"+str(c)
                    training_history.append( pref_calc.getperformance() )         
                     
                    total_fn+=f_fn 
                    total_fp+=f_fp
                    total_tn+=f_tn
                    total_tp+=f_tp
                #Calcualte overall stats. 
                overall_stats = PerformanceCalculation(total_fn,total_fp,total_tn,total_tp)
                overall_stats.run_params = "overall\tnumber_of_folds: "+str(k_folds)+"\t--g:"+str(g)+"\t--j:"+str(j)+"\t--c:"+str(c)

                #Check to see if best. 
                if best_stats == None or overall_stats.getMCC() > best_stats.getMCC(): 
                    best_stats =  overall_stats
                    best_pred_results_list = pred_results_list   
                    
                    of = open(out_base+".best-stats.txt",'w')
                    of.write( best_stats.getperformance()+"\n" )
                    of.close()
                    
                    of = open(out_base+".preds.txt",'w')
                    of.write("\n".join(best_pred_results_list))
                    of.close()
                if verbose: print( overall_stats.getperformance() ) 
                training_history.append( overall_stats.getperformance()  )            

    of = open(out_base+".all-stats.txt",'w')
    of.write( "\n".join(training_history) )
    of.close()            

elif "SIG" == kernal_type or  "sigmoid" == kernal_type:
    print("Comming soon.")
else: 
    print("ERROR kernal type not recognized.")


