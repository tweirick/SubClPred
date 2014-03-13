"""
@author: Tyler Weirick
@date: 2013-08-18
==============================================================================
This program takes sets of vectors pertaining to classes and finds a subset 
of features better correlated to classification. It will output the findings 
as a file with the element descriptions sepatated by spaces. The input 
vectors need to be in description:value format.
==============================================================================

python scale_features.py --comma_delimited_vec_file_names "../data/LO1_proteinlevel.faa.on.trembl.psiblastout.acflat.faa.pruneBJXZ.fasta.AACOMP-DIPEP-SPLITAA-PLPredPhysChem22_nolog.EL_DESC_VAL.vec,../data/4CL_proteinlevel.faa.on.trembl.psiblastout.acflat.faa.pruneBJXZ.fasta.AACOMP-DIPEP-SPLITAA-PLPredPhysChem22_nolog.EL_DESC_VAL.vec" --out_file_name test_red.txt
"""

from sklearn.svm import SVC
from sklearn.datasets import load_digits
from sklearn.feature_selection import RFE
from sklearn.feature_selection import SelectKBest, f_regression, f_classif
from sklearn.cluster import KMeans
from sklearn import datasets
import numpy as np
import pylab as pl
import argparse 
from glob import glob
def unique(a):
    order = np.lexsort(a.T)
    a = a[order]
    diff = np.diff(a, axis=0)
    ui = np.ones(len(a), 'bool')
    ui[1:] = (diff != 0).any(axis=1) 
    return a[ui]

#=============================================================================
#                               Classes
#=============================================================================

class clusterset():
    '''Used for storing and interpreting feature scaling. 
    '''
    def __init__(self):
        self.data       = []
        self.true_preds = []
        self.el_titles  = []

    def makefeaturereductionarray(self,class_file_list):
        '''
        Vectorize a set of classes. This creates two lists, the first 
        containing a 2D array of values obtained from the input vectors
        the second vector contains one entry for each row in the 2D 
        vector which corresponds to the class it belongs to.
        output: 
        self.data = np.array([ [0.93423,0.25662,...]...])
        self.true_preds = np.array([ 0,0,...1,1,...])
        '''
        out_vecs   = []
        tpred_list = []
        first_round = True 
        class_number = 0
        line_dict = {}
        vecs_seen = set()
        last_set = set()
        disjoint_set = set()
        for file_name in class_file_list:
            tmp_set = set()
            for file_name in class_file_list:
                 nfile = open(file_name,'r')
                 line = nfile.readline()
                 for vec_el in line.strip().split()[1:]:
                     #print(vec_el)
                     tmp_set.add( "".join( vec_el.split(":")[:-1] ) )
                     #el_val   = vec_el.split(":")[-1].strip()
                 if last_set == set(): 
                     last_set = tmp_set 
                 for e in last_set ^ tmp_set:
                      #print(line)
                      print(e)
                      disjoint_set.add(e)
        
        for file_name in class_file_list:
            print(file_name)
            domain_els_normal = set()
            domain_els_extra  = set()

            for line in open(file_name,'r'):
                tmp_vals = []
                #For a vec line get all values and store as an list
                if not line in vecs_seen:
                    vecs_seen.add(line)
                    
                    line_len = len( line.strip().split() )
                    for vec_el in line.strip().split()[1:]:
                        el_title = "".join( vec_el.split(":")[:-1] ) 
                        el_val   = vec_el.split(":")[-1].strip()
                        #if len(vec_el.split(":") ) != 2: print(line.strip().split())
                        if first_round:  self.el_titles.append(el_title)      
                        if not el_title in disjoint_set:
                            tmp_vals.append(float(el_val))
                    first_round = False
                    #Add the row to another list. 
                    #print( len(tmp_vals))
                    out_vecs.append( tmp_vals )
                    #np.array(tmp_vals).astype(np.float32)
                    #For each row there must be a corresponding class number. 
                    tpred_list.append(class_number)
            first_round = False
            #Increment to indicate a new class. 
            class_number+=1
 
        self.data       = np.array(out_vecs).astype(np.float32)
        self.true_preds = np.array(tpred_list).astype(np.float32)    
#=============================================================================
#                               Functions
#=============================================================================

def getargs():
    parser = argparse.ArgumentParser()
    parser.add_argument('--vec_file_glob',
                       help='''Should be in EL_DESC_VAL format.''',
                       default=None)

    parser.add_argument('--comma_delimited_vec_file_names',
                       help='''Should be in EL_DESC_VAL format, 
                       takes presidence over vec_file_set if given.''',
                       default=None)

    parser.add_argument('--out_file_name',
                       help='''Should be in EL_DESC_VAL format, 
                       takes presidence over vec_file_set if given.''',
                       default=None)

    args = parser.parse_args()
    out_file_name = args.out_file_name
    vec_file_glob = args.vec_file_glob
    comma_delimited_vec_file_names = args.comma_delimited_vec_file_names

    if vec_file_glob != None: 
        file_name_set =  glob(vec_file_glob)

    if comma_delimited_vec_file_names != None:
        file_name_set = comma_delimited_vec_file_names.split(",")

    return file_name_set,out_file_name


#=============================================================================
#                               Main Program
#=============================================================================

file_name_glob,out_file_name = getargs()

#Make a cluster class
c = clusterset()
#Add generate vectors in the class. 
c.makefeaturereductionarray(file_name_glob)

#print(len( c.el_titles  ))
#Generate and array of boolean values, the T|F values depend on wether the 
#standard deviation of a given column is zero, indicating a column of all 
#the same values. This can cause errors durring the anova_filter fitting. 
bool_arr = (np.std(c.data, axis=0) == 0)
for i in range(len(bool_arr)-1,0,-1):
   if bool_arr[i]:
       #Remove the column from the data
       c.data = np.delete(c.data,i,1)
       #Also remove the column title from the element titles list. 
       #Important as equivalents is determined by order. 
       c.el_titles.pop(i)
#print(len( c.el_titles  ))

#print(c.el_titles)

#Can also use chi2, f_classif, or f_regression.
anova_filter = SelectKBest(f_regression, k='all')

z = anova_filter.fit(c.data,c.true_preds)

#Calculate stats for the run.  
score_array = z.scores_

>>>>>>>>>>>>>>>>>>>> File 1
#for i in range(len(z.scores_)):
>>>>>>>>>>>>>>>>>>>> File 2
#for i in range(len(z.scores_)):
>>>>>>>>>>>>>>>>>>>> File 3
#for i in sorted(list(z.scores_)):
<<<<<<<<<<<<<<<<<<<<
#    print(z.scores_[i],c.el_titles[i])


f_score_average = np.mean(z.scores_)
f_score_10th    = np.percentile(z.scores_,10)
f_score_25th    = np.percentile(z.scores_,25)
f_score_75th    = np.percentile(z.scores_,75.0)
f_score_90th    = np.percentile(z.scores_,90.0)
mean            = np.mean(z.scores_)
#mode            = np.mode(score_array)
median          = np.median(z.scores_)
min_val         = np.amin(z.scores_)
max_val         = np.amax(z.scores_)
#std_dev         = np.std(z.scores_)


#Get data for desired conditions.
 
score_above_mean_dict = {}
score_above_10th_dict = {}
score_above_25th_dict = {}
score_above_75th_dict = {}
score_above_90th_dict = {}

scores_array = []
for i in z.get_support(indices=True):

    scores_array.append( [ c.el_titles[i],z.scores_[i] ] )

    if z.scores_[i] > f_score_10th:
        score_above_10th_dict.update( {c.el_titles[i]:z.scores_[i]}  )    
    if z.scores_[i] > f_score_25th:
        score_above_25th_dict.update( {c.el_titles[i]:z.scores_[i]} )
    if z.scores_[i] > f_score_average:
        score_above_mean_dict.update( {c.el_titles[i]:z.scores_[i]} )
    if z.scores_[i] > f_score_75th:
        score_above_75th_dict.update({c.el_titles[i]:z.scores_[i]})
    if z.scores_[i] > f_score_90th:
        score_above_90th_dict.update({c.el_titles[i]:z.scores_[i]})

sorted_scores_array = sorted(scores_array, key=lambda x: x[1],reverse=True)
#print( sorted_scores_array )

stat_out_list = []
for el in sorted_scores_array:
    stat_out_list.append( el[0]+"\t"+str(el[1])   )
out_file = open(out_file_name + ".feature-vals.txt",'w')
out_file.write("\n".join(stat_out_list))
out_file.close()

#Output the
out_list = []
for e in sorted(score_above_10th_dict, key=score_above_mean_dict.get):
    out_list.append(e)
out_file = open(out_file_name+".reducedto10th.txt",'w')
out_file.write(" ".join(out_list))
out_file.close()
out_list = []

for e in sorted(score_above_25th_dict, key=score_above_75th_dict.get):
    out_list.append(e)
out_file = open(out_file_name+".reducedto25th.txt",'w')
out_file.write(" ".join(out_list))
out_file.close()
 
out_list = []
for e in sorted(score_above_mean_dict, key=score_above_mean_dict.get):
    out_list.append(e)
out_file = open(out_file_name+".reducedtomean.txt",'w')
out_file.write(" ".join(out_list))
out_file.close()
out_list = []
for e in sorted(score_above_75th_dict, key=score_above_75th_dict.get):
    out_list.append(e)
out_file = open(out_file_name+".reducedto75th.txt",'w')
out_file.write(" ".join(out_list))
out_file.close()
out_list = []
for e in sorted(score_above_90th_dict, key=score_above_90th_dict.get):
    out_list.append(e)
out_file = open(out_file_name+".reducedto90th.txt",'w')
out_file.write(" ".join(out_list))
out_file.close()


#sorted_scores_array
#out_file = open(out_file_name+".top500.txt",'w')
#out_file.write(" ".join( [str(i[0]) for i in sorted_scores_array[:500]]  ))
#out_file.close()

#sorted_scores_array
#out_file = open(out_file_name+".top400.txt",'w')
#out_file.write(" ".join( [str(i[0]) for i in sorted_scores_array[:400]]  ))
#out_file.close()

#sorted_scores_array
#out_file = open(out_file_name+".top300.txt",'w')
#out_file.write(" ".join( [str(i[0]) for i in sorted_scores_array[:300]]  ))
#out_file.close()


#sorted_scores_array
out_file = open(out_file_name+".top200.txt",'w')
out_file.write(" ".join( [str(i[0]) for i in sorted_scores_array[:200]]  )) 
out_file.close()

#sorted_scores_array
out_file = open(out_file_name+".top100.txt",'w')
out_file.write(" ".join( [ str(i[0]) for i in sorted_scores_array[:100]]  )) 
out_file.close()

#sorted_scores_array
out_file = open(out_file_name+".top50.txt",'w')
out_file.write(" ".join( [ str(i[0]) for i in sorted_scores_array[:50]]  )) 
out_file.close()


#sorted_scores_array
out_file = open(out_file_name+".4.txt",'w')
out_file.write(" ".join( [ str(i[0]) for i in sorted_scores_array[:4]]  ))
out_file.close()

#sorted_scores_array
out_file = open(out_file_name+".3.txt",'w')
out_file.write(" ".join( [ str(i[0]) for i in sorted_scores_array[:3]]  ))
out_file.close()

#sorted_scores_array
out_file = open(out_file_name+".2.txt",'w')
out_file.write(" ".join( [ str(i[0]) for i in sorted_scores_array[:2]]  ))
out_file.close()




out_list = []
out_list.append("Total number of elements: "        +str(len(z.scores_)))
out_list.append("Elements above average: "+str(len(score_above_mean_dict)))
out_list.append("Elements above 75th percentile: "+str(len(score_above_75th_dict)))
out_list.append("Elements above 90th percentile: "+str(len(score_above_90th_dict)))
out_list.append("MIN: "+str(min_val))
out_list.append("25th_PER: "+str(f_score_25th))
out_list.append("MEAN: "+str(mean))
out_list.append("MEDIAN: "+str(median))
out_list.append("75th_PER: "+str(f_score_75th))
out_list.append("MAX: "+str(max_val))
#out_list.append("STANDARD_DEV: "+str(std_dev))
#out_list.append(" ".join(out_list))
out_file = open(out_file_name+".reduced.vec.log",'w')
out_file.write("\n".join(out_list))
out_file.close()

#Output full list of weights 
#Output elements over threshold as .vecids
