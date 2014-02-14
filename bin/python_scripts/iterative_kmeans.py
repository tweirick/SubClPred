import fileinput, random, argparse
import numpy as np
from math import pow
from random  import random as rnd
import Orange
from math import sqrt
from sklearn.cluster import KMeans
from sklearn.decomposition import PCA
from sklearn.preprocessing import scale
import pylab as pl
seed = 42
seed = 0
seed = 10
seed = 5

def getargs():
    parser = argparse.ArgumentParser()
    parser.add_argument('--input_file',
                       help='''A vector file in EL_DESC_VAL format.''',required=True)    
    parser.add_argument('--print_variance',
                       help='''True will print cluster numbers and variances,
                            false will print the clusters and their IDS ''',
                       default=True)
    args = parser.parse_args()
    print_variance = str(args.print_variance)
    in_file = args.input_file
    if print_variance.lower() == "f" or print_variance.lower() == "false":
        print_variance = False        
    elif print_variance.lower() == "t" or print_variance.lower() == "true":
        print_variance = True
    else:        
        print("error must be bool")
        exit()
    return in_file,print_variance

def getvectors(file_name):
    dl = []
    out_list = []
    for line in fileinput.input():
        spl = line.split()
        out_list.append( spl )
    return out_list      

def intraclusterdistance(Ci,q=2):
    #Ci is a cluster of vectors. 
    #Xj is an n-dimensional feature.
    #Ai is the centroid of Ci  
    #Ti is the size of the cluster ?
    Ti = len(Ci)
    Ai = np.mean(Ci,axis=0)
    CiAi = Ci - Ai 
    #Zi = np.sqrt( np.linalg.norm(CiAi-CiAi)  )**float(q)
    Zi = np.sqrt( np.einsum('ij,ij->i', CiAi,CiAi )  )**float(q) 
    #Zi = np.sqrt( np.einsum('i,i->i',CiAi,CiAi)  )**float(q)
    Si = ( (1.0/Ti) * np.sum(Zi)   )**(1.0/q)
    return Si

def daviesbouldin(C,q=2):
     
     len_c = len(C)
     intracluster_dists =  []
     for Ci in C: 
         intracluster_dists.append( intraclusterdistance(Ci,q)  )
     S = np.array(intracluster_dists)
     M = np.empty(shape=(len(C),len(C)))     
     for i in range(0,len_c): 
         max_Rj = None
         Ai = np.mean(C[i],axis=0)
         for j in range(0,len_c): 
             Aj = np.mean(C[j],axis=0)
             M[i][j] = np.sqrt( np.linalg.norm(Ai-Aj)  )**float(q)

     D = np.empty(shape=(len_c ))
     Di = 0
     for i in range(0,len_c):
         max_Rj = None
         for j in range(0, len_c ):
             if i != j:
                 Rij = (S[i] + S[j])/M[i][j] 
                 if max_Rj > Rij or max_Rj == None:
                     max_Rj = Rij
         assert max_Rj != None
         D[i] = max_Rj
         Di+=max_Rj

     DB = Di / float( len_c  ) 
     return DB
  

def dokmeans(np_vecs,in_file,k_clusters):

    k  = k_clusters
    km              = KMeans(n_clusters=k,max_iter=10000,init="k-means++")
    kmf             = km.fit(np_vecs)

    dict_cnt        = {}
    membership_dict = {}
    #Convert into other data configs 
    for i in range( len(km.labels_) ):
        ln = str(km.labels_[i])
        if ln in dict_cnt:   
            dict_cnt[ln].append(         np_vecs[i]  )
            membership_dict[ln].append( line_ids[i] )
        else:
            dict_cnt[ln]       =[np_vecs[i] ]
            membership_dict[ln]=[line_ids[i]]
    
    fu = []
    for e in dict_cnt:
        fu.append(dict_cnt[e])
    fu = np.array(fu)
     
    kmeansDB = daviesbouldin(  fu  ) 

    in_file = in_file.split("/")[-1]

    intercluster_lengths = [] 
    if show_clustering_performance:
        #Make arrays of the centers of the clusters, 
        #Make arrays fo of the clusters   
        intra_cluster_dist  = []
        inter_cluster_dist  = []
        cluster_centroids   = []
        
        for cluster_el in sorted(dict_cnt):    
            cluster_centroid = np.mean( np.array( dict_cnt[cluster_el]),axis=0 )
            el_diffs             = cluster_centroid - dict_cnt[cluster_el] 
            sqrd_el_diff = (el_diffs)**2
            euc_dist = np.sqrt( np.sum(sqrd_el_diff,axis=0)  )
            for ed in euc_dist:
                intra_cluster_dist.append( ed  )
            intercluster_lengths.append( cluster_centroid )        

        total_centroid = np.mean( np.array( intercluster_lengths  ),axis=0 )
        el_diffs      = total_centroid - np.array( intercluster_lengths )
        sqrd_el_diff  = (el_diffs)**2
        inter_cluster_dist = np.sqrt( np.sum( sqrd_el_diff,axis=0  ) )
        intra_cluster_var = np.var( np.array(intra_cluster_dist) )
        inter_cluster_var = np.var( np.array(inter_cluster_dist) )

        #print(in_file+"\t number_of_clusters \t"+str(len(membership_dict))+"\t DBindex \t"+str(kmeansDB)+
        #    "\t inter_cluster \t"+str(inter_cluster_var)+"\t intra_cluster \t"+str(intra_cluster_var) )

    else:
        for w in sorted(membership_dict):
            #print w,membership_dict[w]
            seq_list = []
            for seq_id in membership_dict[w]:
                seq_list.append(seq_id)
            print("cluster-"+str(w)+"\t"+"len:"+str(len(seq_list))+"\t"+" ".join(seq_list))
    return k,kmeansDB

def textvecstolist(in_file): 
    domain_set = set()
    rl         = []
    line_ids   = []
    for line in open(in_file,'r'):
        spl = line.split()
        vl = []
        dl = []
        for domain_el in spl[1:]:
            spel = domain_el.split(":")
            vl.append(float(spel[1]))
            domain_set.add(spel[0])
        rl.append( vl )
        line_ids.append( spl[0] )
    return rl,line_ids,domain_set



#=============================================================================
#                              Main Program                 
#=============================================================================
#Get command line arguments 
in_file,show_clustering_performance = getargs()
#Get vectors from standard input.
rl,line_ids,domain_set = textvecstolist(in_file)
#Scale vecs 
np_vecs     = scale( np.array( rl  ) ,with_std=False,with_mean=False)


dbindex_cnt = {}

for j in range(50):
    min_k        = None
    min_db_index = None
    for i in range(2,30):
        k, kmeansDB = dokmeans(np_vecs,in_file,i)
        if min_db_index == None or kmeansDB < min_db_index:
             min_k        = str(k)
             min_db_index = kmeansDB

    if min_k in dbindex_cnt:
        dbindex_cnt[min_k].append(min_db_index)
    else:  
        dbindex_cnt[min_k]=[min_db_index]        
    


for dict_el in dbindex_cnt:
    print(dict_el, len(dbindex_cnt[dict_el]), np.average(np.array(dbindex_cnt[dict_el])))






