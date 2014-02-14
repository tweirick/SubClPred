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
np.random.seed(seed+10)
random.seed(seed)

#http://www.biomedcentral.com/content/pdf/1471-2105-3-36.pdf
#Clustering of the SOM easily reveals distinct gene expression patterns: results of a reanalysis of lymphoma study
#http://www.biomedcentral.com/1471-2105/3/36/
#http://www.plosone.org/article/info%3Adoi%2F10.1371%2Fjournal.pone.0038817
#http://www.umich.edu/~mycology/resources/Publications/hoegger-et-al.pdf
#Phylogenetic comparison and classification of laccase and related multicopper oxidase protein sequences.
#http://www.ncbi.nlm.nih.gov/pubmed/16650005
#The Laccase Engineering Database: a classification and analysis system for laccases and related multicopper oxidases
#http://www.ncbi.nlm.nih.gov/pmc/articles/PMC3077825/
#Laccase versus Laccase-Like Multi-Copper Oxidase: A Comparative Study of Similar Enzymes with Diverse Substrate Spectra 
#http://www.plosone.org/article/info%3Adoi%2F10.1371%2Fjournal.pone.0065633;jsessionid=277E9CDF36E8C0C72C3D40C505CC86D3
#The interplay of descriptor-based computational analysis with pharmacophore modeling builds the basis for a novel classification scheme for feruloyl esterases
#http://www.sciencedirect.com/science/article/pii/S0734975010001175


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


in_file,show_clustering_performance = getargs()


rl         = []
line_ids   = []
for line in open(in_file,'r'):#fileinput.input():
    spl = line.split()
    vl = []
    dl = []
    for domain_el in spl[1:]:
        spel = domain_el.split(":")
        vl.append(float(spel[1]))
        domain_set.add(spel[0])
    rl.append( vl )
    line_ids.append( spl[0] )

np_vecs     = scale( np.array( rl  ) )

km              = KMeans(n_clusters=k,n_init=k)
kmf             = km.fit(np_vecs)
dict_cnt        = {}
membership_dict = {}

for i in range( len(km.labels_) ):
    ln = str(km.labels_[i])
    if ln in dict_cnt:   
        dict_cnt[ln].append(         np_vecs[i]  )
        membership_dict[ln].append( line_ids[i] )
    else:
        dict_cnt[ln]       =[np_vecs[i] ]
        membership_dict[ln]=[line_ids[i]]


in_file = in_file.split("/")[-1]


if show_clustering_performance:
    #Make arrays of the centers of the clusters, 
    #Make arrays fo of the clusters   
    intra_cluster     = []
    cluster_centroids = []
    for cluster_el in sorted(dict_cnt):
        cluster_centroids.append(  np.mean( np.array( dict_cnt[cluster_el]),axis=0))
        #print( np.mean( np.array( dict_cnt[cluster_el]),axis=0))
        intra_cluster.append(      np.var(  np.array(dict_cnt[cluster_el])))
        #print(np.var(  np.array(dict_cnt[cluster_el])))
    inter_cluster_var = np.var( np.array(cluster_centroids ) )
    intra_cluster_var = np.mean( np.array(intra_cluster ) )
    print(in_file+"\tnumber_of_clusters\t"+str(k)+
        "\tinter_cluster\t"+str(inter_cluster_var)+"\tintra_cluster\t"+
        str(intra_cluster_var) 
    )
else:
    for w in sorted(membership_dict):
        #print w,membership_dict[w]
        seq_list = []
        for seq_id in membership_dict[w]:
            seq_list.append(seq_id)
        print("cluster-"+str(w)+"\t"+"len:"+str(len(seq_list))+"\t"+" ".join(seq_list))