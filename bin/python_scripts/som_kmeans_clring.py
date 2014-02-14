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
    #print("Ti",Ti)
    #print(Ti)
    Ai = np.mean(Ci,axis=0)
    #print(Ai)
    CiAi = Ci - Ai 
    #print(CiAi)
    #Zi = np.sqrt( np.linalg.norm(CiAi-CiAi)  )**float(q)
    Zi = np.sqrt( np.einsum('ij,ij->i', CiAi,CiAi )  )**float(q) 
    #Zi = np.sqrt( np.einsum('i,i->i',CiAi,CiAi)  )**float(q)
    #print(Zi)
    Si = ( (1.0/Ti) * np.sum(Zi)   )**(1.0/q)
    #print(Si)
    return Si

def daviesbouldin(C,q=2):
     
     len_c = len(C)
     intracluster_dists =  []
     for Ci in C: 
         intracluster_dists.append( intraclusterdistance(Ci,q)  )
     S = np.array(intracluster_dists)

     M = np.empty(shape=(len(C),len(C)))     

     for i in range(0,len_c): 
         Ai = np.mean(C[i],axis=0)
         for j in range(0,len_c): 
             Aj = np.mean(C[j],axis=0)
             #M[i][j] = np.sqrt( np.linalg.norm(C[i]-C[j]) )**float(q)
             M[i][j] = np.sqrt( np.linalg.norm(Ai-Aj)  )**float(q) 
             #M[i][j] = np.sqrt( np.einsum('ij,ij->i',Ai,Aj)  )**float(q)
      
     #Find differences between clusters 
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
     #print(Di,len_c)
     DB = Di / float( 100.0  ) 
     return DB
   
def getclusternumber(orange_map,vec_len):
    cluster_cnt = 0
    
    max_neurons = 0
    for n in orange_map:
        neuron_cnt = 0
        for i in n.instances:
            neuron_cnt+=1
        if max_neurons < neuron_cnt: 
           max_neurons = neuron_cnt 

    for neuron_cutoff in range(2,max_neurons):
        all_clusters = []
        for n in orange_map:
             no_lable = []
             for i in n.instances:
                 #Get rid of tags at the end of data.  
                 no_lable.append( i.native(0)[:-1] )

             if len(no_lable) > neuron_cutoff:
                 all_clusters.append(no_lable)
             
        if daviesbouldin( np.array(all_clusters)  ) <= 0.03: 
            break
    cluster_cnt = len(all_clusters)
    return cluster_cnt     

"""
#cluster_members = np.array(no_lable)
all_clusters[i]=np.array(no_lable)
i+=1
#print("daviesbouldin")
#print( daviesbouldin(cluster_members) )
#Ti = len( cluster_members  )
#cluster_centroid = np.mean(cluster_members,axis=0)
"""


def kmeansperformance():
    intarcluster_var,intercluster_var = -1,-1
    return intarcluster_var,intercluster_var

in_file,show_clustering_performance = getargs()

#Get vectors from standard input.
out_list   = []
domain_set = set()
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
d           = Orange.data.Domain( [Orange.feature.Continuous(x) for x in sorted(domain_set) ] )
orange_vecs = Orange.data.Table(d,np_vecs)   
#print(str(len(line_ids))+" vectors of length "+str(len(vl)) +" imported.")
som     = Orange.projection.som.SOMLearner(map_shape=(10, 10),initialize=Orange.projection.som.InitializeRandom,epochs=5000)
som_map = som(orange_vecs)
k       = getclusternumber(som_map,len(np_vecs))
#print(str(k)+" quality clusters found")
#k=10
km              = KMeans(n_clusters=k,n_init=k)
kmf             = km.fit(np_vecs)
#kmf             = km.fit(rl)
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

reduced_data  = PCA(n_components=2).fit_transform(np_vecs)
reduced_data  = scale( reduced_data )
kmeans = KMeans(init='k-means++', n_clusters=k, n_init=k)
kmeans.fit(reduced_data)
# Step size of the mesh. Decrease to increase the quality of the VQ.
h = .01     # point in the mesh [x_min, m_max]x[y_min, y_max].
x_min, x_max = reduced_data[:, 0].min() + 1, reduced_data[:, 0].max() - 1
y_min, y_max = reduced_data[:, 1].min() + 1, reduced_data[:, 1].max() - 1
#print(reduced_data)
#x_min, x_max = -0.15626329669343397, 0.057864973853256174 
#y_min, y_max = -40.0   , 20.0 
xx, yy = np.meshgrid(np.arange(x_min, x_max, h), np.arange(y_min, y_max, h))
# Obtain labels for each point in mesh. Use last trained model.
Z = kmeans.predict(np.c_[xx.ravel(), yy.ravel()])
# Put the result into a color plot
Z = Z.reshape(xx.shape)
pl.figure(1)
pl.clf()
#print(xx,yy)
pl.imshow(Z, interpolation='nearest',
          extent=(
          xx.min(), xx.max(), yy.min(), yy.max()),
          cmap=pl.cm.Paired,
          aspect='auto', origin='lower')
pl.plot(reduced_data[:, 0], reduced_data[:, 1], 'k.', markersize=2)
# Plot the centroids as a white X
centroids = kmeans.cluster_centers_
pl.scatter(centroids[:, 0], centroids[:, 1],
           marker='x', s=169, linewidths=3,
           color='w', zorder=10)

#title_str = 'K-means clustering of (PCA-reduced data) from\n'+in_file+'\nCentroids are marked with white cross'
title_str = 'K-means clustering of (PCA-reduced data) from\n'+in_file+'\nk='+str(k)+' Centroids are marked with white cross'

pl.title( title_str,fontsize=10  )
pl.xlim(x_min, x_max)
pl.ylim(y_min, y_max)
pl.xticks(())
pl.yticks(())
#pl.show()
pl.savefig(in_file+"k-means-clusters.png")
