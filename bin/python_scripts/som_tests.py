import numpy as np
#from mvpa2.suite import *
from random  import random as rnd
import Orange
import random
from sklearn.cluster import KMeans

np.random.seed(42)

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


#def somcnt(vecs,sx=30,sy=20,iters=400,lr=0.05):
#    x = []
#    for i in range(vecs):
#        z = [rnd(),rnd(),rnd()]
#        x.append( z  )
#        #print(z)
#    colors = np.array(x)
#    som = SimpleSOMMapper((sx, sy), iters, learning_rate=lr)
#    som.train(colors)
#    #for i in colors: 
#    #    print(i)
#    print( len(som.K) )
#    pl.imshow(som.K, origin='lower')
#    pl.show()


#for i in range(300):
#    print( str(rnd())+"\t"+str(rnd())+"\t"+str(rnd()) )

def makesomevecs(V):
    x = []
    for i in range(V):
        z = [round(rnd(),3),round(rnd(),3),round(rnd(),3)]
        x.append(z)
    return x

def getvectors(file_name):
    out_list = []
    for line in open(file_name,'r'):
        spl = line.split()
        out_list.append( spl )
    return out_list      


from math import sqrt 
def getclusternumber(orange_map,vec_len):
    cluster_cnt = 0
    for n in orange_map:

         no_lable = []
         for i in n.instances:
             #Get rid of tags at the end of data.  
             no_lable.append( i.native(0)[:-1] )
             
         cluster_members = np.array(no_lable)
         #print(cluster_members)
          
         Ti = len( cluster_members  )
         cluster_centroid = np.mean(cluster_members,axis=0)
         Si = 0
         if Ti > 1 and float(Ti)/float(vec_len) > 0.02:
             #for vec_el in cluster_members:
             #    Si+=np.absolute(vec_el-cluster_members)**2
             #Si = srqt(Si/Ti)
             #print(Si)
             cluster_cnt+=1
    return cluster_cnt


file_name = "feruloyl_esterases.faa.AACOMP.SPACED_VALS.vec"
file_name = "LO1_314.tsv"
v = np.array( getvectors(file_name) )
vec_numb = len(v)
print( len(v) )


random.seed(0)
som = Orange.projection.som.SOMLearner(map_shape=(10, 10),initialize=Orange.projection.som.InitializeRandom,epochs=4000)
som_map = som(Orange.data.Table("fe.tab"))
k = getclusternumber(som_map,vec_numb)


km = KMeans(n_clusters=k)
kmf = km.fit(v)

print(kmf)
print(km.__dict__)

dict_cnt = {}
for i in km.labels_:
    ln = str(i)
    if ln in dict_cnt:   
        dict_cnt[ln]+=1
    else:
        dict_cnt[ln]=1



for w in sorted(dict_cnt, key=dict_cnt.get, reverse=True):
    print w,dict_cnt[w]

#random.seed(0)
#som = Orange.projection.som.SOMLearner(map_shape=(10, 10),initialize=Orange.projection.som.InitializeRandom,epochs=4000)
#som_map = som(Orange.data.Table("fe.tab"))
#k = getclusternumber(som_map,242)
#print(k)
#feruloyl_esterases.faa.AACOMP.SPACED_VALS.vec





