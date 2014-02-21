import fileinput, random, argparse
from   random import random 
import numpy as np
import pylab as pl
import Orange
from   sklearn.cluster import KMeans
from   sklearn.decomposition import PCA
from   sklearn.preprocessing import scale
from   time import time
#np.random.seed(seed+10)
#random.seed(seed)
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
                       required=True,
                       help='''A vector file in EL_DESC_VAL format.''' )   

    parser.add_argument('--som_x',
                       required=True,
                       help='''Number of neurons in the x direction ex: 5 . ''')

    parser.add_argument('--som_y',
                       required=True,
                       help='''Number of neurons in the y direction ex: 5 .''')
    
    parser.add_argument('--out_dir',
                       default="",
                       help='''
                       If true will generate a graph of a PCA on the k-means 
                       results.''')
    
    args = parser.parse_args()

    out_dir = args.out_dir
    in_file  = args.input_file
    x        = int(args.som_x)
    y        = int(args.som_y)

    if make_pca.lower() == "f" or make_pca.lower() == "false":
        make_pca = False        
    elif make_pca.lower() == "t" or make_pca.lower() == "true":
        make_pca = True
    else:        
        print("error must be bool")
        exit()
    return in_file,x,y,out_dir


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
    Zi = np.sqrt( np.einsum('ij,ij->i', CiAi,CiAi )  )**float(q) 
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
  
 
def getclusternumber(orange_map,vec_len):
    """

    """
    neuron_centriods = []
    neuron_clusters  = []
    #Convert orange map to numpy arrays 
    #Need an array of clusters. 
    #For the dendrogram approach.
    vector_clusters = []
    cluster_means  = []
    for neuron in orange_map:
        vector_cluster = []
        if len(neuron.instances) > 0:
            for neuron_el in neuron.instances:
                vector_cluster.append( neuron_el.native(0) )
            vector_clusters.append( np.array(vector_cluster) )
            cluster_means.append( np.mean( np.array(vector_cluster),axis=0  ) )
             
    np_vector_clusters = np.array( vector_clusters  )
    cluster_centroids  = np.array(cluster_means) 
      
    while True:
        distance_dict = {} 
        initial_db = daviesbouldin( np_vector_clusters  ) 
        for i in range(0,len(np_vector_clusters)): 
            for j in range(0, i-1  ):            
                prototype_clusters = np.copy( np_vector_clusters  )             
                combined_cluster      =  np.concatenate( (np_vector_clusters[i],np_vector_clusters[j]) )
                prototype_clusters[j] = combined_cluster
                dbindex               = daviesbouldin( np.delete( prototype_clusters, i, axis=0 )  )
                distance_dict.update( { str(i)+":"+str(j): dbindex } )
           
        for min_val_key in sorted(distance_dict, key=distance_dict.get):
           if distance_dict[min_val_key] >= initial_db:
               print("SOMDB", distance_dict[min_val_key])
               return np_vector_clusters
           i,j = min_val_key.split(":")
           i = int(i)
           j = int(j)
           combined_cluster      =  np.concatenate((np_vector_clusters[i],np_vector_clusters[j]) )
           np_vector_clusters[j] = combined_cluster
           np_vector_clusters    =  np.delete( np_vector_clusters, i, axis=0 )
           break 


def dokmeans(np_vecs,in_file,k_clusters,domain_set,out_dir):

    #This is to keep runs over-writing each others names. 
    #The time will prevent this if running jobs linerarly, 
    #However, running them in parrallel while the change is very
    #small could overlap. Thus the random number at the end. 
     
    tstr = str(time())+"."+str(random.random()) 

    if out_dir != "":
        in_file = in_file.split("/")[-1]

    k       = len(k_clusters)
    vec_len = len(k_clusters[0][0])
    init_vecs = np.ndarray( shape=(k,vec_len)  )
    for cl_i in range(len(k_clusters)):
        cluster_centroid = np.mean( k_clusters[cl_i]  ,axis=0 ) 
        for cl_j in range(len(cluster_centroid)):
            init_vecs[cl_i][cl_j] = cluster_centroid[cl_j]        
    km              = KMeans(n_clusters=k,max_iter=10000,init=init_vecs)
    kmf             = km.fit(np_vecs)
      
    dict_cnt        = {}
    membership_dict = {}
    #Convert into other data configs 
    tsv_cluster_assignments     = []
    for i in range( len(km.labels_) ):
        ln = str(km.labels_[i])
        tsv_cluster_assignments.append( line_ids[i]+"\t"+ln )
        if ln in dict_cnt:   
            dict_cnt[ln].append(         np_vecs[i]  )
            membership_dict[ln].append( line_ids[i] )
        else:
            dict_cnt[ln]       =[np_vecs[i] ]
            membership_dict[ln]=[line_ids[i]]
     
    #Write tsv file. 
    tsv_file = open(out_dir+in_file+"."+tstr+".clstr.tsv",'w')
    tsv_file.write( "\n".join(tsv_cluster_assignments) )
    tsv_file.close()    
    
    fu = []
    for e in dict_cnt:
        fu.append(dict_cnt[e])
    fu = np.array(fu)
    db_index = daviesbouldin(  fu  ) 
    in_file = in_file.split("/")[-1]
    intercluster_lengths = [] 
    
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

    total_centroid     = np.mean( np.array( intercluster_lengths  ),axis=0 )
    el_diffs           = total_centroid - np.array( intercluster_lengths )
    sqrd_el_diff       = (el_diffs)**2
    inter_cluster_dist = np.sqrt( np.sum( sqrd_el_diff,axis=0  ) )
    intra_cluster_var  = np.var( np.array(intra_cluster_dist) )
    inter_cluster_var  = np.var( np.array(inter_cluster_dist) )

    stat_str = (in_file+"\tnumber_of_clusters\t"+str(len(membership_dict))+
                "\t DBINDEX \t"+str(db_index)+"\t inter_cluster \t"+
                str(inter_cluster_var)+"\t intra_cluster \t"+
                str(intra_cluster_var)+"\n" )

    #Write stats
    stat_file = open(out_dir+in_file+"."+tstr+".stat.txt",'w')
    stat_file.write( stat_str )
    stat_file.close()

    #Write cluster tsv and cluster 
    cluster_lengths             = []
    cluster_and_ids_on_one_line = []
    for w in sorted(membership_dict):
        #print w,membership_dict[w]
        seq_list = []
        for seq_id in membership_dict[w]:
           seq_list.append(seq_id)
        cluster_lengths.append( "cluster-"+str(w)+"\t"+str(len( seq_list  ) ) )
        cluster_and_ids_on_one_line.append("cluster-"+str(w)+"\t "+" ".join(seq_list))
    
    of = open(out_dir+in_file+"."+tstr+".cluster-lens.txt",'w')
    of.write("\n".join(cluster_lengths))
    of.close()
        
    of = open(out_dir+in_file+"."+tstr+".cluster-members.txt",'w')
    of.write("\n".join(cluster_and_ids_on_one_line))
    of.close()

    """    
    reduced_data  = PCA(n_components=2).fit_transform(np_vecs)
    reduced_data  = scale( reduced_data )
    kmeans = KMeans(init='k-means++', n_clusters=k, n_init=k)
    kmeans.fit(reduced_data)
    # Step size of the mesh. Decrease to increase the quality of the VQ.
    h = .02     # point in the mesh [x_min, m_max]x[y_min, y_max].
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
    title_str = 'K-means clustering of (PCA-reduced data) from\n'+in_file+'\nk='+str(k)+' Centroids are marked with white cross'
    pl.title( title_str,fontsize=10  )
    pl.xlim(x_min, x_max)
    pl.ylim(y_min, y_max)
    pl.xticks(())
    pl.yticks(())
    pl.savefig(in_file+"k-means-clusters.png")
    """


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


def runOrangeSOM(domain_set,np_vecs,som_size_x,som_size_y):

    d           = Orange.data.Domain( [Orange.feature.Continuous(x) for x in sorted(domain_set) ] )

    orange_vecs = Orange.data.Table(d,np_vecs)

    som = Orange.projection.som.SOMLearner(
        map_shape=(som_size_x,som_size_y),
        initialize=Orange.projection.som.InitializeRandom,
        epochs=10000)

    som_map = som(orange_vecs)
    
    return som_map

#=============================================================================
#                              Main Program                 
#=============================================================================
#Get command line arguments 
in_file,x,y,out_dir = getargs()
#Get vectors from standard input.
rl,line_ids,domain_set = textvecstolist(in_file)

#Create output file into 
in_file = in_file.split("/")[-1]

#Scale vecs 
np_vecs     = scale( np.array( rl  ) ,with_std=False,with_mean=False)
#Run SOM 
som_map =  runOrangeSOM(domain_set,np_vecs,x,y)
#Find k and centroids
k_clusters = getclusternumber( som_map,len(np_vecs) )
#Do kmeans clustering
dokmeans(np_vecs,in_file,k_clusters,domain_set,out_dir)

