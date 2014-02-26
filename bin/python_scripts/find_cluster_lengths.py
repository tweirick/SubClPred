'''
'''

from glob import glob

#file_glob = glob("data/*som_to_kmeans/*PLPredPhysChem22_nolog*.stat.txt")
file_glob = glob("data/*som_to_kmeans/*.stat.txt")

out_dict = {}
for file_name in file_glob:
    for line in open(file_name):
        sp_line      = line.split()
        cluster_numb = sp_line[2].strip()
        seq_type     = "_"+file_name.split("lac_pt.faa.pruneBJXZUO.fasta.")[-1].split(".EL_DESC_VAL.vec.")[0]+"_"
        neurons      = file_name.split("/")[-2].split("_clusters_from_")[-1].split("-som_to_kmeans")[0]
        
        if not neurons in out_dict:
            out_dict.update({neurons:{}})
            #print("\n"+neurons+"\n")
        if not seq_type in out_dict[neurons]:
            out_dict[neurons].update({seq_type:{}})
            #print("\n"+seq_type+"\n")
        if not cluster_numb in out_dict[neurons][seq_type]:
            #print(neurons,seq_type,cluster_numb)
            #print(cluster_numb in out_dict[neurons][seq_type])
            out_dict[neurons][seq_type].update({cluster_numb:0})       
            #print(cluster_numb in out_dict[neurons][seq_type]) 
        out_dict[neurons][seq_type][cluster_numb]+=1
        #print(seq_type,neurons,cluster_numb)

title_str = ["x"]        
matrix_list = []
for neuron_type_el in sorted(out_dict): 
    title_str.append(neuron_type_el)
    row_list = [neuron_type_el]
    title_str = ["x"]
    for vec_type in sorted( out_dict[neuron_type_el] ): 
        title_str.append(vec_type)
        cnt_str = ""
        for el_counts in sorted(out_dict[neuron_type_el][vec_type], key=out_dict[neuron_type_el][vec_type].get,reverse=True):
            cnt_str = el_counts
            #cnt_str = cnt_str+el_counts+":"+str(out_dict[neuron_type_el][vec_type][el_counts])+"/"           
            break 
        row_list.append(cnt_str)
    matrix_list.append("\t".join(row_list))  



print("\t".join(title_str))
print( "\n".join(matrix_list)      )
           
