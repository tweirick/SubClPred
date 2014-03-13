'''
The purpose of this program is to explore the relationships amongst various 
clustering instances. It take cluster assignments in the following format. 

cluster-0        Q12541 Q70KY3 P06811 Q12542 Q02081 Q9Y780 J9PBQ8 J9PBR2
cluster-1        I1SB14 Q5I7J0 H8ZRU2 Q12729 Q96WM9 Q12717 Q9UVY4 Q6A1A1 Q50JG4 
cluster-2        Q72HW2 G8A520 G8A555 D2KZ04 Q9VX11 A1Z6F6  
'''

import argparse
from   glob import glob
import numpy as np


parser = argparse.ArgumentParser(description=__doc__)
parser.add_argument('--file_set',
                    required=True,
                    help='''A file or regex must be in format show in doc string.''' ) 
parser.add_argument('--output_matrix',
                    required=True,
                    help='''Name and path to output matrix to.''' )

args          = parser.parse_args()
file_glob     = glob(args.file_set)
output_matrix = args.output_matrix
name_2d_dict  = {}
#max_cnt       = len(file_glob)


max_cnt=0
for file_name in file_glob: 
   
    
    for line in open(file_name.replace(".cluster-members.txt",".stat.txt"),'r'):
        cluster_number = int( line.split()[2] )

    if cluster_number == 12:
        max_cnt+=1
        for line in open(file_name,'r'):
            cluster_members = line.strip().split("\t")[-1].split()    
          
            while len(cluster_members) > 0:
               match_against = cluster_members.pop(0)

               if match_against in name_2d_dict:
                   name_2d_dict[match_against][match_against]+=1
               else: 
                   name_2d_dict.update({match_against:{match_against:1}})
               for ac_el in cluster_members:
                   #if not ac_el in name_2d_dict:
                   #    name_2d_dict.update({ac_el:{match_against:1}})
                   if ac_el in name_2d_dict[match_against]:
                       name_2d_dict[match_against][ac_el]+=1
                   else: 
                       name_2d_dict[match_against].update({ac_el:1})
                   if not ac_el in name_2d_dict:
                       #name_2d_dict[ac_el][ac_el]+=1
                       #else:
                       name_2d_dict.update({ac_el: {ac_el:0} }) 

                   if match_against in name_2d_dict[ac_el]:
                      name_2d_dict[ac_el][match_against]+=1
                   else: 
                      name_2d_dict[ac_el].update({match_against:1})
               

#Make a matrix from the dict
matrix_out_list = []
#print("x",end="\t")
title_line_list = ["x"]
for ac_id in sorted(name_2d_dict):
    title_line_list.append(ac_id)
matrix_out_list.append( "\t".join(title_line_list) )

numb_always_together = 0
numb_never_together  = 0
numb_above_median    = 0 
numb_below_median    = 0
top_75_percent       = 0
top_90_precent       = 0
top_65_percent       = 0
top_50_percent       = 0

num_matrix = []
for ac_id in sorted(name_2d_dict):
    out_line = [ac_id]
    out_numbs = []
    for ac_id_1 in sorted(name_2d_dict):
        if ac_id_1 in name_2d_dict[ac_id]:
            cell_cnt = name_2d_dict[ac_id][ac_id_1]
            if ac_id == ac_id_1: 
                 #assert cell_cnt == max_cnt
                 numb_always_together+=1
            else: 
                 #assert cell_cnt <= max_cnt
                 if cell_cnt == max_cnt: numb_always_together+=1
                 if cell_cnt > float(max_cnt)*0.9:  top_90_precent+=1                  
                 if cell_cnt > float(max_cnt)*0.75: top_75_percent+=1
                 if cell_cnt > float(max_cnt)*0.65: top_65_percent+=1
                 if cell_cnt > float(max_cnt)*0.50: top_50_percent+=1
        else: 
            cell_cnt = 0
            numb_never_together+=1
        out_line.append(  str(cell_cnt)  )
        out_numbs.append( int(cell_cnt)  )
        
    matrix_out_list.append( "\t".join(out_line) )
    num_matrix.append(out_numbs)

id_list = sorted(name_2d_dict)
np_matrix = np.array(num_matrix)

matrix_out_file = open(output_matrix,'w')
matrix_out_file.write( "\n".join(matrix_out_list) )
matrix_out_file.close()


print("total_els",len(name_2d_dict)*len(name_2d_dict))
print("numb_always_together",numb_always_together)
print("numb_never_together",numb_never_together)
print("top_90_precent",top_90_precent)
print("top_75_percent",top_75_percent)
print("top_65_percent",top_65_percent)
print("top_50_percent",top_50_percent)
print(len(name_2d_dict))
print("Precent no overlap",(numb_always_together+numb_never_together)/float( len(name_2d_dict)*len(name_2d_dict)-1)  )



cluster_set_list = []
for i in range(0,len(num_matrix)):
    for j in range(0,i-1):
        if num_matrix[i][j] >= float(max_cnt)*0.75:
            #Add to a cluster set.
            id_i = id_list[i]
            id_j = id_list[j]
            
            cluster_found = False 
            for c_i in range(0,len(cluster_set_list)):
                if id_i in cluster_set_list[c_i] or id_j in cluster_set_list[c_i]:
                    cluster_found = True
                    cluster_set_list[c_i].add(id_i)
                    cluster_set_list[c_i].add(id_j) 
            if not cluster_found:
                cluster_set_list.append( set([id_i,id_j] ))
 
print(cluster_set_list)
print("len(cluster_set_list)",len(cluster_set_list))
el_len = 0
el_cnt = 0
cnt_dict = {}
for el in cluster_set_list:
    print(len(el),sorted(el))
    el_len+=len(el)
    for id_el in sorted(el):
         if id_el in cnt_dict:cnt_dict[id_el]+=1
         else: cnt_dict[id_el]=1
for el in sorted(cnt_dict,key=cnt_dict.get):
    print(el,cnt_dict[el])
print(len(cluster_set_list))
print(el_len)

exit()


new_set_list = []
cnt=0
no_joins = True
#while no_joins:
no_joins = False
for i in range(0,len(cluster_set_list)):
    new_set = set()
    for j in range(0,len(cluster_set_list)):
        is_same   = cluster_set_list[i] ^ cluster_set_list[j] != set()
        is_subset = cluster_set_list[i].issubset(cluster_set_list[j] ) or cluster_set_list[j].issubset( cluster_set_list[i] )
        if i != j and is_same or is_subset :    
            no_joins = True
            new_set = new_set.union(cluster_set_list[i])
            new_set = new_set.union(cluster_set_list[j])
            cnt+=1    

    if new_set == set():
        new_set_list.append( cluster_set_list[i] )
    else:
        new_set_list.append(new_set)
        #for set_el in new_set_list: 
        #    if set_el - new_set == set():
        #         in_set = True
        #if not new_set_list:         
        #    new_set_list.append(new_set)


for i in range(0,len(new_set_list)):
    non_redundant = []
    for j in range(0,len(new_set_list)):
        if new_set_list[i] == new_set_list[j] and i !=j :
            print("!!!!")
        else: 
            non_redundant.append( new_set_list[j]  ) 
    break

# cluster_set_list = new_set_list
new_set_list = non_redundant


el_len = 0
print("-----------")
cnt_dict = {}
for el in new_set_list:
    print(len(el),sorted(el))
    el_len+=len(el)
    for id_el in sorted(el):
         if id_el in cnt_dict:cnt_dict[id_el]+=1
         else: cnt_dict[id_el]=1
print("++++++++++++")
print("len(cluster_set_list)", len(cluster_set_list) )
print("len(new_set_list)"    , len(new_set_list))
print("cnt",cnt)
print("el_len",el_len)


#for el in sorted(cnt_dict,key=cnt_dict.get):
#    print(el,cnt_dict[el])

