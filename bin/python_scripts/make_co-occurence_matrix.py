'''
The purpose of this program is to explore the relationships amongst various 
clustering instances. It take cluster assignments in the following format. 

cluster-0        Q12541 Q70KY3 P06811 Q12542 Q02081 Q9Y780 J9PBQ8 J9PBR2
cluster-1        I1SB14 Q5I7J0 H8ZRU2 Q12729 Q96WM9 Q12717 Q9UVY4 Q6A1A1 Q50JG4 
cluster-2        Q72HW2 G8A520 G8A555 D2KZ04 Q9VX11 A1Z6F6  
'''

import argparse
from glob import glob

parser = argparse.ArgumentParser(description=__doc__)
parser.add_argument('--file_set',
                    required=True,
                    help='''A file or regex must be in format show in doc string.''' ) 

args      = parser.parse_args()
file_glob = glob(args.file_set)

name_2d_dict = {}

max_cnt = len(file_glob)

for file_name in file_glob: 
    
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
print("x",end="\t")
for ac_id in sorted(name_2d_dict):
    print(ac_id,end="\t")
print() 

numb_always_together = 0
numb_never_together  = 0
numb_above_median    = 0 
numb_below_median    = 0

for ac_id in sorted(name_2d_dict):
    print(ac_id,end="\t")
    for ac_id_1 in sorted(name_2d_dict):
        if ac_id_1 in name_2d_dict[ac_id]:

            cell_cnt = name_2d_dict[ac_id][ac_id_1]

            if ac_id == ac_id_1: 
                 assert cell_cnt == max_cnt
                 numb_always_together+=1
            else: 
                 assert cell_cnt <= max_cnt
                 if cell_cnt == max_cnt: numb_always_together+=1
            print(name_2d_dict[ac_id][ac_id_1],end="\t")
        else: 
            numb_never_together+=1
            print(0,end="\t")
    print()

print("total_els",len(name_2d_dict)*len(name_2d_dict))
print("numb_always_together",numb_always_together)
print("numb_never_together",numb_never_together)
print("Precent no overlap",(numb_always_together+numb_never_together)/float( len(name_2d_dict)*len(name_2d_dict)-1)  ))

