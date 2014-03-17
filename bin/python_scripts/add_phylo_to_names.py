

import string 

tree_file = "la.tree"
def_file  = "id_specs.txt"



id_dict = {}
for line in open(def_file,'r'):
    uid,specs = line.split()
    #print(uid)
    id_dict.update({uid:specs})

for line in open(tree_file,"r"):
    if line[0] in string.ascii_letters:
        length  = line.split(":")[-1]
        print(line.split(":")[0].split("|")[1]+""+id_dict[line.split(":")[0].split("|")[-1]]+":"+length.strip() ) 
    else:
        print(line.strip())

