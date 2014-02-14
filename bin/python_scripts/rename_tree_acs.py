'''
The program renames the names on a phylogenic tree. Using a tsv file of substitutions. 
'''
import argparse

parser = argparse.ArgumentParser(description=__doc__)

parser.add_argument(
            '--rename_file',
            required=True,
            type=str,
            help='A tsv file with current ids in the first column and new names in the second column.')

parser.add_argument(
            '--tree_file',
            required=True,
            type=str,
            help='A tree file in newark format.')

args = parser.parse_args()
rename_file = args.rename_file
file_name   = args.tree_file


rename_dict = {}
for line in open(rename_file,'r'):
    current_id,new_id = line.split()
    rename_dict[current_id]=new_id 


for line in open(file_name,'r'):  

    if not line[0] in "(),\n;":

        if "|" in line:
            line_id = line.split("|")[1] 
        else: 
            line_id = line.strip()

        if  line_id in rename_dict:
            print( rename_dict[ line_id ]  )  
        else: 
            print("ERROR,id",line_id,"not found. Exiting.")
            exit()
    else: 
       print(line.strip())        
