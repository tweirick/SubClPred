
import argparse
from glob import glob
parser = argparse.ArgumentParser()

parser.add_argument('--cluster_file',
                    required=True,
                    help='Takes a file name')

parser.add_argument('--fasta_file',
                    default="",
                    required=True,
                    help='Takes a file name as a string.')

parser.add_argument('--out_dir',
                    default="",
                    help='Takes a file name as a string.')

args            = parser.parse_args()
fasta_file_name = args.fasta_file
ac_file_name    = args.cluster_file
out_dir         = args.out_dir


#out_dir  = args.out_dir
one_ac_per_line = False
fasta_dict = {}
fasta_name = None
fasta_file = open(fasta_file_name,'r')
while True:
    line = fasta_file.readline()
    if len(line) == 0 or line[0] == ">":
        if fasta_name != None:
            fasta_sequence = "".join(fasta_data) 
            fasta_dict.update( {fasta_name:fasta_sequence} )
        if len(line) == 0:break
        fasta_sequence = ""          #Redundancy, just in case.
        fasta_name     = line.split("|")[1]#Get new name
        fasta_data     = []          #Reset fasta data.
    else:
        fasta_data.append(line.strip())
fasta_file.close()


out_file_name = out_dir+ac_file_name.split("/")[-1]

print(out_file_name)

if one_ac_per_line: 
    #@todo: change this to handle 0\tac_number
    cnt=0
    for line in open(ac_file_name,'r'):    
        out_str_list = []
        for id_el in line.split()[1:]: 
            out_str_list.append( ">"+id_el+"\n"+fasta_dict[id_el]  )
        f = open(out_file_name+"."+str(cnt)+".faa",'w')
        f.write("\n".join(out_str_list))
        f.close()
        cnt+=1   
else: 
    cnt=0
    for line in open(ac_file_name,'r'):    
        out_str_list = []
        for id_el in line.split()[1:]: 
            out_str_list.append( ">"+id_el+"\n"+fasta_dict[id_el]  )
 
        
        f = open(out_file_name+".cluster-"+str(cnt)+".faa",'w')
        f.write("\n".join(out_str_list))
        f.close()
        cnt+=1   
