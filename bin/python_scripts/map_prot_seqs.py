
from glob import glob

files_to_map_name = glob("/home/tyler/Desktop/SubClPred/laccasedb/*.faa" )
map_to_file_name  = "data/2013-01-23.1_laccases_seqs_faa/LO1_tp_combined.faa"

map_to_file = open(map_to_file_name,'r')
fasta_name = None
fasta_dict = {}

while True: 
    line = map_to_file.readline()
    if len(line) == 0 or line[0] == ">":
        if fasta_name != None: 
            fasta_txt = "".join(fasta_list)
            assert  not fasta_txt in fasta_dict
            fasta_dict.update( { fasta_txt  :fasta_name} ) 
            assert fasta_txt in fasta_dict
        if len(line) == 0: break
        fasta_list = []
        fasta_name = line.strip()
    else: 
        fasta_list.append( line.strip() )



for file_name in files_to_map_name:
    mfile = open(file_name,'r')
    fasta_name = None
    out_list   = []
    while True:
        line = mfile.readline()
        if len(line) == 0 or line[0] == ">":
            if fasta_name != None:
                fasta_txt = "".join(fasta_list)
                if fasta_txt in fasta_dict:
                    out_list.append(fasta_dict[fasta_txt]+"\n"+fasta_txt)
            if len(line) == 0: break
            fasta_list = []
            fasta_name = line.strip()
        else: 
            fasta_list.append( line.strip() )

    ofile = open(file_name+".mapped-faa",'w')
    ofile.write("\n".join(out_list))
    ofile.close()
