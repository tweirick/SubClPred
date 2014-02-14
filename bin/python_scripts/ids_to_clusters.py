

fasta_file_name = "/home/tyler/Desktop/find_enzyme_subclasses/data/2013-01-23.1_laccases_seqs_faa/LO1_tp_combined.faa.pruneBJXZU.fasta.greaterthan61chars.faa"
ac_file_name    = "/home/tyler/Desktop/find_enzyme_subclasses/tyler_5x5n_8_clusters_shenCTD.txt"
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

if one_ac_per_line: 
    #@todo: change this to handle 0\tac_number
    cnt=0
    for line in open(ac_file_name,'r'):    
        out_str_list = []
        for id_el in line.split()[1:]: 
            out_str_list.append( ">"+id_el+"\n"+fasta_dict[id_el]  )
        f = open(file_name+"."+str(cnt)+".faa",'w')
        f.write("\n".join(out_str_list))
        f.close()
        cnt+=1   
else: 
    cnt=0
    for line in open(ac_file_name,'r'):    
        out_str_list = []
        for id_el in line.split()[1:]: 
            out_str_list.append( ">"+id_el+"\n"+fasta_dict[id_el]  )
        f = open(ac_file_name+".cluster-"+str(cnt)+".faa",'w')
        f.write("\n".join(out_str_list))
        f.close()
        cnt+=1   