'''

'''
from glob import glob
import argparse

def getargs(ver='%prog 0.0'):

    parser = argparse.ArgumentParser(description='Take a blastclust cluster '+
             'file as input and remove all clusters from fasta file set.')    
    parser.add_argument('--file_set', 
                        help='')
    
    parser.add_argument('--file_suffix', 
                        default=".pruneBJXZU.fasta",
                        help='')
    
    parser.add_argument('--remove_over_x_percent', 
                        default=0.0,
                        help='')

    args = parser.parse_args()
    
    return glob(args.file_set),args.file_suffix,args.remove_over_x_percent


def build_fasta_dict(input_file_name):
    fasta_seq_name_list = []
    temp_fasta_seq      = []
    fasta_name          = None 
    file = open(input_file_name,'r')

    while True:
        line = file.readline()
        if line == "" or line[0] == ">":

            if fasta_name != None:
                fasta_seq_name_list.append([fasta_name,''.join(temp_fasta_seq)])
                temp_fasta_seq = []
            
            if line == "": break
            fasta_name = line.strip()
        else:
            temp_fasta_seq.append(line.strip()) 
    
    return fasta_seq_name_list

    
def buildnewfasta(fasta_name_seq_list,max_percent_illegal_chars):
    """
    This will discard function will splite the new_dict into two lists. One
    with the entries to keep and the other to discard. This is based on 
    whether there is a less than or equal to percentage of illegal chars in
    the sequence than remove_over_x_percent
    """
    illegal_char_set = {"B","J","X","Z","U","O","0"}    

    out_list, discard_list = [],[]
    for fasta_name_seq in fasta_name_seq_list:
        fasta_name,fasta_seq = fasta_name_seq[0],fasta_name_seq[1]
        illegal_char_cnt = 0
        
        #This could be done with the count() function but I have encountered 
        #problems with count() and very large strings. It was returning 
        #incorrect results silently.  
        for aa in fasta_seq:
            if aa.upper() in illegal_char_set:
                illegal_char_cnt+=1.0 
        fasta_entry = "\n".join(fasta_name_seq)
        percent_illegal_chars = illegal_char_cnt/float(len(fasta_seq))
        assert type(percent_illegal_chars) == float
        if percent_illegal_chars <= float(max_percent_illegal_chars):    
            out_list.append(fasta_entry)
        else:
            discard_list.append(fasta_entry)

    return out_list,discard_list



#PRE_PSI-BLAST_52412/*catted_52412.faa
#file_set = glob("/home/TWeirick/FEATURE_GENERATING_PROGRAMS/40per_Sets/40per_all_catted.faa")#

file_set,file_suffix,remove_over_x_percent = getargs(ver='%prog 0.0')
print("File_Name","Fastas_In_Input","Fastas_In_Kept_Output","Files_In_Discarded_Output")
for file_name in file_set:
    
    fasta_list = build_fasta_dict(file_name)
   
    out_list,discard_list = buildnewfasta(fasta_list,remove_over_x_percent)

    print(file_name,
          str(len(fasta_list)),
          str(len(out_list)),
          str(len(discard_list)))

    out_file = open(file_name+file_suffix,'w')
    out_file.write('\n'.join(out_list))
    out_file.close()
    
    out_file = open(file_name+file_suffix+".discard",'w')
    out_file.write('\n'.join(discard_list))
    out_file.close()
