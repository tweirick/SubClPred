import argparse
from glob import glob
parser = argparse.ArgumentParser()    
parser.add_argument('--file_set',
                    help='Takes a file name as a string.')
parser.add_argument('--min_chars',
                    help='all seqs under this number will be excluded. ')

args = parser.parse_args()

infile_name = args.file_set

min_chars = int(args.min_chars)
    
for file_name in glob(infile_name):
    file = open(file_name,"r")
    out_list = []
    fasta_name = ""
    larger_than_cnt  = 0
    smaller_than_cnt = 0
    while True:
        line = file.readline()
        #Read through line by line. This is done iteratively to allow for very
        #large files. 
        if line == "" or line[0] == '>': #If the start of a fasta entry enter data 
            
            if fasta_name != "":
                fasta_seq = ''.join(fasta_data)
                if len(fasta_seq) > min_chars:        
                    larger_than_cnt+=1            
                    out_list.append(fasta_name.strip())
                    out_list.append(fasta_seq.strip())
                else:
                    smaller_than_cnt+=1
            
            if line == "":
                break
            
            fasta_data = []
            fasta_name = line.strip()    
            
        else:
            fasta_data.append(line.strip())  
            
    file.close()     
    
    outfile = open(file_name+".greaterthan"+str(min_chars)+"chars.faa","w")
    outfile.write("\n".join(out_list))
    outfile.close()
    
    print(smaller_than_cnt,"removed",larger_than_cnt,"kept.")
    

