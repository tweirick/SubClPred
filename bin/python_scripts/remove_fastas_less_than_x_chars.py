import argparse
from glob import glob
parser = argparse.ArgumentParser()    
parser.add_argument('--file_set',
                    required=True,
                    help='Takes a file name as a string.')

parser.add_argument('--out_dir',
                    default="",
                    help='Takes a file name as a string.')

parser.add_argument('--min_chars',
                    required=True,
                    help='all seqs under this number will be excluded. ')

args = parser.parse_args()

infile_name = args.file_set
out_dir     = args.out_dir 
min_chars   = int(args.min_chars)
    

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
    
    if out_dir != "": 
        file_name = file_name.split("/")[-1]

    outfile = open(out_dir+file_name+".greaterthan"+str(min_chars)+"chars.faa","w")
    outfile.write("\n".join(out_list))
    outfile.close()
    
    outfile = open(out_dir+file_name+".greaterthan"+str(min_chars)+"chars.log","w")
    outfile.write( str(smaller_than_cnt)+" removed "+str(larger_than_cnt)+" kept.\n" )
    outfile.close()
 
