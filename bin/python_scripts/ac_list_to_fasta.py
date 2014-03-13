'''
@author: Tyler Weirick
@date:   3/13/2013
This program accepts a flat file containg some type of fasta ids
such as an AC number and prints fasta from a fasta file which contain 
one of the ids. Prints matches to standard output.  
'''
#The hexversion is a 32-bit number with the following layout:
#Bits (big endian order)        Meaning
#1-8    PY_MAJOR_VERSION (the 2 in 2.1.0a3)
#9-16   PY_MINOR_VERSION (the 1 in 2.1.0a3)
#17-24  PY_MICRO_VERSION (the 0 in 2.1.0a3)
#25-28  PY_RELEASE_LEVEL (0xA for alpha, 0xB for beta, 0xC for release candidate and 0xF for final)
#29-32  PY_RELEASE_SERIAL (the 3 in 2.1.0a3, zero for final releases)
#Thus 2.1.0a3 is hexversion 0x020100a3.

from sys import hexversion
from glob import glob
py_version = hex(hexversion)

#Get the top comments.
desc=__doc__

if py_version > "0x30200f0":
    import argparse
    def getargs(ver='%prog 0.0'):
        '''
        This function handles the command line arguments. 
        '''
        parser = argparse.ArgumentParser(
            #Get the head comments.
            description=desc,
            formatter_class=argparse.RawDescriptionHelpFormatter)
        parser.add_argument('--fasta_file',
            help='''A fasta file containing the IDs you wish to extract.''')
        parser.add_argument('--flat_file',
            help='''A text file with one ID per line.''')
        args = parser.parse_args()
        return args.fasta_file,args.flat_file


fasta_file,flat_file = getargs()

ac_set = set()
for line in open(flat_file,'r'):
    ac = line.strip()
    if not ac in ac_set:
        ac_set.add(ac)

print_line = False
for line in open(fasta_file,'r'):
    if line [0] == ">":
        if line.split("|")[1] in ac_set:
            print_line = True
        else:
            print_line = False 
    if print_line:
        print(line.strip())



