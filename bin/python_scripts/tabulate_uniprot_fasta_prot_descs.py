import argparse
from glob import glob
parser = argparse.ArgumentParser()
parser.add_argument('--file_set',
                    help='Takes a file name or regex as a string.')
args = parser.parse_args()
file_set = glob(args.file_set)

for file_name in file_set: 
    for line in open(file_name,'r'):
        if len(line) > 0 and line[0] == ">":
            print( line.split()[0]+"\t"+" ".join(line.split("OS=")[0].split()[1:])  )

