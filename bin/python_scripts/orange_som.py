import fileinput, random, argparse
import numpy as np
from math import pow
from random  import random as rnd
import Orange
from math import sqrt
from sklearn.cluster import KMeans
from sklearn.decomposition import PCA
from sklearn.preprocessing import scale
import pylab as pl

seed = 42
seed = 0
seed = 10
seed = 5
np.random.seed(seed+10)
random.seed(seed)


def getargs():
    parser = argparse.ArgumentParser()
    parser.add_argument('--input_file',
                       help='''A vector file in EL_DESC_VAL format.''',required=True)
    
    parser.add_argument('--print_variance',
                       help='''True will print cluster numbers and variances,
                            false will print the clusters and their IDS ''',
                       default=True)
    args = parser.parse_args()
    print_variance = str(args.print_variance)
    in_file = args.input_file

    if print_variance.lower() == "f" or print_variance.lower() == "false":
        print_variance = False        
    elif print_variance.lower() == "t" or print_variance.lower() == "true":
        print_variance = True
    else:        
        print("error must be bool")
        exit()
    
    return in_file,print_variance