#!/usr/bin/python
'''
@author: Tyler Weirick
@date:   3/13/2013

This program is meant to make the submission of many PBS scripts easier. It 
requires two inputs a template file and a CSV file or in the case of only 
one value to insert into a template a regular expression can be given. 
 
#!/bin/bash
#PBS -q batch
#PBS -l nodes=1:ppn=12
#PBS -l walltime=100:00:00
#PBS -N psiblast.{0}
#PBS -j oe
module load blast+
cd $PBS_O_WORKDIR
date
psiblast -db {1} -query {2} -out {2}.psiblastout -num_descriptions 5000 -evalue 1e-39 -num_iterations 3  -num_threads 12
date
 
'''
import csv
from glob import glob
from time import time
#The hexversion is a 32-bit number with the following layout:
#Bits (big endian order)	Meaning
#1-8	PY_MAJOR_VERSION (the 2 in 2.1.0a3)
#9-16	PY_MINOR_VERSION (the 1 in 2.1.0a3)
#17-24	PY_MICRO_VERSION (the 0 in 2.1.0a3)
#25-28	PY_RELEASE_LEVEL (0xA for alpha, 0xB for beta, 0xC for release candidate and 0xF for final)
#29-32	PY_RELEASE_SERIAL (the 3 in 2.1.0a3, zero for final releases)
#Thus 2.1.0a3 is hexversion 0x020100a3.


from sys import hexversion
from glob import glob
py_version = hex(hexversion)

#Get the top comments.
desc=open(__file__).read().split("'''")[1]

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
        parser.add_argument('--template_file', 
            help='''A template file to make pbs scripts.''')
        parser.add_argument('--file_set',
            help='''
                 Either a csv file for the template separated by 
                 or a file name or regex if there is only one value to add
                 to the template. 
                 ''')
        parser.add_argument('--is_flat_file',
            help='''
                 Boolean accepting T or F. T if you are passing a tab 
                 delimited flat file. F for file name or regex.
                 ''',default=False)
        parser.add_argument('--output_to_directory',
            help='''
                 The directory to output the qsub files to.
                 ''',default=None)

        parser.add_argument('--shell_script',
            help='''
                 name and location of where to output shell script 
                 if blank will be (template_name).submit.sh
                 ''',default=None)

        args = parser.parse_args()
        return args.template_file,args.file_set,args.is_flat_file,args.output_to_directory,args.shell_script



#=============================================================================
#                         Main Program 
#=============================================================================


temp_file,file_set,is_flat_file,out_dir,qss =  getargs(ver='%prog 0.0')

#This will be a shell script that the arguments are passed to. 
#It will contain the commands to run all of the qsub commands.
#I tried doing by printing and pasting them but things would
#somewhat randomly get skipped. Plus it something goes wrong
#it is easy to resubmit with a script. 
qsub_submission_script_name = qss
#Location and name of the template file. 
template_file_path_and_name     = temp_file
#To use for first part of qsub file names 
base_name              = template_file_path_and_name.split("/")[-1]
#File extension to use for
file_extension         = ".qsub"
sh_txt_list            = ["#!/bin/bash"]


#Get template, required so throw error and exit it there is a problem. 
template_file = open(template_file_path_and_name,'r')
template_txt  = template_file.read()
template_file.close()

file_number = 0
is_csv = False

if type(is_flat_file) != bool:
    
    if is_flat_file.upper() == "T" or  is_flat_file.upper() == "TRUE":
        is_flat_file = True
    elif is_flat_file.upper() == "F" or  is_flat_file.upper() == "FALSE":
        is_flat_file = False
    else:
        print("ERROR: "+str(is_flat_file)+" not recognized as boolean value.")

base_qsub = "sleep 0.1\nqsub "

if is_flat_file:
    print(is_flat_file)
    is_csv = False
    print("is_csv = False")
    if is_csv:
        for csv_file_name in glob(file_set):
            with open(csv_file,'r') as csvfile:
                csv_iter = csv.reader(csvfile, delimiter=',', quotechar='"')
                for row in csv_iter:
                    #Convert iterator to a tuple
                    #this is required for the format function
                    row_tuple = tuple(row)
                    qsub_file_name = out_dir+base_name+str(time())+file_extension
                    qsub_file = open(qsub_file_name,'w')
                    qsub_file.write( template_txt.format( *row_tuple  ) )
                    qsub_file.close()
                    #Make text for submission of the qsub file. These will be written 
                    #to a bash shell script for easy submission.  
                    sh_txt_list.append(base_qsub+out_dir+"/"+qsub_file_name)
                    file_number+=1
    else:
         for tab_file_name in glob(file_set):
             for line in open(tab_file_name,'r'):
                 #sp_line = line.strip().split(",")
                 
                 if len(line) > 1:
                     qsub_file_name = out_dir+tab_file_name+"-"+base_name+"."+str(time())+file_extension
                     qsub_file = open(qsub_file_name,'w')
                     qsub_file.write( template_txt.format(line)   )
                     qsub_file.close()
                     sh_txt_list.append(base_qsub+qsub_file_name)
                     file_number+=1
else:          
    #todo: Need check for only one insert area. 
    for file_name in glob(file_set):
        qsub_file_name = out_dir+base_name+"."+str(time())+file_extension 
        qsub_file = open(qsub_file_name,'w')
        qsub_file.write( template_txt.format(file_name) )
        qsub_file.close()
        sh_txt_list.append(base_qsub+qsub_file_name)
        file_number+=1

#Make the mass submission shell script
print(qsub_submission_script_name)
script_file = open(qsub_submission_script_name,"w")
script_file.write( "\n".join(sh_txt_list) )
script_file.close()


