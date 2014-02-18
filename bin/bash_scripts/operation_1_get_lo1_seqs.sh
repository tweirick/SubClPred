#!/bin/bash
#This is the product of some work on uniprot.org see docs for an explination on how this was made.
wget "http://www.uniprot.org/uniprot/?query=%28%0D%0Aname%3a%22laccase%22+OR%0D%0Aname%3a%22urishiol+oxidase%22+OR%0D%0Aname%3a%22urushiol+oxidase%22+OR%0D%0Aname%3a%22p-diphenol+oxidase%22+OR%0D%0Aec%3a1.10.3.2%0D%0A%29+AND+%0D%0A%28existence%3a%22evidence+at+protein+level%22+OR+%0D%0Aexistence%3a%22evidence+at+transcript+level%22%29+%0D%0AAND+fragment%3ano+%0D%0AAND+NOT+name%3a%22Catechol+oxidase%22+%0D%0AAND+NOT+ec%3a1.10.3.1&format=fasta" -O  data/2013-12-10.1_laccase_faa/lac_pt.faa 