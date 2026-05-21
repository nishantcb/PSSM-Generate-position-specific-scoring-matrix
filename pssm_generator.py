#!/usr/bin/env python3
# -*- coding: utf-8 -*-
"""
Created on Mon Jul 18 11:44:31 2024

@author: pratik & anand
"""

import sys
import glob
import os

dir_path = sys.argv[1]
os.makedirs(dir_path+'/pssm_raw1', exist_ok=True)
os.makedirs(dir_path+'/pssm_raw', exist_ok=True)

listdir=glob.glob(dir_path+'/*.fasta')

for i in listdir:
    filename=i.split('/')[-1].rsplit('.',1)[0]
    cmd = "ncbi-blast-2.16/bin/psiblast -out "+dir_path+"/pssm_raw1/"+filename+".homologs -outfmt 7 -query "+dir_path+"/"+filename+".fasta -db ./swissprot/swissdb -evalue 0.1 -word_size 3 -max_target_seqs 6000 -num_threads 14 -gapopen 11 -gapextend 1 -matrix BLOSUM62 -comp_based_stats D -num_iterations 3 -out_pssm "+dir_path+"/pssm_raw1/"+filename+".cptpssm -out_ascii_pssm "+dir_path+"/pssm_raw/"+filename+".pssm"
    print('\n',str(cmd))
    os.system(cmd)
