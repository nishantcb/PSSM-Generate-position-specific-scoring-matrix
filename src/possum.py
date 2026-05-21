#! /usr/bin/env python
# -*- coding: utf-8 -*-
# Authors: Chris (chris@nohup.cc) & Young (young@nohup.cc)

def usage():
    print("possum.py usage:")
    print("python possum.py <options> <source files> ")
    print("-i,--input: input a file in fasta format.")
    print("-o,--output: output a file of the generated feature.")
    print("-t,--type: specify a feature encoding algorithm.")
    print("-p,--pssmdir: specify the directory of pssm files.")
    print("-h,--help: show the help information.")

import fileinput
import sys, getopt
from os import listdir
from os.path import isfile, join
import re
import numpy as np
from possum_ft import *

# Parse command-line arguments
opts, args = getopt.getopt(sys.argv[1:], 'i:o:t:p:a:b:h', ['input=', 'output=', 'type=', 'pssmdir=', 'argument=', 'veriable=', 'help='])
inputFile = ""
outputFile = ""
algoType = ""
pssmdir = ""
argument = ""
veriable = ""

for opt, arg in opts:
    if opt in ('-i', '--input'):
        inputFile = arg
    elif opt in ('-o', '--output'):
        outputFile = arg
    elif opt in ('-t', '--type'):
        algoType = arg
    elif opt in ('-p', '--pssmdir'):
        pssmdir = arg
    elif opt in ('-a', '--argument'):
        argument = int(arg)
    elif opt in ('-b', '--veriable'):
        veriable = int(arg)
    elif opt in ('-h', '--help'):
        usage()
        sys.exit(2)
    else:
        usage()
        sys.exit(2)

check_head = re.compile(r'>')

smplist = []
smpcnt = 0
for line, strin in enumerate(fileinput.input(inputFile)):
    if not check_head.match(strin):
        smplist.append(strin.strip())
        smpcnt += 1

# Get all PSSM files from the directory
onlyfiles = [f for f in listdir(pssmdir) if isfile(join(pssmdir, f))]

fastaDict = {}

for fi in onlyfiles:
    cntnt = ''
    pssmContentMatrix = readToMatrix(fileinput.input(pssmdir + '/' + fi))
    pssmContentMatrix = np.array(pssmContentMatrix)
    sequence = pssmContentMatrix[:, 0]
    seqLength = len(sequence)
    for i in range(seqLength):
        cntnt += sequence[i]
    if cntnt in fastaDict:
        continue
    fastaDict[cntnt] = fi

# Create the list of files corresponding to the sequences in the input file
finalist = []
for smp in smplist:
    finalist.append(pssmdir + '/' + fastaDict[smp])

# Open the output file
with open(outputFile, 'w') as file_out:
    for fi in finalist:
        # Pass in the original matrix
        input_matrix = fileinput.input(fi)
        # Output a feature vector
        feature_vector = calculateDescriptors(input_matrix, algoType, argument, veriable)
        np.savetxt(file_out, feature_vector, delimiter=",")

