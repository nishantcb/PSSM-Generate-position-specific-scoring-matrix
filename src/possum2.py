#! /usr/bin/env python
# -*- coding: utf-8 -*-
# Authors: Chris (chris@nohup.cc) & Young (young@nohup.cc)

import sys
import getopt
import logging
from os import listdir
from os.path import isfile, join
import numpy as np
from possum_ft import *

# Configure logging
logging.basicConfig(level=logging.INFO, format='%(asctime)s - %(levelname)s - %(message)s')

def usage():
    print("possum.py usage:")
    print("python possum.py <options> <source files>")
    print("-i,--input: input a file in fasta format.")
    print("-o,--output: output a file of the generated feature.")
    print("-t,--type: specify a feature encoding algorithm.")
    print("-p,--pssmdir: specify the directory of PSSM files.")
    print("-h,--help: show the help information.")

# Parse command-line arguments
try:
    opts, args = getopt.getopt(sys.argv[1:], 'i:o:t:p:a:b:h', 
                                ['input=', 'output=', 'type=', 'pssmdir=', 'argument=', 'veriable=', 'help='])
except getopt.GetoptError as err:
    print(str(err))
    usage()
    sys.exit(2)

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

if not inputFile or not outputFile or not algoType or not pssmdir:
    print("Missing required arguments!")
    usage()
    sys.exit(2)

# Read sequences from input FASTA file
check_head = re.compile(r'>')
smplist = []
try:
    with open(inputFile, 'r', encoding='utf-8', errors='replace') as f:
        for line in f:
            if not check_head.match(line):
                smplist.append(line.strip())
except Exception as e:
    logging.error(f"Error reading input file {inputFile}: {e}")
    sys.exit(1)

# Get all PSSM files from the directory
try:
    onlyfiles = [f for f in listdir(pssmdir) if isfile(join(pssmdir, f))]
except Exception as e:
    logging.error(f"Error accessing PSSM directory {pssmdir}: {e}")
    sys.exit(1)

# Map sequences to PSSM files
fastaDict = {}
for fi in onlyfiles:
    try:
        with open(join(pssmdir, fi), 'r', encoding='utf-8', errors='replace') as f:
            pssmContentMatrix = readToMatrix(f)
            pssmContentMatrix = np.array(pssmContentMatrix)
            sequence = pssmContentMatrix[:, 0]
            cntnt = ''.join(sequence)
            if cntnt not in fastaDict:
                fastaDict[cntnt] = fi
    except Exception as e:
        logging.error(f"Error processing file {fi}: {e}")
        continue

# Match input sequences with PSSM files
finalist = []
for smp in smplist:
    if smp in fastaDict:
        finalist.append(join(pssmdir, fastaDict[smp]))
    else:
        logging.warning(f"Sequence not found in PSSM directory: {smp}")

# Generate feature vectors and save to output file
try:
    with open(outputFile, 'w') as file_out:
        for fi in finalist:
            try:
                with open(fi, 'r', encoding='utf-8', errors='replace') as input_matrix:
                    feature_vector = calculateDescriptors(input_matrix, algoType, argument, veriable)
                    np.savetxt(file_out, feature_vector, delimiter=",")
            except Exception as e:
                logging.error(f"Error processing file {fi}: {e}")
except Exception as e:
    logging.error(f"Error writing to output file {outputFile}: {e}")
    sys.exit(1)

logging.info(f"Feature extraction completed successfully. Output saved to {outputFile}.")
