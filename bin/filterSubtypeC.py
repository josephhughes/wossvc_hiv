#! /usr/bin/env python3


import argparse
from Bio import SeqIO

# read in user input

parser = argparse.ArgumentParser()
parser.add_argument('--fasta', required=True, help='Required fasta file containing all samples')
parser.add_argument('--subtypes', required=True, help='tab delimited text file with stanford HIV DB subtype for each sample')
parser.add_argument('--db', required=True, help='fasta file containing existing seqIDs and seqs added to db')
parser.add_argument('--out',required=True, help='output fasta file with the existing db sequences and new sequences appended')
args=parser.parse_args()

# fasta file of new samples
fastaIn = ""
# subtypes of new samples
subtypeFile = ""
#existing db file to be appended to
db=""

if args.fasta:
    fastaIn = args.fasta
else:
    print ("Fasta file required!")
if args.subtypes:
    subtypeFile = args.subtypes
else:
    print ("Text file containing subtypes required!")
if args.db:
    db = args.db
else:
    print("File required - fasta file containing db")


# fill dictionary with subtype (key) to list of sample headers (value)
subtype2sample = {}
subtypesIn = open(subtypeFile, "r")

for line in subtypesIn:
    sample_info = line.split("\t")
    sample = sample_info[0]
    subtype = sample_info[1].split("(")[0].rstrip()
    if subtype in subtype2sample.keys():
        subtype2sample[subtype].append(sample)
    else:
        subtype2sample[subtype]=[sample]

subtypesIn.close()


# open existing db and add the new subtype C samples

f1 = open(args.out, "w")
with open(db) as f:
    for line in f:
        f1.write(line)
            
            
for key, value in sorted(subtype2sample.items()):
    if key =='C':
        fileIn=open(fastaIn, "r")
        for record in SeqIO.parse(fileIn, 'fasta'):
            for item in value:
                if item.strip() == record.id:
                    f1.write(">" + record.id + "\n")
                    f1.write(str(record.seq) + "\n")

