#!/usr/bin/env python
#-*- coding: utf-8 -*-

from pysam import AlignmentFile
import argparse
import os,sys,inspect

currentdir = os.path.dirname(os.path.abspath(inspect.getfile(inspect.currentframe())))
parentdir = os.path.dirname(os.path.dirname(currentdir))
sys.path.insert(0,parentdir)

# self import
from trackcluster import convert

parser=argparse.ArgumentParser()
parser.add_argument("-b", "--bamfile",
                    help="the sorted and indexed bam file")
parser.add_argument("-o", "--out", default="bigg.bed",
                    help="the output file name")
parser.add_argument("-s", "--score", default=21,
                    help="The min mapq score used to keep a read")

args = parser.parse_args()

# make a file using the functions
samfile=AlignmentFile(args.bamfile)

fw=open(args.out, "w")

for n, record in enumerate(samfile):
    # add mapq filter to rm the secondary and supplementary mapping
    if record.mapq>=args.score:
        try:
            bigg=convert.sam_to_bigGenePred(record, samfile)
            fw.write(bigg.to_str())
            fw.write("\n")
        except ValueError:
            pass
    #if n>100:
        #break

fw.close()
samfile.close()
