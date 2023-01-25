#!/usr/bin/env python
# -*- coding: utf-8 -*-
import argparse
import os,sys,inspect

currentdir = os.path.dirname(os.path.abspath(inspect.getfile(inspect.currentframe())))
parentdir = os.path.dirname(os.path.dirname(currentdir))
sys.path.insert(0,parentdir)

from trackcluster.tracklist import read_bigg, write_bigg

parser=argparse.ArgumentParser()
parser.add_argument("-b", "--biggfile",
                    help="the bigg bed file")
parser.add_argument("-o", "--out", default="mapped_len.txt",
                    help="the output file name")

args = parser.parse_args(args=None if sys.argv[1:] else ['--help'])

# make a file using the functions
outfile=args.out

bigg_l=read_bigg(args.biggfile)
with open (outfile, "w") as fw:
    for bigg in bigg_l:
        bigg.get_exon()
        fw.write(bigg.name+"\t"+str(bigg.exonlen)+"\t"+bigg.geneName+"\t"+bigg.ttype+"\n")


