#!/usr/bin/env python
#-*- coding: utf-8 -*-
import argparse
import os,sys,inspect

currentdir = os.path.dirname(os.path.abspath(inspect.getfile(inspect.currentframe())))
parentdir = os.path.dirname(os.path.dirname(currentdir))
sys.path.insert(0,parentdir)

from trackcluster import convert

parser=argparse.ArgumentParser()
parser.add_argument("-i", "--gff",
                    help="the sorted gff file")
parser.add_argument("-o", "--out", default="bigg.bed",
                    help="the output bigGenePred file name")

parser.add_argument("-k", "--key", default="Name",
                    help="The key used as gene name in gff line, like Name=let-1")

args = parser.parse_args()

# make a file using the functions

fw=open(args.out, "w")

gff = convert.GFF(args.gff, args.key)
bigg_list = convert.gff_to_bigGenePred(gff)

for bigg in bigg_list:
    fw.write(bigg.to_str())
    fw.write("\n")


fw.close()
